package PrunedMSMEqualLengths;

import HeuristicsEqualLengths.Greedy;
import HeuristicsUnequalLengths.ItakuraUnequalLengths;
import MSMDistances.*;

import java.util.Arrays;
import java.util.Objects;

public class PrunedAdvancedLBU {

    final double c;

    final double[] X;
    final double[] Y;
    double constant;
    double[] remainingCostX;
    double[] remainingCostY;

    double[] ts1;
    double[] ts2;

    double upperBound;
    double[] upperBoundArray;
    final int m;

    final String UB;


    /**
     * compute MSM distance with pruned version -> prune the dynamic programming table
     * with update UB and Advanced LB
     * @param c  cost for merge and split
     * @param UB Upper Bound Strategy: g = greedy, i = Itakura with d = 3/4
     */
    public PrunedAdvancedLBU(double c, double[] X, double[] Y, double constant, String UB) {
        this.c = c;
        this.X = X;
        this.Y = Y;
        this.m = X.length;
        this.constant = constant;
        if (m != Y.length) {
            throw new RuntimeException("Unequal lengths of time series!");

        }

        this.UB = UB;
        this.setRemainingCost();
        this.setUpperBound();

        // first column and row is inf -> need to extend time series
        this.ts1 = new double[m + 1];
        this.ts2 = new double[m + 1];
        ts1[0] = Double.POSITIVE_INFINITY;
        ts2[0] = Double.POSITIVE_INFINITY;

        System.arraycopy(X, 0, ts1, 1, m + 1 - 1);
        System.arraycopy(Y, 0, ts2, 1, m + 1 - 1);

    }

    public void setUpperBound() {
        if (Objects.equals(this.UB, "g")) {
            Greedy greedy = new Greedy(X, Y, c);
            this.upperBoundArray = greedy.computeGreedyArray();
            this.upperBound = upperBoundArray[0] + 0.0000001;
        }
        if (Objects.equals(this.UB, "i")) {
            ItakuraUnequalLengths itakuraUnequalLengths = new ItakuraUnequalLengths(c, 3 / 4.);
            this.upperBoundArray = itakuraUnequalLengths.itakuraHeuristicArray(X, Y);
            this.upperBound = upperBoundArray[0]+ 0.0000001;
        }
    }

    public void setRemainingCost() {
        MsmConstant msmConstantX = new MsmConstant(X, this.c, this.constant);
        MsmConstant msmConstantY = new MsmConstant(Y, this.c, this.constant);

        this.remainingCostX = msmConstantX.computeRemainingCostPerEntry();
        this.remainingCostY = msmConstantY.computeRemainingCostPerEntry();
    }


    /**
     * compute the Msm distance table with pruned entries
     *
     * @return array: first Double: Distance, second Double: relative amount of pruned cells
     */
    public double[] msmDistPruned() {

        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        double[] tmpArray = new double[m + 1];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0

        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule


        double tmp = 0;

        // index for pruning start and end of row
        int sc = 1;
        int ec = 1;


        int countCutCells = 0;


        // remember if an entry smaller than UB was found -> cannot cut
        boolean smallerFound;
        int ecNext;


        //initialize first row
        // the first entry of this row is inf, this inf is used for computing the second row
        Arrays.fill(tmpArray, Double.POSITIVE_INFINITY);

        //row index
        for (int i = 1; i < tmpArray.length; i++) {
            smallerFound = false;
            final double xi = ts1[i];
            // the index for the pruned end cannot be lower than the diagonal
            // All entries on the diagonal have to be equal or smaller than
            // the upper bound (Euclidean distance = diagonal path)
            ecNext = i;

            final double remainingCosti = remainingCostX[i];

            //column index
            for (int j = sc; j < tmpArray.length; j++) {

                final double yj = ts2[j];
                double d1, d2, d3;
                d1 = tmp + Math.abs(xi - yj);
                //merge
                d2 = tmpArray[j] + C(xi, ts1[i - 1], yj);
                //split
                d3 = tmpArray[j - 1] + C(yj, xi, ts2[j - 1]);

                // store old entry before overwriting
                tmp = tmpArray[j];
                tmpArray[j] = Math.min(d1, Math.min(d2, d3));

                // PruningExperiments strategy
                double lb = getLowerBound(i, j, xi, remainingCosti);
                if ((tmpArray[j] + lb) > this.upperBound) {
                    if (!smallerFound) sc = j + 1;
                    if (j > ec) {
                        tmpArray = this.fillWithInf(j + 1, m + 1, tmpArray);
                        break;
                    }
                } else {
                    smallerFound = true;
                    ecNext = j + 1;
                }
                if(i==j) {
                    this.upperBound = tmpArray[j] + this.upperBoundArray[j] +0.00000001;
                }

            }

            tmpArray = this.fillWithInf(1, sc, tmpArray);


            // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
            tmp = Double.POSITIVE_INFINITY;
            ec = ecNext;
            countCutCells = countCutCells + (ts1.length - (ec)) + (sc - 1);
        }

        double totalSizeTable = (double) m * m;
        double percentageCutCells = countCutCells / (totalSizeTable);

        //     System.out.println("table size: " + totalSizeTable + " counter Simple better: " + counterSimpleBetter);
        return new double[]{tmpArray[m], percentageCutCells};
    }

    public double[] fillWithInf(int start, int end, double[] array) {
        for (int i = start; i < end; i++) {
            array[i] = Double.POSITIVE_INFINITY;
        }
        return array;
    }

    public double getLowerBound(int xCoord, int yCoord, double xi, double remainingCosti) {
        double lowerBound;

        // To apply the triangle inequality the time series are cut at the current index. For a correct estimation, the cost for the potential next
        // move operation has to be subtracted (the first elements of two time series are aligned (move)
        double simpleLB = Math.abs(xCoord - yCoord) * c;

        double dist;
        if (xCoord < m && yCoord < m && (xCoord > 1 || yCoord > 1))
            dist = Math.max(Math.abs(xi - this.constant) + Math.abs(ts2[yCoord] - this.constant) - c, 0);
        else if ((xCoord == m) && (yCoord == m)) {
            dist = 0.;
        } else {
            dist = Math.abs(xi - this.constant) + Math.abs(ts2[yCoord] - this.constant);
        }

        lowerBound = Math.max(Math.abs(remainingCosti - remainingCostY[yCoord] - (((m - xCoord) - (m - yCoord))) * this.c) - dist, 0) - 0.000001;

        return Math.max(lowerBound, simpleLB);

    }

    // compute the Euclidean distance of the time series = the diagonal of the table of the DP
    public double[] upperBoundEuclid(double[] X, double[] Y) {
        double[] diagonalBounds = new double[this.m];

        double sumDiagonalDistance = 0.;
        // compute the array in reverse, first entry of the array is the total euclidean distance
        for (int i = m - 1; i >= 0; i--) {
            sumDiagonalDistance += Math.abs(X[i] - Y[i]);
            diagonalBounds[i] = sumDiagonalDistance;
        }

        return diagonalBounds;
    }

    public double C(double new_point, double x, double y) {

        // c - cost of Split/Merge operation. Change this value to what is more
        // appropriate for your data.

        if (new_point < Math.min(x, y) || new_point > Math.max(x, y)) {
            return this.c + Math.min(Math.abs(new_point - x), Math.abs(new_point - y));
        }

        return this.c;
    }



}
