package PrunedMSMEqualLengths;

import HeuristicsUnequalLengths.*;
import HeuristicsEqualLengths.Greedy;

import java.util.Arrays;
import java.util.Objects;


/**
 * Pruned version with simple LB
 * Band is computed using the remaining cost per entry [i,i] from the greedy algorithm
 */
public class PrunedSimpleLBUB {
    final double c;

    final double[] X;
    final double[] Y;

    final String UB;

    double[] ts1;
    double[] ts2;

    double upperBound;
    double[] upperBoundArray;
    final int m;


    /**
     * compute MSM distance with pruned version -> prune the dynamic programming table
     * with update UB and Simple LB
     * with additional band regarding max merge and split operations
     * @param c  cost for merge and split
     * @param UB Upper Bound Strategy: g = greedy, i = Itakura with d = 3/4
     */
    public PrunedSimpleLBUB(double c, double[] X, double[] Y, String UB ) {
        this.c = c;
        this.X = X;
        this.Y = Y;
        this.m = X.length;
        if (m != Y.length) {
            throw new RuntimeException("Unequal lengths of time series!");

        }
        this.UB = UB;

        this.setUpperBound();
        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule
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

    /**
     * compute the Msm distance table with pruned entries
     *
     * @return pair: first Double: Distance, second Double: relative amount of pruned cells
     */
    public double[] msmDistPruned() {

        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        double[] tmpArray = new double[m + 1];


        //value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0
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

      //  int counterBandwidth =0;
        //row index
        for (int i = 1; i < tmpArray.length; i++) {

            //compute bandwidth regarding the upper bound
            int bandwidth = this.computeBandwidth();
            int start = Math.max(sc, i-bandwidth);

            int end = Math.min(i+bandwidth+1,tmpArray.length);


            final double xi = ts1[i];
            // the index for the pruned end cannot be lower than the diagonal
            // All entries on the diagonal have to be equal or smaller than
            // the upper bound (Euclidean distance = diagonal path)
            ecNext = i;
            smallerFound = false;

            //column index
            for (int j = start; j < end; j++) {



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
                double lb = getLowerBound(i, j);
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
                    this.upperBound = tmpArray[j] + this.upperBoundArray[j] +0.00001;
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

        return new double[]{tmpArray[m], percentageCutCells};
    }


    public int computeBandwidth(){

        return (int) Math.ceil(this.upperBound/ this.c);
    }

    /**
     * @param start inklusiv
     * @param end   exklusiv
     * @param array array
     * @return filled array
     */
    public double[] fillWithInf(int start, int end, double[] array) {
        for (int i = start; i < end; i++) {
            array[i] = Double.POSITIVE_INFINITY;
        }
        return array;
    }


    public double getLowerBound(int xCoord, int yCoord) {

        return Math.abs(xCoord - yCoord) * c;

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
