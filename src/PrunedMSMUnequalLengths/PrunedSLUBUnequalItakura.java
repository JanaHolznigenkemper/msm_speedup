package PrunedMSMUnequalLengths;


import java.util.Arrays;

import HeuristicsUnequalLengths.ItakuraUnequalLengths;

public class PrunedSLUBUnequalItakura {

    final double c;

    final double[] X;
    final double[] Y;


    double[] ts1;
    double[] ts2;

    double upperBound;
    double[] upperBoundArray;
    int remainingCostPointer;

    final int m;
    final int n;

    final int max;
    final int min;
    final int diff;



    // compute MSM distance with pruned version -> prune the dynamic table
    // for unequal lengths
    public PrunedSLUBUnequalItakura(double c, double[] X, double[] Y) {
        this.c = c;
        this.X = X;
        this.Y = Y;
        this.m = X.length;
        this.n = Y.length;
        this.diff = Math.abs(m - n);

        this.max = Math.max(m, n);
        this.min = Math.min(m,n);

        // first column and row is inf -> need to extend time series


        // the longer time series has to be on the y-Axis
        if (m == max) {
            this.ts1 = new double[m + 1];
            this.ts2 = new double[n + 1];
            this.ts1[0] = Double.POSITIVE_INFINITY;
            this.ts2[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(X, 0, this.ts1, 1, m + 1 - 1);
            System.arraycopy(Y, 0, this.ts2, 1, n + 1 - 1);
        } else {
            this.ts1 = new double[n + 1];
            this.ts2 = new double[m + 1];
            this.ts1[0] = Double.POSITIVE_INFINITY;
            this.ts2[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(X, 0, this.ts2, 1, m + 1 - 1);
            System.arraycopy(Y, 0, this.ts1, 1, n + 1 - 1);
        }

        this.setUpperBound();

    }

    public void setUpperBound() {

        ItakuraUnequalLengths itakuraUnequalLengths = new ItakuraUnequalLengths(this.c, 3/4.);

        this.upperBoundArray = itakuraUnequalLengths.itakuraHeuristicArray(this.X, this.Y);
        this.upperBound = this.upperBoundArray[0] + 0.0000001;
        // start at position 1 to update UB; Entry 0 corresponds to total distance
        this.remainingCostPointer = 1;
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
        double[] tmpArray = new double[min + 1];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0
        double tmp = 0;

        // index for pruning start and end of row
        int sc = 1;
        int ec = 1;

        // count cells that are cut during the pruning
        int countCutCells = 0;

        // remember if an entry smaller than UB was found -> cannot cut
        boolean smallerFound;
        int ecNext;



        //initialize first row
        // the first entry of this row is inf, this inf is used for computing the second row
        Arrays.fill(tmpArray, Double.POSITIVE_INFINITY);



        //i row index
        for (int i = 1; i < max + 1; i++) {

            //compute bandwidth regarding the upper bound and the cost for merge and split operations
            int bandwidth = this.computeBandwidth();
            int start = Math.max(sc, i-bandwidth);

            int end = Math.min(i+bandwidth+1,tmpArray.length);



            smallerFound = false;
            final double xi = ts1[i];
            // the index for the pruned end cannot be lower than the diagonal
            // All entries on the diagonal have to be equal or smaller than
            // the upper bound (Euclidean distance = diagonal path)
            ecNext = i;


            //column index
            for (int j = start; j < end; j++) {

                final double yj = ts2[j];
                double d1, d2, d3;
                //split
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
                        tmpArray = this.fillWithInf(j + 1, min + 1, tmpArray);
                        break;
                    }
                } else {
                    smallerFound = true;
                    ecNext = j + 1;
                }

                //update Upper Bound if the slanted diagonal is hit
                if (j == Math.ceil(min / ((double) max) * i)) {
                    this.upperBound = tmpArray[j] + this.upperBoundArray[remainingCostPointer] + 0.00001;
                    remainingCostPointer++;
                }





            }

            tmpArray = this.fillWithInf(1, sc, tmpArray);

            // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
            tmp = Double.POSITIVE_INFINITY;
            ec = ecNext;
            countCutCells = countCutCells + (ts1.length - (ec)) + (sc - 1);

        }


        double totalSizeTable = (double) max * min;
        double percentageCutCells = countCutCells / (totalSizeTable);
        return new double[]{tmpArray[min], percentageCutCells};
    }

    /**
     * @param start inklusiv
     * @param end   exklusiv
     * @param array array to fill with inf
     * @return filled with inf
     */
    public double[] fillWithInf(int start, int end, double[] array) {
        for (int i = start; i < end; i++) {
            array[i] = Double.POSITIVE_INFINITY;
        }
        return array;
    }

    /**
     * Compute simple lower bound: Minimum number of merge or split operations
     *
     * @param xCoord current x Coordinate
     * @param yCoord current y Coordinate
     * @return double simple lower bound
     */
    public double getLowerBound(int xCoord, int yCoord) {
        return Math.abs((max - xCoord) - (min - yCoord)) * this.c;
    }

    public int computeBandwidth(){
        return (int) Math.ceil(this.upperBound/ this.c);
    }

    /**
     * cost of Split/Merge operation
     *
     * @param new_point point to merge/ split to
     * @param x         xcoord
     * @param y         ycoord
     * @return cost for merge/split
     */
    public double C(double new_point, double x, double y) {
        if (new_point < Math.min(x, y) || new_point > Math.max(x, y)) {
            return this.c + Math.min(Math.abs(new_point - x), Math.abs(new_point - y));
        }

        return this.c;
    }


    

}
