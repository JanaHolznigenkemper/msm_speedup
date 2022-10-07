package PrunedMSMUnequalLengths;



import java.util.Arrays;

import HeuristicsUnequalLengths.GreedyUnequalLengths;

public class PrunedSLUBUnequalGreedy {
    final double c;

    final double[] X;
    final double[] Y;


    double[] ts1;
    double[] ts2;

    double upperBound;
    double[] upperBoundArray;

    final int m;
    final int n;

    final int max;
    final int diff;

    //for non-quadratic arrays, we need to precompute the indices, where to update the upper bound
    int[] updateIdxi;
    int[] updateIdxj;

    // compute MSM distance with pruned version -> prune the dynamic table
    // for unequal lengths
    public PrunedSLUBUnequalGreedy(double c, double[] X, double[] Y) {
        this.c = c;
        this.X = X;
        this.Y = Y;
        this.m = X.length;
        this.n = Y.length;
        this.diff = Math.abs(m - n);

        this.max = Math.max(m, n);

        // first column and row is inf -> need to extend time series
        this.ts1 = new double[m + 1];
        this.ts2 = new double[n + 1];
        ts1[0] = Double.POSITIVE_INFINITY;
        ts2[0] = Double.POSITIVE_INFINITY;

        System.arraycopy(X, 0, ts1, 1, m + 1 - 1);
        System.arraycopy(Y, 0, ts2, 1, n + 1 - 1);

        this.setUpperBound();
        computeUpdateIndex();
    }

    public void setUpperBound() {
        GreedyUnequalLengths greedyUnequalLengths = new GreedyUnequalLengths(this.X, this.Y, this.c);
        double[] upperBoundGreedy = greedyUnequalLengths.computeGreedyArrayUnequalLengths();
        this.upperBoundArray = upperBoundGreedy;
        this.upperBound = upperBoundGreedy[0] + 0.0000001;

    }

    public void computeUpdateIndex() {
        this.updateIdxi = new int[this.max];
        this.updateIdxj = new int[this.max];

        int counti;
        int countj;

        if (this.max == this.m) {
            counti = 1;
            countj = 1;
            for (int k = 0; k < max; k++) {
                updateIdxi[k] = counti++;
            }
            for (int k = 0; k < diff; k++) {
                updateIdxj[k] = 1;
            }
            for (int k = diff; k < this.max; k++) {
                updateIdxj[k] = countj++;
            }
        }


        if (this.max == this.n) {
            counti = 1;
            countj = 1;
            for (int k = 0; k < max; k++) {
                updateIdxj[k] = countj++;
            }
            for (int k = 0; k < diff; k++) {
                updateIdxi[k] = 1;
            }
            for (int k = diff; k < this.max; k++) {
                updateIdxi[k] = counti++;
            }
        }

    }

    public int computeBandwidth(){

        return (int) Math.ceil(this.upperBound/ this.c);
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
        double[] tmpArray = new double[n + 1];
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

        int pointerUpdateIndexes = 0;
        int currentUpdatei = this.updateIdxi[pointerUpdateIndexes];
        int currentUpdatej = this.updateIdxj[pointerUpdateIndexes];

        //initialize first row
        // the first entry of this row is inf, this inf is used for computing the second row
        Arrays.fill(tmpArray, Double.POSITIVE_INFINITY);


        //i row index
        for (int i = 1; i < m + 1; i++) {

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
                        tmpArray = this.fillWithInf(j + 1, n + 1, tmpArray);
                        break;
                    }
                } else {
                    smallerFound = true;
                    ecNext = j + 1;
                }

                //update Upper Bound if the respective update coordinates are hit
                if(i==currentUpdatei && j==currentUpdatej){
                    this.upperBound = tmpArray[j] + this.upperBoundArray[pointerUpdateIndexes+1] + 0.00001;
                    pointerUpdateIndexes++;
                    if(pointerUpdateIndexes<updateIdxi.length) {
                        currentUpdatei = this.updateIdxi[pointerUpdateIndexes];
                        currentUpdatej = this.updateIdxj[pointerUpdateIndexes];
                    }


                }






            }

            //  System.out.println(this.upperBound);
            tmpArray = this.fillWithInf(1, sc, tmpArray);

            // set tmp to infinity since the move computation in the next row is not possible and accesses tmp
            tmp = Double.POSITIVE_INFINITY;
            ec = ecNext;
            countCutCells = countCutCells + (ts1.length - (ec)) + (sc - 1);

        }


        double totalSizeTable = (double) m * n;
        double percentageCutCells = countCutCells / (totalSizeTable);

 //       System.out.println("counter Bandwidth useful: " +counterBandwidth);
        return new double[]{tmpArray[n], percentageCutCells};
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
        return Math.abs((m - xCoord) - (n - yCoord)) * this.c;
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


    public static void main(String[] args) {


        double[] ts1 = new double[]{ 14396.0, 14396.0, 14395.0, 14397.0, 14395.0, 14394.0, 14394.0, 14393.0, 14393.0, 14392.0, 14392.0, 14391.0, 14391.0, 14390.0, 14390.0, 14389.0, 14389.0, 14388.0, 14388.0, 14387.0, 14387.0, 14386.0, 14386.0, 14385.0, 14385.0, 14384.0, 14384.0, 14383.0, 14383.0, 14382.0, 14382.0, 14381.0, 14381.0, 14380.0, 14380.0, 14379.0, 14379.0, 14378.0, 14378.0, 14377.0, 14377.0, 14376.0, 14376.0, 14375.0, 14375.0, 14374.0, 14374.0, 14373.0, 14373.0, 14372.0, 14372.0, 14371.0, 14371.0, 14370.0, 14370.0, 14369.0, 14369.0, 14368.0, 14368.0, 14367.0, 14367.0, 14366.0, 14366.0, 14365.0, 14365.0, 14364.0, 14364.0, 14363.0, 14363.0, 14362.0, 14362.0, 14361.0, 14361.0, 14360.0, 14360.0, 14359.0, 14359.0, 14358.0, 14358.0, 14357.0, 14357.0, 14356.0, 14356.0, 14355.0, 14355.0, 14354.0, 14354.0, 14353.0, 14353.0, 14352.0, 14352.0, 14351.0, 14351.0, 14350.0, 14350.0, 14349.0, 14349.0, 14348.0, 14348.0, 14347.0, 14347.0, 14346.0, 14346.0, 14345.0, 14345.0, 14344.0, 14344.0, 14343.0, 14343.0, 14342.0, 14342.0, 14341.0, 14341.0, 14340.0, 14340.0, 14339.0, 14339.0, 14338.0, 14338.0, 14337.0, 14337.0, 14336.0, 14336.0, 14335.0, 14335.0, 14334.0, 14334.0, 14333.0, 14333.0, 14332.0, 14332.0, 14331.0, 14331.0, 14330.0, 14330.0, 14329.0, 14329.0, 14328.0, 14328.0, 14327.0, 14327.0, 14326.0, 14326.0, 14325.0, 14325.0, 14324.0, 14324.0, 14323.0, 14323.0, 14322.0, 14322.0, 14321.0, 14320.0, 14320.0, 14319.0, 14319.0, 14318.0, 14318.0, 14317.0, 14317.0, 14316.0, 14316.0, 14315.0, 14315.0, 14314.0, 14314.0, 14313.0, 14313.0, 14312.0, 14312.0, 14311.0, 14311.0, 14310.0, 14310.0, 14309.0, 14309.0, 14308.0, 14308.0, 14307.0, 14307.0, 14306.0, 14306.0, 14305.0, 14305.0, 14304.0, 14304.0, 14303.0, 14303.0, 14302.0, 14302.0, 14301.0, 14301.0, 14300.0, 14300.0, 14299.0, 14299.0, 14298.0, 14298.0, 14297.0, 14297.0, 14296.0, 14296.0, 14295.0, 14295.0, 14294.0, 14294.0, 14293.0, 14293.0, 14292.0, 14292.0, 14291.0, 14291.0, 14290.0, 14290.0, 14289.0, 14289.0, 14288.0, 14288.0, 14287.0, 14287.0, 14286.0, 14286.0, 14285.0, 14285.0, 14284.0, 14284.0, 14283.0, 14283.0, 14282.0, 14282.0, 14281.0, 14281.0, 14280.0, 14280.0, 14279.0, 14279.0, 14278.0, 14278.0, 14277.0, 14277.0, 14276.0, 14276.0, 14275.0, 14275.0, 14274.0, 14274.0, 14273.0, 14273.0, 14272.0, 14272.0, 14271.0, 14271.0, 14270.0, 14270.0, 14269.0, 14269.0, 14268.0, 14268.0, 14267.0, 14267.0, 14266.0, 14266.0, 14265.0, 14265.0, 14264.0, 14264.0, 14263.0, 14263.0, 14262.0, 14262.0, 14261.0, 14261.0, 14260.0, 14260.0, 14259.0, 14259.0, 14258.0, 14258.0, 14257.0, 14257.0, 14256.0, 14256.0, 14255.0, 14255.0, 14254.0, 14254.0, 14253.0, 14253.0, 14252.0, 14252.0, 14251.0, 14251.0, 14250.0, 14250.0, 14249.0, 14249.0, 14248.0, 14248.0, 14247.0, 14247.0, 14246.0, 14246.0, 14245.0, 14245.0, 14244.0, 14244.0, 14243.0, 14243.0, 14242.0, 14242.0, 14241.0, 14241.0, 14240.0, 14240.0, 14239.0, 14239.0, 14238.0, 14238.0, 14237.0, 14237.0, 14236.0, 14236.0, 14235.0, 14235.0, 14234.0, 14234.0, 14233.0, 14233.0, 14232.0, 14232.0, 14231.0, 14231.0, 14230.0, 14230.0, 14229.0, 14229.0, 14228.0, 14228.0, 14227.0, 14227.0, 14226.0, 14226.0, 14225.0, 14225.0, 14224.0, 14224.0, 14223.0, 14223.0, 14222.0, 14222.0, 14221.0, 14221.0, 14220.0, 14220.0, 14219.0, 14219.0, 14218.0, 14218.0, 14217.0, 14217.0, 14216.0, 14216.0, 14215.0, 14215.0, 14214.0, 14214.0, 14213.0, 14213.0, 14212.0, 14212.0, 14211.0, 14211.0, 14210.0, 14210.0, 14209.0, 14209.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14209.0, 14210.0, 14211.0, 14212.0, 14213.0, 14214.0, 14214.0, 14215.0, 14215.0, 14216.0, 14216.0, 14217.0, 14217.0, 14218.0, 14218.0, 14218.0, 14219.0, 14219.0, 14220.0, 14220.0, 14221.0, 14222.0, 14222.0, 14223.0, 14224.0, 14224.0, 14225.0, 14225.0, 14226.0, 14226.0, 14227.0, 14227.0, 14227.0, 14228.0, 14228.0, 14228.0, 14229.0, 14229.0, 14229.0, 14230.0, 14231.0, 14231.0, 14232.0, 14232.0, 14233.0, 14234.0, 14234.0, 14234.0, 14235.0, 14235.0, 14236.0, 14236.0, 14237.0, 14237.0, 14238.0, 14238.0, 14239.0, 14239.0, 14240.0, 14240.0, 14241.0, 14241.0, 14242.0};
        double[] ts2 = new double[]{14261.0, 14262.0, 14262.0, 14263.0, 14263.0, 14264.0, 14264.0, 14265.0, 14265.0, 14266.0, 14266.0, 14266.0, 14267.0, 14267.0, 14268.0, 14268.0, 14269.0, 14269.0, 14270.0, 14270.0, 14270.0, 14271.0, 14271.0, 14272.0, 14272.0, 14273.0, 14273.0, 14274.0, 14274.0, 14274.0, 14275.0, 14275.0, 14276.0, 14276.0, 14277.0, 14277.0, 14278.0, 14278.0, 14279.0, 14279.0, 14279.0, 14280.0, 14280.0, 14281.0, 14281.0, 14282.0, 14282.0, 14283.0, 14283.0, 14283.0, 14284.0, 14284.0, 14285.0, 14285.0, 14286.0, 14286.0, 14287.0, 14287.0, 14288.0, 14288.0, 14288.0, 14289.0, 14289.0, 14290.0, 14290.0, 14291.0, 14291.0, 14292.0, 14292.0, 14292.0, 14293.0, 14293.0, 14294.0, 14294.0, 14295.0, 14295.0, 14296.0, 14296.0, 14296.0, 14297.0, 14297.0, 14298.0, 14298.0, 14299.0, 14299.0, 14300.0, 14300.0, 14301.0, 14301.0, 14301.0, 14302.0, 14302.0, 14303.0, 14303.0, 14304.0, 14304.0, 14305.0, 14305.0, 14305.0, 14306.0, 14306.0, 14307.0, 14307.0, 14308.0, 14308.0, 14309.0, 14309.0, 14310.0, 14310.0, 14310.0, 14311.0, 14311.0, 14312.0, 14312.0, 14313.0, 14313.0, 14314.0, 14314.0, 14314.0, 14315.0, 14315.0, 14316.0, 14316.0, 14317.0, 14317.0, 14318.0, 14318.0, 14318.0, 14319.0, 14319.0, 14320.0, 14320.0, 14321.0, 14321.0, 14322.0, 14322.0, 14323.0, 14323.0, 14323.0, 14324.0, 14324.0, 14325.0, 14325.0, 14326.0, 14326.0, 14327.0, 14327.0, 14327.0, 14328.0, 14328.0, 14329.0, 14329.0, 14330.0, 14330.0, 14331.0, 14331.0, 14332.0, 14332.0, 14332.0, 14333.0, 14333.0, 14334.0, 14334.0, 14335.0, 14335.0, 14336.0, 14336.0, 14336.0, 14337.0, 14337.0, 14338.0, 14338.0, 14339.0, 14339.0, 14340.0, 14340.0, 14341.0, 14341.0, 14341.0, 14342.0, 14342.0, 14343.0, 14343.0, 14344.0, 14344.0, 14345.0, 14345.0, 14345.0, 14346.0, 14346.0, 14347.0, 14347.0, 14348.0, 14348.0, 14349.0, 14349.0, 14349.0, 14350.0, 14350.0, 14351.0, 14351.0, 14352.0, 14352.0, 14353.0, 14353.0, 14354.0, 14354.0, 14354.0, 14355.0, 14355.0, 14356.0, 14356.0, 14357.0, 14357.0, 14358.0, 14358.0, 14358.0, 14359.0, 14359.0, 14360.0, 14360.0, 14361.0, 14361.0, 14362.0, 14362.0, 14363.0, 14363.0, 14363.0, 14364.0, 14364.0, 14365.0, 14365.0, 14366.0, 14366.0, 14367.0, 14367.0, 14367.0, 14368.0, 14368.0, 14369.0, 14369.0, 14370.0, 14370.0, 14371.0, 14371.0, 14371.0, 14372.0, 14372.0, 14373.0, 14373.0, 14374.0, 14374.0, 14375.0, 14375.0, 14376.0, 14376.0, 14376.0, 14377.0, 14377.0, 14378.0, 14378.0, 14379.0, 14379.0, 14380.0, 14380.0, 14380.0, 14381.0, 14381.0, 14382.0, 14382.0, 14383.0, 14383.0, 14384.0, 14384.0, 14385.0, 14385.0, 14385.0, 14386.0, 14386.0, 14387.0, 14387.0, 14388.0, 14388.0, 14389.0, 14389.0, 14389.0, 14390.0, 14390.0, 14391.0, 14391.0, 14392.0, 14392.0, 14393.0, 14393.0, 14393.0, 14394.0, 14394.0, 14395.0, 14395.0, 14396.0, 14396.0, 14397.0, 14397.0, 14398.0, 14398.0, 14398.0, 14399.0, 14399.0, 14400.0, 14400.0, 14401.0, 14401.0, 14402.0, 14402.0, 14402.0, 14403.0, 14403.0, 14404.0, 14404.0, 14405.0, 14405.0, 14406.0, 14406.0, 14407.0, 14407.0, 14407.0, 14408.0, 14408.0, 14409.0, 14409.0, 14410.0, 14410.0, 14411.0, 14411.0, 14411.0, 14412.0, 14412.0, 14413.0, 14413.0, 14414.0, 14414.0, 14415.0, 14415.0, 14415.0, 14416.0, 14416.0, 14417.0, 14417.0, 14418.0, 14418.0, 14419.0, 14419.0, 14420.0, 14420.0, 14420.0, 14421.0, 14421.0, 14422.0, 14422.0, 14423.0, 14423.0, 14424.0, 14424.0, 14424.0, 14425.0, 14425.0, 14426.0, 14426.0, 14427.0, 14427.0, 14428.0, 14428.0, 14429.0, 14429.0, 14429.0, 14430.0, 14430.0, 14431.0, 14431.0, 14432.0, 14432.0, 14433.0, 14433.0, 14433.0, 14434.0, 14434.0, 14435.0, 14435.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14443.0, 14452.0, 14460.0, 14469.0, 14478.0, 14486.0, 14488.0, 14481.0, 14474.0, 14468.0, 14461.0, 14454.0, 14447.0, 14438.0, 14429.0, 14421.0, 14412.0, 14403.0, 14394.0, 14398.0, 14410.0, 14422.0, 14434.0, 14446.0, 14458.0, 14464.0, 14449.0, 14435.0, 14421.0, 14406.0, 14392.0, 14378.0, 14385.0, 14398.0, 14410.0, 14422.0, 14435.0, 14447.0, 14451.0, 14444.0, 14437.0, 14429.0, 14422.0, 14415.0, 14408.0, 14403.0, 14399.0, 14394.0, 14390.0, 14385.0, 14381.0, 14385.0, 14395.0, 14404.0, 14414.0, 14424.0, 14433.0, 14439.0, 14428.0, 14418.0, 14407.0, 14397.0, 14386.0};

        PrunedSLUBUnequalGreedy prunedSLUBUnequalGreedy = new PrunedSLUBUnequalGreedy(0.1, ts1, ts2);
        System.out.println(prunedSLUBUnequalGreedy.msmDistPruned()[0]);
    }
}
