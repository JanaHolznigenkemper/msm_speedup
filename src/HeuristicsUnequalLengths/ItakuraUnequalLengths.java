package HeuristicsUnequalLengths;
import java.util.Arrays;

public class ItakuraUnequalLengths {
    final double c;
    final double d;

    /**
     * Computing the Itakura Heuristic for time series of unequal lengths
     *
     * @param c Cost for Split and Merge
     * @param d Size of the parallelogram: d \in (0,1], d=1 is only the diagonal/ or slanted diagonal
     */
    public ItakuraUnequalLengths(double c, double d) {
        this.c = c;
        this.d = d;
    }


    /**
     * Compute the Itakura Heuristic for two time series of different lengths
     *
     * @param Xinput time series (y coordinate)
     * @param Yinput time series (x coordinate)
     * @return heuristic distance
     */
    public double itakuraHeuristic(double[] Xinput, double[] Yinput) {
        final int m = Xinput.length;
        final int n = Yinput.length;
        final int max = Math.max(m, n);
        final int min = Math.min(m, n);


        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        int tmpArrayLength = min + 1;
        double[] tmpArray = new double[tmpArrayLength];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0

        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule
        double[] X;
        double[] Y;

        // the longer time series has to be on the y-Axis
        if (m == max) {
            X = new double[m + 1];
            Y = new double[n + 1];
            X[0] = Double.POSITIVE_INFINITY;
            Y[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(Xinput, 0, X, 1, m + 1 - 1);
            System.arraycopy(Yinput, 0, Y, 1, n + 1 - 1);
        } else {
            X = new double[n + 1];
            Y = new double[m + 1];
            X[0] = Double.POSITIVE_INFINITY;
            Y[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(Xinput, 0, Y, 1, m + 1 - 1);
            System.arraycopy(Yinput, 0, X, 1, n + 1 - 1);
        }

        double tmp;

        tmpArray[0] = 0;


        // first row is inf
        for (int k = 1; k < tmpArrayLength; k++) {
            tmpArray[k] = Double.POSITIVE_INFINITY;
        }


        /*first line are only splits
         * for the first row the band is roundceil(2/3i) = 0
         * it is enough to compute the first move entry above
         * The rest of the array has to be set to inf to avoid false computation
         * It is sufficient just to set all right entries of the band to inf, since in
         * the recursion rule there is the special case of the first entry included.
         */

        // i is the row index
        for (int i = 1; i < max + 1; i++) {

            // determine start and end for Y
            // start at min 1, since the first entry is a special case

            double start1 = d * min / ((double) max) * i;
            double start2 = 1 / d * min / ((double) max) * i - (1 - d) / d * min;
            int start = (int) Math.ceil(Math.max(start1, start2));

            double end1 = 1 / d * min / ((double) max) * i;
            double end2 = d * min / ((double) max) * i + (1 - d) * min;
            int end = (int) Math.ceil(Math.min(end1, end2));


            final double xi = X[i];

            tmp = tmpArray[start - 1];

            //first entry is always inf
            tmpArray[0] = Double.POSITIVE_INFINITY;
            for (int k = 0; k < start; k++) {
                tmpArray[k] = Double.POSITIVE_INFINITY;
            }
            for (int k = end + 1; k <= min; k++) {
                tmpArray[k] = Double.POSITIVE_INFINITY;
            }


            for (int j = start; j <= end; j++) {
                double yj = Y[j];
                double d1, d2, d3;

                //move
                d1 = tmp + Math.abs(xi - yj);
                //merge
                d2 = tmpArray[j] + C(xi, X[i - 1], yj);
                //split
                d3 = tmpArray[j - 1] + C(yj, xi, Y[j - 1]);

                tmp = tmpArray[j];
                tmpArray[j] = Math.min(d1, Math.min(d2, d3));

            }


        }


        return tmpArray[min];
    }

    /**
     * Compute the Itakura Heuristic for two time series of different lengths
     *
     * @param Xinput time series (y coordinate)
     * @param Yinput time series (x coordinate)
     * @return heuristic distance
     */
    public double[] itakuraHeuristicArray(double[] Xinput, double[] Yinput) {
        final int m = Xinput.length;
        final int n = Yinput.length;
        final int max = Math.max(m, n);
        final int min = Math.min(m, n);

        // array with remaining cost
        double[] remainingCost = new double[max + 1];
        remainingCost[max] = 0.;
        int pointerRemainingCost = max-1;


        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        int tmpArrayLength = min + 1;
        double[] tmpArray = new double[tmpArrayLength];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0

        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule
        double[] X;
        double[] Y;

        // the longer time series has to be on the y-Axis
        if (m == max) {
            X = new double[m + 1];
            Y = new double[n + 1];
            X[0] = Double.POSITIVE_INFINITY;
            Y[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(Xinput, 0, X, 1, m + 1 - 1);
            System.arraycopy(Yinput, 0, Y, 1, n + 1 - 1);
        } else {
            X = new double[n + 1];
            Y = new double[m + 1];
            X[0] = Double.POSITIVE_INFINITY;
            Y[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(Xinput, 0, Y, 1, m + 1 - 1);
            System.arraycopy(Yinput, 0, X, 1, n + 1 - 1);
        }

        double tmp;

        tmpArray[0] = 0;


        // first row is inf
        for (int k = 1; k < tmpArrayLength; k++) {
            tmpArray[k] = Double.POSITIVE_INFINITY;
        }


        /*first line are only splits
         * for the first row the band is roundceil(2/3i) = 0
         * it is enough to compute the first move entry above
         * The rest of the array has to be set to inf to avoid false computation
         * It is sufficient just to set all right entries of the band to inf, since in
         * the recursion rule there is the special case of the first entry included.
         */

        // i is the row index
        for (int i = 1; i < max + 1; i++) {

            // determine start and end for Y
            // start at min 1, since the first entry is a special case

            double start1 = d * min / ((double) max) * i;
            double start2 = 1 / d * min / ((double) max) * i - (1 - d) / d * min;
            int start = (int) Math.ceil(Math.max(start1, start2));

            double end1 = 1 / d * min / ((double) max) * i;
            double end2 = d * min / ((double) max) * i + (1 - d) * min;
            int end = (int) Math.ceil(Math.min(end1, end2));


            final double xi = X[i];

            tmp = tmpArray[start - 1];

            //first entry is always inf
            tmpArray[0] = Double.POSITIVE_INFINITY;
            for (int k = 0; k < start; k++) {
                tmpArray[k] = Double.POSITIVE_INFINITY;
            }
            for (int k = end + 1; k <= min; k++) {
                tmpArray[k] = Double.POSITIVE_INFINITY;
            }


            for (int j = start; j <= end; j++) {
                double yj = Y[j];
                double d1, d2, d3;

                //move
                d1 = tmp + Math.abs(xi - yj);
                //merge
                d2 = tmpArray[j] + C(xi, X[i - 1], yj);
                //split
                d3 = tmpArray[j - 1] + C(yj, xi, Y[j - 1]);

                tmp = tmpArray[j];
                tmpArray[j] = Math.min(d1, Math.min(d2, d3));



                if (j == Math.ceil(min / ((double) max) * i)) {
                    remainingCost[pointerRemainingCost] = tmpArray[j];
                    if(pointerRemainingCost>0) {
                        pointerRemainingCost--;
                    }
                }

            }


        }


        return remainingCost;
    }


    public double C(double new_point, double x, double y) {

        // c - cost of Split/Merge operation. Change this value to what is more
        // appropriate for your data.

        if (new_point < Math.min(x, y) || new_point > Math.max(x, y)) {
            return this.c + Math.min(Math.abs(new_point - x), Math.abs(new_point - y));
        }

        return this.c;
    }

    public static void main(String[] args) {

        // double[] ts1 = new double[]{14393.0, 14392.0, 14392.0, 14391.0, 14391.0, 14390.0, 14390.0, 14389.0, 14389.0, 14388.0, 14388.0, 14387.0, 14387.0, 14386.0, 14386.0, 14385.0, 14385.0, 14384.0, 14384.0, 14383.0, 14383.0, 14382.0, 14382.0, 14381.0, 14381.0, 14380.0, 14380.0, 14379.0, 14379.0, 14378.0, 14378.0, 14377.0, 14377.0, 14376.0, 14376.0, 14375.0, 14375.0, 14374.0, 14374.0, 14373.0, 14373.0, 14372.0, 14372.0, 14371.0, 14371.0, 14370.0, 14370.0, 14369.0, 14369.0, 14368.0, 14368.0, 14367.0, 14367.0, 14366.0, 14366.0, 14365.0, 14365.0, 14364.0, 14364.0, 14363.0, 14363.0, 14362.0, 14362.0, 14361.0, 14361.0, 14360.0, 14360.0, 14359.0, 14359.0, 14358.0, 14358.0, 14357.0, 14357.0, 14356.0, 14356.0, 14355.0, 14355.0, 14354.0, 14354.0, 14353.0, 14353.0, 14352.0, 14352.0, 14351.0, 14351.0, 14350.0, 14350.0, 14349.0, 14349.0, 14348.0, 14348.0, 14347.0, 14347.0, 14346.0, 14346.0, 14345.0, 14345.0, 14344.0, 14344.0, 14343.0, 14343.0, 14342.0, 14342.0, 14341.0, 14341.0, 14340.0, 14340.0, 14339.0, 14339.0, 14338.0, 14338.0, 14337.0, 14337.0, 14336.0, 14336.0, 14335.0, 14335.0, 14334.0, 14334.0, 14333.0, 14333.0, 14332.0, 14332.0, 14331.0, 14331.0, 14330.0, 14330.0, 14329.0, 14329.0, 14328.0, 14328.0, 14327.0, 14327.0, 14326.0, 14326.0, 14325.0, 14325.0, 14324.0, 14324.0, 14323.0, 14323.0, 14322.0, 14322.0, 14321.0, 14320.0, 14320.0, 14319.0, 14319.0, 14318.0, 14318.0, 14317.0, 14317.0, 14316.0, 14316.0, 14315.0, 14315.0, 14314.0, 14314.0, 14313.0, 14313.0, 14312.0, 14312.0, 14311.0, 14311.0, 14310.0, 14310.0, 14309.0, 14309.0, 14308.0, 14308.0, 14307.0, 14307.0, 14306.0, 14306.0, 14305.0, 14305.0, 14304.0, 14304.0, 14303.0, 14303.0, 14302.0, 14302.0, 14301.0, 14301.0, 14300.0, 14300.0, 14299.0, 14299.0, 14298.0, 14298.0, 14297.0, 14297.0, 14296.0, 14296.0, 14295.0, 14295.0, 14294.0, 14294.0, 14293.0, 14293.0, 14292.0, 14292.0, 14291.0, 14291.0, 14290.0, 14290.0, 14289.0, 14289.0, 14288.0, 14288.0, 14287.0, 14287.0, 14286.0, 14286.0, 14285.0, 14285.0, 14284.0, 14284.0, 14283.0, 14283.0, 14282.0, 14282.0, 14281.0, 14281.0, 14280.0, 14280.0, 14279.0, 14279.0, 14278.0, 14278.0, 14277.0, 14277.0, 14276.0, 14276.0, 14275.0, 14275.0, 14274.0, 14274.0, 14273.0, 14273.0, 14272.0, 14272.0, 14271.0, 14271.0, 14270.0, 14270.0, 14269.0, 14269.0, 14268.0, 14268.0, 14267.0, 14267.0, 14266.0, 14266.0, 14265.0, 14265.0, 14264.0, 14264.0, 14263.0, 14263.0, 14262.0, 14262.0, 14261.0, 14261.0, 14260.0, 14260.0, 14259.0, 14259.0, 14258.0, 14258.0, 14257.0, 14257.0, 14256.0, 14256.0, 14255.0, 14255.0, 14254.0, 14254.0, 14253.0, 14253.0, 14252.0, 14252.0, 14251.0, 14251.0, 14250.0, 14250.0, 14249.0, 14249.0, 14248.0, 14248.0, 14247.0, 14247.0, 14246.0, 14246.0, 14245.0, 14245.0, 14244.0, 14244.0, 14243.0, 14243.0, 14242.0, 14242.0, 14241.0, 14241.0, 14240.0, 14240.0, 14239.0, 14239.0, 14238.0, 14238.0, 14237.0, 14237.0, 14236.0, 14236.0, 14235.0, 14235.0, 14234.0, 14234.0, 14233.0, 14233.0, 14232.0, 14232.0, 14231.0, 14231.0, 14230.0, 14230.0, 14229.0, 14229.0, 14228.0, 14228.0, 14227.0, 14227.0, 14226.0, 14226.0, 14225.0, 14225.0, 14224.0, 14224.0, 14223.0, 14223.0, 14222.0, 14222.0, 14221.0, 14221.0, 14220.0, 14220.0, 14219.0, 14219.0, 14218.0, 14218.0, 14217.0, 14217.0, 14216.0, 14216.0, 14215.0, 14215.0, 14214.0, 14214.0, 14213.0, 14213.0, 14212.0, 14212.0, 14211.0, 14211.0, 14210.0, 14210.0, 14209.0, 14209.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14208.0, 14209.0, 14210.0, 14211.0, 14212.0, 14213.0, 14214.0, 14214.0, 14215.0, 14215.0, 14216.0, 14216.0, 14217.0, 14217.0, 14218.0, 14218.0, 14218.0, 14219.0, 14219.0, 14220.0, 14220.0, 14221.0, 14222.0, 14222.0, 14223.0, 14224.0, 14224.0, 14225.0, 14225.0, 14226.0, 14226.0, 14227.0, 14227.0, 14227.0, 14228.0, 14228.0, 14228.0, 14229.0, 14229.0, 14229.0, 14230.0, 14231.0, 14231.0, 14232.0, 14232.0, 14233.0, 14234.0, 14234.0, 14234.0, 14235.0, 14235.0, 14236.0, 14236.0, 14237.0, 14237.0, 14238.0, 14238.0, 14239.0, 14239.0, 14240.0, 14240.0, 14241.0, 14241.0, 14242.0};
        //double[] ts2 = new double[]{14261.0, 14262.0, 14262.0, 14263.0, 14263.0, 14264.0, 14264.0, 14265.0, 14265.0, 14266.0, 14266.0, 14266.0, 14267.0, 14267.0, 14268.0, 14268.0, 14269.0, 14269.0, 14270.0, 14270.0, 14270.0, 14271.0, 14271.0, 14272.0, 14272.0, 14273.0, 14273.0, 14274.0, 14274.0, 14274.0, 14275.0, 14275.0, 14276.0, 14276.0, 14277.0, 14277.0, 14278.0, 14278.0, 14279.0, 14279.0, 14279.0, 14280.0, 14280.0, 14281.0, 14281.0, 14282.0, 14282.0, 14283.0, 14283.0, 14283.0, 14284.0, 14284.0, 14285.0, 14285.0, 14286.0, 14286.0, 14287.0, 14287.0, 14288.0, 14288.0, 14288.0, 14289.0, 14289.0, 14290.0, 14290.0, 14291.0, 14291.0, 14292.0, 14292.0, 14292.0, 14293.0, 14293.0, 14294.0, 14294.0, 14295.0, 14295.0, 14296.0, 14296.0, 14296.0, 14297.0, 14297.0, 14298.0, 14298.0, 14299.0, 14299.0, 14300.0, 14300.0, 14301.0, 14301.0, 14301.0, 14302.0, 14302.0, 14303.0, 14303.0, 14304.0, 14304.0, 14305.0, 14305.0, 14305.0, 14306.0, 14306.0, 14307.0, 14307.0, 14308.0, 14308.0, 14309.0, 14309.0, 14310.0, 14310.0, 14310.0, 14311.0, 14311.0, 14312.0, 14312.0, 14313.0, 14313.0, 14314.0, 14314.0, 14314.0, 14315.0, 14315.0, 14316.0, 14316.0, 14317.0, 14317.0, 14318.0, 14318.0, 14318.0, 14319.0, 14319.0, 14320.0, 14320.0, 14321.0, 14321.0, 14322.0, 14322.0, 14323.0, 14323.0, 14323.0, 14324.0, 14324.0, 14325.0, 14325.0, 14326.0, 14326.0, 14327.0, 14327.0, 14327.0, 14328.0, 14328.0, 14329.0, 14329.0, 14330.0, 14330.0, 14331.0, 14331.0, 14332.0, 14332.0, 14332.0, 14333.0, 14333.0, 14334.0, 14334.0, 14335.0, 14335.0, 14336.0, 14336.0, 14336.0, 14337.0, 14337.0, 14338.0, 14338.0, 14339.0, 14339.0, 14340.0, 14340.0, 14341.0, 14341.0, 14341.0, 14342.0, 14342.0, 14343.0, 14343.0, 14344.0, 14344.0, 14345.0, 14345.0, 14345.0, 14346.0, 14346.0, 14347.0, 14347.0, 14348.0, 14348.0, 14349.0, 14349.0, 14349.0, 14350.0, 14350.0, 14351.0, 14351.0, 14352.0, 14352.0, 14353.0, 14353.0, 14354.0, 14354.0, 14354.0, 14355.0, 14355.0, 14356.0, 14356.0, 14357.0, 14357.0, 14358.0, 14358.0, 14358.0, 14359.0, 14359.0, 14360.0, 14360.0, 14361.0, 14361.0, 14362.0, 14362.0, 14363.0, 14363.0, 14363.0, 14364.0, 14364.0, 14365.0, 14365.0, 14366.0, 14366.0, 14367.0, 14367.0, 14367.0, 14368.0, 14368.0, 14369.0, 14369.0, 14370.0, 14370.0, 14371.0, 14371.0, 14371.0, 14372.0, 14372.0, 14373.0, 14373.0, 14374.0, 14374.0, 14375.0, 14375.0, 14376.0, 14376.0, 14376.0, 14377.0, 14377.0, 14378.0, 14378.0, 14379.0, 14379.0, 14380.0, 14380.0, 14380.0, 14381.0, 14381.0, 14382.0, 14382.0, 14383.0, 14383.0, 14384.0, 14384.0, 14385.0, 14385.0, 14385.0, 14386.0, 14386.0, 14387.0, 14387.0, 14388.0, 14388.0, 14389.0, 14389.0, 14389.0, 14390.0, 14390.0, 14391.0, 14391.0, 14392.0, 14392.0, 14393.0, 14393.0, 14393.0, 14394.0, 14394.0, 14395.0, 14395.0, 14396.0, 14396.0, 14397.0, 14397.0, 14398.0, 14398.0, 14398.0, 14399.0, 14399.0, 14400.0, 14400.0, 14401.0, 14401.0, 14402.0, 14402.0, 14402.0, 14403.0, 14403.0, 14404.0, 14404.0, 14405.0, 14405.0, 14406.0, 14406.0, 14407.0, 14407.0, 14407.0, 14408.0, 14408.0, 14409.0, 14409.0, 14410.0, 14410.0, 14411.0, 14411.0, 14411.0, 14412.0, 14412.0, 14413.0, 14413.0, 14414.0, 14414.0, 14415.0, 14415.0, 14415.0, 14416.0, 14416.0, 14417.0, 14417.0, 14418.0, 14418.0, 14419.0, 14419.0, 14420.0, 14420.0, 14420.0, 14421.0, 14421.0, 14422.0, 14422.0, 14423.0, 14423.0, 14424.0, 14424.0, 14424.0, 14425.0, 14425.0, 14426.0, 14426.0, 14427.0, 14427.0, 14428.0, 14428.0, 14429.0, 14429.0, 14429.0, 14430.0, 14430.0, 14431.0, 14431.0, 14432.0, 14432.0, 14433.0, 14433.0, 14433.0, 14434.0, 14434.0, 14435.0, 14435.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14436.0, 14443.0, 14452.0, 14460.0, 14469.0, 14478.0, 14486.0, 14488.0, 14481.0, 14474.0, 14468.0, 14461.0, 14454.0, 14447.0, 14438.0, 14429.0, 14421.0, 14412.0, 14403.0, 14394.0, 14398.0, 14410.0, 14422.0, 14434.0, 14446.0, 14458.0, 14464.0, 14449.0, 14435.0, 14421.0, 14406.0, 14392.0, 14378.0, 14385.0, 14398.0, 14410.0, 14422.0, 14435.0, 14447.0, 14451.0, 14444.0, 14437.0, 14429.0, 14422.0, 14415.0, 14408.0, 14403.0, 14399.0, 14394.0, 14390.0, 14385.0, 14381.0, 14385.0, 14395.0, 14404.0, 14414.0, 14424.0, 14433.0, 14439.0, 14428.0, 14418.0, 14407.0, 14397.0, 14386.0};

        double[] ts1 = new double[]{4, 5, 6, 5};
        double[] ts2 = new double[]{3, 2, 1, 0, 1, 7, 8, 9};
        ItakuraUnequalLengths itakuraUnequalLengths = new ItakuraUnequalLengths(1, 1.);
        System.out.println(Arrays.toString(itakuraUnequalLengths.itakuraHeuristicArray(ts1, ts2)));


    }
}
