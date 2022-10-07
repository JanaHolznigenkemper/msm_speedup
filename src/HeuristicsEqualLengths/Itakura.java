package HeuristicsEqualLengths;


public class Itakura {

    final double c;
    final double d;

    /**
     * Computing the Itakura Heuristic for time series of unequal lengths
     *
     * @param c Cost for Split and Merge
     * @param d Size of the parallelogram: d \in (0,1], d=1 is only the diagonal/ or slanted diagonal
     */
    public Itakura(double c, double d) {
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


        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        int tmpArrayLength = m + 1;
        double[] tmpArray = new double[tmpArrayLength];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0

        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule
        double[] X;
        double[] Y;

        X = new double[m + 1];
        Y = new double[m + 1];
        X[0] = Double.POSITIVE_INFINITY;
        Y[0] = Double.POSITIVE_INFINITY;
        System.arraycopy(Xinput, 0, X, 1, m + 1 - 1);
        System.arraycopy(Yinput, 0, Y, 1, m + 1 - 1);

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
        for (int i = 1; i < m + 1; i++) {

            // determine start and end for Y
            // start at min 1, since the first entry is a special case

            double start1 = d * i;
            double start2 = 1 / d * i - (1 - d) / d * m;
            int start = (int) Math.ceil(Math.max(start1, start2));

            double end1 = 1 / d * i;
            double end2 = d * i + (1 - d) * m;
            int end = (int) Math.ceil(Math.min(end1, end2));


            final double xi = X[i];

            tmp = tmpArray[start - 1];

            //first entry is always inf
            tmpArray[0] = Double.POSITIVE_INFINITY;
            for (int k = 0; k < start; k++) {
                tmpArray[k] = Double.POSITIVE_INFINITY;
            }
            for (int k = end + 1; k <= m; k++) {
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


        return tmpArray[m];
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


        // array with remaining cost
        double[] remainingCost = new double[m + 1];
        remainingCost[m] = 0.;


        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        int tmpArrayLength = m + 1;
        double[] tmpArray = new double[tmpArrayLength];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0

        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule
        double[] X;
        double[] Y;

        X = new double[m + 1];
        Y = new double[m + 1];
        X[0] = Double.POSITIVE_INFINITY;
        Y[0] = Double.POSITIVE_INFINITY;
        System.arraycopy(Xinput, 0, X, 1, m + 1 - 1);
        System.arraycopy(Yinput, 0, Y, 1, m + 1 - 1);

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
        for (int i = 1; i < m + 1; i++) {

            // determine start and end for Y
            // start at min 1, since the first entry is a special case

            double start1 = d * i;
            double start2 = 1 / d * i - (1 - d) / d * m;
            int start = (int) Math.ceil(Math.max(start1, start2));

            double end1 = 1 / d * i;
            double end2 = d * i + (1 - d) * m;
            int end = (int) Math.ceil(Math.min(end1, end2));


            final double xi = X[i];

            tmp = tmpArray[start - 1];

            //first entry is always inf
            tmpArray[0] = Double.POSITIVE_INFINITY;
            for (int k = 0; k < start; k++) {
                tmpArray[k] = Double.POSITIVE_INFINITY;
            }
            for (int k = end + 1; k <= m; k++) {
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


                if (j == i) {
                    remainingCost[m - j] = tmpArray[j];
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


}
