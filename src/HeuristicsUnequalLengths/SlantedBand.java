package HeuristicsUnequalLengths;

public class SlantedBand {

    final double c;
    final int bandwidth;

    public SlantedBand(double c, int bandwidth) {
        this.c = c;
        this.bandwidth = bandwidth;
    }

    // Implementation of the msm distance with Slanted Band and efficient storage usage.
    public double slantedBandHeuristic(double[] Xinput, double[] Yinput) {
        final int m = Xinput.length;
        final int n = Yinput.length;
        final int max = Math.max(m,n);
        final int min = Math.min(m,n);
        // Create an array with one extra entry, regarding the whole matrix we initialize
        // MsmDistAStar.Entry [0,0] is set to 0
        // the first row and the first column with inf --> Every entry follows the same computational rules
        int tmpArrayLength = min+1;
        double[] tmpArray = new double[tmpArrayLength];
        // value storing the first "real value" of the array before overwriting it
        // the first value of the first row has to be 0

        //need to extend the time series adding a first entry set to infinity
        // to facilitate the computation rule

        double[] X;
        double[] Y;

        // the longer time series has to be on the y-Axis
        if(m==max){
            X = new double[m + 1];
            Y = new double[n + 1];
            X[0] = Double.POSITIVE_INFINITY;
            Y[0] = Double.POSITIVE_INFINITY;
            System.arraycopy(Xinput, 0, X, 1, m + 1 - 1);
            System.arraycopy(Yinput, 0, Y, 1, n + 1 - 1);
        }else {
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

            //for non-quadratic arrays the diagonal does not always correspond to the [i,i] entries
            // determine shifted diagonal entry
            int diagj = (int) Math.ceil( ((double)min/max)*i);


            // determine start and end for Y
            // start at min 1, since the first entry is a special case
            int start = Math.max(1, diagj - bandwidth);

            int end = Math.min(min, diagj + bandwidth);


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

    public double C(double new_point, double x, double y) {

        // c - cost of Split/Merge operation. Change this value to what is more
        // appropriate for your data.

        if (new_point < Math.min(x, y) || new_point > Math.max(x, y)) {
            return this.c + Math.min(Math.abs(new_point - x), Math.abs(new_point - y));
        }

        return this.c;
    }

    public static void main(String[] args) {

        double[] ts2 = new double[]{4, 5, 6, 7};
        double[] ts1 = new double[]{3, 2, 1, 0, 1, 7, 5};
        SlantedBand slantedBand = new SlantedBand(1, 3);
        System.out.println(slantedBand.slantedBandHeuristic(ts1, ts2));

    }
}
