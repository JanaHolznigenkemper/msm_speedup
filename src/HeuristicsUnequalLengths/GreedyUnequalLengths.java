package HeuristicsUnequalLengths;



/**
 * Compute GreedyUnequalLengths Algorithm for two time series of different length
 */
public class GreedyUnequalLengths {

    final double[] X;
    final double[] Y;
    final double c;
    final int m;
    final int n;
    //determine which time series is the longest.
    final int min;
    final int max;
    final int diff;
    final double[] shortTs;
    final double[] longTs;

    public GreedyUnequalLengths(double[] x, double[] y, double c) {
        X = x;
        Y = y;
        this.c = c;
        this.m = x.length;
        this.n = Y.length;
        if (this.m > this.n) {
            this.max = this.m;
            this.min = this.n;
            this.shortTs = y;
            this.longTs = x;
        } else {
            this.max = this.n;
            this.min = this.m;
            this.shortTs = x;
            this.longTs = y;
        }
        this.diff = max - min;

    }

    /**
     * Compute the greedy heuristic backwards to update the upper bound in the pruned MSM Version
     *
     * @return array of intermediate upper bounds
     */
    public double[] computeGreedyArrayUnequalLengths() {
        double[] greedyArray = new double[this.max + 1];

        // compute upper Bounds for every diagonal entry
        // the upper bound is computed from right to left: Possibility to update the upper Bound when a diagonal entry is computed
        // upper Bound for the last entry is 0, because there is nothing left to compute
        greedyArray[max] = 0.;

        // assume that the time series are aligned at the end
        double distCurrent = Math.abs(X[m - 1] - Y[n - 1]);
        double distTmp = distCurrent;
        double shortCurrent;
        double longCurrent;
        double shortTmp = this.shortTs[min - 1];
        double longTmp = this.longTs[max - 1];
        int rel = (shortTmp > longTmp) ? 1 : 2;


        greedyArray[max - 1] = distCurrent;

// first loop for suitable alignment
        for (int i = 2; i <= this.max - this.diff; i++) {
            shortCurrent = this.shortTs[this.min - i];
            longCurrent = this.longTs[this.max - i];
            distCurrent = shortCurrent - longCurrent;

            if ((rel == 1) && distCurrent > 0) {
                if (distCurrent > 2 * c && distTmp > 2 * c) {
                    greedyArray[this.max - i] = 2 * c + Math.abs(shortCurrent - shortTmp) + Math.abs(longCurrent - longTmp) + greedyArray[this.max - i + 1];
                } else {
                    greedyArray[this.max - i] = distCurrent + greedyArray[this.max - i + 1];
                }
            } else if ((rel == 1) && distCurrent <= 0) {
                greedyArray[this.max - i] = -1 * distCurrent + greedyArray[this.max - i + 1];
                rel = 2;
            } else if ((rel == 2) && distCurrent <= 0) {
                if (Math.abs(distCurrent) > 2 * c && Math.abs(distTmp) > 2 * c) {
                    greedyArray[this.max - i] = 2 * c + Math.abs(shortCurrent - shortTmp) + Math.abs(longCurrent - longTmp) + greedyArray[this.max - i + 1];
                } else {
                    greedyArray[this.max - i] = -1 * distCurrent + greedyArray[this.max - i + 1];
                }
            } else {
                greedyArray[this.max - i] = distCurrent + greedyArray[this.max - i + 1];
                rel = 1;
            }

            distTmp = distCurrent;
            shortTmp = shortCurrent;
            longTmp = longCurrent;
        }

        shortCurrent = this.shortTs[0];
        //Aligning the rest of the longer time series to the shorter on
        //Cost for merges and moves to the last point in the shorter TS
        for (int i = max - diff + 1; i <= max; i++) {
            longCurrent = this.longTs[max - i];

            distCurrent = shortCurrent - longCurrent;

            //the algorithm changes in that way, that the 2c intervals changes to a c interval since no split needs to be done
            // for move operations also c has to be added
            if ((rel == 1) && distCurrent > 0) {
                if (distCurrent > c && distTmp > c) {
                    greedyArray[this.max - i] = c + Math.abs(shortCurrent - shortTmp) + Math.abs(longCurrent - longTmp) + greedyArray[this.max - i + 1];
                } else {
                    greedyArray[this.max - i] = c + distCurrent + greedyArray[this.max - i + 1];
                }
            } else if ((rel == 1) && distCurrent <= 0) {
                greedyArray[this.max - i] = -1 * distCurrent + c + greedyArray[this.max - i + 1];
                rel = 2;
            } else if ((rel == 2) && distCurrent <= 0) {
                if (Math.abs(distCurrent) > c && Math.abs(distTmp) > c) {
                    greedyArray[this.max - i] = c + Math.abs(shortCurrent - shortTmp) + Math.abs(longCurrent - longTmp) + greedyArray[this.max - i + 1];
                } else {
                    greedyArray[this.max - i] = -1 * distCurrent + c + greedyArray[this.max - i + 1];
                }
            } else {
                greedyArray[this.max - i] = distCurrent + c + greedyArray[this.max - i + 1];
                rel = 1;
            }

            distTmp = distCurrent;
            shortTmp = shortCurrent;
            longTmp = longCurrent;
        }

        return greedyArray;
    }

    public double upperBoundOnlyMove() {
        double upperBound = 0;
        int n = Y.length;
        int m = X.length;

        if (n >= m) {

            for (int i = 0; i < m; i++) {
                upperBound += Math.abs(X[i] - Y[i]);

            }
            for (int i = m; i < n; i++) {
                upperBound += Math.abs(X[m - 1] - Y[i]) + this.c;

            }
        } else {

            for (int i = 0; i < n; i++) {
                upperBound += Math.abs(X[i] - Y[i]);

            }
            for (int i = n; i < m; i++) {
                upperBound += Math.abs(X[i] - Y[n - 1]) + this.c;

            }

        }

        return upperBound;
    }


}
