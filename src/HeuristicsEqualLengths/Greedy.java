package HeuristicsEqualLengths;


/**
 * Computing the greedy heuristic remaining cost array for time series of equal lengths
 */
public class Greedy {

    final double[] X;
    final double[] Y;
    final double c;
    final int m;


    public Greedy(double[] x, double[] y, double c) {
        X = x;
        Y = y;
        this.c = c;
        this.m = x.length;
    }

    /**
     * Compute the greedy heuristic backwards to update the upper bound in the pruned MSM Version
     *
     * @return array of intermediate upper bounds
     */
    public double[] computeGreedyArray() {
        double[] greedyArray = new double[this.m + 1];

        // compute upper Bounds for every diagonal entry
        // the upper bound is computed from right to left: Possibility to update the upper Bound when a diagonal entry is computed
        // upper Bound for the last entry is 0, because there is nothing left to compute
        greedyArray[m] = 0.;

        // assume that the time series are aligned at the end
        double distCurrent = Math.abs(X[m - 1] - Y[m - 1]);
        double distTmp = distCurrent;
        double xCurrent;
        double yCurrent;
        double xTmp = this.X[m - 1];
        double yTmp = this.Y[m - 1];
        int rel = (xTmp > yTmp) ? 1 : 2;


        greedyArray[m - 1] = distCurrent;

// first loop for suitable alignment
        for (int i = 2; i <= this.m; i++) {
            xCurrent = this.X[this.m - i];
            yCurrent = this.Y[this.m - i];
            distCurrent = xCurrent - yCurrent;

            if ((rel == 1) && distCurrent > 0) {
                if (distCurrent > 2 * c && distTmp > 2 * c) {
                    greedyArray[this.m - i] = 2 * c + Math.abs(xCurrent - xTmp) + Math.abs(yCurrent - yTmp) + greedyArray[this.m - i + 1];
                } else {
                    greedyArray[this.m - i] = distCurrent + greedyArray[this.m - i + 1];
                }
            } else if ((rel == 1) && distCurrent <= 0) {
                greedyArray[this.m - i] = -1 * distCurrent + greedyArray[this.m - i + 1];
                rel = 2;
            } else if ((rel == 2) && distCurrent <= 0) {
                if (Math.abs(distCurrent) > 2 * c && Math.abs(distTmp) > 2 * c) {
                    greedyArray[this.m - i] = 2 * c + Math.abs(xCurrent - xTmp) + Math.abs(yCurrent - yTmp) + greedyArray[this.m - i + 1];
                } else {
                    greedyArray[this.m - i] = -1 * distCurrent + greedyArray[this.m - i + 1];
                }
            } else {
                greedyArray[this.m - i] = distCurrent + greedyArray[this.m - i + 1];
                rel = 1;
            }

            distTmp = distCurrent;
            xTmp = xCurrent;
            yTmp = yCurrent;
        }

        return greedyArray;
    }



}

