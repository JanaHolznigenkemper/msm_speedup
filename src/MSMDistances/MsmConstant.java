package MSMDistances;


//Compute the distance from a given time series to a constant time series in a simpler manner than the MSM algorithm 
public class MsmConstant {

    double c;
    double constantTs;
    double lowerTx;
    double upperTx;
    double[] ts;
    double[] remainingCost;

    public MsmConstant(double[] ts, double c, double constantTs) {
        this.c = c;
        this.constantTs = constantTs;

        // the thresholds are constants representing a value in y-axis (The 2c area around the constant)
        this.lowerTx = constantTs - 2 * c;
        this.upperTx = constantTs + 2 * c;
        this.ts = ts;
    }

    public double[] computeRemainingCostPerEntry() {
        double[] rCostsAll = new double[this.ts.length+1];
        // the last entry in the array is 0, because the remaining cost for the last entry are 0
        rCostsAll[this.ts.length] = 0.;
        // the case pointer indicates the case of the last regarded time series point:
        // Case 0: point is within thresholds
        // Case 1: point is above the upper threshold
        // Case 2: point is below the lower threshold
        int casePointer = -1;

        for (int i = ts.length - 1; i >= 0; i--) {
            double tsi = ts[i];
            // Case 1
            if (tsi >= this.upperTx) {
                if (casePointer == 0 || casePointer == 2 || casePointer == -1) {
                    // if the case switches, the cost for the entry are the move cost
                    rCostsAll[i] = this.computeMoveCost(tsi) + rCostsAll[i+1];
                } else {// the case don't change -> cost are merge + split cost + remaining move cost depending on the position of the points
                    // this case will never appear for i=ts.length-1 (last entry)
                    rCostsAll[i] = this.computeMergeCost(tsi, ts[i + 1]) + rCostsAll[i+1];
                }
                casePointer = 1;
            } else if (tsi <= this.lowerTx) {
                // same as before -> case switches
                if (casePointer == 0 || casePointer == 1 || casePointer == -1) {
                    // if the case switches, the cost for the entry are the move cost
                    rCostsAll[i] = this.computeMoveCost(tsi) + rCostsAll[i+1];
                } else { // the case don't change -> cost are merge + split cost + remaining move cost depending on the position of the points
                    // this case will never appear for i=ts.length-1 (last entry)
                    rCostsAll[i] = this.computeMergeCost(tsi, ts[i + 1]) + rCostsAll[i+1];
                }
                casePointer = 2;
            } else { //  point is within the thresholds -> just compute move cost
                rCostsAll[i] = this.computeMoveCost(tsi) + rCostsAll[i+1];
                casePointer = 0;
            }

            rCostsAll[i] = rCostsAll[i];

        }
        this.remainingCost = rCostsAll;
        return rCostsAll;
    }


// the remaining cost per point are stored in i+1 entry of the remaining cost entry
    public double getRemainingCost(int idx) {
        return this.remainingCost[idx+1];
    }

    private double computeMoveCost(double currentPoint) {
        return Math.abs(currentPoint - this.constantTs);
    }

    private double computeMergeCost(double currentPoint, double predecessor) {
        // two times c are added for merge and split operation. (the time series are of equal length)
        // depending on how the predeceasing point is located, some extra move cost need to be included
        // die move kosten eines Peaks werden stückweise bezahlt von rechts nach links
        return Math.max(0, Math.abs(this.constantTs-currentPoint) -Math.abs(this.constantTs-predecessor) ) + 2 * this.c;
    }

    public double getTotalCost() {
        return remainingCost[0];
    }

}
