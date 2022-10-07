package HeuristicsUnequalLengths;

import MSMDistances.*;


public class ConstantUB {

    public static double computeUB(double[] ts1, double[] ts2, double c) {

        int m = ts1.length;
        int n = ts2.length;

        MsmConstant msmConstantTs1 = new MsmConstant(ts1, c, 0.);
        MsmConstant msmConstantTs2 = new MsmConstant(ts2, c, 0.);

        double distTs1 = msmConstantTs1.computeRemainingCostPerEntry()[0];
        double distTs2 = msmConstantTs2.computeRemainingCostPerEntry()[0];


        double additionalMergeCost = Math.abs(m-n)*c;

        return distTs1+distTs2+additionalMergeCost;

    }

    public static void main(String[] args) {
        double[] ts1 = new double[]{4, 5, 6, 7, 8, 1, 3, 5, 4, 1, 2, 1, 5, 8};
        double[] ts2 = new double[]{3, 2, 1, 0, 1, 7, 5};
        double c = 0.1;
        System.out.println(computeUB(ts1, ts2, c));
    }
}
