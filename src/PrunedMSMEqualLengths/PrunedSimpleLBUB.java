package MsmDistArray;

import DifferentLengths.GreedyUnequalLengths;
import DifferentLengths.ItakuraUnequalLengths;
import Heuristics.Greedy;

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

    public static void main(String[] args) {
        double[] ts1 = new double[]{-0.36292, -0.36533, -0.38781, -0.40559, -0.399, -0.36658, -0.32454, -0.31733, -0.32653, -0.30488, -0.27717, -0.27173, -0.27581, -0.24862, -0.1957, -0.17907, -0.15147, -0.10054, -0.084118, -0.065608, -0.052431, -0.026914, 0.0054007, 0.033009, 0.076304, 0.093246, 0.10506, 0.13602, 0.18768, 0.22596, 0.24258, 0.26308, 0.27772, 0.29006, 0.30805, 0.37268, 0.42465, 0.40928, 0.41294, 0.45655, 0.49336, 0.51103, 0.51428, 0.54094, 0.60233, 0.64291, 0.66612, 0.71046, 0.76087, 0.8163, 0.87716, 0.88793, 0.89985, 0.93091, 0.93008, 0.92788, 0.90655, 0.85907, 0.80804, 0.76913, 0.73797, 0.69771, 0.65587, 0.63161, 0.6541, 0.66006, 0.63329, 0.62136, 0.57943, 0.50622, 0.40865, 0.35417, 0.30272, 0.24196, 0.27155, 0.34256, 0.43689, 0.58298, 0.73326, 0.8232, 0.88302, 0.95737, 1.0084, 1.0608, 1.1154, 1.1657, 1.226, 1.2691, 1.3209, 1.3457, 1.3501, 1.3594, 1.3367, 1.3083, 1.3274, 1.3673, 1.3051, 1.1814, 1.0646, 0.93708, 0.82142, 0.72218, 0.63747, 0.5445, 0.48614, 0.48196, 0.48165, 0.46303, 0.4166, 0.36243, 0.33754, 0.32625, 0.3138, 0.30292, 0.31265, 0.35124, 0.37247, 0.37968, 0.3914, 0.41284, 0.41733, 0.39966, 0.4142, 0.44431, 0.4555, 0.44651, 0.45895, 0.4943, 0.51386, 0.50947, 0.49211, 0.43731, 0.39433, 0.36379, 0.28034, 0.22146, 0.17492, 0.13246, 0.11939, 0.10663, 0.086135, 0.054239, 0.025793, 0.0039366, -0.030783, -0.086942, -0.11675, -0.1135, -0.11863, -0.13839, -0.14415, -0.15816, -0.1865, -0.22017, -0.2573, -0.26692, -0.26567, -0.27497, -0.29599, -0.31555, -0.30112, -0.29076, -0.29076, -0.25991, -0.23502, -0.24308, -0.26211, -0.2734, -0.28585, -0.29777, -0.31639, -0.3372, -0.3464, -0.36742, -0.39534, -0.42598, -0.43205, -0.4036, -0.42034, -0.41929, -0.38227, -0.39335, -0.38844, -0.36721, -0.36198, -0.34504, -0.31795, -0.30666, -0.31691, -0.32747, -0.31544, -0.29913, -0.29882, -0.2894, -0.25646, -0.24778, -0.22896, -0.198, -0.18566, -0.14184, -0.068327, -0.023986, -0.016979, -0.076484, -0.134, -0.13745, -0.15199, -0.15251, -0.12773, -0.093425, -0.055045, -0.05149, -0.063934, -0.089556, -0.11737, -0.15816, -0.20208, -0.2321, -0.24433, -0.20511, -0.1568, -0.17238, -0.21191, -0.26859, -0.31137, -0.3236, -0.3395, -0.3372, -0.31555, -0.27675, -0.23429, -0.22028, -0.1774, -0.12187, -0.080458, -0.032038, 0.029035, 0.078814, 0.13414, 0.22271, 0.3207, 0.44306, 0.53958, 0.57263, 0.6564, 0.77049, 0.86304, 0.88124, 0.74874, 0.48081, 0.16886, -0.20334, -0.41082, -0.56842, -0.63001, -0.67205, -0.69224, -0.70102, -0.69935, -0.71075, -0.76084, -0.74233, -0.73647, -0.76722, -0.76555, -0.75174, -0.72738, -0.67279, -0.59393, -0.64183, -0.59948, -0.52784, -0.48381, -0.46363, -0.4493, -0.41469, -0.39566, -0.33897, -0.24893, -0.18462, -0.16266, -0.11016, -0.0056846, 0.011257, -0.0062075, 0.053193, 0.19877, 0.22857, 0.22878, 0.28138, 0.30073, 0.3253, 0.33032, 0.30857, 0.2979, 0.24363, 0.33555, 0.36086, 0.35846, 0.35804, 0.35741, 0.33242, 0.32436, 0.33942, 0.3069, 0.29979, 0.33095, 0.33001, 0.33451, 0.36609, 0.39192, 0.44118, 0.47611, 0.48771, 0.49378, 0.52275, 0.59595, 0.66299, 0.71883, 0.79005, 0.87266, 0.98153, 1.0767, 1.2143, 1.3547, 1.4473, 1.5495, 1.6294, 1.6838, 1.7758, 1.8382, 1.8637, 1.9083, 1.9632, 1.9752, 1.9048, 1.9058, 1.9231, 1.9071, 1.8778, 1.8578, 1.7831, 1.7087, 1.6883, 1.5238, 1.3502, 1.2929, 1.2024, 1.123, 1.0236, 0.90979, 0.83941, 0.74445, 0.62105, 0.54575, 0.49127, 0.34465, 0.2499, 0.16248, 0.082893, -0.015306, -0.1112, -0.16255, -0.24925, -0.30561, -0.34943, -0.38248, -0.42138, -0.42745, -0.45014, -0.48036, -0.48831, -0.49302, -0.50923, -0.50567, -0.4858, -0.45924, -0.42368, -0.3419, -0.23847, -0.12218, 0.10454, 0.38136, 0.58989, 0.71151, 0.87517, 1.1659, 1.8489, 2.3128, 3.0017, 3.5318, 3.7205, 3.5666, 3.1838, 2.5993, 1.9607, 1.2488, 0.75292, 0.094919, -0.39942, -0.71451, -1.0295, -1.2201, -1.3328, -1.4206, -1.4831, -1.5349, -1.5767, -1.6347, -1.6788, -1.6789, -1.6754, -1.699, -1.7108, -1.7157, -1.726, -1.7346, -1.7503, -1.7653, -1.7538, -1.7688, -1.783, -1.7817, -1.7911, -1.8002, -1.8008, -1.8129, -1.8128, -1.7942, -1.7945, -1.8035, -1.8127, -1.815, -1.8244, -1.8314, -1.8304, -1.8421, -1.8386, -1.8469, -1.8425, -1.8345, -1.8489, -1.833, -1.8386, -1.8701, -1.8554, -1.849, -1.8578, -1.864, -1.8539, -1.8455, -1.8449, -1.8539, -1.8684, -1.8715, -1.8777, -1.8615, -1.8437, -1.8562, -1.8538, -1.8526, -1.8493, -1.8494, -1.8763, -1.8805, -1.8582, -1.8447, -1.8362, -1.8401, -1.8611};
        double[] ts2 = new double[]{-0.11628, -0.11643, -0.10543, -0.089785, -0.065172, -0.050817, -0.039526, -0.026389, -0.022588, -0.013289, 0.0024313, 0.028041, 0.061732, 0.091954, 0.11195, 0.13856, 0.17269, 0.19753, 0.21635, 0.24266, 0.26543, 0.28654, 0.30425, 0.3086, 0.32502, 0.3482, 0.36296, 0.37451, 0.373, 0.38705, 0.40975, 0.41001, 0.41436, 0.42152, 0.42905, 0.4461, 0.45595, 0.46019, 0.47145, 0.49108, 0.50905, 0.5237, 0.55462, 0.57728, 0.59477, 0.62356, 0.65655, 0.67105, 0.67765, 0.67924, 0.66086, 0.64607, 0.62791, 0.60237, 0.58569, 0.57466, 0.56053, 0.53706, 0.50953, 0.48839, 0.46957, 0.46274, 0.46695, 0.46488, 0.46466, 0.47831, 0.49104, 0.48975, 0.47499, 0.44392, 0.41123, 0.39303, 0.37742, 0.35506, 0.35473, 0.37971, 0.41104, 0.44392, 0.48388, 0.52562, 0.55857, 0.584, 0.61097, 0.64714, 0.67581, 0.69109, 0.69625, 0.70282, 0.71016, 0.70064, 0.68385, 0.66684, 0.65035, 0.63599, 0.61798, 0.5968, 0.56145, 0.52853, 0.49359, 0.44838, 0.40399, 0.37196, 0.35133, 0.3306, 0.32506, 0.32222, 0.32668, 0.33222, 0.32167, 0.32063, 0.32111, 0.32561, 0.3427, 0.3599, 0.37049, 0.38141, 0.40366, 0.42931, 0.45477, 0.47455, 0.48994, 0.50182, 0.51484, 0.53049, 0.54448, 0.55555, 0.57031, 0.58341, 0.57953, 0.5778, 0.57806, 0.56983, 0.54935, 0.5175, 0.49499, 0.46222, 0.43104, 0.41687, 0.40126, 0.38621, 0.36809, 0.34927, 0.33384, 0.32417, 0.31414, 0.30119, 0.2862, 0.27713, 0.28285, 0.28713, 0.28144, 0.2762, 0.27233, 0.26469, 0.26163, 0.26011, 0.26059, 0.26823, 0.27059, 0.2682, 0.28314, 0.30314, 0.31012, 0.32716, 0.34192, 0.34211, 0.35078, 0.35558, 0.34964, 0.33879, 0.31956, 0.30668, 0.30639, 0.30152, 0.27709, 0.26495, 0.26823, 0.26369, 0.27026, 0.28259, 0.28672, 0.28063, 0.28761, 0.31167, 0.32447, 0.32975, 0.33178, 0.32938, 0.33465, 0.34879, 0.35167, 0.35665, 0.37687, 0.38455, 0.39521, 0.41827, 0.43894, 0.45289, 0.47078, 0.49197, 0.52695, 0.56673, 0.60016, 0.62769, 0.65444, 0.68496, 0.70175, 0.73042, 0.77072, 0.79256, 0.78227, 0.77965, 0.78386, 0.7767, 0.78304, 0.76109, 0.71673, 0.69009, 0.66584, 0.64304, 0.62226, 0.59488, 0.57193, 0.51828, 0.49602, 0.49012, 0.46868, 0.46669, 0.46746, 0.47451, 0.49495, 0.52352, 0.53746, 0.55994, 0.59677, 0.62595, 0.6495, 0.66629, 0.68972, 0.71352, 0.73581, 0.77541, 0.78131, 0.70308, 0.71289, 0.74448, 0.74456, 0.72828, 0.70669, 0.64038, 0.5251, 0.35724, 0.31779, 0.28856, 0.2762, 0.26605, 0.261, 0.26938, 0.28344, 0.30174, 0.33657, 0.383, 0.41503, 0.43193, 0.46444, 0.50307, 0.56178, 0.62363, 0.66053, 0.72688, 0.82209, 0.90891, 0.98504, 1.0242, 1.0553, 1.1008, 1.1106, 1.158, 1.2015, 1.2296, 1.2569, 1.2572, 1.218, 1.1677, 1.1557, 1.1186, 1.0952, 1.1323, 1.1718, 1.2029, 1.1774, 1.1096, 1.0705, 0.9808, 0.90818, 0.79935, 0.59592, 0.58278, 0.51787, 0.4065, 0.33277, 0.19616, 0.094537, 0.029701, -0.082922, -0.23839, -0.32378, -0.42869, -0.53817, -0.64102, -0.74932, -0.86232, -0.98391, -1.101, -1.21, -1.3114, -1.4173, -1.5403, -1.6547, -1.7535, -1.8384, -1.9189, -1.9884, -2.0558, -2.1074, -2.1991, -2.3491, -2.4473, -2.5454, -2.6899, -2.7983, -2.9016, -3.011, -3.0901, -3.1282, -3.1463, -3.2208, -3.2375, -3.2557, -3.2594, -3.2272, -3.1939, -3.1735, -3.1332, -3.0362, -3.078, -3.0839, -3.0376, -2.9873, -2.9525, -2.9014, -2.8431, -2.7628, -2.7013, -2.6313, -2.5531, -2.4453, -2.3673, -2.2904, -2.2055, -2.1081, -1.9985, -1.9415, -1.862, -1.7784, -1.6948, -1.6395, -1.5921, -1.5319, -1.4389, -1.3797, -1.3041, -1.2161, -1.1503, -1.0714, -0.99863, -0.91914, -0.83184, -0.73058, -0.63714, -0.55246, -0.44603, -0.32581, -0.22499, -0.1427, -0.036869, 0.081732, 0.31668, 0.47661, 0.65817, 0.78777, 0.84777, 0.83736, 0.77323, 0.65086, 0.50385, 0.34609, 0.25491, 0.11144, 0.011214, -0.059046, -0.12148, -0.13333, -0.14174, -0.15174, -0.15879, -0.16473, -0.15204, -0.13621, -0.11912, -0.11536, -0.11562, -0.10484, -0.093918, -0.085542, -0.081003, -0.07842, -0.075209, -0.076944, -0.076944, -0.071335, -0.058456, -0.043216, -0.039636, -0.033658, -0.019968, -0.017459, -0.0068309, -0.0014433, -0.0015171, -0.00088979, 0.0016195, 0.005937, 0.0038336, 0.0016195, -0.0013326, -0.0074582, 0.0027634, 0.01424, 0.015162, 0.015827, 0.017044, 0.01435, 0.012579, 0.016232, 0.01221, 0.0084094, 0.015199, 0.018299, 0.014277, 0.0070809, 0.011694, 0.007044, -0.003805, -0.001185, -0.0063881, -0.018086, -0.013621, -0.01148, -0.013658, -0.013953, -0.015909, -0.023473, -0.03126, -0.032699, -0.037053, -0.047939, -0.058936, -0.052773, -0.057091};
        PrunedSimpleLBUB prunedSimpleLBUB = new PrunedSimpleLBUB(0.1,ts1, ts2, "g");
        double distPrunedWithBand = prunedSimpleLBUB.msmDistPruned()[0];
        System.out.println("pruned with band: " +distPrunedWithBand);
        MsmOneArray msmOneArray = new MsmOneArray(0.1);
        System.out.println(msmOneArray.MSM_Distance(ts1,ts2));
    }
}
