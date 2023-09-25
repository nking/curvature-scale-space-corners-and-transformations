package algorithms;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class Tmp2Test extends TestCase {
    
    public void test0() {

        // number of throws needed for a fair coin to get 3 consecutive heads?
        int pOutlier = 50; // 50% chance heads, 50% chance tails
        int n95 = RANSACAlgorithmIterations.numberOfSubsamples(95., pOutlier, 3);
        int n99 = RANSACAlgorithmIterations.numberOfSubsamples(99., pOutlier, 3);

        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int nTests = 100; // number of tests to run
        int maxN = 10000; // maximum number of coin tosses in a test
        int t;
        int nC;
        TIntList nThrows = new TIntArrayList(); // for each test, has the number of tosses until 3 consecutive heads, else -1
        for (int nT = 0; nT < nTests; ++nT) {
            // count consecutive number of heads
            int nH = 0;
            for (nC = 0; nC < maxN; ++nC) {
                // t=0 is heads, t=1 is tails.
                t = rand.nextInt(2);
                if (t == 0) {
                    ++nH;
                } else {
                    nH = 0;
                }
                if (nH == 3) {
                    break;
                }
            }
            if (nH == 3) {
                nThrows.add(nC);
            } else {
                nThrows.add(-1);
            }
        }
        System.out.printf("number of throws to get 3 consecutive heads=\n%s\n",
                Arrays.toString(nThrows.toArray()));

        double[] count = MiscMath0.convertIntToDouble(nThrows.toArray());
        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(count);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;
        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(count);
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(count);
        System.out.printf("median of absolute deviation of x, the median, the min, and the max = \n  %s with r0=%.3f r1=%.3f\n",
                FormatArray.toString(mADMinMax, "%.3f"), r0, r1);
        System.out.printf("medianAndIQR = %s\n",
                FormatArray.toString(medianAndIQR, "%.3f"));
        System.out.printf("avgAndStDev = %s\n",
                FormatArray.toString(avgAndStDev, "%.3f"));

        System.out.printf("expected for 95 percent certainty=%d\n", n95);
        System.out.printf("expected for 99 percent certainty=%d\n", n99);
    }
}
