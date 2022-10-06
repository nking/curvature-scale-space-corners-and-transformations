package algorithms.imageProcessing.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RANSACAlgorithmIterationsTest extends TestCase{
    
    public RANSACAlgorithmIterationsTest() {
    }

    /**
     * Test of numberOfSubsamplesOfSize7For95PercentInliers method, of class RANSACAlgorithmIterations.
     */
    public void testNumberOfSubsamplesOfSize7For95PercentInliers() {
        //p |  5%  10%  20%  25%  30%  40%  50%
        //--------------------------------------
        //7 |  3   5    13   21   35   106  382
        int[] outlierPercents = new int[]{5, 10, 20, 25, 30, 40, 50};
        int[] expectedM = new int[]{3, 5, 13, 21, 35, 106, 382};
        int m;
        int i;
        for (i = 0; i < outlierPercents.length; ++i) {
            m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(outlierPercents[i]);
            assertEquals(expectedM[i], m);
            System.out.printf("p=8, outlierPercent=%d, m=%d\n", outlierPercents[i],
                    RANSACAlgorithmIterations.numberOfSubsamplesFor95PercentInliers(outlierPercents[i], 8));
        }

    }
    
}
