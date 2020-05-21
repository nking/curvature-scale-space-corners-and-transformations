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
        int m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(5);
        assertEquals(3, m);
        
        m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(10);
        assertEquals(5, m);
        
        m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(20);
        assertEquals(13, m);
        
        m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(25);
        assertEquals(21, m);
        
        m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(30);
        assertEquals(35, m);
        
        m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(40);
        assertEquals(106, m);
        
        m = RANSACAlgorithmIterations.numberOfSubsamplesOfSize7For95PercentInliers(50);
        assertEquals(382, m);
    }
    
}
