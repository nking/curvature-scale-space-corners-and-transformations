package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LowPassFilterTest extends TestCase {
    
    public LowPassFilterTest() {
    }

    public void testLowpassfilter() {
       
        int nCols = 8;
        int nRows = 4;
        
        LowPassFilter lpFilter = new LowPassFilter();
        
        double[][] lp = lpFilter.lowpassfilter(nRows, nCols, 0.45f, 15);
        
        double[][] expected = new double[4][];
        expected[0] = new double[]{1.00000000e+00, 1.00000000e+00, 
            9.99999978e-01, 9.95804952e-01, 4.06672274e-02, 9.95804952e-01,
            9.99999978e-01, 1.00000000e+00};
        expected[1] = new double[]{9.99999978e-01, 9.99999376e-01,
            9.99280614e-01, 4.88445808e-01, 1.48928501e-03, 4.88445808e-01,
            9.99280614e-01, 9.99999376e-01};
        expected[2] = new double[]{4.06672274e-02, 1.67875976e-02,
            1.48928501e-03, 5.24749584e-05, 1.29367381e-06, 5.24749584e-05,
            1.48928501e-03, 1.67875976e-02};
        expected[3] = new double[]{9.99999978e-01, 9.99999376e-01, 
            9.99280614e-01, 4.88445808e-01, 1.48928501e-03, 4.88445808e-01,
            9.99280614e-01, 9.99999376e-01};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertTrue(Math.abs(expected[i][j] - lp[i][j]) < 0.1);
            }
        }
    }
    
}
