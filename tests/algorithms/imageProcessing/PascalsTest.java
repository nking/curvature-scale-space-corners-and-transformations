package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PascalsTest extends TestCase {
    
    public PascalsTest() {
    }
    
    public void testPascals() throws Exception {
        /*
        // fsecond deriv:
        int[] seed0 = new int[]{1, 1, -2, -2, 1, 1};
        int seed0N = 3;
        int finalN = 256;
        GaussianHelperForTests.printPascalsTriangle2(seed0, seed0N, finalN,
            true);
        */
        
        //first deriv:
        int[] seed0 = new int[]{1, 2, 0, -2, -1};
        int seed0N = 4;
        int finalN = 32;
        GaussianHelperForTests.printPascalsTriangle2(seed0, seed0N, finalN,
            true);
    }

}
