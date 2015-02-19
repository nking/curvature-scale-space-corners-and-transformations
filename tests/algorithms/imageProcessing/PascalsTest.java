package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class PascalsTest {
    
    public PascalsTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
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
