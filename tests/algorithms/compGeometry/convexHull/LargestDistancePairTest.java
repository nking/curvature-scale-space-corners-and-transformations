package algorithms.compGeometry.convexHull;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LargestDistancePairTest extends TestCase {
    
    public LargestDistancePairTest() {
    }

    public void testFindLargestDistancePair() throws Exception {
        
        /*            2,6
         *
         *
         *     0,2   2,2 3,2
         *                      7,1
         *            2,0
         *         
         */
        int n = 6;
        long[] x = new long[]{0, 2, 7, 2, 2, 3};
        long[] y = new long[]{2, 2, 1, 6, 0, 2};
        
        long[] expResult = new long[]{7, 1, 2, 6};
        
        long[] result = LargestDistancePair.findLargestDistancePair(x, y);
        
        System.out.printf("result=%s\n", Arrays.toString(result));
    }
    
}
