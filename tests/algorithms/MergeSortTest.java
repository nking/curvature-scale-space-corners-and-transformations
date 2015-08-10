package algorithms;

import java.util.Arrays;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MergeSortTest extends TestCase {
    
    public MergeSortTest() {
    }
    
    public void testSortByDecr() throws Exception {
        
        int[] a = new int[]{1, 2, 3, 4, 5, 6};

    	MergeSort.sortByDecr(a);

    	int[] expectedA = new int[]{6, 5, 4, 3, 2, 1};
        
        assertTrue(Arrays.equals(expectedA, a));
    }
}
