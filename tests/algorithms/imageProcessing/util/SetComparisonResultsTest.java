package algorithms.imageProcessing.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SetComparisonResultsTest extends TestCase {
    
    public SetComparisonResultsTest() {
    }
    
    public void test0() throws Exception {
        
        SetComparisonResults s0 = new SetComparisonResults(100, 1, 90);
        
        SetComparisonResults s1 = new SetComparisonResults(100, 1, 90);
        
        assertTrue(s0.equals(s1));
        
        assertTrue(s0.equals(s0));
        
        assertTrue(s1.equals(s1));
        
        assertTrue(s0.compareTo(s1) == 0);
        
        //----------------------
        s0 = new SetComparisonResults(100, 0, 90);
        
        s1 = new SetComparisonResults(100, 1, 90);
        
        assertFalse(s0.equals(s1));
        
        assertTrue(s0.compareTo(s1) == -1);
        
        assertFalse(s1.equals(s0));
        
        assertTrue(s1.compareTo(s0) == 1);
        
        //----------------------
        s0 = new SetComparisonResults(100, 0, 90);
        
        s1 = new SetComparisonResults(100, 0, 95);
        
        assertFalse(s0.equals(s1));
        
        assertTrue(s0.compareTo(s1) == 1);
        
    }
}
