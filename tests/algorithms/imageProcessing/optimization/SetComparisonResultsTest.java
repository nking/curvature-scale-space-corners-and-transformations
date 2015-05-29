package algorithms.imageProcessing.optimization;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SetComparisonResultsTest extends TestCase {
    
    public SetComparisonResultsTest() {
    }
    
    public void test0() throws Exception {
        
        /*
        SetComparisonResults(int nExpected, int nOverrun, int nMatched,
            int nExpectedBorder, int nMatchedBorder)
        */
        
        /*
        -- test that overrun always produces !best
        -- test that no overrun, and border match chooses the largest border
           match
        -- test that no overrun, and no border match chooses the largest
           expected sky match
        -- test that no pixels matching produces equal
        */
        
        int nExpected = 100;
        int nOverrun = 1;
        int nMatched = 90;
        int nExpectedBorder = 10;
        int nMatchedBorder = 0;
        
        SetComparisonResults s0 = new SetComparisonResults(nExpected, nOverrun, 
            nMatched, nExpectedBorder, nMatchedBorder);
        
        SetComparisonResults s1 = new SetComparisonResults(nExpected, nOverrun, 
            nMatched, nExpectedBorder, nMatchedBorder);
        
        assertTrue(s0.equals(s1));
        
        assertTrue(s0.equals(s0));
        
        assertTrue(s1.equals(s1));
        
        assertTrue(s0.compareTo(s1) == 0);
        
        //----------------------
        // the difference in overrun is smaller than eps
        s0 = new SetComparisonResults(nExpected, 0, 
            nMatched, nExpectedBorder, nMatchedBorder);
        
        s1 = new SetComparisonResults(nExpected, nOverrun, 
            nMatched, nExpectedBorder, nMatchedBorder);
        
        assertTrue(s0.equals(s1));
        
        assertTrue(s0.compareTo(s1) == 0);
        
        assertTrue(s1.equals(s0));
        
        assertTrue(s1.compareTo(s0) == 0);
        
        //----------------------
        nExpected = 100;
        nOverrun = 0;
        nMatched = 90;
        nExpectedBorder = 10;
        nMatchedBorder = 10;
        
        s0 = new SetComparisonResults(nExpected, 0, 
            nMatched, nExpectedBorder, nMatchedBorder);
        
        s1 = new SetComparisonResults(nExpected, nOverrun, 
            nMatched, nExpectedBorder, nMatchedBorder/2);
        
        assertFalse(s0.equals(s1));
        
        assertTrue(s0.compareTo(s1) == -1);
        
        assertFalse(s1.equals(s0));
        
        assertTrue(s1.compareTo(s0) == 1);
        
        //----------------------
        nMatched = 90; 
        
        s0 = new SetComparisonResults(
            nExpected, nOverrun, 
            nMatched, nExpectedBorder, nMatchedBorder);
        
        s1 = new SetComparisonResults(nExpected, nOverrun, 
            nMatched + 5, nExpectedBorder, nMatchedBorder);
        
        assertFalse(s0.equals(s1));
        
        assertTrue(s0.compareTo(s1) == 1);
        
    }
}
