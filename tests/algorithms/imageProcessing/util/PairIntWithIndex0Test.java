package algorithms.imageProcessing.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PairIntWithIndex0Test extends TestCase {
    
    public PairIntWithIndex0Test() {
    }

    public void testEquals() {
       
        PairIntWithIndex0 p = new PairIntWithIndex0(2, 4, 0);
        
        PairIntWithIndex0 p2 = new PairIntWithIndex0(2, 4, 1);
        
        assertTrue(p.equals(p2));
    }

    public void testHashCode() {
        
        PairIntWithIndex0 p = new PairIntWithIndex0(2, 4, 0);
        
        PairIntWithIndex0 p2 = new PairIntWithIndex0(2, 4, 1);
        
        assertEquals(p.hashCode(), p2.hashCode());
    }

    public void testCopy() {
        
        PairIntWithIndex0 p = new PairIntWithIndex0(2, 4, 0);
        
        PairIntWithIndex0 p2 = (PairIntWithIndex0) p.copy();
        
        assertTrue(p.equals(p2));
        assertEquals(p.hashCode(), p2.hashCode());
    }

}
