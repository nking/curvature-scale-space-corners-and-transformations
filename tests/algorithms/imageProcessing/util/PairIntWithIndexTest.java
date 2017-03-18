package algorithms.imageProcessing.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PairIntWithIndexTest extends TestCase {
    
    public PairIntWithIndexTest() {
    }

    public void testEquals() {
       
        PairIntWithIndex p = new PairIntWithIndex(2, 4, 0);
        
        PairIntWithIndex p2 = new PairIntWithIndex(2, 4, 1);
        
        assertTrue(p.equals(p2));
    }

    public void testHashCode() {
        
        PairIntWithIndex p = new PairIntWithIndex(2, 4, 0);
        
        PairIntWithIndex p2 = new PairIntWithIndex(2, 4, 1);
        
        assertEquals(p.hashCode(), p2.hashCode());
    }

    public void testCopy() {
        
        PairIntWithIndex p = new PairIntWithIndex(2, 4, 0);
        
        PairIntWithIndex p2 = (PairIntWithIndex) p.copy();
        
        assertTrue(p.equals(p2));
        assertEquals(p.hashCode(), p2.hashCode());
    }

}
