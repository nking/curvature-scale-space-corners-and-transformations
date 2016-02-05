package algorithms.util;

import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PairIntPairTest extends TestCase {
    
    public PairIntPairTest() {
    }
    
    public void testEquals() throws Exception {
        
        PairIntPair p1 = new PairIntPair(1, 2, 10, 12);
        
        PairIntPair p2 = new PairIntPair(1, 2, 10, 12);
        
        assertEquals(1, p1.getX1());
        assertEquals(2, p1.getY1());
        assertEquals(10, p1.getX2());
        assertEquals(12, p1.getY2());
        
        assertTrue(p1.equals(p2));
        
        p2 = new PairIntPair(10, 12, 1, 2);
        
        assertTrue(p1.equals(p2));
        
        p2 = new PairIntPair(10, 12, 2, 2);
        
        assertFalse(p1.equals(p2));        
    }
    
    public void testSetContains() throws Exception {
        
        PairIntPair p1 = new PairIntPair(1, 2, 10, 12);
        
        PairIntPair p2 = new PairIntPair(10, 12, 1, 2);
        
        Set<PairIntPair> set = new HashSet<PairIntPair>();
        
        set.add(p1);
        
        assertTrue(set.contains(p1));
        assertTrue(set.contains(p2));
        
        p2 = new PairIntPair(10, 12, 2, 2);
        
        assertFalse(set.contains(p2));
    }
}
