package algorithms.util;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TrioIntTest extends TestCase {
    
    public TrioIntTest() {
    }
    
    public void test0() {
        
        TrioInt t = new TrioInt(1, 2, 3);
        
        assertEquals(1, t.getX());
        assertEquals(2, t.getY());
        assertEquals(3, t.getZ());
    }
}
