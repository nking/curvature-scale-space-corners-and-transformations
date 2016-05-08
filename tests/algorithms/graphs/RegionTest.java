package algorithms.graphs;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RegionTest extends TestCase {
    
    public RegionTest() {
    }
    
    public void test0() throws Exception {
        
        /*
        3           
        2  @  @    
        1  @  @    
        0
           0  1  2  3  4  5  6
        */
        
        Set<PairInt> set0 = new HashSet<PairInt>();
        set0.add(new PairInt(0, 1));
        set0.add(new PairInt(0, 2));
        set0.add(new PairInt(1, 1));
        set0.add(new PairInt(1, 2));
        
        Region region0 = new Region(set0);
        assertTrue(region0.contains(0, 1));
        assertTrue(region0.contains(0, 2));
        assertTrue(region0.contains(new PairInt(1, 1)));
        assertTrue(region0.contains(1, 2));
        
        assertTrue(region0.perimeterContains(0, 1));
        assertTrue(region0.perimeterContains(0, 2));
        assertTrue(region0.perimeterContains(new PairInt(1, 1)));
        assertTrue(region0.perimeterContains(1, 2));
        
        assertEquals(4, region0.size());
        
        
        Set<PairInt> set1 = new HashSet<PairInt>();
        set1.add(new PairInt(0, 3));
        set1.add(new PairInt(1, 3));
        set1.add(new PairInt(2, 1));
        set1.add(new PairInt(2, 2));
        set1.add(new PairInt(2, 3));
        Region region1 = new Region(set1);
        assertEquals(5, region1.size());
        
        region0.mergeIntoThis(region1);
        
        assertEquals(0, region1.size());
        assertEquals(9, region0.size());
        
        /*
        3  #  #  #
        2  @  @  #  
        1  @  @  #  
        0
           0  1  2  3  4  5  6
        */
        assertFalse(region0.perimeterContains(new PairInt(1, 2)));
        assertTrue(region0.perimeterContains(0, 1));
        assertTrue(region0.perimeterContains(0, 2));
        assertTrue(region0.perimeterContains(0, 3));
        assertTrue(region0.perimeterContains(1, 3));
        assertTrue(region0.perimeterContains(2, 3));
        assertTrue(region0.perimeterContains(2, 2));
        assertTrue(region0.perimeterContains(2, 1));
        assertTrue(region0.perimeterContains(1, 1));
        
        assertTrue(region0.contains(new PairInt(1, 2)));
        assertTrue(region0.contains(0, 1));
        assertTrue(region0.contains(0, 2));
        assertTrue(region0.contains(0, 3));
        assertTrue(region0.contains(1, 3));
        assertTrue(region0.contains(2, 3));
        assertTrue(region0.contains(2, 2));
        assertTrue(region0.contains(2, 1));
        assertTrue(region0.contains(1, 1));
        
        assertFalse(region0.contains(new PairInt(10, 2)));
    }
    
}
