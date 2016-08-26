package algorithms.search;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NearestNeighbor1DTest extends TestCase {
    
    public NearestNeighbor1DTest() {
    }
    
    public void test0() {
        
        // 1  3  5  23  67 89  99
        NearestNeighbor1D nn = new NearestNeighbor1D(100);
        nn.insert(1); nn.insert(3); nn.insert(5);
        nn.insert(23); nn.insert(67);
        nn.insert(89); nn.insert(00);
        
        TIntSet set = nn.findClosest(10);
        assertEquals(1, set.size());
        assertEquals(5, set.iterator().next());
    
        TIntSet expected = new TIntHashSet();
        expected.add(3);
        expected.add(5);
        set = nn.findClosest(4);
        assertEquals(expected.size(), set.size());
        TIntIterator iter = set.iterator();
        while (iter.hasNext()) {
            int t = iter.next();
            assertTrue(expected.remove(t));
        }
        assertEquals(0, expected.size());
    }
}
