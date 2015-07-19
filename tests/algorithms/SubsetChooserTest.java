package algorithms;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SubsetChooserTest extends TestCase {
    
    public SubsetChooserTest() {
    }
    
    public void testGetNextSubset() throws Exception {
        
        int n = 3;
        int k = 2;
        
        SubsetChooser chooser = new SubsetChooser(n, k);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(0, 1));
        expected.add(new PairInt(0, 2));
        expected.add(new PairInt(1, 2));
        
        int[] selectedIndexes = new int[k];
        
        while (chooser.getNextSubset(selectedIndexes) != -1) {
            int idx0 = selectedIndexes[0];
            int idx1 = selectedIndexes[1];
            if (idx0 > idx1) {
                int swap = idx0;
                idx0 = idx1;
                idx1 = swap;
            }
            PairInt p = new PairInt(idx0, idx1);
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }
}
