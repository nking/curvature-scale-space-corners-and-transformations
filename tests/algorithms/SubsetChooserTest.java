package algorithms;

import algorithms.misc.MiscMath;
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
    
    public void testGetNextSubset2() throws Exception {
        
        // k=2
        
        int n = 68;
        long np = (MiscMath.computeNDivNMinusK(n, 2)) >> 1;
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int i0 = 0; i0 < n; ++i0) {
            for (int i1 = (i0 + 1); i1 < n; ++i1) {
                PairInt p = new PairInt(i0, i1);
                expected.add(p);
            }
        }
        
        System.out.println("expected.size=" + expected.size() + " np=" + np);
        
        assertTrue(expected.size() == np);
        
        int[] selectedIndexes = new int[2];
        
        SubsetChooser chooser = new SubsetChooser(n, 2);
        
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
