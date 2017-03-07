package algorithms.imageProcessing;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ConnectedPointsFinderTest extends TestCase {
    
    public ConnectedPointsFinderTest() {
    }
    
    public void test0() {
        
        /*
        4
        3       #  #
        2       #  #
        1 #  # 
        0 #  #
          0  1  2  3  4
        */
        int w = 5;
        int h = 5;
        TIntSet pixIdxs = new TIntHashSet();
        pixIdxs.add(0);
        pixIdxs.add((1 * w) + 0);
        pixIdxs.add((0 * w) + 1);
        pixIdxs.add((1 * w) + 1);
        pixIdxs.add((2 * w) + 2);
        pixIdxs.add((3 * w) + 2);
        pixIdxs.add((2 * w) + 3);
        pixIdxs.add((3 * w) + 3);
        
        ConnectedPointsFinder finder 
            = new ConnectedPointsFinder(w, h);
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(pixIdxs);
        
        assertEquals(2, finder.getNumberOfGroups());
        assertEquals(4, finder.getNumberofGroupMembers(0));
        assertEquals(4, finder.getNumberofGroupMembers(1));
        
        finder = new ConnectedPointsFinder(w, h);
        finder.setMinimumNumberInCluster(1);
        finder.setToUse8Neighbors();
        finder.findConnectedPointGroups(pixIdxs);
        
        assertEquals(1, finder.getNumberOfGroups());
        assertEquals(8, finder.getNumberofGroupMembers(0));
    }
}
