package algorithms.compGeometry.clustering;

import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class FixedDistanceGroupFinderTest extends TestCase {
    
    public FixedDistanceGroupFinderTest() {
    }

    public void testFindGroupsOfPoints() {
        
        /*
        5
        4             *
        3
        2    *
        1 *  *        *
        0
          0  1  2  3  4  5
        */
        float[] x = new float[]{0, 1, 1, 4, 4};
        float[] y = new float[]{1, 1, 2, 1, 4};
        
        FixedDistanceGroupFinder groupFinder = new FixedDistanceGroupFinder(x, y);
        
        groupFinder.findGroupsOfPoints(1);
        
        int[] gn = groupFinder.getGroupNumbers();
        for (int i = 0; i < gn.length; ++i) {
            assertTrue(gn[i] > -1);
        }
        assertTrue(gn[0] == gn[1]);
        assertTrue(gn[0] == gn[2]);
        
        List<Set<Integer>> sortedGroupList = groupFinder.getDescendingSortGroupList();
        assertTrue(sortedGroupList.size() == 3);
        assertTrue(sortedGroupList.get(0).size() == 3);
        assertTrue(sortedGroupList.get(1).size() == 1);
        assertTrue(sortedGroupList.get(2).size() == 1);
    }

}
