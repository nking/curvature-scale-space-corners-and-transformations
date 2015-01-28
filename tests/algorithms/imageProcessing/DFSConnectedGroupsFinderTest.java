package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class DFSConnectedGroupsFinderTest extends TestCase {

    public DFSConnectedGroupsFinderTest(String testName) {
        super(testName);
    }
    
    @Test
    public void testFindConnectedPointGroups() {
        
        System.out.println("findConnectedPointGroups");
        
        /*
        0  1  2  3  4
        1  @  @
        2        @
        3
        4  @  @
        */
        Set<PairInt> points = new HashSet<PairInt>();
        
        PairInt p = new PairInt(1, 1);
        points.add(p);
        p = new PairInt(2, 1);
        points.add(p);
        p = new PairInt(3, 2);
        points.add(p);
        
        p = new PairInt(1, 4);
        points.add(p);
        p = new PairInt(2, 4);
        points.add(p);
        
        int imageWidth = 5;
        int imageHeight = 5;
        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(2);
        finder.findConnectedPointGroups(points, imageWidth, imageHeight);
        
        assertTrue(finder.getNumberOfGroups() == 2);
        
        List<Set<PairInt>> result = finder.getGroupMembershipList();
        
        Set<PairInt> expected = new HashSet<PairInt>(points);
        
        for (Set<PairInt> set : result) {
            for (PairInt p2 : set) {
                assertTrue(expected.contains(p2));
                expected.remove(p2);
            }
        }
        
        assertTrue(expected.isEmpty());
        
        boolean firstIs3 = false;
        int n0 = result.get(0).size();
        int n1 = result.get(1).size();
        if (n0 == 3) {
            firstIs3 = true;
        }
        assertTrue((firstIs3 && (n0 == 3) && (n1 == 2)) || (!firstIs3 && (n0 == 2) && (n1 == 3)));
        
        assertNotNull(finder.getXY(0));
        
        assertNotNull(finder.getXY(1));
        
        if (firstIs3) {
            assertTrue(finder.getXY(0).size() == 3);
            assertTrue(finder.getXY(1).size() == 2);
        } else {
            assertTrue(finder.getXY(1).size() == 3);
            assertTrue(finder.getXY(0).size() == 2);
        }
        
        //----
        
        finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(3);
        finder.findConnectedPointGroups(points, imageWidth, imageHeight);
        
        assertTrue(finder.getNumberOfGroups() == 1);
        
    }

}
