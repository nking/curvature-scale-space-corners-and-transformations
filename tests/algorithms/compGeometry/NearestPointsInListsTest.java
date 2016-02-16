package algorithms.compGeometry;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class NearestPointsInListsTest extends TestCase {
    
    public NearestPointsInListsTest() {
    }
    
    public void testFindNeighbors() throws Exception {
        
        /*
          5
          4
          3      *
          2   *  *  *           *
          1      *
          0
           0  1  2  3  4  5  6  7  8  9
        */
        List<Set<PairInt>> pointsList = new ArrayList<Set<PairInt>>();
        Set<PairInt> set = new HashSet<PairInt>();
        set.add(new PairInt(1, 2));
        pointsList.add(set);
        set = new HashSet<PairInt>();
        set.add(new PairInt(2, 1));
        pointsList.add(set);
        set = new HashSet<PairInt>();
        set.add(new PairInt(2, 2));
        pointsList.add(set);
        set = new HashSet<PairInt>();
        set.add(new PairInt(2, 3));
        pointsList.add(set);
        set = new HashSet<PairInt>();
        set.add(new PairInt(3, 2));
        pointsList.add(set);
        set = new HashSet<PairInt>();
        set.add(new PairInt(7, 2));
        pointsList.add(set);
        
        NearestPointsInLists np = new NearestPointsInLists(pointsList);
        
        Map<Integer, PairInt> resultsMap = np.findNeighbors(2, 2, 2);
        
        Collection<PairInt> results = resultsMap.values();
        
        assertTrue(results.contains(new PairInt(1, 2)));
        assertTrue(results.contains(new PairInt(2, 1)));
        assertTrue(results.contains(new PairInt(2, 2)));
        assertTrue(results.contains(new PairInt(2, 3)));
        assertTrue(results.contains(new PairInt(3, 2)));
        
        assertFalse(results.contains(new PairInt(7, 2)));
    }
}
