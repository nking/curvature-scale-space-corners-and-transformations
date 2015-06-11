package algorithms.compGeometry;

import algorithms.util.PairInt;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class NearestPointsTest extends TestCase {
    
    public NearestPointsTest() {
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
        int[] x = new int[]{1, 2, 2, 2, 3, 7};
        int[] y = new int[]{2, 1, 2, 3, 2, 2};
        
        NearestPoints nearestPoints = new NearestPoints(x, y);
         
        Set<PairInt> results = nearestPoints.findNeighbors(2, 2, 2);
        
        assertTrue(results.contains(new PairInt(1, 2)));
        assertTrue(results.contains(new PairInt(2, 1)));
        assertTrue(results.contains(new PairInt(2, 2)));
        assertTrue(results.contains(new PairInt(2, 3)));
        assertTrue(results.contains(new PairInt(3, 2)));
        
        assertFalse(results.contains(new PairInt(7, 2)));
    }
}
