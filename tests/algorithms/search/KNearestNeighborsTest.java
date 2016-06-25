package algorithms.search;

import algorithms.util.PairFloat;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

/**
 *
 * @author nichole
 */
public class KNearestNeighborsTest extends TestCase {
    
    public KNearestNeighborsTest() {
    }

    public void testFindNearest() {
       
        /*
        0 2 4 6 8
        2
        4   
        6
        8
        */
     
        float[] x = new float[25];
        float[] y = new float[25];
        int count = 0;
        
        for (int i = 0; i < 10; i += 2) {
            for (int j = 0; j < 10; j += 2) {
                x[count] = i;
                y[count] = j;
                count++;
            }
        }
        int maxX = 10;
        int maxY = 10;
        int k = 5;
       
        KNearestNeighbors kNN = new KNearestNeighbors(x, y);
        
        List<PairFloat> nearest = kNN.findNearest(k, 4, 4);
        
        assertEquals(5, nearest.size());
        
        Set<PairFloat> expected = new HashSet<PairFloat>();
        expected.add(new PairFloat(4, 4));
        expected.add(new PairFloat(4, 2));
        expected.add(new PairFloat(4, 6));
        expected.add(new PairFloat(2, 4));
        expected.add(new PairFloat(6, 4));
    
        for (PairFloat p : nearest) {
            assertTrue(expected.remove(p));
        }
       
        assertEquals(0, expected.size());
    }    
}
