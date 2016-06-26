package algorithms.search;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class NearestNeighbor2DTest extends TestCase {
    
    public NearestNeighbor2DTest() {
    }
    
    public void test0() {
        
        // simple 10 x 10 grid with gaps of 1
        Set<PairInt> points = getTestData();
        
        /*
        0 2 4 6 8
        2
        4   
        6
        8
        */
        
        int maxX = 10;
        int maxY = 10;
        int k = 5;

        points.remove(new PairInt(4, 4));
        
        NearestNeighbor2D knn2D = new
            NearestNeighbor2D(points, maxX, maxY);
        
        Set<PairInt> nearest = knn2D.findClosest(4, 4);
        
        assertEquals(4, nearest.size());
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
        // -----------------------------------
        nearest = knn2D.findClosest(4, 4, 2);
        
        assertEquals(4, nearest.size());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
        
        // -----------------------------------
        points.add(new PairInt(4, 4));
        
        knn2D = new NearestNeighbor2D(points, maxX, maxY);
        
        nearest = knn2D.findClosest(4, 4);
        
        assertEquals(1, nearest.size());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 4));
    
        for (PairInt p2 : nearest) {
            assertTrue(expected.remove(p2));
        }
       
        assertEquals(0, expected.size());
    }
    
    private Set<PairInt> getTestData() {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = 0; i < 10; i += 2) {
            for (int j = 0; j < 10; j += 2) {
                PairInt p = new PairInt(i, j);
                points.add(p);
            }
        }
        
        return points;
    }
    
}
