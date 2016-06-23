package algorithms.search;

import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class KNearestNeighbors2DTest extends TestCase {
    
    public KNearestNeighbors2DTest() {
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
        
        KNearestNeighbors2D knn2D = new
            KNearestNeighbors2D(points, maxX, maxY);
        
        // reuse the vector
        FixedSizeSortedVector<PairDistance> vec =
            new FixedSizeSortedVector<PairDistance>(
                2*k, PairDistance.class);
        
        knn2D.findClosest(4, 4, 2, vec);
        
        assertEquals(5, vec.getNumberOfItems());
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 4));
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        PairDistance[] r = vec.getArray();
        for (int i = 0; i < k; ++i) {
            PairDistance pd = r[i];
            assertTrue(expected.remove(pd.p2));
        }
       
        assertEquals(0, expected.size());

        PairDistance[] output = new PairDistance[k];
        
        knn2D.findClosestUsingGreedy(4, 4, 2, k, output);
    
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(4, 4));
        expected.add(new PairInt(4, 2));
        expected.add(new PairInt(4, 6));
        expected.add(new PairInt(2, 4));
        expected.add(new PairInt(6, 4));
    
        r = output;
        for (int i = 0; i < k; ++i) {
            PairDistance pd = r[i];
            assertTrue(expected.remove(pd.p2));
        }
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
