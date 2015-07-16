package algorithms.compGeometry;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PointPartitionerTest extends TestCase {
    
    public PointPartitionerTest() {
    }

    public void testRandomSubsets() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
            
        int n = 100;
        
        Set<PairInt> allPoints = new HashSet<PairInt>();
        
        PairIntArray set1 = new PairIntArray();
        
        for (int i = 0; i < n; ++i) {
            int x = sr.nextInt(1000);
            int y = sr.nextInt(1000);
            PairInt p = new PairInt(x, y);
            while (allPoints.contains(p)) {
                x = sr.nextInt(1000);
                y = sr.nextInt(1000);
                p = new PairInt(x, y);
            }
            set1.add(x, y);
            allPoints.add(p);
        }
        
        PointPartitioner partitioner = new PointPartitioner();
        List<PairIntArray> subsets = partitioner.randomSubsets(set1, 30);
        
        assertTrue(subsets.size() == 4);
        
        int nPoints = 0;
        for (PairIntArray subset : subsets) {
            nPoints += subset.getN();
        }
        assertTrue(nPoints == allPoints.size());
        
        for (PairIntArray subset : subsets) {
            
            assertTrue(subset.getN() <= 30);
            
            for (int i = 0; i < subset.getN(); ++i) {
                PairInt p = new PairInt(subset.getX(i), subset.getY(i));
                boolean removed = allPoints.remove(p);
                assertTrue(removed);
            }
        }
        
        assertTrue(allPoints.isEmpty());
    }
    
}
