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
    
    public void testPartition() throws Exception {
        
        int width = 1000;
        int height = 1000;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
            
        int n = 1000;
        
        Set<PairInt> allPoints = new HashSet<PairInt>();
        
        PairIntArray set1 = new PairIntArray();
        
        for (int i = 0; i < n; ++i) {
            int x = sr.nextInt(width);
            int y = sr.nextInt(height);
            PairInt p = new PairInt(x, y);
            while (allPoints.contains(p)) {
                x = sr.nextInt(width);
                y = sr.nextInt(height);
                p = new PairInt(x, y);
            }
            set1.add(x, y);
            allPoints.add(p);
        }
        
        PointPartitioner partitioner = new PointPartitioner();
        PairIntArray[] partitions = partitioner.partition(set1, 2);
        
        assertTrue(partitions.length == 4);
        
        int nPoints = 0;
        for (PairIntArray subset : partitions) {
            nPoints += subset.getN();
        }
        assertTrue(nPoints == allPoints.size());
        
        double nExpectedPerBin = (float)nPoints/4.;
        // 3 times random error:
        double nEps = 3*Math.sqrt(nPoints);
                
        for (PairIntArray subset : partitions) {
             
            assertTrue(subset.getN() >= (nExpectedPerBin - nEps));
            
            for (int i = 0; i < subset.getN(); ++i) {
                PairInt p = new PairInt(subset.getX(i), subset.getY(i));
                boolean removed = allPoints.remove(p);
                assertTrue(removed);
            }
        }
        
        assertTrue(allPoints.isEmpty());
    }
    
    public void testReduceByBinSampling() throws Exception {
        
        int width = 1000;
        int height = 1000;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
            
        int n = 1000;
        
        Set<PairInt> allPoints = new HashSet<PairInt>();
        
        PairIntArray set1 = new PairIntArray();
        
        for (int i = 0; i < n; ++i) {
            int x = sr.nextInt(width);
            int y = sr.nextInt(height);
            PairInt p = new PairInt(x, y);
            while (allPoints.contains(p)) {
                x = sr.nextInt(width);
                y = sr.nextInt(height);
                p = new PairInt(x, y);
            }
            set1.add(x, y);
            allPoints.add(p);
        }
        
        PointPartitioner partitioner = new PointPartitioner();
        PairIntArray[] partitions = partitioner.reduceByBinSampling(set1, 2, 10);
        
        assertTrue(partitions.length == 2);
        
        int nPoints = 0;
        for (PairIntArray subset : partitions) {
            nPoints += subset.getN();
        }
        assertTrue(nPoints == allPoints.size());
        
        PairIntArray primaryPartition = partitions[0];
        assertTrue(primaryPartition.getN() == (4 * 10));
        
        for (PairIntArray subset : partitions) {
                         
            for (int i = 0; i < subset.getN(); ++i) {
                PairInt p = new PairInt(subset.getX(i), subset.getY(i));
                boolean removed = allPoints.remove(p);
                assertTrue(removed);
            }
        }
        
        assertTrue(allPoints.isEmpty());
    }
    
}
