package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.PairInt;
import java.security.SecureRandom;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PISortTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     *
     * @param testName
     */
    public PISortTest(String testName) {
        super(testName);
    }

    /**
     *
     * @throws Exception
     */
    public void testMergeSort0() throws Exception {
        
        int nPoints = 3;
        
        PairInt[] points = new PairInt[nPoints];
        points[0] = new PairInt(10, 100);
        points[1] = new PairInt(10, 90);
        points[2] = new PairInt(1, 1);
        
        PISort.mergeSortByXThenY(points);
        
        assertOrdered(points);
    }
    
    /**
     *
     * @throws Exception
     */
    public void testQuickSort0() throws Exception {
        
        int nPoints = 3;
        
        PairInt[] points = new PairInt[nPoints];
        points[0] = new PairInt(10, 100);
        points[1] = new PairInt(10, 90);
        points[2] = new PairInt(1, 1);
        
        PISort.quickSortByXThenY(points);
        
        assertOrdered(points);
    }
    
    /**
     *
     * @throws Exception
     */
    public void testMergeSort1() throws Exception {
        
        PairInt[] points = generateRandomPoints(100);
        
        PISort.mergeSortByXThenY(points);
        
        assertOrdered(points);
    }
    
    /**
     *
     * @throws Exception
     */
    public void testQuickSort1() throws Exception {
        
        PairInt[] points = generateRandomPoints(100);
        
        PISort.quickSortByXThenY(points);
        
        assertOrdered(points);
    }
    
    private PairInt[] generateRandomPoints(int nPoints) throws Exception {
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        log.info("SEED=" + seed);
        
        PairInt[] points = new PairInt[nPoints];
        
        for (int i = 0; i < nPoints; ++i) {
            int x = sr.nextInt(3000);
            int y = sr.nextInt(3000);
            points[i] = new PairInt(x, y);
        }
        
        return points;
    }

    private void assertOrdered(PairInt[] points) {
        
        PairInt lastP = points[0];
        
        for (int i = 1; i < points.length; ++i) {
            
            PairInt p = points[i];
            
            boolean ord = (p.getX() > lastP.getX()) ||
                ((p.getX() == lastP.getX()) && (p.getY() > lastP.getY()));
            
            assertTrue(ord);
        }
    }
    
}
