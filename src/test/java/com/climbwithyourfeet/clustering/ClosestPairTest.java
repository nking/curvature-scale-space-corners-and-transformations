package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.PairFloat;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ClosestPairTest extends TestCase {
    
    public ClosestPairTest(String testName) {
        super(testName);
    }

    public void testBruteForceMinDistance() {
        
        List<Float> xPoints = new ArrayList<Float>();
        List<Float> yPoints = new ArrayList<Float>();        
        xPoints.add(10.f); yPoints.add(10.f);
        xPoints.add(11.f); yPoints.add(11.f);
        
        List<PairFloat> p = new ArrayList<PairFloat>();
        p.add(new PairFloat(10.f, 10.f));
        p.add(new PairFloat(11.f, 11.f));
        
        ClosestPair instance = new ClosestPair();
        ClosestPair.ClosestPairFloat result = instance.bruteForceMinDistance(p);
        
        assertNotNull(result);
        
        float dist = (float)Math.sqrt(2.);
        
        assertTrue(Math.abs(result.separation - dist) < 0.1);
        PairFloat r0 = result.point0;
        PairFloat r1 = result.point1;
        
        boolean found0 = false;
        boolean found1 = false;
        if ((r0.getX() == 10.f) && (r0.getY() == 10.f)) {
            found0 = true;
        } else if ((r1.getX() == 10.f) && (r1.getY() == 10.f)) {
            found0 = true;
        }
        
        if ((r0.getX() == 11.f) && (r0.getY() == 11.f)) {
            found1 = true;
        } else if ((r1.getX() == 11.f) && (r1.getY() == 11.f)) {
            found1 = true;
        }
        
        assertTrue(found0);
        assertTrue(found1);
    }
    
    public void testFindClosestPair() {
        
        List<Float> xPoints = new ArrayList<Float>();
        List<Float> yPoints = new ArrayList<Float>();
        
        xPoints.add(100.f); yPoints.add(100.f);
        xPoints.add(110.f); yPoints.add(200.f);
        xPoints.add(10.f); yPoints.add(10.f);
        xPoints.add(11.f); yPoints.add(11.f);
        xPoints.add(300.f); yPoints.add(100.f);
        xPoints.add(410.f); yPoints.add(200.f);
        xPoints.add(500.f); yPoints.add(100.f);
        xPoints.add(510.f); yPoints.add(200.f);
        xPoints.add(210.f); yPoints.add(300.f);
        
        ClosestPair instance = new ClosestPair();
                
        ClosestPair.ClosestPairFloat result = instance.findClosestPair(xPoints, 
            yPoints);
        
        assertNotNull(result);
        
        float dist = (float)Math.sqrt(2.);
        
        assertTrue(Math.abs(result.separation - dist) < 0.1);
        PairFloat r0 = result.point0;
        PairFloat r1 = result.point1;
        
        boolean found0 = false;
        boolean found1 = false;
        if ((r0.getX() == 10.f) && (r0.getY() == 10.f)) {
            found0 = true;
        } else if ((r1.getX() == 10.f) && (r1.getY() == 10.f)) {
            found0 = true;
        }
        
        if ((r0.getX() == 11.f) && (r0.getY() == 11.f)) {
            found1 = true;
        } else if ((r1.getX() == 11.f) && (r1.getY() == 11.f)) {
            found1 = true;
        }
        
        assertTrue(found0);
        assertTrue(found1);
    }
    
}
