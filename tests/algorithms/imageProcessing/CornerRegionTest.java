package algorithms.imageProcessing;

import algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CornerRegionTest extends TestCase {
    
    public CornerRegionTest() {
    }
    
    public void test0() throws Exception {
        
        CornerRegion cr = new CornerRegion(0, 3, 1);
        
        cr.set(0, 0.24f, 267, 84);
        cr.set(1, 0.26f, 268, 85);
        cr.set(2, 0.21f, 268, 86);
        
        float orientation = cr.getRelativeOrientationInDegrees();
        
        assertTrue(Math.abs(orientation - 337.5) < 2);
    }
    
    public void test1() throws Exception {
        
        CornerRegion cr = new CornerRegion(0, 3, 1);
        
        cr.set(0, 0.34f, 143, 255);
        cr.set(1, 0.43f, 142, 255);
        cr.set(2, 0.35f, 141, 255);
        
        boolean caughtException = false;
        
        try {
            double orientation = cr.getRelativeOrientation();
        } catch (CornerRegionDegneracyException e) {
            caughtException = true;
        }
        
        assertTrue(caughtException);
    }
    
    public void test2() throws Exception {
        
        CornerRegion cr = new CornerRegion(0, 5, 2);
        
        cr.set(0, 0.13f, 143, 254);
        cr.set(1, 0.34f, 143, 255);
        cr.set(2, 0.43f, 142, 255);
        cr.set(3, 0.35f, 141, 255);
        cr.set(4, 0.28f, 140, 255);
        
        boolean caughtException = false;
        
        try {
            double orientation = cr.getRelativeOrientation();
            assertTrue(Math.abs(orientation - 1) < 1);
        } catch (CornerRegionDegneracyException e) {
            caughtException = true;
        }
        
        assertFalse(caughtException);
    }
}
