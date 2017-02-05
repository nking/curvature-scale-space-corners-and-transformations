package algorithms.imageProcessing.features.mser;

import algorithms.util.PairIntArray;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class EllipseHelperTest extends TestCase {
    
    public EllipseHelperTest() {
    }
    
    public void testWithinBounds() {
        
        PairIntArray ellipse = new PairIntArray();
        
        double a = 4;
        double b = 10;
        int xc = 50;
        int yc = 30;
        
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            double x = a * Math.cos(t) + xc;
            double y = b * Math.sin(t) + yc;
            ellipse.add((int)Math.round(x), (int)Math.round(y));
        }
        
        EllipseHelper eh = new EllipseHelper(xc, yc, ellipse);
        
        assertTrue(eh.isWithin(xc, yc));
        assertFalse(eh.isWithin(0, 0));
        assertFalse(eh.isWithin(100, 100));
    }
    
    public void testIntersects() {
        
        PairIntArray ellipse = new PairIntArray();
        
        double a = 4;
        double b = 10;
        int xc = 50;
        int yc = 30;
        
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            double x = a * Math.cos(t) + xc;
            double y = b * Math.sin(t) + yc;
            ellipse.add((int)Math.round(x), (int)Math.round(y));
        }
        
        EllipseHelper eh = new EllipseHelper(xc, yc, ellipse);
        
        PairIntArray ellipse2 = new PairIntArray();
        
        int xc2 = 150;
        int yc2 = 130;
        
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            double x = a * Math.cos(t) + xc2;
            double y = b * Math.sin(t) + yc2;
            ellipse2.add((int)Math.round(x), (int)Math.round(y));
        }
        
        EllipseHelper eh2 = new EllipseHelper(xc2, yc2, ellipse2);
        
        assertFalse(eh.intersects(eh2));
        
        PairIntArray ellipse3 = new PairIntArray();
        
        int xc3 = xc + 2;
        int yc3 = yc - 3;
        
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            double x = a * Math.cos(t) + xc3;
            double y = b * Math.sin(t) + yc3;
            ellipse3.add((int)Math.round(x), (int)Math.round(y));
        }
        
        EllipseHelper eh3 = new EllipseHelper(xc3, yc3, ellipse3);
        
        assertTrue(eh.intersects(eh3));
    }
}
