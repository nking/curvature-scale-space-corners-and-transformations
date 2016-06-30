package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ColorHistogramTest extends TestCase {
    
    public ColorHistogramTest() {
    }
    
    public void test0() {
        
        Image img0 = new Image(3, 3);
        Image img1 = new Image(3, 3);
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                img1.setRGB( i, j, 127, 127, 127);
                if (j == 0) {
                    img0.setRGB( i, j, 127, 127, 127);
                } else if (j == 1) {
                    img0.setRGB( i, j, 0, 255, 0);
                } else {
                    img0.setRGB( i, j, 127, 127, 127);
                }
                points.add(new PairInt(i, j));
            }
        }
        
        ColorHistogram ch = new ColorHistogram();
        int[][] hist0 = ch.histogram(img0, points);
        int[][] hist1 = ch.histogram(img1, points);
        
        assertEquals(2, hist0.length);
        assertEquals(16, hist0[0].length);
        
        assertEquals(3, hist0[0][0]);
        assertEquals(6, hist0[0][5]);
        
        assertEquals(3, hist0[1][15]);
        assertEquals(6, hist0[1][5]);
        
        assertEquals(9, hist1[0][5]);
        assertEquals(9, hist1[1][5]);
        
        float sim = ch.intersection(hist0, hist1);
                
        float expected = (6.f/9.f);
                
        assertTrue(Math.abs(sim - expected) < 0.1);
    }
}
