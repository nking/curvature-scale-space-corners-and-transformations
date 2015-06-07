package algorithms.imageProcessing;

import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PostLineThinnerCorrectionsTest extends TestCase {
    
    public PostLineThinnerCorrectionsTest() {
    }

    public void testCorrectForArtifacts() {
        
    }

    public void testSumOver8Neighborhood() {
        
    }

     public void testAddOnePixelBorder() throws Exception {
        
        int w = 10;
        int h = 10;
        GreyscaleImage img = new GreyscaleImage(w, h);
        img.fill(200);
        
        PostLineThinnerCorrections pc = new PostLineThinnerCorrections();
        
        GreyscaleImage img2 = pc.addOnePixelBorders(img);
        
        assertTrue(pc.hasAtLeastOneBorderPixel(img));
        
        assertFalse(pc.hasAtLeastOneBorderPixel(img2));
        
        assertTrue(img.getWidth() == w);
        assertTrue(img.getHeight() == h);
        
        assertTrue(img2.getWidth() == (w + 2));
        assertTrue(img2.getHeight() == (h + 2));
        
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                int v = img.getValue(col, row);
                int v2 = img2.getValue(col + 1, row + 1);
                assertTrue(v == v2);
            }
        }
        
        GreyscaleImage img3 = pc.removeOnePixelBorders(img2);
        
        assertTrue(pc.hasAtLeastOneBorderPixel(img));
        assertTrue(pc.hasAtLeastOneBorderPixel(img3));
        assertFalse(pc.hasAtLeastOneBorderPixel(img2));
        
        assertTrue(img3.getWidth() == w);
        assertTrue(img3.getHeight() == h);
        
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                int v = img.getValue(col, row);
                int v2 = img2.getValue(col + 1, row + 1);
                int v3 = img3.getValue(col, row);
                assertTrue(v == v2);
                assertTrue(v == v3);
            }
        }
    }

    public void testHasAtLeastOneBorderPixel() {
        
    }
    
}
