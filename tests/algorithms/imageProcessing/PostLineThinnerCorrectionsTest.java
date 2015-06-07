package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PostLineThinnerCorrectionsTest extends TestCase {
    
    public PostLineThinnerCorrectionsTest() {
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
        
        int w = 10;
        int h = 10;
        GreyscaleImage img = new GreyscaleImage(w, h);
        img.fill(200);
        
        PostLineThinnerCorrections pc = new PostLineThinnerCorrections();
        
        assertTrue(pc.hasAtLeastOneBorderPixel(img));
        
        img = new GreyscaleImage(w, h);
        assertFalse(pc.hasAtLeastOneBorderPixel(img));
    }
    
    public void testCorrectForArtifacts() {
        
    }
    
    public void testCorrectForHoleArtifacts1() throws Exception {
        
        /* 
                      1                   7
                      1*              2   6
                 1    0    1          1   5 
                      1               0   4
                      1              -1   3
        
           -2   -1    0    1    2
            2   3     4    5    6
        */ 
        
        int w = 10;
        int h = 10;
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int col = 3; col <= 5; col++) {
            if (col != 4) {
                points.add(new PairInt(col, 5));
            }
        }
        for (int row = 3; row <= 7; row++) {
            if (row != 5) {
                points.add(new PairInt(4, row));
            }
        }
        
        assertTrue(points.size() == 6);
                
        PostLineThinnerCorrections pc = new PostLineThinnerCorrections();
        
        pc.correctForHoleArtifacts1(points, w, h);
        
        /* 
                      1                   7
                      1*              2   6
                 1    0    1          1   5 
                      1               0   4
                      1              -1   3
        
           -2   -1    0    1    2
            2   3     4    5    6
        
        
                      1                   7
                      1               2   6
                      1               1   5 
                      1*              0   4
                      1              -1   3
        
           -2   -1    0    1    2
            2    3    4    5    6
        */
        for (int row = 3; row <= 7; row++) {
            assertTrue(points.contains(new PairInt(4, row)));
        }
        
        assertFalse(points.contains(new PairInt(3, 5)));
        assertFalse(points.contains(new PairInt(5, 5)));
    }

}
