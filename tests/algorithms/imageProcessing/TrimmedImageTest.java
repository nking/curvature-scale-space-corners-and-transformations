package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TrimmedImageTest extends TestCase {
    
    public void testTrimmed() {
        
        int w = 10;
        int h = 10;
        
        int rOffsetX = 2;
        int rOffsetY = 3;
        int rW = 3;
        int rH = 4;
        int w2 = 4;
        int h2 = 4;
      
        Image img = new Image(w, h);
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            img.setRGB(i, 0, 0, 255);
        }
        
        for (int i = rOffsetX; i < (rOffsetX + rW); ++i) {
            for (int j = rOffsetY; j < (rOffsetY + rH); ++j) {
                img.setRGB(i, j, 255, 0, 0);
            }
        }
        
        TrimmedImage tImg = new TrimmedImage(img, rOffsetX, 
            rOffsetX + w2, rOffsetY, rOffsetY + h2);
        
        assertEquals(w2, tImg.getTrimmed().getWidth());
        assertEquals(h2, tImg.getTrimmed().getHeight());
        
        assertEquals(w, tImg.getOrigWidth());
        assertEquals(h, tImg.getOrigHeight());
        
        assertEquals(rOffsetX, tImg.getXOffset());
        assertEquals(rOffsetY, tImg.getYOffset());
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 4; ++j) {
                assertEquals(255, tImg.getTrimmed().getR(i, j));
                assertEquals(0, tImg.getTrimmed().getG(i, j));
                assertEquals(0, tImg.getTrimmed().getB(i, j));
            }
        }
         
        for (int i = 3; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                assertEquals(0, tImg.getTrimmed().getR(i, j));
                assertEquals(0, tImg.getTrimmed().getG(i, j));
                assertEquals(255, tImg.getTrimmed().getB(i, j));
            }
        }
        
        Image fullFrameTrimmed = tImg.copyTrimmedToFullFrame();
        
        assertEquals(w, fullFrameTrimmed.getWidth());
        assertEquals(h, fullFrameTrimmed.getHeight());
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (i < 2 || i > 5 || j <= 2 || j >= 7) {
                    assertEquals(0, fullFrameTrimmed.getR(i, j));
                    assertEquals(0, fullFrameTrimmed.getG(i, j));
                    assertEquals(0, fullFrameTrimmed.getB(i, j));
                } else if (i == 5) {
                    assertEquals(0, fullFrameTrimmed.getR(i, j));
                    assertEquals(0, fullFrameTrimmed.getG(i, j));
                    assertEquals(255, fullFrameTrimmed.getB(i, j));
                } else {
                    assertEquals(255, fullFrameTrimmed.getR(i, j));
                    assertEquals(0, fullFrameTrimmed.getG(i, j));
                    assertEquals(0, fullFrameTrimmed.getB(i, j));
                }
            }
        }
        
        /*
        9
        8
        7
        6       r  r  r  b
        5       r  r  r  b
        4       r  r  r  b
        3       r  r  r  b
        2
        1
        0       
          0  1  2  3  4  5  6  7  8  9
        */
        
    }
    
}
