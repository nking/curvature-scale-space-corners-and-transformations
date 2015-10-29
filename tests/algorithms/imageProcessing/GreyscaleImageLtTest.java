package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GreyscaleImageLtTest extends TestCase {
    
    public GreyscaleImageLtTest() {
    }
    
    public void test0() throws Exception {
        
        // simple read and write tests
        
        for (int i = 0; i < 3; ++i) {
        
            int w = 3;
            int h = 2;
            GreyscaleImageLt img;
            
            if (i == 0) {
                img = new GreyscaleImageLt(w, h);
            } else if (i == 1) {
                img = new GreyscaleImageLt(w, h, true);
            } else {
                img = new GreyscaleImageLt(w, h, false);
            }

            assertEquals(w*h, img.getNPixels());
            assertEquals(w, img.getWidth());
            assertEquals(h, img.getHeight());
            img.debugPrint();

            int v = 10;

            img.setValue(1, 1, v);        
            assertEquals(v, img.getValue(1, 1));

            int pixIdx = img.getInternalIndex(1, 1);
            assertEquals(v, img.getValue(pixIdx));

            GreyscaleImageLt img2 = img.copyImage();
            assertEquals(img.getNPixels(), img2.getNPixels());
            assertEquals(w, img2.getWidth());
            assertEquals(h, img2.getHeight());
            assertEquals(v, img2.getValue(1, 1));
            assertEquals(v, img2.getValue(pixIdx));

            v = 200;
            img2.setValue(1, 1, v);
            assertEquals(v, img2.getValue(1, 1));
            assertEquals(v, img2.getValue(pixIdx));

            img.resetTo(img2);
            assertEquals(v, img.getValue(1, 1));
            assertEquals(v, img.getValue(pixIdx));
        }
        
    }
    
}
