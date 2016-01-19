package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageTest extends TestCase {
    
    public void test0() throws Exception {
        
        // simple read and write tests
        
        for (int i = 0; i < 3; ++i) {
        
            int w = 3;
            int h = 2;
            Image img;
            
            if (i == 0) {
                img = new Image(w, h);
            } else if (i == 1) {
                img = new Image(w, h, true);
            } else {
                img = new Image(w, h, false);
            }

            assertEquals(w*h, img.getNPixels());
            assertEquals(w, img.getWidth());
            assertEquals(h, img.getHeight());
            img.debugPrint();

            int r = 10;
            int g = 200;
            int b = 34;
            int rgb = (((r & 0x0ff) << 16) | ((g & 0x0ff) << 8) | 
                (b & 0x0ff));

            img.setRGB(1, 1, rgb);        
            assertEquals(rgb, img.getRGB(1, 1));
            assertEquals(r, img.getR(1, 1));
            assertEquals(g, img.getG(1, 1));
            assertEquals(b, img.getB(1, 1));

            int pixIdx = img.getInternalIndex(1, 1);
            assertEquals(r, img.getR(pixIdx));
            assertEquals(g, img.getG(pixIdx));
            assertEquals(b, img.getB(pixIdx));

            Image img2 = img.copyImage();
            assertEquals(img.getNPixels(), img2.getNPixels());
            assertEquals(w, img2.getWidth());
            assertEquals(h, img2.getHeight());
            assertEquals(rgb, img2.getRGB(1, 1));
            assertEquals(r, img2.getR(1, 1));
            assertEquals(g, img2.getG(1, 1));
            assertEquals(b, img2.getB(1, 1));        
            assertEquals(r, img2.getR(pixIdx));
            assertEquals(g, img2.getG(pixIdx));
            assertEquals(b, img2.getB(pixIdx));


            ImageExt img3 = img.copyToImageExt();
            assertEquals(img.getNPixels(), img3.getNPixels());
            assertEquals(w, img3.getWidth());
            assertEquals(h, img3.getHeight());
            assertEquals(rgb, img3.getRGB(1, 1));
            assertEquals(r, img3.getR(1, 1));
            assertEquals(g, img3.getG(1, 1));
            assertEquals(b, img3.getB(1, 1));        
            assertEquals(r, img3.getR(pixIdx));
            assertEquals(g, img3.getG(pixIdx));
            assertEquals(b, img3.getB(pixIdx));

            GreyscaleImage img4 = img.copyToGreyscale();
            assertEquals(img.getNPixels(), img4.getNPixels());
            assertEquals(w, img4.getWidth());
            assertEquals(h, img4.getHeight());
            System.out.println("img4:");
            img4.debugPrint();
            assertTrue(img4.getValue(1, 1) > 0);
            assertTrue(img4.getValue(pixIdx) > 0);

            GreyscaleImage img5 = img.copyToGreyscale();
            assertEquals(img.getNPixels(), img5.getNPixels());
            assertEquals(w, img5.getWidth());
            assertEquals(h, img5.getHeight());
            assertTrue(img5.getValue(1, 1) > 0);
            assertTrue(img5.getValue(pixIdx) > 0);
            assertEquals(img.getNPixels(), img5.getNPixels());

            r = 200;
            img2.setRGB(1, 1, r, g, b);
            assertEquals(r, img2.getR(1, 1));
            assertEquals(g, img2.getG(1, 1));
            assertEquals(b, img2.getB(1, 1));        
            assertEquals(r, img2.getR(pixIdx));

            img.resetTo(img2);
            assertEquals(r, img.getR(1, 1));
            assertEquals(g, img.getG(1, 1));
            assertEquals(b, img.getB(1, 1));        
            assertEquals(r, img.getR(pixIdx));
        }
        
    }
    
    public void testCopy() throws Exception {
        
        String fileName2 = "two_circles_color2.png";
            
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        GreyscaleImage imgR2 = ImageIOHelper.readImage(filePath2).copyRedToGreyscale();
        GreyscaleImage imgG2 = ImageIOHelper.readImage(filePath2).copyGreenToGreyscale();
        GreyscaleImage imgB2 = ImageIOHelper.readImage(filePath2).copyBlueToGreyscale();

        assertEquals(100, imgR2.getWidth());
        assertEquals(100, imgG2.getWidth());
        assertEquals(100, imgB2.getWidth());
        assertEquals(100, imgR2.getHeight());
        assertEquals(100, imgG2.getHeight());
        assertEquals(100, imgB2.getHeight());
        
        // image has a white circle
        
        int nr = 0;
        int ng = 0;
        int nb = 0;
        for (int i = 0; i < imgR2.getWidth(); ++i) {
            for (int j = 0; j < imgR2.getHeight(); ++j) {
                int v = imgR2.getValue(i, j);
                if (v > 0) {
                    nr++;
                }
                v = imgG2.getValue(i, j);
                if (v > 0) {
                    ng++;
                }
                v = imgB2.getValue(i, j);
                if (v > 0) {
                    nb++;
                }
            }
        }   
        assertTrue(((float)nr/(float)imgR2.getNPixels()) > 0.8);
        assertTrue((nr > ng) && (ng > 0));
        assertTrue((ng > nb) && (nb > 0));
    }
    
}
