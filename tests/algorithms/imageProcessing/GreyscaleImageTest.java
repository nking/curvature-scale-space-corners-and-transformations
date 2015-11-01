package algorithms.imageProcessing;

import java.security.SecureRandom;
import java.util.HashMap;
import java.util.Map;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GreyscaleImageTest extends TestCase {
    
    public GreyscaleImageTest() {
    }
    
    public void test0() throws Exception {
        
        // simple read and write tests
        
        for (int i = 0; i < 3; ++i) {
        
            int w = 3;
            int h = 2;
            GreyscaleImage img;
            
            if (i == 0) {
                img = new GreyscaleImage(w, h);
            } else if (i == 1) {
                img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32);
            } else {
                img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits64);
            }

            assertEquals(w*h, img.getNPixels());
            assertEquals(w, img.getWidth());
            assertEquals(h, img.getHeight());
            img.debugPrint();

            int v11 = 11;
            int v21 = 21;

            img.setValue(1, 1, v11);  
            img.setValue(2, 1, v21);  
            assertEquals(v11, img.getValue(1, 1));
            assertEquals(v21, img.getValue(2, 1));

            int pix11Idx = img.getInternalIndex(1, 1);
            int pix21Idx = img.getInternalIndex(2, 1);
            assertEquals(v11, img.getValue(pix11Idx));
            assertEquals(v21, img.getValue(pix21Idx));

            GreyscaleImage img2 = img.copyImage();
            assertEquals(img.getNPixels(), img2.getNPixels());
            assertEquals(w, img2.getWidth());
            assertEquals(h, img2.getHeight());
            assertEquals(v11, img2.getValue(1, 1));
            assertEquals(v11, img2.getValue(pix11Idx));
            assertEquals(v21, img2.getValue(2, 1));
            assertEquals(v21, img2.getValue(pix21Idx));

            v11 = 200;
            img2.setValue(1, 1, v11);
            assertEquals(v11, img2.getValue(1, 1));
            assertEquals(v11, img2.getValue(pix11Idx));
            assertEquals(v21, img2.getValue(2, 1));
            assertEquals(v21, img2.getValue(pix21Idx));

            img.resetTo(img2);
            assertEquals(v11, img.getValue(1, 1));
            assertEquals(v11, img.getValue(pix11Idx));
            assertEquals(v21, img.getValue(2, 1));
            assertEquals(v21, img.getValue(pix21Idx));
        }
        
    }
    
    public void testFill0() {
        int w = 10;
        int h = 10;
        int v = 50;
        GreyscaleImage img = new GreyscaleImage(w, h);
        img.fill(v);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                assertEquals(v, img.getValue(i, j));
            }
        }
        
        int f = 2;
        img.multiply(f);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                assertEquals(v*f, img.getValue(i, j));
            }
        }
    }
    
    public void testFill1() {
        int w = 10;
        int h = 10;
        int v = -50;
        GreyscaleImage img = new GreyscaleImage(w, h, 
            GreyscaleImage.Type.Bits32Signed);
        img.fill(v);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v2 = img.getValue(i, j);
                assertEquals(v, v2);
            }
        }
        
        int f = 2;
        img.multiply(f);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                assertEquals(v*f, img.getValue(i, j));
            }
        }
    }
    
    public void testFill2() {
        int w = 10;
        int h = 10;
        int v = -50;
        GreyscaleImage img = new GreyscaleImage(w, h, 
            GreyscaleImage.Type.Bits64Signed);
        img.fill(v);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                assertEquals(v, img.getValue(i, j));
            }
        }
        
        int f = 2;
        img.multiply(f);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                assertEquals(v*f, img.getValue(i, j));
            }
        }
    }
    
    public void testImageType0() throws Exception {
       
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        final int w = 100;
        final int h = 100;
        
        Map<Integer, Integer> expected = new HashMap<Integer, Integer>();
        
        for (int i = 0; i < (w * h); ++i) {
            int v = sr.nextInt(256);
            expected.put(Integer.valueOf(i), Integer.valueOf(v));
        }
        
        for (int i = 0; i < 2; ++i) {
            GreyscaleImage img = null;
            if (i == 0) {
                img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32);
            } else if (i == 1) {
                img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits64);
            }
            for (int ii = 0; ii < img.getNPixels(); ++ii) {
                img.setValue(ii, expected.get(Integer.valueOf(ii)).intValue());
            }
            
            for (int ii = 0; ii < img.getNPixels(); ++ii) {
                int v = img.getValue(ii);
                int ve = expected.get(Integer.valueOf(ii)).intValue();
                assertEquals(ve, v);
            }
        }
    }
    
    public void testImageType1() throws Exception {
       
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        final int w = 100;
        final int h = 100;
        
        Map<Integer, Integer> expected = new HashMap<Integer, Integer>();
        
        for (int i = 0; i < (w * h); ++i) {
            int v = sr.nextInt(256) - 255;
            expected.put(Integer.valueOf(i), Integer.valueOf(v));
        }
        
        for (int i = 0; i < 2; ++i) {
            GreyscaleImage img = null;
            if (i == 0) {
                img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32Signed);
            } else if (i == 1) {
                img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits64Signed);
            }
            for (int ii = 0; ii < img.getNPixels(); ++ii) {
                img.setValue(ii, expected.get(Integer.valueOf(ii)).intValue());
            }
            
            for (int ii = 0; ii < img.getNPixels(); ++ii) {
                int v = img.getValue(ii);
                int ve = expected.get(Integer.valueOf(ii)).intValue();
                assertEquals(ve, v);
            }
        }
    }
    
    public void testImageType2() throws Exception {
       
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        final int w = 100;
        final int h = 100;
        
        Map<Integer, Integer> expected = new HashMap<Integer, Integer>();
        
        for (int i = 0; i < (w * h); ++i) {
            int v = sr.nextInt();
            expected.put(Integer.valueOf(i), Integer.valueOf(v));
        }
        
        GreyscaleImage img = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32FullRangeInt);

        for (int ii = 0; ii < img.getNPixels(); ++ii) {
            int v = expected.get(Integer.valueOf(ii)).intValue();
            img.setValue(ii, v);
        }

        for (int ii = 0; ii < img.getNPixels(); ++ii) {
            int v = img.getValue(ii);
            int ve = expected.get(Integer.valueOf(ii)).intValue();
            assertEquals(ve, v);
        }
    }
}
