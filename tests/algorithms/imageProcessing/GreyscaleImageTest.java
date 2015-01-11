package algorithms.imageProcessing;

import java.security.SecureRandom;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class GreyscaleImageTest {
    
    public GreyscaleImageTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSetValue() throws Exception {
        
        int nPoints = 1000;
        int width = (int)Math.sqrt(nPoints);
        int height = width;
        
        GreyscaleImage img = new GreyscaleImage(width, height);
         
        assertTrue(img.getWidth() == width);
        assertTrue(img.getHeight() == height);
        
        int[] x = new int[nPoints];
        int[] y = new int[nPoints];

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);

        for (int i = 0; i < nPoints; i++) {
            x[i] = sr.nextInt(width);
            y[i] = sr.nextInt(height);
            img.setValue(x[i], y[i], 1);
            assertTrue(img.getValue(x[i], y[i]) == 1);
        }
        
        
        boolean boundaryExceptionCaught = false;
        try {
            img.setValue(width + 1, 0, 1);
        } catch(Throwable t) {
            boundaryExceptionCaught = true;
        }
        assertTrue(boundaryExceptionCaught);
        
        boundaryExceptionCaught = false;
        try {
            img.setValue(0, height + 1, 1);
        } catch(Throwable t) {
            boundaryExceptionCaught = true;
        }
        assertTrue(boundaryExceptionCaught);
       
        GreyscaleImage imgCopy = img.copyImage();
        assertTrue(imgCopy.getWidth() == width);
        assertTrue(imgCopy.getHeight() == height);
        for (int col = 0; col < imgCopy.getWidth(); col++) {
            for (int row = 0; row < imgCopy.getHeight(); row++) {
                assertTrue(img.getValue(col, row) == imgCopy.getValue(col, row));
            }
        }
        assertTrue(imgCopy.getNPixels() == img.getNPixels());
        
        Image imageCopy = img.copyImageToGreen();
        for (int col = 0; col < imgCopy.getWidth(); col++) {
            for (int row = 0; row < imgCopy.getHeight(); row++) {
                assertTrue(imageCopy.getG(col, row) == imgCopy.getValue(col, row));
            }
        }
        
        GreyscaleImage emptyImage = new GreyscaleImage(width, height);
        img.resetTo(emptyImage);
        for (int col = 0; col < imgCopy.getWidth(); col++) {
            for (int row = 0; row < imgCopy.getHeight(); row++) {
                assertTrue(img.getValue(col, row) == 0);
            }
        }
    }

    public void testInternalIndexUse() throws Exception {
        
        GreyscaleImage img = new GreyscaleImage(10, 19);
        
        for (int i = 0; i < img.getNPixels(); i++) {
            img.setValue(i, i);
        }
        
        for (int i = 0; i < img.getNPixels(); i++) {
            int v = img.getValue(i);
            assertTrue(v == i);
        }
    }
}
