package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
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
    
    public void testGetXY() throws Exception {
        
        GreyscaleImage img = new GreyscaleImage(2, 3);

        PairIntArray xy = new PairIntArray();
        
        img.getXY(xy, 0);
        assertTrue(xy.getN() == 1);
        assertTrue(xy.getX(0) == 0 && xy.getY(0) == 0);
        
        img.getXY(xy, 1);
        assertTrue(xy.getN() == 2);
        assertTrue(xy.getX(1) == 1 && xy.getY(1) == 0);
        
        img.getXY(xy, 2);
        assertTrue(xy.getN() == 3);
        assertTrue(xy.getX(2) == 0 && xy.getY(2) == 1);
        
        img.getXY(xy, 3);
        assertTrue(xy.getN() == 4);
        assertTrue(xy.getX(3) == 1 && xy.getY(3) == 1);
        
        img.getXY(xy, 4);
        assertTrue(xy.getN() == 5);
        assertTrue(xy.getX(4) == 0 && xy.getY(4) == 2);
        
        img.getXY(xy, 5);
        assertTrue(xy.getN() == 6);
        assertTrue(xy.getX(5) == 1 && xy.getY(5) == 2);
        
        /*
        w   h
        2 x 3 image
        
        idx  col  row
        0     0    0
        1     1    0
        2     0    1
        3     1    1
        4     0    2
        5     1    2
        */
    }
    
    public void testFill() throws Exception {
        
        GreyscaleImage img = new GreyscaleImage(2, 2);
        
        img.fill(10);
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            assertTrue(img.getValue(idx) == 10);
        }
    }
}
