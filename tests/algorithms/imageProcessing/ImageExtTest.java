package algorithms.imageProcessing;

import java.awt.Color;
import java.security.SecureRandom;
import junit.framework.TestCase;
import org.junit.Test;


/**
 *
 * @author nichole
 */
public class ImageExtTest extends TestCase {

    public void testInit() {
        
        int w = 10;
        int h = 10;
        ImageExt img = new ImageExt(w, h);
        
        assertNotNull(img);
        
        int n = w * h;
        
        assertTrue(img.nPixels == n);
        assertTrue(img.r.length == n);
        assertTrue(img.g.length == n);
        assertTrue(img.b.length == n);
        assertTrue(img.cieX.length == n);
        assertTrue(img.cieY.length == n);
        assertTrue(img.hue.length == n);
        assertTrue(img.saturation.length == n);
        assertTrue(img.brightness.length == n);
        
    }

    public void testSetRadiusForPopulateOnDemand() {

        // this also implicitly tests calculateColorIncludingNeighbors
        
        int w = 10;
        int h = 10;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
        
        int count = 0;
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                count++;
                img.setRGB(col, row, 100 + count, 110 + count, 120 + count);
            }
        }
        
        // assert extra colors are all zeros at first
        for (int idx = 0; idx < n; idx++) {
            assertTrue(img.cieX[idx] == 0.);
            assertTrue(img.cieY[idx] == 0.);
            assertTrue(img.hue[idx] == 0.);
            assertTrue(img.saturation[idx] == 0.);
            assertTrue(img.brightness[idx] == 0.);
        }
        
        img.setRadiusForPopulateOnDemand(5);
        img.getCIEX(w/2, h/2);
        for (int idx = 0; idx < n; idx++) {
            assertTrue(img.cieX[idx] > 0.);
            assertTrue(img.cieY[idx] > 0.);
            assertTrue(img.hue[idx] > 0.);
            assertTrue(img.saturation[idx] > 0.);
            assertTrue(img.brightness[idx] > 0.);
        }
        
        // ------ same test but w/ radius = 0 ------
        img = new ImageExt(w, h);
        
        count = 0;
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                count++;
                img.setRGB(col, row, 100 + count, 110 + count, 120 + count);
            }
        }
        
        img.setRadiusForPopulateOnDemand(0);
        int x = w/2;
        int y = h/2;
        img.getCIEX(x, y);
        
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                
                int idx = img.getInternalIndex(col, row);
                
                if ((col == x) && (row == y)) {
                    assertTrue(img.cieX[idx] > 0.);
                    assertTrue(img.cieY[idx] > 0.);
                    assertTrue(img.hue[idx] > 0.);
                    assertTrue(img.saturation[idx] > 0.);
                    assertTrue(img.brightness[idx] > 0.);
                } else {
                    assertTrue(img.cieX[idx] == 0.);
                    assertTrue(img.cieY[idx] == 0.);
                    assertTrue(img.hue[idx] == 0.);
                    assertTrue(img.saturation[idx] == 0.);
                    assertTrue(img.brightness[idx] == 0.);
                }
            }
        }
    }
    
    public void testIndexes() throws Exception {
                
        int w = 10;
        int h = 10;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
        img.setRadiusForPopulateOnDemand(n);
        
        int nTests = 100;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
                
        for (int i = 0; i < nTests; i++) {
            
            int x = sr.nextInt(10);
            int y = sr.nextInt(10);
            
            int idx = img.getInternalIndex(x, y);
            
            int col = img.getCol(idx);
            int row = img.getRow(idx);
            
            assertTrue(x == col);
            assertTrue(y == row);
        }
    }
    
    public void testGetters() throws Exception {
                
        int w = 10;
        int h = 10;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
        img.setRadiusForPopulateOnDemand(n);
        
        int nTests = 10;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < nTests; i++) {
            
            int x = sr.nextInt(10);
            int y = sr.nextInt(10);
            
            int r = sr.nextInt(256);
            int g = sr.nextInt(256);
            int b = sr.nextInt(256);
            img.setRGB(x, y, r, g, b);
            
            assertTrue(img.getR(x, y) == r);
            assertTrue(img.getG(x, y) == g);
            assertTrue(img.getB(x, y) == b);
        }
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            int r = img.getR(x, y);
            int g = img.getG(x, y);
            int b = img.getB(x, y);
            
            float[] cieXY = cieC.rgbToXYChromaticity(r, g, b);
            float[] hsb = new float[3];
            Color.RGBtoHSB(r, g, b, hsb);

            assertTrue(Math.abs(img.getCIEX(x, y) - cieXY[0]) < 0.01);
            assertTrue(Math.abs(img.getCIEX(idx) - cieXY[0]) < 0.01);
            
            assertTrue(Math.abs(img.getCIEY(x, y) - cieXY[1]) < 0.01);
            assertTrue(Math.abs(img.getCIEY(idx) - cieXY[1]) < 0.01);
            
            assertTrue(Math.abs(img.getHue(x, y) - hsb[0]) < 0.01);
            assertTrue(Math.abs(img.getHue(idx) - hsb[0]) < 0.01);
            
            assertTrue(Math.abs(img.getSaturation(x, y) - hsb[1]) < 0.01);
            assertTrue(Math.abs(img.getSaturation(idx) - hsb[1]) < 0.01);
            
            assertTrue(Math.abs(img.getBrightness(x, y) - hsb[2]) < 0.01);
            assertTrue(Math.abs(img.getBrightness(idx) - hsb[2]) < 0.01);
            
        }
        
    }

    public void testCalculateColorIncludingNeighbors() throws Exception {
        
        int w = 10;
        int h = 10;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int idx = 0; idx < n; idx++) {
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            int r = idx*1;
            int g = idx*2;
            int b = idx*3;
            
            img.setRGB(x, y, r, g, b);
        }
        
        int neighborRadius = 1;
        
        for (int idx = 0; idx < n; idx++) {
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            img.calculateColorIncludingNeighbors(idx, neighborRadius);
            
            for (int col = (x - neighborRadius); col <= (x + neighborRadius); 
                col++) {
            
                if ((col < 0) || (col > (w - 1))) {
                    continue;
                }

                for (int row = (y - neighborRadius); row <= 
                    (y + neighborRadius); row++) {

                    if ((row < 0) || (row > (h - 1))) {
                        continue;
                    }
                
                    int index = img.getInternalIndex(col, row);
                    
                    int expectedR = index*1;
                    int expectedG = index*2;
                    int expectedB = index*3;
                    
                    float[] expectedCIEXY = cieC.rgbToXYChromaticity(
                        expectedR, expectedG, expectedB);
                    
                    float[] expectedHSB = new float[3];
                    Color.RGBtoHSB(expectedR, expectedG, expectedB, 
                        expectedHSB);
                        
                    assertTrue(Math.abs(img.cieX[index] - expectedCIEXY[0]) < 0.01);
                    assertTrue(Math.abs(img.cieY[index] - expectedCIEXY[1]) < 0.01);
            
                    assertTrue(Math.abs(img.hue[index] - expectedHSB[0]) < 0.01);
                    assertTrue(Math.abs(img.saturation[index] - expectedHSB[1]) < 0.01);
                    assertTrue(Math.abs(img.brightness[index] - expectedHSB[2]) < 0.01);
                }
            }
        }
    }
    
    private ImageExt getImageExt0() {
        int w = 10;
        int h = 10;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
                        
        for (int idx = 0; idx < n; idx++) {
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            int r = idx*1;
            int g = idx*2;
            int b = idx*3;
            
            img.setRGB(x, y, r, g, b);
        }
        for (int idx = 0; idx < n; idx++) {
            img.calculateColor(idx);
        }
        
        return img;
    }

    public void testCopyImage() {
        
        ImageExt img = getImageExt0();
        
        Image img2 = img.copyImage();
        
        assertTrue(img2 instanceof ImageExt);
        
        ImageExt image2 = (ImageExt)img2;
        
        assertTrue(image2.getNPixels() == img.getNPixels());
        assertTrue(image2.getWidth() == img.getWidth());
        assertTrue(image2.getHeight() == img.getHeight());
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            
            assertTrue(image2.r[idx] == img.r[idx]);
            assertTrue(image2.g[idx] == img.g[idx]);
            assertTrue(image2.b[idx] == img.b[idx]);
            assertTrue(Math.abs(image2.cieX[idx] - img.cieX[idx]) < 0.01);
            assertTrue(Math.abs(image2.cieY[idx] - img.cieY[idx]) < 0.01);
            assertTrue(Math.abs(image2.hue[idx] - img.hue[idx]) < 0.01);
            assertTrue(Math.abs(image2.saturation[idx] - img.saturation[idx]) < 0.01);
            assertTrue(Math.abs(image2.brightness[idx] - img.brightness[idx]) < 0.01);
        }
    }

    public void testResetTo() {
        
        ImageExt img = getImageExt0();
        
        ImageExt image2 = new ImageExt(img.getWidth(), img.getHeight());
        
        image2.resetTo(img);
        
        assertTrue(image2.getNPixels() == img.getNPixels());
        assertTrue(image2.getWidth() == img.getWidth());
        assertTrue(image2.getHeight() == img.getHeight());
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            
            assertTrue(image2.r[idx] == img.r[idx]);
            assertTrue(image2.g[idx] == img.g[idx]);
            assertTrue(image2.b[idx] == img.b[idx]);
            assertTrue(Math.abs(image2.cieX[idx] - img.cieX[idx]) < 0.01);
            assertTrue(Math.abs(image2.cieY[idx] - img.cieY[idx]) < 0.01);
            assertTrue(Math.abs(image2.hue[idx] - img.hue[idx]) < 0.01);
            assertTrue(Math.abs(image2.saturation[idx] - img.saturation[idx]) < 0.01);
            assertTrue(Math.abs(image2.brightness[idx] - img.brightness[idx]) < 0.01);
        }
        
    }
    
    
}
