package algorithms.imageProcessing;

import algorithms.matrix.MatrixUtil;
import java.awt.Color;
import java.security.SecureRandom;
import junit.framework.TestCase;


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
        
        assertTrue(img.getNPixels() == n);
        
    }

    public void testSetRadiusForPopulateOnDemand() {

        // this also implicitly tests calculateColorIncludingNeighbors
        
        int w = 5;
        int h = 5;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
        
        int count = 0;
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                count++;
                img.setRGB(col, row, 100 + count, 110 + count, 120 + count);
            }
        }
        
        img.setRadiusForPopulateOnDemand(5);
        img.getCIEX(w/2, h/2);
        for (int idx = 0; idx < n; idx++) {
            assertTrue(img.getCIEX(idx) > 0.);
            assertTrue(img.getCIEY(idx) > 0.);
            assertTrue(img.getHue(idx) > 0.);
            assertTrue(img.getSaturation(idx) > 0.);
            assertTrue(img.getBrightness(idx) > 0.);
            assertTrue(img.getLuma(idx) > 0.);
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
                
                assertTrue(img.getCIEX(idx) > 0.);
                assertTrue(img.getCIEY(idx) > 0.);
                assertTrue(img.getHue(idx) > 0.);
                assertTrue(img.getSaturation(idx) > 0.);
                assertTrue(img.getBrightness(idx) > 0.);
                assertTrue(img.getLuma(idx) > 0.);
            }
        }
    }
    
    public void testIndexes() throws Exception {
                
        int w = 5;
        int h = 5;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
        img.setRadiusForPopulateOnDemand(n);
        
        int nTests = 100;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
                
        for (int i = 0; i < nTests; i++) {
            
            int x = sr.nextInt(w);
            int y = sr.nextInt(h);
            
            int idx = img.getInternalIndex(x, y);
            
            int col = img.getCol(idx);
            int row = img.getRow(idx);
            
            assertTrue(x == col);
            assertTrue(y == row);
        }
    }
    
    public void testGetters() throws Exception {
                
        int w = 5;
        int h = 5;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
        img.setRadiusForPopulateOnDemand(n);
                
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < 10; i++) {
            
            int idx = sr.nextInt(n);
            
            int r = sr.nextInt(255);
            int g = sr.nextInt(255);
            int b = sr.nextInt(255);
            img.setRGB(idx, r, g, b);
            
            assertTrue(img.getR(idx) == r);
            assertTrue(img.getG(idx) == g);
            assertTrue(img.getB(idx) == b);
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
            
            double[][] rgbToLumaMatrix = new double[3][];
            rgbToLumaMatrix[0] = new double[]{0.256, 0.504, 0.098};
            rgbToLumaMatrix[1] = new double[]{-0.148, -0.291, 0.439};
            rgbToLumaMatrix[2] = new double[]{0.439, -0.368, -0.072};
            double[] yuv = MatrixUtil.multiply(rgbToLumaMatrix, 
                new double[]{r, g, b});

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
            
            assertTrue(Math.abs(img.getLuma(x, y) - yuv[0]) < 0.01);
            assertTrue(Math.abs(img.getLuma(idx) - yuv[0]) < 0.01);
        }
        
    }

    public void testCalculateColorIncludingNeighbors() throws Exception {
        
        int w = 5;
        int h = 5;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int idx = 0; idx < n; idx++) {
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            int r = idx*1;
            int g = idx*2;
            int b = idx*1;
            
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
                    int expectedB = index*1;
                    
                    float[] expectedCIEXY = cieC.rgbToXYChromaticity(
                        expectedR, expectedG, expectedB);
                    
                    float[] expectedHSB = new float[3];
                    Color.RGBtoHSB(expectedR, expectedG, expectedB, 
                        expectedHSB);
                   
                    double[][] rgbToLumaMatrix = new double[3][];
                    rgbToLumaMatrix[0] = new double[]{0.256, 0.504, 0.098};
                    rgbToLumaMatrix[1] = new double[]{-0.148, -0.291, 0.439};
                    rgbToLumaMatrix[2] = new double[]{0.439, -0.368, -0.072};
                    double[] expectedYUV = MatrixUtil.multiply(rgbToLumaMatrix,
                        new double[]{expectedR, expectedG, expectedB});

                    assertTrue(Math.abs(img.getCIEX(index) - expectedCIEXY[0]) < 0.01);
                    assertTrue(Math.abs(img.getCIEY(index) - expectedCIEXY[1]) < 0.01);
            
                    assertTrue(Math.abs(img.getHue(index) - expectedHSB[0]) < 0.01);
                    assertTrue(Math.abs(img.getSaturation(index) - expectedHSB[1]) < 0.01);
                    assertTrue(Math.abs(img.getBrightness(index) - expectedHSB[2]) < 0.01);
                    assertTrue(Math.abs(img.getLuma(index) - expectedYUV[0]) < 0.01);
                }
            }
        }
    }
    
    private ImageExt getImageExt0() {
        int w = 5;
        int h = 5;
        int n = w * h;
        
        ImageExt img = new ImageExt(w, h);
                        
        for (int idx = 0; idx < n; idx++) {
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            int r = idx*1;
            int g = idx*2;
            int b = idx*1;
            
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
            
            assertTrue(image2.getR(idx) == img.getR(idx));
            assertTrue(image2.getG(idx) == img.getG(idx));
            assertTrue(image2.getB(idx) == img.getB(idx));
            assertTrue(Math.abs(image2.getCIEX(idx) - img.getCIEX(idx)) < 0.01);
            assertTrue(Math.abs(image2.getCIEY(idx) - img.getCIEY(idx)) < 0.01);
            assertTrue(Math.abs(image2.getHue(idx) - img.getHue(idx)) < 0.01);
            assertTrue(Math.abs(image2.getSaturation(idx) - img.getSaturation(idx)) < 0.01);
            assertTrue(Math.abs(image2.getBrightness(idx) - img.getBrightness(idx)) < 0.01);
            assertTrue(Math.abs(image2.getLuma(idx) - img.getLuma(idx)) < 0.01);
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
            
            assertTrue(image2.getR(idx) == img.getR(idx));
            assertTrue(image2.getG(idx) == img.getG(idx));
            assertTrue(image2.getB(idx) == img.getB(idx));
            assertTrue(Math.abs(image2.getCIEX(idx) - img.getCIEX(idx)) < 0.01);
            assertTrue(Math.abs(image2.getCIEX(idx) - img.getCIEX(idx)) < 0.01);
            assertTrue(Math.abs(image2.getHue(idx) - img.getHue(idx)) < 0.01);
            assertTrue(Math.abs(image2.getSaturation(idx) - img.getSaturation(idx)) 
                < 0.01);
            assertTrue(Math.abs(image2.getBrightness(idx) - img.getBrightness(idx)) 
                < 0.01);
            assertTrue(Math.abs(image2.getLuma(idx) - img.getLuma(idx)) < 0.01);
        }
        
    }
    
    public void testCopySubImage() {
        
        int w = 10;
        int h = 10;
        
        ImageExt img = new ImageExt(w, h);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                img.setRGB(pixIdx, i, i, j);
                img.calculateColorIncludingNeighbors(pixIdx, 0);
            }
        }
        
        double tol = 0.01;
        
        ImageExt img2 = (ImageExt)img.copySubImage(5, w, 5, h);
        for (int i = 5; i < w; ++i) {
            for (int j = 5; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                int pixIdx2 = img2.getInternalIndex(i - 5, j - 5);
                int r = img.getR(pixIdx);
                int g = img.getG(pixIdx);
                int b = img.getB(pixIdx);
                int r2 = img2.getR(pixIdx2);
                int g2 = img2.getG(pixIdx2);
                int b2 = img2.getB(pixIdx2);
                assertEquals(r, r2);
                assertEquals(g, g2);
                assertEquals(b, b2);
                assertTrue(Math.abs(img.cieX[pixIdx] -
                    img2.cieX[pixIdx2]) < tol);
                assertTrue(Math.abs(img.cieY[pixIdx] -
                    img2.cieY[pixIdx2]) < tol);
                assertTrue(Math.abs(img.hue[pixIdx] -
                    img2.hue[pixIdx2]) < tol);
                assertTrue(Math.abs(img.saturation[pixIdx] -
                    img2.saturation[pixIdx2]) < tol);
                assertTrue(Math.abs(img.brightness[pixIdx] -
                    img2.brightness[pixIdx2]) < tol);
                assertTrue(Math.abs(img.luma[pixIdx] -
                    img2.luma[pixIdx2]) < tol);
            }
        }
        
    }
    
}
