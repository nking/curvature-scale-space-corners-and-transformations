package algorithms.imageProcessing;

import algorithms.misc.Complex;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageProcessorTest extends TestCase {

    public ImageProcessorTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {

        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testCreateEmptyImg() {
                
        Image input = new Image(10, 10);
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                input.setRGB(i, j, i, i, i);
            }
        }
        
        Image result = new Image(input.getWidth(), input.getWidth());
        
        assertNotNull(result);
        assertTrue(result.getNPixels() == input.getNPixels());
        assertTrue(result.r.length == input.r.length);
        assertTrue(result.g.length == input.g.length);
        assertTrue(result.b.length == input.b.length);
    }
    
    public void testCopyImage() throws Exception {
                    
        String filePath = ResourceFinder.findFileInTestResources("tajmahal.png");
        
        Image img = ImageIOHelper.readImage(filePath);
        
        ImageProcessor ImageProcessor = new ImageProcessor();
                
        Image img2 = img.copyImage();
        
        assertNotNull(img2);
                
        assertTrue(Arrays.equals(img2.r, img.r));
        assertTrue(Arrays.equals(img2.g, img.g));
        assertTrue(Arrays.equals(img2.b, img.b));
    }
   
    @Test   
    public void testApplySobelKernel() throws Exception {
            
        String filePath = ResourceFinder.findFileInTestResources("tajmahal.png");
        
        Image img = ImageIOHelper.readImage(filePath);
                        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        ImageDisplayer.displayImage("", img);
        
        ImageProcessor.applySobelKernel(img);
        
        ImageDisplayer.displayImage("sobel", img);
                
    }
    
    public void testApplySobelKernel2() throws Exception {
                        
        //valve_gaussian.png
        
        String filePath = ResourceFinder.findFileInTestResources("lena.jpg");
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        Image imgX = ImageIOHelper.readImage(filePath);
        ImageDisplayer.displayImage("", imgX);
        SobelX s0 = new SobelX();
        Kernel kernelX = s0.getKernel();
        float norm0 = s0.getNormalizationFactor();
        ImageProcessor.applyKernel(imgX, kernelX, norm0);
        ImageDisplayer.displayImage("Sobel X", imgX);
        
        Image imgY = ImageIOHelper.readImage(filePath);
        SobelY s1 = new SobelY();
        Kernel kernelY = s1.getKernel();
        norm0 = s1.getNormalizationFactor();
        ImageProcessor.applyKernel(imgY, kernelY, norm0);
        ImageDisplayer.displayImage("Sobel Y", imgY);
        
        Image img = ImageIOHelper.readImage(filePath);
        ImageProcessor.applySobelKernel(img);
                
        ImageDisplayer.displayImage("Sobel", img);
        
    }
    
    /*     -45    90    45          y/x
                -  |  +
            0 -----|----- 0
                +  |  -   
            45    90    -45
        
           when X is 0: if Y > 0, theta is 90
           when Y is 0: if X > 0, theta is 0
    */
    public void testCalculateTheta() throws Exception {
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        assertTrue(ImageProcessor.calculateTheta(0, 1) == 90);
        assertTrue(ImageProcessor.calculateTheta(0, -1) == 90);
        assertTrue(ImageProcessor.calculateTheta(1, 0) == 0);
        assertTrue(ImageProcessor.calculateTheta(-1, 0) == 0);
        
        
        assertTrue(ImageProcessor.calculateTheta(1, 1) == 45);
        assertTrue(ImageProcessor.calculateTheta(-1, 1) == -45);
        assertTrue(ImageProcessor.calculateTheta(-1, -1) == 45);
        assertTrue(ImageProcessor.calculateTheta(1, -1) == -45);
    }

     @Test
    public void testCombineConvolvedImages_Image_Image() {
        
        int w = 10;
        int h = 10;
        Image imageX = new Image(w, h);
        Image imageY = new Image(w, h);
        
        for (int col = 0; col < w; col++) {
             for (int row = 0; row < h; row++) {
                 if ((col & 1) == 0) {
                     imageX.setRGB(col, row, 1, 1, 1);
                 } else {
                     imageY.setRGB(col, row, 1, 1, 1);
                 }
             }
        }
        
        ImageProcessor instance = new ImageProcessor();
        Image result = instance.combineConvolvedImages(imageX, imageY);
        
        for (int col = 0; col < w; col++) {
             for (int row = 0; row < h; row++) {
                 int r = result.getR(col, row);
                 int g = result.getG(col, row);
                 int b = result.getB(col, row);
                 assertTrue(r == 1);
                 assertTrue(g == 1);
                 assertTrue(b == 1);
             }
        }
    }

    @Test
    public void testCombineConvolvedImages_GreyscaleImage_GreyscaleImage() {
        
        int w = 10;
        int h = 10;
        GreyscaleImage imageX = new GreyscaleImage(w, h);
        GreyscaleImage imageY = new GreyscaleImage(w, h);
        
        for (int col = 0; col < w; col++) {
             for (int row = 0; row < h; row++) {
                 if ((col & 1) == 0) {
                     imageX.setValue(col, row, 1);
                 } else {
                     imageY.setValue(col, row, 1);
                 }
             }
        }
        
        ImageProcessor instance = new ImageProcessor();
        GreyscaleImage result = instance.combineConvolvedImages(imageX, imageY);
        
        for (int col = 0; col < w; col++) {
             for (int row = 0; row < h; row++) {
                 int v = result.getValue(col, row);
                 assertTrue(v == 1);
             }
        }
    }
    
    @Test
    public void testComputeTheta_Image_Image() {
        
        int w = 10;
        int h = 10;
        Image imageX = new Image(w, h);
        Image imageY = new Image(w, h);
        
        // for 90 degrees
        imageX.setRGB(1, 1, 0, 0, 0);
        imageY.setRGB(1, 1, 1, 1, 1);
        
        // 90
        imageX.setRGB(2, 2, 0, 0, 0);
        imageY.setRGB(2, 2, -1, -1, -1);
        
        //0
        imageX.setRGB(3, 3, 1, 1, 1);
        imageY.setRGB(3, 3, 0, 0, 0);
        
        //0
        imageX.setRGB(4, 4, -1, -1, -1);
        imageY.setRGB(4, 4, 0, 0, 0);
        
        //45
        imageX.setRGB(5, 5, 1, 1, 1);
        imageY.setRGB(5, 5, 1, 1, 1);
        
        //-45
        imageX.setRGB(6, 6, -1, -1, -1);
        imageY.setRGB(6, 6, 1, 1, 1);
        
        //45
        imageX.setRGB(7, 7, -1, -1, -1);
        imageY.setRGB(7, 7, -1, -1, -1);
        
        //-45
        imageX.setRGB(8, 8, 1, 1, 1);
        imageY.setRGB(8, 8, -1, -1, -1);
        
        ImageProcessor instance = new ImageProcessor();
        Image result = instance.computeTheta(imageX, imageY);
                
        // 1,1  90
        // 2,2  90
        // 3,3  0
        // 4,4  0
        // 5,5  45
        // 6,6  -45
        // 7,7  45
        // 8,8  -45

        assertTrue(result.getR(1, 1) == 90);
        assertTrue(result.getG(1, 1) == 90);
        assertTrue(result.getB(1, 1) == 90);
        
        assertTrue(result.getR(2, 2) == 90);
        assertTrue(result.getG(2, 2) == 90);
        assertTrue(result.getB(2, 2) == 90);
        
        assertTrue(result.getR(3, 3) == 0);
        assertTrue(result.getG(3, 3) == 0);
        assertTrue(result.getB(3, 3) == 0);
        
        assertTrue(result.getR(4, 4) == 0);
        assertTrue(result.getG(4, 4) == 0);
        assertTrue(result.getB(4, 4) == 0);
        
        assertTrue(result.getR(5, 5) == 45);
        assertTrue(result.getG(5, 5) == 45);
        assertTrue(result.getB(5, 5) == 45);
        
        assertTrue(result.getR(6, 6) == -45);
        assertTrue(result.getG(6, 6) == -45);
        assertTrue(result.getB(6, 6) == -45);
        
        assertTrue(result.getR(7, 7) == 45);
        assertTrue(result.getG(7, 7) == 45);
        assertTrue(result.getB(7, 7) == 45);
        
        assertTrue(result.getR(8, 8) == -45);
        assertTrue(result.getG(8, 8) == -45);
        assertTrue(result.getB(8, 8) == -45);
    }
    
    @Test
    public void testComputeTheta_GreyscaleImage_GreyscaleImage() {
        
        int w = 10;
        int h = 10;
        GreyscaleImage imageX = new GreyscaleImage(w, h);
        GreyscaleImage imageY = new GreyscaleImage(w, h);
        
        // for 90 degrees
        imageX.setValue(1, 1, 0);
        imageY.setValue(1, 1, 1);
        
        // 90
        imageX.setValue(2, 2, 0);
        imageY.setValue(2, 2, -1);
        
        //0
        imageX.setValue(3, 3, 1);
        imageY.setValue(3, 3, 0);
        
        //0
        imageX.setValue(4, 4, -1);
        imageY.setValue(4, 4, 0);
        
        //45
        imageX.setValue(5, 5, 1);
        imageY.setValue(5, 5, 1);
        
        //-45
        imageX.setValue(6, 6, -1);
        imageY.setValue(6, 6, 1);
        
        //45
        imageX.setValue(7, 7, -1);
        imageY.setValue(7, 7, -1);
        
        //-45
        imageX.setValue(8, 8, 1);
        imageY.setValue(8, 8, -1);
        
        ImageProcessor instance = new ImageProcessor();
        GreyscaleImage result = instance.computeTheta(imageX, imageY);
        
        // 1,1  90
        // 2,2  90
        // 3,3  0
        // 4,4  0
        // 5,5  45
        // 6,6  -45
        // 7,7  45
        // 8,8  -45
        assertTrue(result.getValue(1, 1) == 90);
        
        assertTrue(result.getValue(2, 2) == 90);
        
        assertTrue(result.getValue(3, 3) == 0);
        
        assertTrue(result.getValue(4, 4) == 0);
        
        assertTrue(result.getValue(5, 5) == 45);
        
        assertTrue(result.getValue(6, 6) == -45);
        
        assertTrue(result.getValue(7, 7) == 45);
        
        assertTrue(result.getValue(8, 8) == -45);
    }
    
    @Test
    public void testSubtractImages() {
        
        int w = 10;
        int h = 10;
        GreyscaleImage image4 = new GreyscaleImage(w, h);
        GreyscaleImage image3 = new GreyscaleImage(w, h);
        
        for (int col = 0; col < w; col++) {
             for (int row = 0; row < h; row++) {
                 image4.setValue(col, row, 4);
                 image3.setValue(col, row, 3);
             }
        }
        
        ImageProcessor instance = new ImageProcessor();
        
        GreyscaleImage result = instance.subtractImages(image4, image3);
       
        for (int col = 0; col < w; col++) {
             for (int row = 0; row < h; row++) {
                 int v = result.getValue(col, row);
                 assertTrue(v == 1);
             }
        }
    }
      
    @Test
    public void testShrinkImageToFirstNonZeros() {
        
        int w = 10;
        int h = 10;
        int border = 1;
        GreyscaleImage image = new GreyscaleImage(w, h);
        
        for (int col = 0 + border; col < (w - border); col++) {
             for (int row = 0 + border; row < (h - border); row++) {
                 image.setValue(col, row, 1);
             }
        }
        
        ImageProcessor instance = new ImageProcessor();
        int[] offsetXY = instance.shrinkImageToFirstNonZeros(image);
        
        // added a buffer of 1 around each border
        assertTrue(image.getWidth() == (w - 2*border + 2));
        assertTrue(image.getHeight() == (h - 2*border + 2));
    }

    public void testBinImageToKeepZeros() throws Exception {
        
        int w0 = 4;
        int h0 = 6;
        int xoff = 2;
        int yoff = 4;
        
        int binFactor = 2;
        
        GreyscaleImage img = new GreyscaleImage(w0, h0);
        img.setXRelativeOffset(xoff);
        img.setYRelativeOffset(yoff);
        
        img.fill(4);
        img.setValue(0, 0, 0);
        img.setValue(2, 4, 0);
        /*
        @ 1 2 3 4
        1
        2
        3
        4   @
        5
        */
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        GreyscaleImage out = ImageProcessor.binImageToKeepZeros(img, binFactor);
        
        assertTrue(out.getWidth() == w0/binFactor);
        assertTrue(out.getHeight() == h0/binFactor);
        
        assertTrue(out.getXRelativeOffset() == xoff/binFactor);
        assertTrue(out.getYRelativeOffset() == yoff/binFactor);
        
        for (int col = 0; col < out.getWidth(); col++) {
            for (int row = 0; row < out.getHeight(); row++) {
                int v = out.getValue(col, row);
                if ((col == 0) && (row == 0)) {
                    assertTrue(v == 0);
                } else if ((col == 1) && (row == 2)) {
                    assertTrue(v == 0);
                } else {
                    assertTrue(v == 4);
                }
            }
        }
    }
    
    public void testBinImage() throws Exception {
                
        int w0 = 4;
        int h0 = 6;
        
        int binFactor = 2;
        
        Image img = new Image(w0, h0);
        
        /*
        @ 1 2 3 4
        1
        2
        3
        4   @
        5
        */
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                if ((col == 0) && (row == 0)) {
                    img.setRGB(col, row, 0, 0, 0);
                } else if ((col == 2) && (row == 4)) {
                    img.setRGB(col, row, 0, 0, 0);
                } else {
                    img.setRGB(col, row, 4, 0, 0);
                }
            }
        }
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        Image out = ImageProcessor.binImage(img, binFactor);
        
        assertTrue(out.getWidth() == w0/binFactor);
        assertTrue(out.getHeight() == h0/binFactor);
        
        for (int col = 0; col < out.getWidth(); col++) {
            for (int row = 0; row < out.getHeight(); row++) {
                int r = out.getR(col, row);
                if ((col == 0) && (row == 0)) {
                    assertTrue(r == 3);
                } else if ((col == 1) && (row == 2)) {
                    assertTrue(r == 3);
                } else {
                    assertTrue(r == 4);
                }
                assertTrue(out.getG(col, row) == 0);
                assertTrue(out.getB(col, row) == 0);
            }
        }
    }
    
    public void testBinImageExt() throws Exception {
                
        int w0 = 4;
        int h0 = 6;
        
        int binFactor = 2;
        
        ImageExt img = new ImageExt(w0, h0);
        
        /*
        @ 1 2 3 4
        1
        2
        3
        4   @
        5
        */
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                if ((col == 0) && (row == 0)) {
                    img.setRGB(col, row, 0, 0, 0);
                } else if ((col == 2) && (row == 4)) {
                    img.setRGB(col, row, 0, 0, 0);
                } else {
                    img.setRGB(col, row, 4, 0, 0);
                }
            }
        }
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        ImageExt out = ImageProcessor.binImage(img, binFactor);
        
        assertTrue(out.getWidth() == w0/binFactor);
        assertTrue(out.getHeight() == h0/binFactor);
        
        for (int col = 0; col < out.getWidth(); col++) {
            for (int row = 0; row < out.getHeight(); row++) {
                int r = out.getR(col, row);
                if ((col == 0) && (row == 0)) {
                    assertTrue(r == 3);
                } else if ((col == 1) && (row == 2)) {
                    assertTrue(r == 3);
                } else {
                    assertTrue(r == 4);
                }
                assertTrue(out.getG(col, row) == 0);
                assertTrue(out.getB(col, row) == 0);
            }
        }
    }
    
    public void testUnbinMask() throws Exception {
       
        int w0 = 4;
        int h0 = 6;
        int xOff = 2;
        int yOff = 10;
        
        int binFactor = 2;
        
        int w1 = w0/binFactor;
        int h1 = h0/binFactor;
        
        GreyscaleImage originalTheta = new GreyscaleImage(w0, h0);
        originalTheta.setXRelativeOffset(xOff);
        originalTheta.setYRelativeOffset(yOff);
        
        GreyscaleImage mask = new GreyscaleImage(w1, h1);
        mask.setXRelativeOffset(xOff/binFactor);
        mask.setYRelativeOffset(yOff/binFactor);
        
        /*
        @ @ 
        1  
        2
        */
        for (int col = 0; col < mask.getWidth(); col++) {
            for (int row = 0; row < mask.getHeight(); row++) {
                if (row == 0) {
                    mask.setValue(col, row, 0);
                } else {
                    mask.setValue(col, row, 4);
                }
            }
        }
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        GreyscaleImage out = ImageProcessor.unbinMask(mask, binFactor, originalTheta);
        
        assertTrue(out.getWidth() == originalTheta.getWidth());
        assertTrue(out.getHeight() == originalTheta.getHeight());
        assertTrue(out.getXRelativeOffset() == originalTheta.getXRelativeOffset());
        assertTrue(out.getYRelativeOffset() == originalTheta.getYRelativeOffset());
        
        for (int col = 0; col < out.getWidth(); col++) {
            for (int row = 0; row < out.getHeight(); row++) {
                if ((row == 0) || (row == 1)) {
                    assertTrue(out.getValue(col, row) == 0);
                } else {
                    assertTrue(out.getValue(col, row) == 4);
                }
            }
        }
    }
    
    public void testUnbinMask2() throws Exception {
       
        int w0 = 5;
        int h0 = 7;
        int xOff = 2;
        int yOff = 10;
        
        int binFactor = 2;
        
        int w1 = w0/binFactor;
        int h1 = h0/binFactor;
        
        GreyscaleImage originalTheta = new GreyscaleImage(w0, h0);
        originalTheta.setXRelativeOffset(xOff);
        originalTheta.setYRelativeOffset(yOff);
        
        GreyscaleImage mask = new GreyscaleImage(w1, h1);
        mask.setXRelativeOffset(xOff/binFactor);
        mask.setYRelativeOffset(yOff/binFactor);
        
        /*
        @ @ 
        1  
        2
        */
        for (int col = 0; col < mask.getWidth(); col++) {
            for (int row = 0; row < mask.getHeight(); row++) {
                if (row == 0) {
                    mask.setValue(col, row, 0);
                } else {
                    mask.setValue(col, row, 4);
                }
            }
        }
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        GreyscaleImage out = ImageProcessor.unbinMask(mask, binFactor, originalTheta);
        
        assertTrue(out.getWidth() == originalTheta.getWidth());
        assertTrue(out.getHeight() == originalTheta.getHeight());
        assertTrue(out.getXRelativeOffset() == originalTheta.getXRelativeOffset());
        assertTrue(out.getYRelativeOffset() == originalTheta.getYRelativeOffset());
        
        for (int col = 0; col < out.getWidth(); col++) {
            for (int row = 0; row < out.getHeight(); row++) {
                if ((row == 0) || (row == 1)) {
                    assertTrue(out.getValue(col, row) == 0);
                } else {
                    assertTrue(out.getValue(col, row) == 4);
                }
            }
        }
    }
    
    public void testPrintImageColorContrastStats() throws Exception {
        
        //SKY avg: 161 min=35 max=124
        String filePath = ResourceFinder.findFileInTestResources(
            "venturi_mountain_j6_0001.png");
        
        Image img = ImageIOHelper.readImage(filePath);
                        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        ImageProcessor.printImageColorContrastStats(img, 161, 501);
        
    }
    
    public void testFFT2D() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage checkerboard = getCheckboard(8);
        int result;
        
        /*
        ImageDisplayer.displayImage("checkerboard", checkerboard);
        
        result = imageProcessor.FFT2D(checkerboard, 1);
        
        assertTrue(result == 1);
        
        ImageDisplayer.displayImage("FFT of checkerboard", checkerboard);
        
        result = imageProcessor.FFT2D(checkerboard, -1);
        
        assertTrue(result == 1);
        
        ImageDisplayer.displayImage("FFT^-1 of FFT checkerboard", checkerboard);
        
        GreyscaleImage expected = getCheckboard(8);
        
        for (int col = 0; col < checkerboard.getWidth(); col++) {
            for (int row = 0; row < checkerboard.getHeight(); row++) {
                int v = checkerboard.getValue(col, row);
                int vExpected = expected.getValue(col, row);
                assertTrue(v == vExpected);
            }
        }
        */
        
        String filePath = ResourceFinder.findFileInTestResources("lena.jpg");        
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        int w = img2.getWidth();
        int h = img2.getHeight();
        
        ImageDisplayer.displayImage("lena", img2);
        
        Complex[][] cc = imageProcessor.convertImage(img2);
                
        Complex[][] ccOut = imageProcessor.apply2DFFT(cc, true);
        
        GreyscaleImage img3 = img2.createWithDimensions();
        
        imageProcessor.writeToImage(img3, ccOut);
                           
        ImageDisplayer.displayImage("FFT of lena", img3);
        
        ccOut = imageProcessor.apply2DFFT(ccOut, false);
        
        imageProcessor.writeToImage(img3, ccOut);
        
        ImageDisplayer.displayImage("FFT^-1 of FFT lena", img3);

        int z = 1;
        
    }
    
    public void testDeconvolve() throws Exception {
        
        // NOT READY FOR USE YET
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        String filePath = ResourceFinder.findFileInTestResources("test_for_deconvolve2.png");        
        GreyscaleImage img0 = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        imageProcessor.applyDeconvolution(img0);
        imageProcessor.applyDeconvolution(img0);
        imageProcessor.applyDeconvolution(img0);
                
        String dirPath = ResourceFinder.findDirectory("bin");
        String filePath2 = dirPath + "/tmpWiener.png";
            
        ImageIOHelper.writeOutputImage(filePath2, img0);
        
        int z = 1;
        
    }
    
    private GreyscaleImage getCheckboard(int width) {
        
        GreyscaleImage checkerboard = new GreyscaleImage(width, width);
        
        for (int i = 0; i < checkerboard.getWidth(); i++) {
            if ((i & 1) == 1) {
                // odd
                for (int j = 1; j < checkerboard.getHeight(); j+=2) {
                    checkerboard.setValue(i, j, 250);
                }
            } else {
                // even
                for (int j = 0; j < checkerboard.getHeight(); j+=2) {
                    checkerboard.setValue(i, j, 250);
                }
            }
        }
        
        return checkerboard;
    }
}
