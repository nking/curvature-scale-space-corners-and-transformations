package algorithms.imageProcessing;

import algorithms.misc.Complex;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageProcessorTest extends TestCase {

    public ImageProcessorTest(String testName) {
        super(testName);
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
    
    public void testReadNonZeroPixels() throws Exception {
        
        int w = 10;
        int h = 10;
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        
        for (int i = 0; i < w; i++) {
            img.setValue(i, i, 1);
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(img);
        
        int count = 0;
        for (int i = 0; i < w; i++) {
            assertTrue(points.contains(new PairInt(i, i)));
            count++;
        }
        
        assertTrue(points.size() == count);
        
        GreyscaleImage img2 = img.createWithDimensions();
        
        imageProcessor.writeAsBinaryToImage(img2, points);
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int v = img.getValue(i, j);
                int v2 = img2.getValue(i, j);
                assertTrue(v == v2);
            }
        }
    }
   
    public void testBiLinearInterpolation0() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Image img = new Image(10, 10);
        Arrays.fill(img.r, 100);
        Arrays.fill(img.g, 200);
        Arrays.fill(img.b, 50);
        
        float x = 2f;
        float y = 2f;
        double[] rgb = imageProcessor.biLinearInterpolation(img, x, y);
        assertTrue(Math.abs(rgb[0] - 100) < 0.1);
        assertTrue(Math.abs(rgb[1] - 200) < 0.1);
        assertTrue(Math.abs(rgb[2] - 50) < 0.1);
        
        x = 2.2f;
        y = 2.7f;
        rgb = imageProcessor.biLinearInterpolation(img, x, y);
        assertTrue(Math.abs(rgb[0] - 100) < 0.1);
        assertTrue(Math.abs(rgb[1] - 200) < 0.1);
        assertTrue(Math.abs(rgb[2] - 50) < 0.1);
        
    }
    
    public void testBiLinearInterpolation1() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        Arrays.fill(img.getValues(), 100);
        
        float x = 2f;
        float y = 2f;
        double v = imageProcessor.biLinearInterpolation(img, x, y);
        assertTrue(Math.abs(v - 100) < 0.1);
        
        x = 2.2f;
        y = 2.7f;
        v = imageProcessor.biLinearInterpolation(img, x, y);
        assertTrue(Math.abs(v - 100) < 0.1);
        
        //-------
        img = new GreyscaleImage(10, 10);
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                img.setValue(i, j, i);
            }
        }
        x = 4.3f;
        y = 6f;
        v = imageProcessor.biLinearInterpolation(img, x, y);
        assertTrue(Math.abs(v - 4.3) < 0.01);
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
    
    public void testApplyAdaptiveMean() throws Exception {
        
        int w = 6;
        int h = 6;
        int dimension = 4;
        /*
        [10] [10] [10] [10] [40] [50]  0
        [10] [10] [10] [10] [40] [50]  1
        [10] [10] [10] [10] [40] [50]  2
        [10] [10] [10] [10] [40] [50]  3
        [20] [20] [20] [20] [40] [50]  4
        [30] [30] [30] [30] [40] [50]  5
         0    1    2    3     4    5
        */
        
        GreyscaleImage out = new GreyscaleImage(w, h);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                out.setValue(i, j, 10);
            }
        }
        for (int i = 0; i < 4; ++i) {
            for (int j = 4; j < 5; ++j) {
                out.setValue(i, j, 20);
            }
        }
        for (int i = 0; i < 4; ++i) {
            for (int j = 5; j < 6; ++j) {
                out.setValue(i, j, 30);
            }
        }
        for (int j = 0; j < 6; ++j) {
            out.setValue(4, j, 40);
        }
        for (int j = 0; j < 6; ++j) {
            out.setValue(5, j, 50);
        }
        /*
        [10] [10] [10] [10] [40] [50]  0
        [10] [10] [10] [10] [40] [50]  1
        [10] [10] [10] [10] [40] [50]  2
        [10] [10] [10] [10] [40] [50]  3
        [20] [20] [20] [20] [40] [50]  4
        [30] [30] [30] [30] [40] [50]  5
         0    1    2    3     4    5
        */
        /*
        expected:
        
        [10] [18] [27] [33] [45] [50]  0
        [13] [13] [28] [34] [45] [50]  1
        [18] [18] [31] [36] [45] [50]  2
        [20] [20] [32] [37] [45] [50]  3
        [25] [25] [35] [38] [45] [50]  4
        [30] [35] [37] [40] [45] [50]  5
         0    1    2    3     4    5
        */
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyBoxcarMean(out, dimension);
                
        assertEquals(out.getValue(0, 0), 10);
        assertEquals(20, out.getValue(0, 3));
        assertEquals(25, out.getValue(0, 4));
        assertEquals(30, out.getValue(0, 5));
        
        assertEquals(31, out.getValue(2, 2));
        
        assertEquals(38, out.getValue(3, 4));
        assertEquals(40, out.getValue(3, 5));
    }
    
    public void testApplyCenteredMean() throws Exception {
        
        int w = 6;
        int h = 6;
        int halfDimension = 1;
        /*
        [10] [10] [10] [10] [40] [50]  0
        [10] [10] [10] [10] [40] [50]  1
        [10] [10] [10] [10] [40] [50]  2
        [10] [10] [10] [10] [40] [50]  3
        [20] [20] [20] [20] [40] [50]  4
        [30] [30] [30] [30] [40] [50]  5
         0    1    2    3     4    5
        */
        
        GreyscaleImage out = new GreyscaleImage(w, h);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                out.setValue(i, j, 10);
            }
        }
        for (int i = 0; i < 4; ++i) {
            for (int j = 4; j < 5; ++j) {
                out.setValue(i, j, 20);
            }
        }
        for (int i = 0; i < 4; ++i) {
            for (int j = 5; j < 6; ++j) {
                out.setValue(i, j, 30);
            }
        }
        for (int j = 0; j < 6; ++j) {
            out.setValue(4, j, 40);
        }
        for (int j = 0; j < 6; ++j) {
            out.setValue(5, j, 50);
        }
        /*  
        [10] [10] [10] [10] [40] [50]  0
        [10] [10] [10] [10] [40] [50]  1
        [10] [10] [10] [10] [40] [50]  2
        [10] [10] [10] [10] [40] [50]  3
        [20] [20] [20] [20] [40] [50]  4
        [30] [30] [30] [30] [40] [50]  5
         0    1    2    3     4    5
        
        halfDimension = 1, sum along rows
        [30] [30] [30] [30] [120] [150]  0
        [30] [30] [30] [30] [120] [150]  1
        [30] [30] [30] [30] [120] [150]  2
        [40] [40] [40] [40] [120] [150]  3
        [60] [60] [60] [60] [120] [150]  4
        [75] [75] [75] [75] [120] [150]  5 
         0    1    2    3     4    5
        
        halfDimension = 1, sum along cols
         [90]  [90]  [90] [180] [300] [405]  0
         [90]  [90]  [90] [180] [300] [405]  1
         [90]  [90]  [90] [180] [300] [405]  2
        [120] [120] [120] [200] [310] [405]  3
        [180] [180] [180] [240] [330] [405]  4
        [225] [225] [225] [270] [345] [405]  5
         0      1      2    3     4     5
        
        0:   10   10   10   20   33   45
        1:   10   10   10   20   33   45
        2:   10   10   10   20   33   45
        3:   13   13   13   22   34   45
        4:   20   20   20   27   37   45
        5:   25   25   25   30   38   45
        */
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyCenteredMean(out, halfDimension);
                
        assertEquals(out.getValue(0, 0), 10);
        assertEquals(13, out.getValue(0, 3));
        assertEquals(20, out.getValue(0, 4));
        assertEquals(25, out.getValue(0, 5));
        
        assertEquals(10, out.getValue(2, 2));
        
        assertEquals(27, out.getValue(3, 4));
        assertEquals(30, out.getValue(3, 5));
        
    }
    
    public void testApplyAdaptiveMeanThresholding() throws Exception {
        
        String fileName = "books_illum3_v0_695x555.png";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        Image img = ImageIOHelper.readImageAsGrayScale(filePath);
        GreyscaleImage img0 = img.copyToGreyscale();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyAdaptiveMeanThresholding(img0);
        
        String bin = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(bin + "/books_thresh.png", img0);
    }
    
    public void testApplyColorSegmentation() throws Exception {
        
        //String fileName = "two_circles_color2.png";
        //String fileName = "two_circles_color.png";
        //String fileName = "cloudy_san_jose.jpg";
        //String fileName = "middlebury_cones_im2.png"; // a limitFrac of 0.1 works well
        //String fileName = "brown_lowe_2003_image1.jpg";
        //String fileName = "venturi_mountain_j6_0010.png";
        String fileName = "books_illum3_v6_695x555.png";
        //String fileName = "brown_lowe_2003_image1.jpg";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        //HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(img);
        //hEq.applyFilter();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        List<Set<PairInt>> clusterSets = 
            imageProcessor.calculateColorSegmentation3(img, /*0.1f,*/ true);
            //imageProcessor.calculateColorSegmentation2(img, true);
        
        int nPoints = count(clusterSets);
        
        int nExtraForDot = 0;
        
        Image img2 = new Image(img.getWidth(), img.getHeight());
        
        for (int i = 0; i < clusterSets.size(); ++i) {
            
            int[] rgb = ImageIOHelper.getNextRGB(i);
            
            Set<PairInt> set = clusterSets.get(i);
            
            ImageIOHelper.addToImage(set, 0, 0, img2, nExtraForDot, rgb[0], 
                rgb[1], rgb[2]);
        }
        String bin = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(bin + "/cluster.png", img2);
        
        
        //assertTrue(nPoints == img.getNPixels());
        
        boolean[] present = new boolean[img.getNPixels()];
        for (Set<PairInt> set : clusterSets) {
            for (PairInt p : set) {
                int pixIdx = img.getInternalIndex(p.getX(), p.getY());
                assertFalse(present[pixIdx]);
                present[pixIdx] = true;
            }
        }
        
        int z = 1;
    }

    private int count(List<Set<PairInt>> setList) {
        
        int c = 0;
        for (Set<PairInt> set : setList) {
            c += set.size();
        }
        
        return c;
    }
    
}
