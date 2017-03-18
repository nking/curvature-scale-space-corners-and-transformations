package algorithms.imageProcessing;

import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MedianSmooth;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessorTest extends TestCase {

    boolean displayImages = false;
    
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
        assertTrue(result.getWidth() == input.getWidth());
        assertTrue(result.getHeight() == input.getHeight());
    }
    
    public void testCopyImage() throws Exception {
                    
        String filePath = ResourceFinder.findFileInTestResources("tajmahal.png");
        
        Image img = ImageIOHelper.readImage(filePath);
                        
        Image img2 = img.copyImage();
        
        assertNotNull(img2);
                
        assertTrue(Arrays.equals(img2.r, img.r));
        assertTrue(Arrays.equals(img2.g, img.g));
        assertTrue(Arrays.equals(img2.b, img.b));
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
        
        //TODO: revisit this test
        
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
    
    public void testComputeTheta_GreyscaleImage_GreyscaleImage() {
        
        int w = 10;
        int h = 10;
        GreyscaleImage imageX = new GreyscaleImage(w, h, 
            GreyscaleImage.Type.Bits64Signed);
        GreyscaleImage imageY = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32Signed);
        
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
        
        //TODO: revisit this test
        //assertTrue(result.getValue(6, 6) == -45);
        
        //assertTrue(result.getValue(7, 7) == 45);
        
        //assertTrue(result.getValue(8, 8) == -45);
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
    
    public void testFFT2D() throws Exception {
        
        /*
        MedianTransform mt = new MedianTransform();
        
        List<GreyscaleImage> transformed = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs = new ArrayList<GreyscaleImage>();
        mt.multiscaleMedianTransform(img30, transformed, coeffs);        
        GreyscaleImage r = mt.reconstructMultiscaleMedianTransform(
            transformed.get(transformed.size() - 1), coeffs);
        //mt.multiscalePyramidalMedianTransform(getGrid(0), transformed, coeffs);        
        //GreyscaleImage r = mt.reconstructPyramidalMultiscaleMedianTransform(
        //    transformed.get(transformed.size() - 1), coeffs);
        
        if (displayImages) {
            for (int i = 0; i < transformed.size(); ++i) {
                ImageDisplayer.displayImage("transformed " + i, transformed.get(i));
            }
            for (int i = 0; i < coeffs.size(); ++i) {
                ImageDisplayer.displayImage("coeffs " + i, coeffs.get(i));
            }
            ImageDisplayer.displayImage("reconstructed ", r);
        }
        
        int z = 1;
        */
        /*
        
        GreyscaleImage checkerboard = getCheckboard(8);
        int result;
        
        if (displayImages)
        ImageDisplayer.displayImage("checkerboard", checkerboard);
        
        result = imageProcessor.FFT2D(checkerboard, 1);
        
        assertTrue(result == 1);
        
        if (displayImages)
        ImageDisplayer.displayImage("FFT of checkerboard", checkerboard);
        
        result = imageProcessor.FFT2D(checkerboard, -1);
        
        assertTrue(result == 1);
        
        if (displayImages)
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
        
        
    }
   
    public void testBiLinearInterpolation0() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Image img = new Image(10, 10);
        for (int i = 0; i < img.getNPixels(); ++i) {
            img.setRGB(i, 100, 200, 50);
        }
        
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
        
        // ------
        String fileName = "lab.gif";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        int idx = fileName.lastIndexOf(".");
        String fileNameRoot = fileName.substring(0, idx);

        GreyscaleImage gsImg = ImageIOHelper.readImage(filePath).copyToGreyscale();
                
        B3SplineFunction spline = new B3SplineFunction();
        gsImg = spline.calculate(gsImg);
        
        MiscDebug.writeImage(gsImg, "_spline_interp_" + fileNameRoot);
    }
    
    public void testBiLinearInterpolation1() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        for (int i = 0; i < img.getNPixels(); ++i) {
            img.setValue(i, 100);
        }
        
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
    
    private GreyscaleImage getGrid(double rotationInDegrees) {
        
        int width = 300;
        int dim = 50;
        
        GreyscaleImage img = new GreyscaleImage(width, width, GreyscaleImage.Type.Bits32FullRangeInt);
        
        Transformer transformer = new Transformer();
        
        /*
           50   100   150   200   250   |
        */
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                int x = (i + 1)*dim;
                int y = (j + 1)*dim;
                float[] xyT = transformer.rotate((float)rotationInDegrees, x, y);
                if (xyT[0] < 0 || (xyT[0] > (width - 1)) || (xyT[1] < 0) ||
                    (xyT[1] > (width - 1))) {
                    continue;
                }
                img.setValue(Math.round(xyT[0]), Math.round(xyT[1]), 250);
            }
        }
        
        return img;
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
        
        String[] fileNames = new String[]{"books_illum3_v0_695x555.png",
            "books_illum3_v6_695x555.png"};
        
        for (String fileName : fileNames) {
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            Image img = ImageIOHelper.readImageAsGrayScale(filePath);
            GreyscaleImage img0 = img.copyToGreyscale();

            ImageProcessor imageProcessor = new ImageProcessor();
            //imageProcessor.applyAdaptiveMeanThresholding(img0);
            imageProcessor.applyAdaptiveMeanThresholding(img0, 11);

            String bin = ResourceFinder.findDirectory("bin");
            if (fileName.contains("v6")) {
                ImageIOHelper.writeOutputImage(bin + "/books_thresh_v6.png", 
                    img0);
            } else {
                ImageIOHelper.writeOutputImage(bin + "/books_thresh_v0.png", 
                    img0);
            }
        }
    }
    
    public void testBlur0() throws Exception {
        
        // color circles rotated by 90, blurred compared to not rotated and blurred
        String cPath = ResourceFinder.findFileInTestResources("two_circles_color2.png");
        Image t = ImageIOHelper.readImage(cPath);
        
        // trimming to 99 by 99
        Image cTrImg = new Image(99, 99);
        for (int i = 0; i < 99; ++i) {
            for (int j = 0; j < 99; ++j) {
                cTrImg.setRGB(i, j, t.getRGB(i, j));
            }
        }
        int m = 99/2;
        
        Image cImg = cTrImg.copyImage();
        
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(m);
        params.setOriginY(m);
        params.setRotationInDegrees(90);
        params.setScale(1.0f);
        params.setTranslationX(0);
        params.setTranslationY(0);
        
        Transformer transformer = new Transformer();
        cTrImg = transformer.applyTransformation(cTrImg, params, 
            cTrImg.getHeight(), cTrImg.getWidth());        
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(cTrImg, 1.0f);
        
        imageProcessor.blur(cImg, 1.0f);
        
        //MiscDebug.writeImage(cImg, "blur_0");
        //MiscDebug.writeImage(cTrImg, "blur_0_tr");
        
        TransformationParameters params2 = new TransformationParameters();
        params2.setOriginX(m);
        params2.setOriginY(m);
        params2.setRotationInDegrees(-90);
        params2.setScale(1.0f);
        params2.setTranslationX(0);
        params2.setTranslationY(0);
        
        cTrImg = transformer.applyTransformation(cTrImg, params2, 
            cTrImg.getWidth(), cTrImg.getHeight());
        
        //MiscDebug.writeImage(cTrImg, "blur_0_tr_revtr");
        
        int w = cImg.getWidth();
        int h = cImg.getHeight();
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int r0 = cImg.getR(i, j);
                int g0 = cImg.getG(i, j);
                int b0 = cImg.getB(i, j);
                
                int r1 = cTrImg.getR(i, j);
                int g1 = cTrImg.getG(i, j);
                int b1 = cTrImg.getB(i, j);
                
                int diffR = Math.abs(r0 - r1);
                int diffG = Math.abs(g0 - g1);
                int diffB = Math.abs(b0 - b1);
              
                assertTrue(diffR <= 1);
                assertTrue(diffG <= 1);
                assertTrue(diffB <= 1);
            }
        }
    }
    
    public void testBlur1() throws Exception {
        
        // color circles rotated by 90, blurred compared to not rotated and blurred
        String cPath = ResourceFinder.findFileInTestResources("closed_curve_square.png");
        Image t = ImageIOHelper.readImage(cPath);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(t, 1.0f);
                
        MiscDebug.writeImage(t, "blur_1");        
        
    }
    
    public void testMedianFilter() throws Exception {
            
        String filePath = ResourceFinder.findFileInTestResources("tajmahal.png");
        
        Image img = ImageIOHelper.readImage(filePath);
                                        
        MedianSmooth med = new MedianSmooth();
        
        GreyscaleImage out = med.calculate(img.copyToGreyscale(), 10, 10);
        
        if (displayImages)
        ImageDisplayer.displayImage("median filtered", out);                
    }
    
    
    public void testIFFTShift() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[][] a = new double[4][8];
        a[0] = new double[]{-0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375};
        a[1] = new double[]{-0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375};
        a[2] = new double[]{-0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375};
        a[3] = new double[]{-0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375};
        
        double[][] aShifted = imageProcessor.ifftShift(a);
        
        double[][] expected = new double[4][8];
        expected[0] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        expected[1] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        expected[2] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        expected[3] = new double[]{0., 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertEquals(expected[i][j], aShifted[i][j]);
            }
        }
        
    }
    
    public void testDilate() {
        int w = 10;
        int h = 7;
        GreyscaleImage img = new GreyscaleImage(w, h);
        img.fill(1);
        
        /*
        6
        5
        4        _  _  _
        3        _  _  _
        2        _  _  _
        1
        0
           0  1  2  3  4  5  6  7  8  9
        */
        for (int i = 2; i <= 4; ++i) {
            for (int j = 2; j <= 4; ++j) {
                img.setValue(i, j, 0);
            }
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage out = imageProcessor.dilate(img);
        for (int i = 2; i <= 4; ++i) {
            for (int j = 2; j <= 4; ++j) {
                int v = out.getValue(i, j);
                if (i == 3 && j == 3) {
                    assertEquals(0, v);
                } else {
                    assertEquals(1, v);
                }
            }
        }
        
        GreyscaleImage out2 = imageProcessor.erode(img);
        for (int i = 2; i <= 4; ++i) {
            for (int j = 2; j <= 4; ++j) {
                int v = out2.getValue(i, j);
                assertEquals(0, v);
            }
        }
    }
   
    public void testCreateTextureTransforms() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "susan-in_plus.png");
        
        Image img = ImageIOHelper.readImage(filePath);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Map<String, GreyscaleImage> tMap = imageProcessor.createTextureTransforms(
            img.copyToGreyscale(), 2);
        
        for (Entry<String, GreyscaleImage> entry : tMap.entrySet()) {
            String label = entry.getKey();
            GreyscaleImage img3 = entry.getValue();
            MiscMath.rescale(img3, 0, 255);
            if (displayImages)
                ImageDisplayer.displayImage(label, img3);
            imageProcessor.applyAdaptiveMeanThresholding(img3, 1);
            if (displayImages)
                ImageDisplayer.displayImage(label, img3);
        }
        
        /*
        L5E5/E5L5, L5S5/S5L5, L5R5/R5L5, E5E5,
            E5S5/S5E5, E5R5/R5E5, S5S5, S5R5/R5S5,
            R5R5
        */
    }
    
    public void testRecursiveBlur() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "closed_curve.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        imageProcessor.applyThinning2(img);
        
        if (displayImages)
            ImageDisplayer.displayImage("_thinned_", img);
        
        float[][] image = new float[img.getWidth()][];
        for (int i = 0; i < img.getWidth(); ++i) {
            image[i] = new float[img.getHeight()];
            for (int j = 0; j < img.getHeight(); ++j) {
                image[i][j] = img.getValue(i, j);
            }
        }
        
        imageProcessor.recursiveBlur(image, SIGMA.ZEROPOINTSEVENONE);
        
        GreyscaleImage out = img.createWithDimensions();
        int w = out.getWidth();
        int h = out.getHeight();
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = Math.round(image[i][j]);
                if (v < 0) {
                    v = 0;
                } else if (v > 255) {
                    v = 255;
                }
                out.setValue(i, j, v);
            }
        }
        
        if (displayImages)
            ImageDisplayer.displayImage("_recursive_blur_", out);
        
    }
    
    public void testCreateZeroCrossingsCurvature() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "closed_curve.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int w = img.getWidth();
        int h = img.getHeight();
                
        float[][] image = new float[h][];
        for (int i = 0; i < h; ++i) {
            image[i] = new float[w];
            for (int j = 0; j < w; ++j) {
                // make it a binary image to see effects of curvature
                if (img.getValue(j, i) > 0) {
                    image[i][j] = 255;
                }
            }
        }
                
        GreyscaleImage out = imageProcessor
            .createZeroCrossingsCurvature(image, 
            SIGMA.getValue(SIGMA.ONE), 2);
        
        if (displayImages)
            ImageDisplayer.displayImage("_zero_crossings_", out);
        
        imageProcessor.applyThinning(out);
        
        if (displayImages)
            ImageDisplayer.displayImage("_zero_crossings_adap_mean_",
                out);            
    }
    
    public void testExploreTextures() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.exploreTextures();
        
    }
    
    public void testBilinearDownSampling() {
        
        int w1 = 10;
        int h1 = 20;
        int w2 = 7;
        int h2 = 14;
        
        GreyscaleImage a = new GreyscaleImage(w1, h1);
        a.fill(10);        
        
        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage aa = imageProcessor.downSample(a, w2, h2, 0, 255);
        float[][] b = imageProcessor.copyToRowMajor(aa);
        
        assertEquals(w2, b[0].length);
        assertEquals(h2, b.length);
        
        for (int i = 0; i < w2; ++i) {
            for (int j = 0; j < h2; ++j) {
                assertEquals(10.0f, b[j][i]);
            }
        }
    }
    
    public void testBilinearDownSampling2() {
        
        int w1 = 10;
        int h1 = 4;
        int w2 = 5;
        int h2 = 4;
        
        GreyscaleImage a = new GreyscaleImage(w1, h1);
        
        for (int i = 0; i < w1; i += 4) {
            for (int j = 0; j < h1; j += 4) {
                a.setValue(i, j, 10);
                if (j + 1 < h1) {
                    a.setValue(i, j + 1, 10);
                }
                if (i + 1 < w1) {
                    a.setValue(i + 1, j, 10);
                }
                if (j + 1 < h1 && i + 1 < w1) {
                    a.setValue(i + 1, j + 1, 10);
                }
            }
        }
        
        /*
        for (int i = 0; i < w1; ++i) {
            StringBuilder sb = new StringBuilder("a col=").append(i);
            for (int j = 0; j < h1; ++j) {
                sb.append(String.format(" [%d]=%d", j, a.getValue(i, j)));
            }
            System.out.println(sb.toString());
        }*/
        
        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage b = imageProcessor.downSample(a, w2, h2, 0, 255);
                
        /*for (int i = 0; i < w2; i += 1) {
            StringBuilder sb = new StringBuilder("CHECK b col=").append(i);
            for (int j = 0; j < h2; j += 1) {
                sb.append(String.format(" [%d]=%d", j, b.getValue(i, j)));
            }
            System.out.println(sb.toString());
        }*/
        assertEquals(10, b.getValue(0, 0));
        assertEquals(10, b.getValue(0, 1));
        assertEquals(10, b.getValue(2, 0));
        assertEquals(10, b.getValue(2, 1));
        assertEquals(10, b.getValue(4, 0));
        assertEquals(10, b.getValue(4, 1));
        assertEquals(0, b.getValue(0, 2));
        assertEquals(0, b.getValue(0, 3));
        assertEquals(0, b.getValue(1, 2));
        assertEquals(0, b.getValue(1, 3));
        assertEquals(0, b.getValue(2, 2));
        assertEquals(0, b.getValue(2, 3));
        assertEquals(0, b.getValue(3, 2));
        assertEquals(0, b.getValue(3, 3));
        assertEquals(0, b.getValue(4, 2));
        assertEquals(0, b.getValue(4, 3));
    }
    
    public void testBilinearDownSampling3() {
        
        int w1 = 10;
        int h1 = 10;
        int w2 = 5;
        int h2 = 5;
        
        GreyscaleImage a = new GreyscaleImage(w1, h1);
        
        for (int i = 0; i < w1; i += 4) {
            for (int j = 0; j < h1; j += 4) {
                a.setValue(i, j, 10);
                if (j + 1 < h1) {
                    a.setValue(i, j + 1, 10);
                }
                if (i + 1 < w1) {
                    a.setValue(i + 1, j, 10);
                }
                if (j + 1 < h1 && i + 1 < w1) {
                    a.setValue(i + 1, j + 1, 10);
                }
            }
        }
        
        /*
        for (int i = 0; i < w1; ++i) {
            StringBuilder sb = new StringBuilder("a col=").append(i);
            for (int j = 0; j < h1; ++j) {
                sb.append(String.format(" [%d]=%d", j, a.getValue(i, j)));
            }
            System.out.println(sb.toString());
        }*/
     
        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage b = imageProcessor.downSample(a, w2, h2, 0, 255);
        
        assertEquals(w2, b.getWidth());
        assertEquals(h2, b.getHeight());
        
        
        for (int i = 0; i < w2; i += 2) {
            for (int j = 0; j < h2; j += 2) {
                assertEquals(10, b.getValue(i, j));
            }
        }
        
        /*
        for (int i = 0; i < w2; i += 1) {
            StringBuilder sb = new StringBuilder("CHECK b col=").append(i);
            for (int j = 0; j < h2; j += 1) {
                sb.append(String.format(" [%d]=%d", j, b.getValue(i, j)));
            }
            System.out.println(sb.toString());
        }*/
    }
    
    public void testBuildPyramid2() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "susan-in_plus.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(
            filePath);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        List<GreyscaleImage> pyramid = 
            imageProcessor.buildPyramid2(img, 32);
        
        //TODO: add more than a visual inspection here
        /*
        for (int i = 0; i < pyramid.size(); ++i) {
            GreyscaleImage pImg = pyramid.get(i);
            String str = Integer.toString(i);
            if (str.length() < 2) {
                str = "0" + str;
            }
            MiscDebug.writeImage(pImg, "_pyramid_" + str);
        }*/
        
    }
    
}
