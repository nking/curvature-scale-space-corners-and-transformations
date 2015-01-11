package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.Arrays;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageProcesserTest extends TestCase {

    public ImageProcesserTest(String testName) {
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
        
        ImageProcesser imageProcesser = new ImageProcesser();
                
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
                        
        ImageProcesser imageProcesser = new ImageProcesser();
        
        ImageDisplayer.displayImage("", img);
        
        imageProcesser.applySobelKernel(img);
        
        ImageDisplayer.displayImage("sobel", img);
                
    }
    
    public void testApplySobelKernel2() throws Exception {
                        
        //valve_gaussian.png
        
        String filePath = ResourceFinder.findFileInTestResources("lena.jpg");
        
        ImageProcesser imageProcesser = new ImageProcesser();
        
        Image imgX = ImageIOHelper.readImage(filePath);
        ImageDisplayer.displayImage("", imgX);
        SobelX s0 = new SobelX();
        Kernel kernelX = s0.getKernel();
        float norm0 = s0.getNormalizationFactor();
        imageProcesser.applyKernel(imgX, kernelX, norm0);
        ImageDisplayer.displayImage("Sobel X", imgX);
        
        Image imgY = ImageIOHelper.readImage(filePath);
        SobelY s1 = new SobelY();
        Kernel kernelY = s1.getKernel();
        norm0 = s1.getNormalizationFactor();
        imageProcesser.applyKernel(imgY, kernelY, norm0);
        ImageDisplayer.displayImage("Sobel Y", imgY);
        
        Image img = ImageIOHelper.readImage(filePath);
        imageProcesser.applySobelKernel(img);
                
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
        
        ImageProcesser imageProcesser = new ImageProcesser();
        
        assertTrue(imageProcesser.calculateTheta(0, 1) == 90);
        assertTrue(imageProcesser.calculateTheta(0, -1) == 90);
        assertTrue(imageProcesser.calculateTheta(1, 0) == 0);
        assertTrue(imageProcesser.calculateTheta(-1, 0) == 0);
        
        
        assertTrue(imageProcesser.calculateTheta(1, 1) == 45);
        assertTrue(imageProcesser.calculateTheta(-1, 1) == -45);
        assertTrue(imageProcesser.calculateTheta(-1, -1) == 45);
        assertTrue(imageProcesser.calculateTheta(1, -1) == -45);
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
        
        ImageProcesser instance = new ImageProcesser();
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
        
        ImageProcesser instance = new ImageProcesser();
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
        
        ImageProcesser instance = new ImageProcesser();
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
        
        ImageProcesser instance = new ImageProcesser();
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
        
        ImageProcesser instance = new ImageProcesser();
        
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
        
        ImageProcesser instance = new ImageProcesser();
        int[] offsetXY = instance.shrinkImageToFirstNonZeros(image);
        
        // added a buffer of 1 around each border
        assertTrue(image.getWidth() == (w - 2*border + 2));
        assertTrue(image.getHeight() == (h - 2*border + 2));
    }

    public void estSegmentationForSky() throws Exception {

        String fileName1 = "brown_lowe_2003_image2.jpg";
        String fileName2 = "brown_lowe_2003_image2_theta.jpg";

        /*
         String fileName1 = "venturi_mountain_j6_0001.png";
         */
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage input = ImageIOHelper.readImageAsGrayScaleB(filePath1);
        GreyscaleImage theta = ImageIOHelper.readImageAsGrayScaleB(
            ResourceFinder.findFileInTestResources(fileName2));
        int image1Width = input.getWidth();
        int image1Height = input.getHeight();

        ImageProcesser imageProcesser = new ImageProcesser();

        imageProcesser.applyImageSegmentationForSky(input, theta);

        ImageDisplayer.displayImage("segmented for sky", input);
        
        int z = 1;
    }
}
