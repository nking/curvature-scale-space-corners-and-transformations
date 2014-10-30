package algorithms.imageProcessing;

import algorithms.ResourceFinder;
import java.io.File;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageIOHelperTest {
    
    public ImageIOHelperTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testReadImage() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "small_test.png");
        
        Image result = ImageIOHelper.readImage(filePath);
        
        assertTrue(result.getR(0, 0) == 255);
        assertTrue(result.getG(1, 0) == 255);
        assertTrue(result.getB(1, 1) == 255);
        
        assertTrue(result.getR(0, 1) == 255);
        assertTrue(result.getG(0, 1) == 255);
        assertTrue(result.getB(0, 1) == 255);
        
    }
    
    @Test
    public void testReadImage_indexed() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "small_test_indexed.png");
        
        Image result = ImageIOHelper.readImage(filePath);
        
        assertTrue(result.getR(0, 0) == 255);
        assertTrue(result.getG(1, 0) == 255);
        assertTrue(result.getB(1, 1) == 255);
        
        assertTrue(result.getR(0, 1) == 255);
        assertTrue(result.getG(0, 1) == 255);
        assertTrue(result.getB(0, 1) == 255);
        
    }
    
    @Test
    public void testReadImage_grayscale() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "small_test_grayscale.png");
        
        Image result = ImageIOHelper.readImage(filePath);
        
        assertTrue(result.getR(1, 1) > 0);
        assertTrue(result.getG(1, 1) > 0);
        assertTrue(result.getB(1, 1) > 0);
        
        assertTrue((result.getR(0, 0) > 0) && (result.getG(0, 0) > 0) &&
            (result.getB(0, 0) > 0));
        
        assertTrue((result.getR(0, 1) > 0) && (result.getG(0, 1) > 0) &&
            (result.getB(0, 1) > 0));
        
        assertTrue((result.getR(1, 0) > 0) && (result.getG(1, 0) > 0) &&
            (result.getB(1, 0) > 0));
    }

    @Test
    public void testReadImageAsGrayScaleG() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "small_test.png");
        
        GreyscaleImage result = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        // reads the G vector only
        assertTrue(result.getValue(0, 0) > 0);
        assertTrue(result.getValue(1, 0) > 0);
        assertTrue(result.getValue(1, 1) > 0);
        
        assertTrue(result.getValue(0, 1) > 0);
        
    }

    @Test
    public void testReadImageAsGrayScaleAvgRGB() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "small_test.png");
        
        GreyscaleImage result = ImageIOHelper.readImageAsGrayScaleAvgRGB(
            filePath);
        
        // the final vector is averaged from R, G and B
        
        assertTrue(result.getValue(0, 0) == 131);
        assertTrue(result.getValue(1, 0) > 0);
        assertTrue(result.getValue(1, 1) > 0);        
        assertTrue(result.getValue(0, 1) > 0);
    }

    @Test
    public void testWriteOutputImage_String_Image() throws Exception {
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        Image img = new Image(1, 1);
        img.setRGB(0, 0, 255, 0, 0);
        
        String filePath = dirPath + "/temp.png";
        
        ImageIOHelper.writeOutputImage(filePath, img);
        
        File file = new File(filePath);
        
        assertTrue(file.exists());
        
        Image result = ImageIOHelper.readImage(filePath);
        
        assertTrue(result.getR(0, 0) == 255);
        assertTrue(result.getG(0, 0) == 0);
        assertTrue(result.getB(0, 0) == 0);
    }

    @Test
    public void testWriteOutputImage_String_GreyscaleImage() throws Exception {
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        GreyscaleImage img = new GreyscaleImage(1, 1);
        img.setValue(0, 0, 255);
        
        String filePath = dirPath + "/temp.png";
        
        ImageIOHelper.writeOutputImage(filePath, img);
        
        File file = new File(filePath);
        
        assertTrue(file.exists());
        
        GreyscaleImage result = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        assertTrue(result.getValue(0, 0) == 255);
        
    }
    
}
