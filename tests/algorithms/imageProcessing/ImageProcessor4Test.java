package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessor4Test extends TestCase {

    public ImageProcessor4Test(String testName) {
        super(testName);
    }
    
    public void test0() throws Exception {
        
        int maxDimension = 256;
        
        String fileName1 = "android_statues_04.jpg";
        
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img = ImageIOHelper.readImageExt(filePath1);

        ImageProcessor imageProcessor = new ImageProcessor();

        int w1 = img.getWidth();
        int h1 = img.getHeight();

        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));

        img = imageProcessor.binImage(img, binFactor1);

        imageProcessor = new ImageProcessor();
        
        
    }
    
}
