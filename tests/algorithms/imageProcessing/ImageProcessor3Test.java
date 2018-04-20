package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessor3Test extends TestCase {

    public ImageProcessor3Test(String testName) {
        super(testName);
    }
    
    public void test0() throws Exception {
        
        int maxDimension = 256;
        
        String fileName1 = "android_statues_01.jpg";
        
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img = ImageIOHelper.readImageExt(filePath1);

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        int w1 = img.getWidth();
        int h1 = img.getHeight();

        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));

        img = imageProcessor.binImage(img, binFactor1);

        imageProcessor = new ImageProcessor();
        GreyscaleImage theta1 = imageProcessor
            .createCIELUVTheta_WideRangeLightness(img, 255);
        //MiscDebug.writeImage(theta1,  "_theta_"  + fileName1Root);
        
        imageProcessor.singlePixelFilter(theta1);
        imageProcessor.singlePixelFilter(theta1);
        MiscDebug.writeImage(theta1,  "_theta_1_"  + fileName1Root);
        
        GreyscaleImage tmp = theta1.copyImage();
        //unsharp mask to bring up the cupcake and icecream
        imageProcessor.applyUnsharpMask(tmp, 0.5f, 5, 0.f);
        MiscDebug.writeImage(tmp,  "_unsharp__1_"  + fileName1Root);
        
    }
    
}
