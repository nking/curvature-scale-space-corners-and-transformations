package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HCPTTest extends TestCase {
    
    public HCPTTest() {
    }

    public void test0() throws IOException {
        
        int nBins = 12;
        int maxDimension = 256;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        String filePath0 = ResourceFinder.findFileInTestResources(
            "android_statues_04.jpg");
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();
        int binFactor0 = (int) Math.ceil(Math.max(
            (float) w0 / maxDimension,
            (float) h0 / maxDimension));
        img0 = imageProcessor.binImage(img0, binFactor0);
        GreyscaleImage pt0 = imageProcessor.createCIELUVTheta(img0, 255);
        
        //MiscDebug.writeImage(pt0, "_lab0_");
        
        String filePath1 = ResourceFinder.findFileInTestResources(
            "android_statues_02.jpg");
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));
        img1 = imageProcessor.binImage(img1, binFactor1);
        GreyscaleImage pt1 = imageProcessor.createCIELUVTheta(img1, 255);
        
        //MiscDebug.writeImage(pt1, "_lab1_");
        
        HCPT hgs0 = new HCPT(pt0, 1, nBins);
        
        HCPT hgs1 = new HCPT(pt1, 1, nBins);
         
        int[] feature0 = new int[nBins];
        hgs0.extractFeature(93, 110, feature0);
        
        int[] feature1 = new int[nBins];
        hgs1.extractFeature(57, 65, feature1);
        
        System.out.println("0=" + Arrays.toString(feature0));
        System.out.println("1=" + Arrays.toString(feature1));
                
        float intersection01 = hgs0.intersection(feature0, feature1);
       
        System.out.println("intersection 0:1=" + intersection01);
                
        // 0 and 1 should be similar
        assertTrue(intersection01 > 0.5);
   
        
        feature1 = new int[nBins];
        hgs1.extractFeature(177, 59, feature1);
        intersection01 = hgs0.intersection(feature0, feature1);
        System.out.println("intersection 0:1=" + intersection01);
        assertTrue(intersection01 < 0.3);
    }
    
}
