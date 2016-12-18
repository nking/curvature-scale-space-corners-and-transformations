package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MedianTransformTest extends TestCase {

    public MedianTransformTest(String testName) {
        super(testName);
    }

    public void testMedianTransform0() throws Exception {
            
        String filePath = ResourceFinder.findFileInTestResources(
            "house.gif");
        //String filePath = ResourceFinder.findFileInTestResources(
        //    "android_statues_03.jpg");
        //String filePath = ResourceFinder.findFileInTestResources(
        //    "new-mexico-sunrise_w725_h490.jpg");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGreyscaleFullRange(filePath);
        
        List<GreyscaleImage> transformed = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs = new ArrayList<GreyscaleImage>();
        
        MedianTransform mt = new MedianTransform();
        
        mt.multiscaleMedianTransform(img, transformed, coeffs);
        
        GreyscaleImage r = mt.reconstructMultiscaleMedianTransform(
            transformed.get(transformed.size() - 1), coeffs);
        
        ImageDisplayer.displayImage("reconstructed ", r);
        
        // --------
        
        List<GreyscaleImage> transformed1 = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs1 = new ArrayList<GreyscaleImage>();
        img = ImageIOHelper.readImageAsGreyscaleFullRange(filePath);
        
        mt.multiscalePyramidalMedianTransform2(img, transformed1, coeffs1);
        
        GreyscaleImage r1 = mt.reconstructPyramidalMultiscaleMedianTransform(
            transformed1.get(transformed1.size() - 1), coeffs1);
        
        ImageDisplayer.displayImage("inexact pyramidal med trans reconstructed ", r1);
        
        // --------
        
        List<GreyscaleImage> transformed2 = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs2 = new ArrayList<GreyscaleImage>();
        img = ImageIOHelper.readImageAsGreyscaleFullRange(filePath);
        
        mt.multiscalePyramidalMedianTransform(img, transformed2, coeffs2);
        
        GreyscaleImage r2 = mt.reconstructPyramidalMultiscaleMedianTransform(
            transformed2.get(transformed2.size() - 1), coeffs2);
        
        ImageDisplayer.displayImage("exact pyramidal med trans reconstructed ", r2);
        
    }
    
}
