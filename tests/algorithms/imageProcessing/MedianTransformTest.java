package algorithms.imageProcessing;

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
            
        String filePath = ResourceFinder.findFileInTestResources("house.gif");
        //String filePath = ResourceFinder.findFileInTestResources(
        //    "new-mexico-sunrise_w725_h490.jpg");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGreyscaleFullRange(filePath);
        
        List<GreyscaleImage> transformed = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs = new ArrayList<GreyscaleImage>();
        
        MedianTransform mt = new MedianTransform();
        
        mt.multiscaleMedianTransform(img, transformed, coeffs);
        
        ImageDisplayer.displayImage("transformed ", 
            transformed.get(transformed.size() - 1));
        
        for (int i = 0; i < coeffs.size(); ++i) {
            ImageDisplayer.displayImage("median transform " + i, coeffs.get(i));
        }
        
        GreyscaleImage r = mt.reconstructMultiscaleMedianTransform(
            transformed.get(transformed.size() - 1), coeffs);
        
        ImageDisplayer.displayImage("reconstructed ", r);
        
        List<GreyscaleImage> transformed2 = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs2 = new ArrayList<GreyscaleImage>();
        img = ImageIOHelper.readImageAsGreyscaleFullRange(filePath);
        
        mt.multiscalePyramidalMedianTransform(img, transformed2, coeffs2);
        
        ImageDisplayer.displayImage("pyramidal transformed ", 
            transformed2.get(transformed2.size() - 1));
        
        for (int i = 0; i < coeffs2.size(); ++i) {
            ImageDisplayer.displayImage("pyramidal median transform " + i, 
                coeffs2.get(i));
        }
        
        GreyscaleImage r2 = mt.reconstructPyramidalMultiscaleMedianTransform(
            transformed2.get(transformed2.size() - 1), coeffs2);
        
        ImageDisplayer.displayImage("reconstructed2 ", r2);
        
        int z = 1;
    }
}
