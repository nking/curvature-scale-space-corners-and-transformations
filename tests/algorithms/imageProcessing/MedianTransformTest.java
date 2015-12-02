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
        ImageProcessor imageProcessor = new ImageProcessor();
        img = imageProcessor.binImage(img, 2);
        
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
        
        ImageDisplayer.displayImage("reconsructed ", r);
        
        int z = 1;
    }
}
