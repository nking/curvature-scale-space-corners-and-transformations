package algorithms.imageProcessing;

import java.util.Arrays;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class FeaturesTest extends TestCase {
    
    public FeaturesTest() {
    }

    public void testIsWithinXBounds0() {
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        GreyscaleImage thetaImg = img.createWithDimensions();
        
        Features features = new Features(img, thetaImg, 2, false);
        
        assertTrue(features.isWithinXBounds(2));
        assertFalse(features.isWithinXBounds(-2));
        
        assertTrue(features.isWithinYBounds(2));
        assertFalse(features.isWithinYBounds(-2));
        
    }
    
    public void testIsWithinXBounds1() {
        
        Image img = new Image(10, 10);
        GreyscaleImage thetaImg = new GreyscaleImage(10, 10);
        
        Features features = new Features(img, thetaImg, 2, false);
        
        assertTrue(features.isWithinXBounds(2));
        assertFalse(features.isWithinXBounds(-2));
        
        assertTrue(features.isWithinYBounds(2));
        assertFalse(features.isWithinYBounds(-2));
        
    }

    public void testExtractIntensity() {
        
        boolean normalize = false;
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        Arrays.fill(img.getValues(), 100);
        
        GreyscaleImage thetaImg = img.createWithDimensions();
        
        Features features = new Features(img, thetaImg, 2, normalize);
        
        IntensityDescriptor desc = features.extractIntensity(5, 5, 0);
        assertNotNull(desc);
        assertTrue(desc.isNormalized() == normalize);
        
        float errSq = desc.sumSquaredError();
        assertTrue(errSq == 0);
        assertTrue(desc.calculateSSD(desc) == 0);
        
    }

    public void testCalculateIntensityStat() {
        
    }
    
    public void testExtractGradient() {
        // not implemented yet
    }
    
}
