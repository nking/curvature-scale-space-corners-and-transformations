package algorithms.imageProcessing;

import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class FeaturesTest extends TestCase {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
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

    public void testCalculateIntensityStat() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1437335464716L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        for (int col = 0; col < img.getWidth(); ++col) {
            for (int row = 0; row < img.getHeight(); ++row) {                
                int v = sr.nextInt(10);
                if (sr.nextBoolean()) {
                    v = 100 + v;
                } else {
                    v = 100 - v;
                }
                img.setValue(col, row, v);
            }
        }
        
        GreyscaleImage img2 = new GreyscaleImage(10, 10);
        for (int col = 0; col < img.getWidth(); ++col) {
            for (int row = 0; row < img.getHeight(); ++row) {                
                int v = sr.nextInt(12);
                if (sr.nextBoolean()) {
                    v = 140 + v;
                } else {
                    v = 140 - v;
                }
                img2.setValue(col, row, v);
            }
        }
        
        GreyscaleImage thetaImg = img.createWithDimensions();
        
        boolean normalize = false;
        
        Features features = new Features(img, thetaImg, 2, normalize);
                
        int x1 = 5;
        int y1 = 5;
        IntensityDescriptor desc1 = features.extractIntensity(x1, y1, 0);        
        
        Features features2 = new Features(img2, thetaImg, 2, normalize);
        int x2 = 6;
        int y2 = 4;
        IntensityDescriptor desc2 = features2.extractIntensity(x2, y2, 0);
        
        FeatureComparisonStat stat = Features.calculateIntensityStat(
            desc1, x1, y1, desc2, x2, y2);
        
        float sqSumIntensityDiff = stat.getSumIntensitySqDiff();
        
        float sqIntensityErr = stat.getImg2PointIntensityErr();
                
        float expSSD = (40 * 40);
        
        double expSqErr = Math.sqrt(12) * 25.;
        
        // 3 * sigma for large confidence.  stdev=12
        assertTrue(Math.abs(sqSumIntensityDiff - expSSD) < 3 * expSqErr);
        
        //assertTrue(Math.abs(sqErr - expSqErr) < 3 * Math.sqrt(expSqErr));
    }
    
    public void testExtractGradient() {
        // not implemented yet
    }
    
}
