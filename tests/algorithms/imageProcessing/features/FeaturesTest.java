package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.Gaussian1DFirstDeriv;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.util.ResourceFinder;
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
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
        GreyscaleImage gradientImg = img.createWithDimensions();
        GreyscaleImage thetaImg = img.createWithDimensions();
        
        Features features = new Features(img, gradientImg, thetaImg, 2, false,
            rOffsets);
        
        assertTrue(features.isWithinXBounds(2));
        assertFalse(features.isWithinXBounds(-2));
        
        assertTrue(features.isWithinYBounds(2));
        assertFalse(features.isWithinYBounds(-2));
        
    }
    
    public void testIsWithinXBounds1() {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        Image img = new Image(10, 10);
        GreyscaleImage gradientImg = new GreyscaleImage(10, 10);
        GreyscaleImage thetaImg = new GreyscaleImage(10, 10);
        
        Features features = new Features(img, gradientImg, thetaImg, 2, false,
            rOffsets);
        
        assertTrue(features.isWithinXBounds(2));
        assertFalse(features.isWithinXBounds(-2));
        
        assertTrue(features.isWithinYBounds(2));
        assertFalse(features.isWithinYBounds(-2));
        
    }

    public void testExtractIntensity() {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        boolean normalize = false;
        
        GreyscaleImage img = new GreyscaleImage(20, 20);
        for (int i = 0; i < img.getNPixels(); ++i) {
            img.setValue(i, 100);
        }
        
        GreyscaleImage gradientImg = img.createWithDimensions();
        GreyscaleImage thetaImg = img.createWithDimensions();
        
        Features features = new Features(img, gradientImg, thetaImg, 2, 
            normalize, rOffsets);
        
        IntensityDescriptor desc = features.extractIntensity(10, 10, 0);
        assertNotNull(desc);
        assertTrue(desc.isNormalized() == normalize);
        
        float errSq = desc.sumSquaredError();
        assertTrue(errSq == 0);
        assertTrue(desc.calculateSSD(desc) == 0);
        
        
    }

    public void testCalculateIntensityStat() throws Exception {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1437335464716L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        GreyscaleImage img = new GreyscaleImage(20, 20);
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
        
        GreyscaleImage img2 = new GreyscaleImage(20, 20);
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
        
        GreyscaleImage gradientImg = img.createWithDimensions();
        Features features = new Features(img, gradientImg, thetaImg, 2,
            normalize, rOffsets);
                
        int x1 = 10;
        int y1 = 10;
        IntensityDescriptor desc1 = features.extractIntensity(x1, y1, 0);        
        
        Features features2 = new Features(img2, gradientImg, thetaImg, 2, 
            normalize, rOffsets);
        int x2 = 11;
        int y2 = 9;
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
    
    public void estExtractIntensity2() throws Exception {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        String filePath = ResourceFinder.findFileInTestResources("blox.gif");
        
        Image img = ImageIOHelper.readImage(filePath);
        
        GreyscaleImage gsImg = img.copyToGreyscale();
                
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // 1      0   -1
        float[] kernel = Gaussian1DFirstDeriv.getKernel(SIGMA.ZEROPOINTSEVENONE);
                
        GreyscaleImage gX = gsImg.copyImage();
        GreyscaleImage gY = gsImg.copyImage();        
        imageProcessor.applyKernel1D(gX, kernel, true);
        imageProcessor.applyKernel1D(gY, kernel, false);
        GreyscaleImage gXY = imageProcessor.combineConvolvedImages(gX, gY);
        
        GreyscaleImage thetaImg = imageProcessor.computeTheta360(gX, gY);

        int blockHalfWidths = 2;
        boolean useNormalizedIntensities = false;
                
        Features features = new Features(img, gXY, thetaImg, 
            blockHalfWidths, useNormalizedIntensities, rOffsets);
        
        int xCenter = 167;
        int yCenter = 109;
        int rotationInDegrees = 0;
        ClrIntensityDescriptor desc = (ClrIntensityDescriptor)features.extractIntensity(xCenter, 
            yCenter, rotationInDegrees);
        
        /*
        [227] [214] [169] [164] [159]     107
        [237] [225] [199] [170] [154]     108
        [239] [231]*[224]*[175] [154]     109
        [241] [226] [206] [174] [154]     110
        [188] [159] [157] [163] [149]     111
         165   166   167   168   169
        */
        
        int[] expectedInternal = new int[]{
            227, 237, 239, 241, 188,
            214, 225, 231, 226, 159,
            169, 199, 224, 206, 157,
            164, 170, 175, 174, 163,
            159, 154, 154, 154, 149
        };
        
        assertTrue(Arrays.equals(expectedInternal, desc.getInternalRedArrayCopy()));
        
        rotationInDegrees = 90;
        desc = (ClrIntensityDescriptor)features.extractIntensity(xCenter, 
            yCenter, rotationInDegrees);
        /*
        [227] [214] [169] [164] [159]     107
        [237] [225] [199] [170] [154]     108
        [239] [231]*[224]*[175] [154]     109
        [241] [226] [206] [174] [154]     110
        [188] [159] [157] [163] [149]     111
         165   166   167   168   169
        
        [159] [   ] [   ] [   ] [   ]     107
        [164] [   ] [   ] [   ] [   ]     108
        [169] [   ]*[   ]*[   ] [   ]     109
        [214] [   ] [   ] [   ] [   ]     110
        [227] [237] [239] [241] [188]     111
         165   166   167   168   169
        
        */
        
        expectedInternal = new int[]{
            159, 164, 169, 214, 227,
            154, 170, 199, 225, 237, 
            154, 175, 224, 231, 239, 
            154, 174, 206, 226, 241, 
            149, 163, 157, 159, 188
        };
        
        assertTrue(Arrays.equals(expectedInternal, desc.getInternalRedArrayCopy()));
        
        
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(thetaImg.getWidth() >> 1);
        params.setOriginY(thetaImg.getHeight() >> 1);
        params.setScale(1);
        params.setTranslationX(0);
        params.setTranslationY(0);
        params.setRotationInDegrees(rotationInDegrees);
                
        Transformer transformer = new Transformer();

        /*
        public IntensityDescriptor extractIntensity(final int xCenter, 
        final int yCenter, final int rotation) {
        */
    }
    
    public void testExtractGradient() {
        // not implemented yet
    }
    
}
