package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class MatchedPointsTransformationCalculatorTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public MatchedPointsTransformationCalculatorTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testCalulateEuclideanGivenScale() throws Exception {
        
        /*
        String fileName1 = "closed_curve.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

        String fileName2 = "closed_curve_translate_scale_rotate20.png";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();
        
        mapper.initialize();
        mapper.createMatchedPointArraysFromContourPeaks();
        
        PairIntArray matchedXY1 = mapper.getMatchedXY1();
        if (matchedXY1.getN() < 3) {
            throw new IllegalStateException("need at least 3 points");
        }
        
        PairIntArray matchedXY2 = mapper.getMatchedXY2();
        
        double centroidX1 = mapper.image1OriginalWidth >> 1;
        double centroidY1 = mapper.image1OriginalHeight >> 1;
        double matchedScale = mapper.getMatchedScale();
        float[] matchedXY1Weights = mapper.getMatchedXY1Weights();
        float[] matchedXY2Weights = mapper.getMatchedXY2Weights();
        
        log.info("matchedScale=" + matchedScale);
        log.info("matchedXY1=" + matchedXY1.toString());
        log.info("matchedXY2=" + matchedXY2.toString());
        log.info("centroidX1=" + centroidX1);
        log.info("centroidY1=" + centroidY1);
        log.info("matchedXY1Weights=" + Arrays.toString(matchedXY1Weights));
        log.info("matchedXY2Weights=" + Arrays.toString(matchedXY2Weights));
        */
        
        double matchedScale = 1.2968395948410034;
        double centroidX1 = 58.0;
        double centroidY1 = 72.0;        
        float[] matchedXY1Weights = new float[] {
            0.18644243f, 0.18644243f, 0.15677878f, 0.15677878f, 0.15677878f, 
            0.15677878f
        };
        float[] matchedXY2Weights = new float[] {
            0.19142407f, 0.19142407f, 0.14760813f, 0.14760813f, 0.1609678f, 
            0.1609678f
        };
        PairIntArray matchedXY1 = new PairIntArray();
        matchedXY1.add(34, 78); 
        matchedXY1.add(35, 72); 
        matchedXY1.add(70, 93); 
        matchedXY1.add(61, 99); 
        matchedXY1.add(60, 45); 
        matchedXY1.add(69, 54);
        
        PairIntArray matchedXY2 = new PairIntArray();
        matchedXY2.add(150, 117); 
        matchedXY2.add(159, 91); 
        matchedXY2.add(196, 147); 
        matchedXY2.add(176, 156); 
        matchedXY2.add(200, 85);
        matchedXY2.add(209, 97);
     
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();
        
        tc.useDebugMode();
        
        TransformationParameters params = tc.calulateEuclideanGivenScale(
            matchedScale,
            matchedXY1, matchedXY1Weights, matchedXY2, matchedXY2Weights, 
            centroidX1, centroidY1);
        
        assertTrue(Math.abs(params.getRotationInDegrees() - 340) < 10.f);
             
        assertTrue(Math.abs(params.getScale() - 1.3) < 0.1);
        
        assertTrue(Math.abs(params.getTranslationX() - 110.9) < (centroidX1*0.02));
        
        assertTrue(Math.abs(params.getTranslationY() - 21.3) < (centroidY1*0.02));
    }

    public static void main(String[] args) {
        
        try {
            
            MatchedPointsTransformationCalculatorTest test = 
                new MatchedPointsTransformationCalculatorTest();
            
            test.testCalulateEuclideanGivenScale();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
    
}
