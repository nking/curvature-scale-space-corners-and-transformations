package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
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
        matchedXY2.add(157, 107); 
        matchedXY2.add(159, 101); 
        matchedXY2.add(190, 143); 
        matchedXY2.add(177, 150); 
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
        
        assertTrue(Math.abs(params.getTranslationX() - 110.6) < (centroidX1*0.02));
        
        assertTrue(Math.abs(params.getTranslationY() - 19.8) < (centroidY1*0.02));
    }

    /*
    public double[] applyTransformation(TransformationParameters params, 
        int centroidX1, int centroidY1,
        double x1, double y1) {
    
    TransformationParameters swapReferenceFrames(TransformationParameters 
        params, int centroidX2, int centroidY2,
        double x1, double y1, double x2, double y2) {
    */
    
    @Test
    public void testApplyTransformation() throws Exception {
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        
        int transX = 125;
        int transY = 14;
        double transXTol = 10.3;
        double transYTol = 5.9;
        double rotation = 25 * Math.PI/180.;
        double scale = 4;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)rotation);
        params.setScale((float)scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        
        for (int i = 0; i < set1.getN(); i++) {
            
            double[] x1y1 = tc.applyTransformation(params, 
                centroidX1, centroidY1,
                (double)set1.getX(i), (double)set1.getY(i));
            
            assertTrue(Math.abs(x1y1[0] - set2.getX(i)) <= 1);
            assertTrue(Math.abs(x1y1[1] - set2.getY(i)) <= 1);
        }
    }
    
    @Test
    public void testSwapReferenceFrames() throws Exception {
        
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        
        int transX = 125;
        int transY = 14;
        double rotation = 25 * Math.PI/180.;
        double scale = 4;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        int centroidX2 = 400;
        int centroidY2 = 400;
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)rotation);
        params.setScale((float)scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
          
        TransformationParameters revParams = tc.swapReferenceFrames(params, 
            centroidX2, centroidY2, 
            set1.getX(0), set1.getX(0), set2.getX(0), set2.getY(0));
      
        for (int i = 0; i < set1.getN(); i++) {
            
            double[] x1y1 = tc.applyTransformation(revParams, 
                centroidX2, centroidY2,
                (double)set2.getX(i), (double)set2.getY(i));
            
            assertTrue(Math.abs(x1y1[0] - set1.getX(i)) <= 1);
            assertTrue(Math.abs(x1y1[1] - set1.getY(i)) <= 1);
        }
        
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
