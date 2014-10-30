package algorithms.imageProcessing;

import algorithms.ResourceFinder;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class TransformationCalculatorTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public TransformationCalculatorTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testSortFromMinToMax() throws Exception {
        
        System.out.println("sortFromMinToMax");
        
        int n = 100;
        
        TransformationParameters mockParams = new TransformationParameters();
        PairIntArray[] mockEdges = new PairIntArray[1];
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());

        TransformationFit[] fits = new TransformationFit[n];
        
        for (int i = 0; i < n; i++) {
            
            int num = sr.nextInt(100);
            
            double chiSqSum = num * sr.nextDouble();
            
            TransformationFit tf = new TransformationFit(mockParams, mockEdges);
           
            tf.setChiSqSum(chiSqSum);
            
            fits[i] = tf;
        }
                
        TransformationRefiner instance = new TransformationRefiner();
        instance.sortFromMinToMax(fits, 0, fits.length - 1);
        
        double lastCSM = fits[0].getChiSqSum();
        for (int i = 1; i < fits.length; i++) {
            double csm = fits[i].getChiSqSum();
            assertTrue(csm >= lastCSM);
            lastCSM = csm;
        }
    }
  
    @Test
    public void testDifferenceOfCurves_PairIntArrayArr_PairIntArrayArr() {
        
        System.out.println("differenceOfCurves");
        
        PairIntArray edge1 = new PairIntArray();
        PairIntArray edge2 = new PairIntArray();
        for (int i = 0; i < 10; i++) {
            edge1.add(0, i);
            edge2.add(2, i);
        }
        PairIntArray[] edges1 = new PairIntArray[]{edge1};
        SearchableCurve[] edges2 = new SearchableCurve[]{
            new SearchableCurve(edge2)};
        
        TransformationRefiner instance = new TransformationRefiner();
        
        double expResult = Math.sqrt(10 * (2*2));
        double result = instance.differenceOfCurves(edges1, edges2);
        assertTrue(expResult == result);
    
        result = instance.differenceOfCurves(edge1, new SearchableCurve(edge2));
        assertTrue(expResult == result);
    }
      
    @Test
    public void testRefineTransformation() throws Exception {
        
        System.out.println("refineTransformation");
       
        String fileName1 = "closed_curve.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

        String fileName2 = "closed_curve_translate_scale_rotate20.png";

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();

        TransformationParameters transformationParams = 
            mapper.createEuclideanTransformation();
     
        PairInt[] indexes = mapper.getMatchedEdgesIndexes();
        
        PairIntArray[] edges1 = new PairIntArray[indexes.length];
        PairIntArray[] edges2 = new PairIntArray[indexes.length];
        
        for (int i = 0; i < indexes.length; i++) {
            PairInt idxes = indexes[i];
            edges1[i] = mapper.getEdges1().get(idxes.getX());
            edges2[i] = mapper.getEdges2().get(idxes.getY());
        }
                
        TransformationRefiner instance = new TransformationRefiner();
        TransformationParameters refinedParameters = 
            instance.refineTransformation(edges1, edges2, transformationParams);
        
        log.info("input params:\n" + transformationParams.toString());
        
        log.info("output params:\n" + refinedParameters.toString());
    }

    private void debugDisplay(PairIntArray[] transformedEdges,
        Image img, String rotDegrees) throws IOException {
        
        for (int i = 0; i < transformedEdges.length; i++) {
            PairIntArray edge = transformedEdges[i];
            debugAddCurveToImage(edge, img, 0, 255, 255, 255);
        }     
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/transformed_edges_" + rotDegrees + ".png", img);
    }
    
    private void debugAddCurveToImage(PairIntArray edge, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                        dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }
    }

}
