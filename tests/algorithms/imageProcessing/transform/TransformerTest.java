package algorithms.imageProcessing.transform;

import java.security.SecureRandom;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class TransformerTest extends TestCase {
    
    public TransformerTest() {
    }
    
    public void testSwapReferenceFramesWRTOrigin() throws Exception {
        
        Transformer transformer = new Transformer();
        
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();
        
        TransformationParameters params0 = new TransformationParameters();
        params0.setRotationInDegrees(0);
        params0.setScale(1);
        params0.setTranslationX(0);
        params0.setTranslationY(0);
        params0.setOriginX(0);
        params0.setOriginY(0);
        
        int x0 = 10;
        int y0 = 20;
        double[] trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - x0) < 0.01);
        assertTrue(Math.abs(trXY[1] - y0) < 0.01);
        
        TransformationParameters revParams0 = tc.swapReferenceFrames(params0);
        double[] revTrXY = transformer.applyTransformation(revParams0, 
            trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
            
        
        params0.setRotationInDegrees(90);
        trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - y0) < 0.01);
        assertTrue(Math.abs(trXY[1] - -x0) < 0.01);
        
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
        
        
        params0.setRotationInDegrees(180);
        trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - -x0) < 0.01);
        assertTrue(Math.abs(trXY[1] - -y0) < 0.01);
        
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
        
        
        params0.setRotationInDegrees(270);
        trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - -y0) < 0.01);
        assertTrue(Math.abs(trXY[1] - x0) < 0.01);

        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
        
        
        params0.setRotationInDegrees(0);
        params0.setTranslationX(10);
        params0.setTranslationY(10);
        trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - (x0 + params0.getTranslationX())) < 0.01);
        assertTrue(Math.abs(trXY[1] - (y0 + params0.getTranslationY())) < 0.01);
        
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
        
        
        params0.setRotationInDegrees(90);
        params0.setTranslationX(10);
        params0.setTranslationY(10);
        trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - (y0 + params0.getTranslationX())) < 0.01);
        assertTrue(Math.abs(trXY[1] - (-x0 + params0.getTranslationY())) < 0.01);
        
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
        
        //------------------
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int n = 100;
        for (int i = 0; i < n; ++i) {
            int x = sr.nextInt(1000);
            int y = sr.nextInt(1000);
            float scale = sr.nextFloat();
            if (sr.nextBoolean()) {
                scale = 1.f/scale;
            }
            int rotDeg = sr.nextInt(360);
            float transX = sr.nextInt(1000);
            float transY = sr.nextInt(1000);
            
            TransformationParameters params = new TransformationParameters();
            params.setRotationInDegrees(rotDeg);
            params.setScale(scale);
            params.setTranslationX(transX);
            params.setTranslationY(transY);
            
            TransformationParameters revParams = tc.swapReferenceFrames(params);
            
            double[] xyTransformed = transformer.applyTransformation(params, 
                x, y);
            
            double[] xyTransformedRevTransformed = 
                transformer.applyTransformation(
                revParams, (float)xyTransformed[0], (float)xyTransformed[1]);
            
            double diffX = Math.abs(xyTransformedRevTransformed[0] - x);
            double diffY = Math.abs(xyTransformedRevTransformed[1] - y);
            
            if ((diffX >= 1) || (diffY >= 1)) {
                System.out.println("diffX=" + diffX +  " diffY=" + diffY);
            }
            assertTrue(diffX < 1);
            
            assertTrue(diffY < 1);
        }

    }

    public void testSwapReferenceFrames() throws Exception {
        
        Transformer transformer = new Transformer();
        
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();
        
        TransformationParameters params0 = new TransformationParameters();
        params0.setRotationInDegrees(0);
        params0.setScale(1);
        params0.setTranslationX(0);
        params0.setTranslationY(0);
        params0.setOriginX(50);
        params0.setOriginY(60);
        
        int x0 = 10;
        int y0 = 20;
        double[] trXY = transformer.applyTransformation(params0, x0, y0);
        assertTrue(Math.abs(trXY[0] - x0) < 0.01);
        assertTrue(Math.abs(trXY[1] - y0) < 0.01);
        
        TransformationParameters revParams0 = tc.swapReferenceFrames(params0);
        double[] revTrXY = transformer.applyTransformation(revParams0, 
            trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - x0) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - y0) < 0.01);
            
        float xa = 0;
        float ya = 0;
        float xb = 50;
        float yb = 0;
        params0.setRotationInDegrees(90);
        params0.setOriginX(50);
        params0.setOriginY(60);
        trXY = transformer.applyTransformation(params0, xa, ya);
        assertTrue(Math.abs(trXY[0] - -10) < 0.01);
        assertTrue(Math.abs(trXY[1] - 110) < 0.01);
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - xa) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - ya) < 0.01);
        
        trXY = transformer.applyTransformation(params0, xb, yb);
        assertTrue(Math.abs(trXY[0] - -10) < 0.01);
        assertTrue(Math.abs(trXY[1] - 60) < 0.01);
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - xb) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - yb) < 0.01);
        
        params0.setRotationInDegrees(180);
        params0.setOriginX(50);
        params0.setOriginY(60);
        trXY = transformer.applyTransformation(params0, xa, ya);
        assertTrue(Math.abs(trXY[0] - 100) < 0.01);
        assertTrue(Math.abs(trXY[1] - 120) < 0.01);
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - xa) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - ya) < 0.01);
        
        trXY = transformer.applyTransformation(params0, xb, yb);
        assertTrue(Math.abs(trXY[0] - 50) < 0.01);
        assertTrue(Math.abs(trXY[1] - 120) < 0.01);
        revParams0 = tc.swapReferenceFrames(params0);
        revTrXY = transformer.applyTransformation(revParams0, trXY[0], trXY[1]);
        assertTrue(Math.abs(revTrXY[0] - xb) < 0.01);
        assertTrue(Math.abs(revTrXY[1] - yb) < 0.01);
        
    }
    
    public void testScale() throws Exception {
        
        Transformer transformer = new Transformer();
        
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(0);
        params.setOriginY(0);
        params.setTranslationX(2);
        params.setTranslationY(-3);
        params.setScale(1);
        params.setRotationInDegrees(45);
        
        /*
        transform the parameters to scale by 4.
        transform the parameters to scale by 1/4
        assert the final params are same as original
        */
        
        TransformationParameters paramsS4 = 
            transformer.applyScaleTransformation(params, 4);
        
        TransformationParameters paramsSDiv4 = 
            transformer.applyScaleTransformation(paramsS4, 0.25f);
        
        assertTrue(Math.abs(paramsSDiv4.getOriginX() - params.getOriginX()) < 0.01);
        assertTrue(Math.abs(paramsSDiv4.getScale() - params.getScale()) < 0.01);
        assertTrue(Math.abs(paramsSDiv4.getTranslationX()- params.getTranslationX()) < 0.01);
        assertTrue(Math.abs(paramsSDiv4.getTranslationY() - params.getTranslationY()) < 0.01);
        assertTrue(Math.abs(paramsSDiv4.getOriginY() - params.getOriginY()) < 0.01);
        assertTrue(Math.abs(paramsSDiv4.getRotationInDegrees() - params.getRotationInDegrees()) < 0.01);
    }
    
    public void testTransformXY() throws Exception {
            
        float rotationInDegrees = 90;
        
        float[][] xy = new float[4][];
        xy[0] = new float[]{0, 2};
        xy[1] = new float[]{2, 0};
        xy[2] = new float[]{0, -3};
        xy[3] = new float[]{-2, 0};
        
        float[][] expected = new float[4][];
        expected[0] = new float[]{2, 0};
        expected[1] = new float[]{0, -2};
        expected[2] = new float[]{-3, 0};
        expected[3] = new float[]{0, 2};
        
        Transformer transformer = new Transformer();
        
        float[][] results = transformer.transformXY(rotationInDegrees, xy);
        
        assertTrue(results.length == expected.length);
        
        for (int i = 0; i < results.length; ++i) {
            assertTrue(results[i].length == expected[i].length);
            
            float x = results[i][0];
            float y = results[i][1];
            
            float xExp = expected[i][0];
            float yExp = expected[i][1];
            
            assertTrue(Math.abs(x - xExp) < 0.1);
            assertTrue(Math.abs(y - yExp) < 0.1);
        }
        
    }
}
