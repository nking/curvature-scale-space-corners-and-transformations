package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.security.SecureRandom;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class EdgeMatcherTest extends TestCase {
    
    private Logger log = Logger.getLogger(EdgeMatcherTest.class.getName());

    public EdgeMatcherTest() {
        
    }
    
    public void testCountMaxMatchable() throws Exception {
        
        PairIntArray[] edges1 = new PairIntArray[2];
        PairIntArray[] edges2 = new PairIntArray[2];
        
        edges1[0] = new PairIntArray();
        edges1[0].add(1, 2);
        edges1[0].add(1, 3);
        
        edges1[1] = new PairIntArray();
        edges1[1].add(10, 2);
        edges1[1].add(10, 3);
        
        edges2[0] = new PairIntArray();
        edges2[0].add(5, 2);
        edges2[0].add(5, 3);
        edges2[0].add(5, 4);
        
        edges2[1] = new PairIntArray();
        edges2[1].add(14, 2);
        edges2[1].add(14, 3);
        edges2[1].add(14, 4);
        
        EdgeMatcher matcher = new EdgeMatcher();
        int nmxsum = matcher.countMaxMatchable(edges1, edges2);
        
        assertTrue(nmxsum == 4);
    }

    public void testEvaluate() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1437080777162L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        int imageWidth = 650;
        int imageHeight = 400;
        float scale = 1.0f;

        float rotInDegrees = (int)(sr.nextFloat() * 360);
        
        int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
        int transY = (int) (0.05f * sr.nextFloat() * imageHeight);
        if (sr.nextBoolean()) {
            transX *= -1;
        }
        if (sr.nextBoolean()) {
            transY *= -1;
        }

        TransformationParameters params = new TransformationParameters();
        params.setRotationInDegrees(rotInDegrees);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        params.setOriginX(0);
        params.setOriginY(0);

        PairIntArray[] edges1 = new PairIntArray[2];
        
        int nPoints = 100;
        
        for (int i = 0; i < edges1.length; ++i) {
            edges1[i] = new PairIntArray();
            for (int j = 0; j < nPoints; ++j) {
                int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                edges1[i].add(x, y);
            }
        }

        // ===  transform the right points  ======
        
        Transformer transformer = new Transformer();

        PairIntArray[] edges2 = new PairIntArray[2];
        
        for (int i = 0; i < edges1.length; ++i) {
            edges2[i] = transformer.applyTransformation(
                params, edges1[i]);
        }
        
        float tolTransX = 8;
        float tolTransY = 8;
        
        EdgeMatcher matcher = new EdgeMatcher();
        
        boolean useGreedyMatching = true;
        
        TransformationPointFit fit = matcher.evaluate(params, edges1, edges2, 
            tolTransX, tolTransY, useGreedyMatching);
       
        int nExpected = edges1.length * nPoints;
        int nMaxMatchable = nExpected;
        int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);
                
        TransformationParameters fitParams = fit.getParameters();
        int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
        float diffRotDeg = getAngleDifference(
            fitParams.getRotationInDegrees(), rotInDegrees);
        float diffScale = Math.abs(fitParams.getScale() - scale);
        float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
        float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

        log.info("nMaxMatchable=" + nMaxMatchable);
        
        log.info("diff =" + String.format(
            " dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d dNPoints=%d meanDiffModel=%f",
            diffRotDeg, diffScale, diffTransX,
            diffTransY, nEps, diffN, (float) fit.getMeanDistFromModel()));
        
        float epsRot = 1;
        float epsTrans = 2;
        boolean converged = false;
        
        if ((diffN <= nEps) && (diffRotDeg <= epsRot)
            && (diffScale < 0.2) && (diffTransX <= epsTrans)
            && (diffTransY <= epsTrans)) {
            converged = true;
        }
        
        assertTrue(converged);
    }
    
    public void testRefine() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1437080777162L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        int imageWidth = 650;
        int imageHeight = 400;
        float scale = 1.0f;

        float rotInDegrees = (int)(sr.nextFloat() * 360);
        
        int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
        int transY = (int) (0.05f * sr.nextFloat() * imageHeight);
        if (sr.nextBoolean()) {
            transX *= -1;
        }
        if (sr.nextBoolean()) {
            transY *= -1;
        }

        TransformationParameters params = new TransformationParameters();
        params.setRotationInDegrees(rotInDegrees);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        params.setOriginX(0);
        params.setOriginY(0);

        PairIntArray[] edges1 = new PairIntArray[2];
        
        int nPoints = 100;
        
        for (int i = 0; i < edges1.length; ++i) {
            edges1[i] = new PairIntArray();
            for (int j = 0; j < nPoints; ++j) {
                int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                edges1[i].add(x, y);
            }
        }

        // ===  transform the right points  ======

        Transformer transformer = new Transformer();

        PairIntArray[] edges2 = new PairIntArray[2];
        
        for (int i = 0; i < edges1.length; ++i) {
            edges2[i] = transformer.applyTransformation(
                params, edges1[i]);
        }
        
        EdgeMatcher matcher = new EdgeMatcher();
        
        TransformationPointFit fit = matcher.refineTransformation(edges1, 
            edges2, params);
       
        int nExpected = edges1.length * nPoints;
        int nMaxMatchable = nExpected;
        int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);
                
        TransformationParameters fitParams = fit.getParameters();
        int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
        float diffRotDeg = getAngleDifference(
            fitParams.getRotationInDegrees(), rotInDegrees);
        float diffScale = Math.abs(fitParams.getScale() - scale);
        float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
        float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

        log.info("nMaxMatchable=" + nMaxMatchable);
        
        log.info("diff =" + String.format(
            " dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d dNPoints=%d meanDiffModel=%f",
            diffRotDeg, diffScale, diffTransX,
            diffTransY, nEps, diffN, (float) fit.getMeanDistFromModel()));
        
        float epsRot = 1;
        float epsTrans = 2;
        boolean converged = false;
        
        if ((diffN <= nEps) && (diffRotDeg <= epsRot)
            && (diffScale < 0.2) && (diffTransX <= epsTrans)
            && (diffTransY <= epsTrans)) {
            converged = true;
        }
        
        assertTrue(converged);
        
        // make small changes in params and see if a better solution is found
        TransformationParameters params2 = new TransformationParameters();
        float rot2 = params.getRotationInDegrees() + 1;
        if (rot2 > 359) {
            rot2 = rot2 - 360;
        }
        params2.setRotationInDegrees(rot2);
        params2.setScale(params.getScale());
        params2.setTranslationX(params.getTranslationX() - 2);
        params2.setTranslationY(params.getTranslationY() + 2);
        
        fit = matcher.refineTransformation(edges1, edges2, params2);
       
        fitParams = fit.getParameters();
        diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
        diffRotDeg = getAngleDifference(
            fitParams.getRotationInDegrees(), rotInDegrees);
        diffScale = Math.abs(fitParams.getScale() - scale);
        diffTransX = Math.abs(fitParams.getTranslationX() - transX);
        diffTransY = Math.abs(fitParams.getTranslationY() - transY);

        log.info("nMaxMatchable=" + nMaxMatchable);
        
        log.info("diff =" + String.format(
            " dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d dNPoints=%d meanDiffModel=%f",
            diffRotDeg, diffScale, diffTransX,
            diffTransY, nEps, diffN, (float) fit.getMeanDistFromModel()));
        
        converged = false;
        
        if ((diffN <= nEps) && (diffRotDeg <= epsRot)
            && (diffScale < 0.2) && (diffTransX <= epsTrans)
            && (diffTransY <= epsTrans)) {
            converged = true;
        }
        
        assertTrue(converged);
        
    }
    
    private float getAngleDifference(float rotDegrees0, float rotDegrees1) {
         /*
         I  |  0
        ---------
         II | III
        */
        int q0 = 0;
        if (rotDegrees0 >= 270) {
            q0 = 3;
        } else if (rotDegrees0 >= 180) {
            q0 = 2;
        } else if (rotDegrees0 >= 90) {
            q0 = 1;
        }
        int q1 = 0;
        if (rotDegrees1 >= 270) {
            q1 = 3;
        } else if (rotDegrees1 >= 180) {
            q1 = 2;
        } else if (rotDegrees1 >= 90) {
            q1 = 1;
        }

        /*
         I  |  0
        ---------
         II | III
        */
        float angleDiff = -1;
        if (q0 == 0){
            if (q1 == 0) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else if (q1 == 1) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 2) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else {
                angleDiff = Math.abs(360 - rotDegrees1 + rotDegrees0);
            }
        } else if (q0 == 1) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 1) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else if (q1 == 2) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            }
        } else if (q0 == 2) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 1) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else if (q1 == 2) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else {
                angleDiff = (rotDegrees1 - rotDegrees0);
            }
        } else if (q0 == 3) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (360 - rotDegrees0 + rotDegrees1);
            } else if (q1 == 1) {
                float diff = (rotDegrees0 - rotDegrees1);
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 2) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            }
        }

        if (angleDiff > 359) {
            angleDiff = angleDiff - 360;
        }

        return angleDiff;
    }

}
