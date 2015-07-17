package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * class to match the edges extracted from two images.  The use of corners and
 * PointMatcher should be preferred over this.
 
 * @author nichole
 */
public final class EdgeMatcher extends AbstractPointMatcher {

    private final Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = true;

    // ======= code that needs testing and revision the most
    /**
     * refine the transformation params to make a better match of edges1 to
     * edges2 where the points within edges in both sets are not necessarily
     * 1 to 1 matches (that is, they are lines with uneven point intervals
     * and densities).
     *
     * TODO: need to add change of scale to this.
     * 
     * @param edges1
     * @param edges2
     * @param params
     * @return
     */
    public TransformationPointFit refineTransformation(PairIntArray[] edges1,
        PairIntArray[] edges2, final TransformationParameters params) {

        if (edges1 == null || edges1.length == 0) {
            throw new IllegalArgumentException("edges1 cannot be null or empty");
        }
        if (edges2 == null || edges2.length == 0) {
            throw new IllegalArgumentException("edges2 cannot be null or empty");
        }
        if (edges1.length != edges2.length) {
            throw new IllegalArgumentException(
                "edges1 and edges2 must be the same length");
        }

        int n = edges1.length;
        
        log.info("refining coordinates with " + n  + " curves");
        
        float rotHalfRangeInDegrees = 6;
        float rotDeltaInDegrees = 2.f;
        float scaleHalfRange = 0.4f; 
        float scaleDelta = 0.2f;
        float transXHalfRange = 10; 
        float transXDelta = 1.f;
        float transYHalfRange = transXHalfRange; 
        float transYDelta = transXDelta;
        
        float tolTransX = 4;
        float tolTransY = 4;
        
        boolean useGreedyMatching = true;
        
        if (!useGreedyMatching) {
            tolTransX = 2;
            tolTransY = 2;
        }
        
        // start with original fit to make results are an improvement
        TransformationPointFit bestFit = evaluate(params, edges1, edges2,
            tolTransX, tolTransY, useGreedyMatching);
        
        PointMatcher pointMatcher = new PointMatcher();        
                
        int nMaxMatchableSum = countMaxMatchable(edges1, edges2);
        
        float scaleWeighted = 0;
        float rotationWeighted = 0;
        float transXWeighted = 0;
        float transYWeighted = 0;
        
        for (int i = 0; i < n; ++i) {
            
            PairIntArray set1 = edges1[i];
            PairIntArray set2 = edges2[i];
                        
            TransformationPointFit fit = pointMatcher.refineTheTransformation(
                params, set1, set2, 
                scaleHalfRange, scaleDelta,
                rotHalfRangeInDegrees, rotDeltaInDegrees, 
                transXHalfRange, transXDelta, 
                transYHalfRange, transYDelta,
                useGreedyMatching);
            
            if (fit == null) {
                continue;
            }
            
            int nMaxMatchable = Math.min(set1.getN(), set2.getN());
            
            float weight = (float)nMaxMatchable/(float)nMaxMatchableSum;
            
            scaleWeighted += (fit.getScale() * weight);
            rotationWeighted += 
                ((fit.getParameters().getRotationInDegrees()/360.f)  * weight);
            transXWeighted += (fit.getTranslationX() * weight);
            transYWeighted += (fit.getTranslationY() * weight);
        }
        
        rotationWeighted *= Math.PI*2;
        
        TransformationParameters paramsWeighted = new TransformationParameters();
        paramsWeighted.setScale(scaleWeighted);
        paramsWeighted.setRotationInRadians(rotationWeighted);
        paramsWeighted.setTranslationX(transXWeighted);
        paramsWeighted.setTranslationY(transYWeighted);
        paramsWeighted.setOriginX(params.getOriginX());
        paramsWeighted.setOriginY(params.getOriginY());
        
        TransformationPointFit weightedFit = evaluate(paramsWeighted, 
            edges1, edges2,tolTransX, tolTransY, useGreedyMatching);
        
        if (pointMatcher.fitIsBetter(bestFit, weightedFit)) {
            bestFit = weightedFit;
log.info("REFINEMENT ACCEPTED");
        }
        
        return bestFit;
    }

    protected int countMaxMatchable(PairIntArray[] edges1, PairIntArray[] edges2) {
        
        int sum = 0;
        
        for (int i = 0; i < edges1.length; ++i) {
            PairIntArray set1 = edges1[i];
            PairIntArray set2 = edges2[i];
            int nMaxMatchable = Math.min(set1.getN(), set2.getN());
            
            sum += nMaxMatchable;
        }
        
        return sum;
    }

    public TransformationPointFit evaluate(TransformationParameters 
        params, PairIntArray[] edges1, PairIntArray[] edges2,
        float tolTransX, float tolTransY,
        boolean useGreedyMatching) {
        
        if (edges1 == null || edges1.length == 0) {
            throw new IllegalArgumentException("edges1 cannot be null or empty");
        }
        if (edges2 == null || edges2.length == 0) {
            throw new IllegalArgumentException("edges2 cannot be null or empty");
        }
        if (edges1.length != edges2.length) {
            throw new IllegalArgumentException(
                "edges1 and edges2 must be the same length");
        }
        
        PointMatcher pointMatcher = new PointMatcher();
        
        Transformer transformer = new Transformer();
        
        int nMaxMatchableSum = countMaxMatchable(edges1, edges2);
        
        /*
        either need a high tolerance with greedy matching or
        a low tolerance and optimal matching.
        The first is N^2 and the 3nd is approx N^4, respectively.
        */
        
        int nMatched = 0;
        float difFromModelWeighted = 0;
        double stDevFromModelWeightedSqSum = 0;
                
        for (int i = 0; i < edges1.length; ++i) {
            PairIntArray set1 = edges1[i];
            PairIntArray set2 = edges2[i];
            
            PairFloatArray tr = transformer.applyTransformation2(
               params, set1);
            
            int nMaxMatchable = Math.min(set1.getN(), set2.getN());
            float weight = (float)nMaxMatchable/(float)nMaxMatchableSum;
         
            TransformationPointFit fit;
            if (useGreedyMatching) {
                fit = pointMatcher.evaluateFitForUnMatchedTransformedGreedy(
                params, tr, set2, tolTransX, tolTransY);
            } else {
                fit = pointMatcher.evaluateFitForUnMatchedTransformedOptimal(
                params, tr, set2, tolTransX, tolTransY);
            }
            
            nMatched += fit.getNumberOfMatchedPoints();
            difFromModelWeighted += (fit.getMeanDistFromModel() * weight);
            stDevFromModelWeightedSqSum += 
                (fit.getStDevFromMean() * fit.getStDevFromMean() * weight);
        }
        
        float stDevFromModelWeighted = (float)Math.sqrt(stDevFromModelWeightedSqSum);
        
        TransformationPointFit fit = new TransformationPointFit(params, nMatched,
            difFromModelWeighted, stDevFromModelWeighted, tolTransX, tolTransY);

        fit.setMaximumNumberMatchable(nMaxMatchableSum);
        
        return fit;
    }

}
