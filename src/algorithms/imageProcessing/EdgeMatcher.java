package algorithms.imageProcessing;

import static algorithms.imageProcessing.EdgeMatcher.minTolerance;
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

    protected static int minTolerance = 5;

    //TODO: this has to be a high number for sets with projection.
    // the solution is sensitive to this value.
    private final float generalTolerance = 8;

    public static float toleranceGridFactor = 4.f;

    protected boolean debug = true;

    /**
     * the maximum number of iterations that the refinement of translation
     * will use in the downhill simplex.  The default value is 50.
     */
    private int dsNMaxIter = 50;
    protected void setDsNMaxIter(int n) {
        dsNMaxIter = n;
    }
    protected int getDsNMaxIter() {
        return dsNMaxIter;
    }
    float nEpsFactor = 2.0f;
    protected void setNEpsFactor(float f) {
        nEpsFactor = f;
    }
    protected float getNEpsFactor() {
        return nEpsFactor;
    }
    protected float cellFactor = 1.25f;
    protected float tolFactor = 0.5f;
    protected void setCellFactor(float f) {
        cellFactor = f;
    }
    protected void setTolFactor(float f) {
        tolFactor = f;
    }
    protected float getCellFactor() {
        return cellFactor;
    }
    protected float getTolFactor() {
        return tolFactor;
    }
  
    // ======= code that needs testing and revision the most
    /**
     * refine the transformation params to make a better match of edges1 to
     * edges2 where the points within edges in both sets are not necessarily
     * 1 to 1 matches (that is, the input is not expected to be matched
     * already).
     *
     * TODO: improve transformEdges to find translation for all edges
     * via a search method rather than trying all pairs of points.
     *
     * @param edges1
     * @param edges2
     * @param params
     * @return
     */
    public TransformationParameters refineTransformation(PairIntArray[] edges1,
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
        
        /*
        REPLACE WITH A SIMPLEX OF SMALL RANGE
        */
        
        float rotHalfRangeInDegrees = 5; 
        float rotDeltaInDegrees = 1;
        float transXHalfRange = 10; 
        float transXDelta = 1;
        float transYHalfRange = transXHalfRange; 
        float transYDelta = transXDelta;
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit[] fits = new TransformationPointFit[n];
        
        TransformationPointFit bestFit = null;
        
        boolean useGreedyMatching = true;
        
        for (int i = 0; i < n; ++i) {
            
            PairIntArray set1 = edges1[i];
            PairIntArray set2 = edges2[i];
            
            TransformationPointFit fit = pointMatcher.refineTheTransformation(
                params, set1, set2, 
                rotHalfRangeInDegrees, rotDeltaInDegrees, 
                transXHalfRange, transXDelta, 
                transYHalfRange, transYDelta,
                useGreedyMatching);
            
            fits[i] = fit;
            
            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
        }
        
        //TODO: consider applying the fits to all points and keeping the
        // best
        
        return (bestFit != null) ? bestFit.getParameters() : null;
    }

}
