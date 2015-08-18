package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.List;

/**
 * class to map contours from edges from one image to another and return the
 * matched inflection points and a transformation matrix that can
 * be applied to the first to put it into the frame of the second.
 * 
 * The algorithm used for matching scale space image contours is documented in 
 * CSSContourMatcherWrapper
 * @see algorithms.imageProcessing.CSSContourMatcherWrapper
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionEdgeMapper extends 
    AbstractCurvatureScaleSpaceInflectionMapper {
        
    /*
    protected final ImageExt image1;
    protected final ImageExt image2;
    // for debugging, keeping a reference of originals
    protected final Image originalImage1;
    protected final Image originalImage2;
    
    public final int image1OriginalWidth;
    public final int image1OriginalHeight;
    protected final int image2OriginalWidth;
    protected final int image2OriginalHeight;
    */
    
    public CurvatureScaleSpaceInflectionEdgeMapper(
        List<PairIntArray> edges1, List<PairIntArray> edges2,
        int offsetImageX1, int offsetImageY1, int offsetImageX2, int offsetImageY2) {
        
        this.edges1 = edges1;
        this.edges2 = edges2;
        
        this.offsetImageX1 = offsetImageX1;
        this.offsetImageY1 = offsetImageY1;
        this.offsetImageX2 = offsetImageX2;
        this.offsetImageY2 = offsetImageY2;
    }
   
    public TransformationParameters createEuclideanTransformationImpl() {
        
        if (bestFittingParameters == null) {
            return null;
        }
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        
        if (debug) {
            tc.useDebugMode();
        }
                
        //TODO: temporarily disabling the refinement while fixing PointMatcher
        if (doRefineTransformations) {
            
            boolean reverseDatasetOrder = bestFittingParameters.getScale() < 1.0;
            
            log.info("BEFORE REFINEMENT:\n" + bestFittingParameters.toString());
            
            PairIntArray[] set1 = getMatchedEdges1InOriginalReferenceFrameArray();
            PairIntArray[] set2 = getMatchedEdges2InOriginalReferenceFrameArray();
            EdgeMatcher matcher = new EdgeMatcher();
            TransformationPointFit fit2 = null;
            if (reverseDatasetOrder) {
                fit2 = matcher.refineTransformation(set2, set1, bestFittingParameters);
            } else {
                fit2 = matcher.refineTransformation(set1, set2, bestFittingParameters);
            }
            
            if (reverseDatasetOrder) {
                bestFittingParameters = tc.swapReferenceFrames(bestFittingParameters);            
            }
            
            if (fit2 != null) {
                log.info("FINAL:\n" + fit2.toString());
                bestFittingParameters = fit2.getParameters();
            }
        }
        
        return bestFittingParameters;
    }

   
    @Override
    protected List<PairIntArray> getEdges(CurvatureScaleSpaceImageMaker imgMaker) {
        
        return imgMaker.getClosedCurves();
    }

    @Override
    protected void createEdges1() {
        
        // already given in constructor
       
    }

    @Override
    protected void createEdges2() {
        
        // already given in constructor
    }

}
