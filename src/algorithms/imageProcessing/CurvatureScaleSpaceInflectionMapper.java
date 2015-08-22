package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.List;

/**
 * class to map contours from one image to another and return the
 * matched inflection points and a transformation matrix that can
 * be applied to the first to put it into the frame of the second.
 * 
 * The algorithm used for matching scale space image contours is documented in 
 * CSSContourMatcherWrapper
 * @see algorithms.imageProcessing.CSSContourMatcherWrapper
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionMapper extends 
    AbstractCurvatureScaleSpaceInflectionMapper {
        
    protected final ImageExt image1;
    protected final ImageExt image2;
    // for debugging, keeping a reference of originals
    protected final Image originalImage1;
    protected final Image originalImage2;
    
    public final int image1OriginalWidth;
    public final int image1OriginalHeight;
    protected final int image2OriginalWidth;
    protected final int image2OriginalHeight;
    
    public CurvatureScaleSpaceInflectionMapper(ImageExt image1, ImageExt image2) {
        
        this.image1 = image1;
        this.image2 = image2;
        
        originalImage1 = (ImageExt)image1.copyImage();
        originalImage2 = (ImageExt)image2.copyImage();
        
        image1OriginalWidth = image1.getWidth();
        image1OriginalHeight = image1.getHeight();
        image2OriginalWidth = image2.getWidth();
        image2OriginalHeight = image2.getHeight();
        
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
        if (false && doRefineTransformations) {
            
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

    protected Image getImage1() {
        return image1;
    }

    protected Image getImage2() {
        return image2;
    }

    Image getOriginalImage1() {
        return originalImage1;
    }

    Image getOriginalImage2() {
        return originalImage2;
    }

    @Override
    protected void createEdges1() {
        
        // note that if the orientation of image2 with respect to image1 is
        // more than 180 degrees, the code in the edge extractor will be forming
        // the edges by reading in the reverse direction in x and y,
        // so the edges will have opposite orientation.
        // The code also reverses edges to append curves read in other directions,
        // so in general, one cannot assure that the points are ordered in
        // a clockwise or counterclockwise manner at this point.
        // the curves tend to be ordered counter clockwise.
        // because the closed curve shapes are not simple convex or concave
        // shapes sometimes, it's not as easy to tell whether points are
        // ordered clockwise in a curve.
        // Tests so far show that testing the contour peak points
        // (left and right) for clockwise order gives the right answer
        // and is less work computationally
        CurvatureScaleSpaceImageMaker imgMaker = new CurvatureScaleSpaceImageMaker(image1);
        if (useLineDrawingMode) {
            imgMaker.useLineDrawingMode();
        }
        if (useOutdoorMode) {
            imgMaker.useOutdoorMode();
        }
        imgMaker.initialize();
        
        edges1 = getEdges(imgMaker);
        offsetImageX1 = imgMaker.getTrimmedXOffset();
        offsetImageY1 = imgMaker.getTrimmedYOffset();
                
    }

    @Override
    protected void createEdges2() {
        
        CurvatureScaleSpaceImageMaker imgMaker = new CurvatureScaleSpaceImageMaker(image2);
        if (useLineDrawingMode) {
            imgMaker.useLineDrawingMode();
        }
        if (useOutdoorMode) {
            imgMaker.useOutdoorMode();
        }
        imgMaker.initialize();
        
        edges2 = getEdges(imgMaker);
        offsetImageX2 = imgMaker.getTrimmedXOffset();
        offsetImageY2 = imgMaker.getTrimmedYOffset();
        
    }

}
