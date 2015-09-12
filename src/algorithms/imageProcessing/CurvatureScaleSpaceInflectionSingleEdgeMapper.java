package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * class to map contours from edge1 from one image to another and return the
 * matched inflection points and a transformation matrix that can
 * be applied to the first to put it into the frame of the second.
 * 
 * The algorithm used for matching scale space image contours is documented in 
 * CSSContourMatcherWrapper
 * @see algorithms.imageProcessing.CSSContourMatcherWrapper
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionSingleEdgeMapper {

    private final int xRelativeOffset1; 
    private final int yRelativeOffset1;
    private final int xRelativeOffset2; 
    private final int yRelativeOffset2;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private CSSContourMatcherWrapper matcherResults = null;
        
    public CurvatureScaleSpaceInflectionSingleEdgeMapper(
        int offsetImageX1, int offsetImageY1, int offsetImageX2, int offsetImageY2) {
        
        this.xRelativeOffset1 = offsetImageX1;
        this.yRelativeOffset1 = offsetImageY1;
        this.xRelativeOffset2 = offsetImageX2;
        this.yRelativeOffset2 = offsetImageY2;
    }
    
    public TransformationParameters matchContours(
        List<CurvatureScaleSpaceContour> contours1,
        List<CurvatureScaleSpaceContour> contours2) {
        
        boolean alreadySorted = true;
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        
        CSSContourMatcherWrapper matcher = new CSSContourMatcherWrapper(
            contours1, contours2, alreadySorted);
        
        boolean didMatch = matcher.matchContours();
        
        if (!didMatch) {
            return null;
        }
        
        this.matcherResults = matcher;
                
        List<CurvatureScaleSpaceContour> m1 = matcher.getSolutionMatchedContours1();
        List<CurvatureScaleSpaceContour> m2 = matcher.getSolutionMatchedContours2();
        if (m1 == null || m2 == null || m1.isEmpty() || m2.isEmpty()) {
            return null;
        }
        assert(m1.size() == m2.size());
                                
        /*
        There may be insignificant low cost matches for very small
        curves, so will only keep a solution when there are as few
        as 2 contours in the match if there are no other matches.
        */
        if (m1.size() == 2) {
        //    return null;
        }

        double cost = matcher.getSolvedCost();
        
        double scale = matcher.getSolvedScale();
        
        int centroidX1 = 0;
        int centroidY1 = 0;
        int centroidX2 = 0;
        int centroidY2 = 0;
                
        TransformationParameters params = null;
                                   
        PairIntArray xy1 = new PairIntArray(m1.size());
        PairIntArray xy2 = new PairIntArray(m2.size());

        List<Float> weights1 = new ArrayList<Float>();
        List<Float> weights2 = new ArrayList<Float>();

        //xy1 and xy2 have the image offsets added
        extract(m1, xy1, weights1, xRelativeOffset1, yRelativeOffset1);
        extract(m2, xy2, weights2, xRelativeOffset2, yRelativeOffset2);
            
        // if scale < 1, we have to swap the order of datasets to avoid
        // numerical errors in some of the methods that are the result of
        // dividing by a small number
        boolean reverseDatasetOrder = scale < 1.0;
        if (reverseDatasetOrder) {
            params = tc.calulateEuclideanGivenScale(1. / scale, 
                xy2, xy1, centroidX2, centroidY2);
        } else {
            params = tc.calulateEuclideanGivenScale(scale, 
                xy1, xy2, centroidX1, centroidY1);
        }
        if (reverseDatasetOrder && (params != null)) {
            params = tc.swapReferenceFrames(params);            
        }
        
        return params;
    }

    public static ScaleSpaceCurveImage createScaleSpaceImage(PairIntArray edge, 
        int edgeIndex) {
        
        CurvatureScaleSpaceCurvesMaker csscMaker = new CurvatureScaleSpaceCurvesMaker();
        
        // if use 2^(1/8) as a sigma factor should result in an error less than 10%
        // in determing the peak of a contour.  smaller factors have smaller
        // errors than that.
        float factor = (float)Math.pow(2, 1./32.);
        
        Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
            csscMaker.createScaleSpaceMetricsForEdge(edge, factor,
            SIGMA.ONE, SIGMA.TWOHUNDREDANDFIFTYSIX);
                       
        ScaleSpaceCurveImage scaleSpaceImage = 
            csscMaker.convertScaleSpaceMapToSparseImage(
            scaleSpaceMap, edgeIndex, edge.getN());
                       
 try {
 String fileSuffix = "edge_" + edgeIndex + "_" + MiscDebug.getCurrentTimeFormatted();
 MiscDebug.printScaleSpaceCurve(scaleSpaceImage, fileSuffix);
 int z = 1;
 } catch (IOException ex) {
 }

        return scaleSpaceImage;
    }
    
    public static List<CurvatureScaleSpaceContour> populateContours(
        ScaleSpaceCurveImage scaleSpaceImage, int edgeIndex, 
        boolean setToExtractWeakCurvesTooIfNeeded) {
        
 try {
 String fileSuffix = "edge_" + edgeIndex + "_" + MiscDebug.getCurrentTimeFormatted();
 MiscDebug.printScaleSpaceCurve(scaleSpaceImage, fileSuffix);
 int z = 1;
 } catch (IOException ex) {
 }
            
        ContourFinder contourFinder = new ContourFinder();
        
        if (setToExtractWeakCurvesTooIfNeeded) {
            contourFinder.overrideTheLowSigmaLimit(1.25f);
        }

        List<CurvatureScaleSpaceContour> result = 
            contourFinder.findContours(scaleSpaceImage, edgeIndex);

        boolean reversed = contourFinder.reverseIfClockwise(result);

        if (reversed) {
            //log.info("EDGES1: contour isCW=true");

            // these are extracted from contourFinder in order of decreasing
            // sigma already, so only need to be sorted if the list was
            // reversed
            Collections.sort(result, new DescendingSigmaComparator());
        }
        
        correctPeaks(result);
        
        /*
        for curves formed via blob boundaries rather than canny edge detector,
        can see a zig zag structure for the points near the top of a sigma=5 to 7
        scale space curve.  this next corrects for some of that.
        */
        removeRedundant(result);

        return result;
    }

    /**
     * when peak details has more than one point, this averages them and
     * replaces the details with a single point.  Note, that using the
     * peak detail to estimate a local orientation for a descriptor subsequently
     * needs to account for this when can see that the average does not
     * equal the curve point at the detail idx.
     * 
     * <pre>
     * For example, for a section of the curve with indexes:
     *  7  8  9  10  11  
     * if idx=9, 10 are the peak details that have been averaged to one value 
     * reported at idx=9, the surrounding points used for determining the
     * tangent to the peak should be idx=8 and idx=11.
      * </pre>
     * @param contours 
     */
    protected static void correctPeaks(List<CurvatureScaleSpaceContour> contours) {
        
        if (contours == null) {
            throw new IllegalArgumentException("contours cannot be null");
        }
        
        if (contours.isEmpty()) {
            return;
        }
        
        // the contours extracted from scale space images using a factor of
        // 2^(1/8) for recursive convolution tend to not have a single peak,
        // so the correction here for the single peak case is not usually
        // needed.  for that rare case, the avg of the other peak is stored
        // instead of both points
        
        for (int i = 0; i < contours.size(); i++) {
            
            CurvatureScaleSpaceContour c1 = contours.get(i);
            
            if (c1.getPeakDetails().length > 1) {                
                CurvatureScaleSpaceImagePoint p0 = c1.getPeakDetails()[0];
                CurvatureScaleSpaceImagePoint p1 = c1.getPeakDetails()[1];
                float t = p0.getScaleFreeLength();
                float s = p0.getSigma();
                int xAvg = Math.round((p0.getXCoord() + p1.getXCoord()) / 2.f);
                int yAvg = Math.round((p0.getYCoord() + p1.getYCoord()) / 2.f);
                CurvatureScaleSpaceImagePoint pAvg =
                    new CurvatureScaleSpaceImagePoint(s, t, xAvg, yAvg,
                    p0.getCoordIdx());
                CurvatureScaleSpaceImagePoint[] p =
                    new CurvatureScaleSpaceImagePoint[]{pAvg};
                c1.setPeakDetails(p);
                contours.set(i, c1);
            }
        }        
    }
    
    private void extract(List<CurvatureScaleSpaceContour> contours, 
        PairIntArray outputXY, List<Float> outputSigmaWeights,
        final int imageOffsetX, final int imageOffsetY) {
        
        float sumSigma = 0;
        
        for (int i = 0; i < contours.size(); i++) {
    
            CurvatureScaleSpaceContour c = contours.get(i);
        
            for (int j = 0; j < c.getPeakDetails().length; j++) {
                
                CurvatureScaleSpaceImagePoint spaceImagePoint = 
                    c.getPeakDetails()[j];
                
                int x = spaceImagePoint.getXCoord() + imageOffsetX;
                int y = spaceImagePoint.getYCoord() + imageOffsetY;
                
                outputXY.add(x, y);
                outputSigmaWeights.add(Float.valueOf(c.getPeakSigma()));
                
                sumSigma += c.getPeakSigma();
            }
        }
        
        
        for (int i = 0; i < outputSigmaWeights.size(); ++i) {
            float w = outputSigmaWeights.get(i)/sumSigma;
            outputSigmaWeights.set(i, w);
        }
    }
    
    private static void removeRedundant(List<CurvatureScaleSpaceContour> contours) {
        
        Set<Integer> indexes = new HashSet<Integer>();
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < contours.size(); ++i) {
            CurvatureScaleSpaceContour contour = contours.get(i);
            for (CurvatureScaleSpaceImagePoint ip : contour.getPeakDetails()) {
                Integer idx = Integer.valueOf(ip.getCoordIdx());
                if (indexes.contains(idx)) {
                    remove.add(i);
                } else {
                    indexes.add(idx);
                }
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i);
            contours.remove(idx);
        }
    }
    
    public CSSContourMatcherWrapper getMatcher() {
        return matcherResults;
    }
}
