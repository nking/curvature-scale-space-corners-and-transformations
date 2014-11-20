package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * class to map contours from one image to another and return the
 * matched inflection points and a transformation matrix that can
 * be applied to the first to put it into the frame of the second.
 * 
 * The algorithm used for matching scale space image contours is documented in 
 * CurvatureScaleSpaceContourMatcher
 * @see algorithms.imageProcessing.CurvatureScaleSpaceContourMatcher
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionMapper {
    
    private final Logger log = Logger.getLogger(this.getClass().getName());
    
    private boolean debug = false;
    
    private boolean useLineDrawingMode = false;
        
    private boolean doNotRefineTransformations = false;
    
    private boolean initialized = false;
    
    private final GreyscaleImage image1; 
    private final GreyscaleImage image2;
    
    private List<PairIntArray> edges1 = null;
    private List<PairIntArray> edges2 = null;
    
    private List<CurvatureScaleSpaceContour> contours1 = new 
        ArrayList<CurvatureScaleSpaceContour>();
    
    private List<CurvatureScaleSpaceContour> contours2 = new 
        ArrayList<CurvatureScaleSpaceContour>();
        
    private List<CurvatureScaleSpaceContour> matchedContours1 = new 
        ArrayList<CurvatureScaleSpaceContour>();
    
    private List<CurvatureScaleSpaceContour> matchedContours2 = new 
        ArrayList<CurvatureScaleSpaceContour>();
    
    private PairIntArray[] matchedXY1 = null;
    
    private PairIntArray[] matchedXY2 = null;
    
    private int offsetImageX1 = 0;
    
    private int offsetImageY1 = 0;
    
    private int offsetImageX2 = 0;
    
    private int offsetImageY2 = 0;
    
    private boolean useOutdoorMode = false;
        
    public void useOutdoorMode() {
        useOutdoorMode = true;
    }
    
    public void useLineDrawingLineMode() {
        this.useLineDrawingMode = true;
    }
    
    public void doNotRefineTransformations() {
        doNotRefineTransformations = true;
    }
    
    public void useDebugMode() {
        debug = true;
    }
    
    public void initialize() {
        
        if (initialized) {
            return;
        }
        
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
        
        CurvatureScaleSpaceImageMaker imgMaker = new 
            CurvatureScaleSpaceImageMaker(image1);
            
        if (useLineDrawingMode) {
            imgMaker.useLineDrawingMode();
        }
        
        if (useOutdoorMode) {
            imgMaker.useOutdoorMode();
        }
        
        imgMaker.initialize();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        edges1 = imgMaker.getClosedCurves();
        
        boolean didReverse = false;
        
        for (int i = 0; i < edges1.size(); i++) {

            PairIntArray curve = edges1.get(i);

            Map<Float, ScaleSpaceCurve> scaleSpaceMap
                = imgMaker.createScaleSpaceMetricsForEdge2(curve);

            ScaleSpaceCurveImage scaleSpaceImage
                = imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i);
                 
            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
                scaleSpaceImage, i);
                    
            PairIntArray testContour = new PairIntArray();
            for (int j = 0; j < result.size(); j++) {
                CurvatureScaleSpaceContour c = result.get(j);
                CurvatureScaleSpaceImagePoint[] points = c.getPeakDetails();
                for (int jj = 0; jj < points.length; jj++) {
                    testContour.add((int)points[jj].getXCoord(), 
                        (int)points[jj].getYCoord());
                }
            }
            
            boolean isCW = curveHelper.curveIsOrderedClockwise(testContour);
            log.fine("EDGES1: contour isCW=" + isCW);
            
            if (isCW) {
                
                didReverse = true;
                
                for (int j = 0; j < result.size(); j++) {
                    
                    CurvatureScaleSpaceContour contour = result.get(j);
                                        
                    CurvatureScaleSpaceContour reversed = 
                        new CurvatureScaleSpaceContour(contour.getPeakSigma(), 
                        1.0f - contour.getPeakScaleFreeLength());
                    
                    CurvatureScaleSpaceImagePoint[] points = 
                        contour.getPeakDetails();
                    if (points.length > 1) {
                        CurvatureScaleSpaceImagePoint tmp = points[0];
                        points[0] = points[1];
                        points[1] = tmp;
                    }
                    for (int jj = 0; jj < points.length; jj++) {
                        points[jj].setScaleFreeLength(1.0f - 
                            points[jj].getScaleFreeLength());
                    }
                    reversed.setPeakDetails(points);
                    reversed.setEdgeNumber(contour.getEdgeNumber());
                                        
                    result.set(j, reversed);
                }
            }
            
            contours1.addAll(result);
        }
        
        if (contours1.isEmpty()) {
            log.info("no contours found in image 1");
            return;
        }
        
        if ((edges1.size() > 1) || didReverse) {
            Collections.sort(contours1, new DescendingSigmaComparator());
        }
        
        float highestPeak1 = contours1.get(0).getPeakSigma();
        
        float lowThresh1 = 0.15f * highestPeak1;
        
        for (int i = (contours1.size() - 1); i > -1; i--) {
            if (contours1.get(i).getPeakSigma() < lowThresh1) {
                contours1.remove(i);
            }
        }
        
        /*
        note that when modifying the contour lists in any way, one has to
        maintain decreasing order by sigma and when sigma is equal, the
        order must be by increasing scale free parameter.
        two of the search methods in the matcher depend upon those properties.
        */
        
        offsetImageX1 = imgMaker.getTrimmedXOffset();
        
        offsetImageY1 = imgMaker.getTrimmedYOffset();
                    
        imgMaker = new CurvatureScaleSpaceImageMaker(image2);
        
        if (useLineDrawingMode) {
            imgMaker.useLineDrawingMode();
        }
        
        if (useOutdoorMode) {
            imgMaker.useOutdoorMode();
        }
        
        imgMaker.initialize();
        
        edges2 = imgMaker.getClosedCurves();
        
        didReverse = false;
        
        for (int i = 0; i < edges2.size(); i++) {

            PairIntArray curve = edges2.get(i);
            
            Map<Float, ScaleSpaceCurve> scaleSpaceMap
                = imgMaker.createScaleSpaceMetricsForEdge2(curve);

            ScaleSpaceCurveImage scaleSpaceImage
                = imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i);
                 
            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
                scaleSpaceImage, i);
            
            PairIntArray testContour = new PairIntArray();
            for (int j = 0; j < result.size(); j++) {
                CurvatureScaleSpaceContour c = result.get(j);
                CurvatureScaleSpaceImagePoint[] points = c.getPeakDetails();
                for (int jj = 0; jj < points.length; jj++) {
                    testContour.add((int)points[jj].getXCoord(), 
                        (int)points[jj].getYCoord());
                }
            }
            
            boolean isCW = curveHelper.curveIsOrderedClockwise(testContour);
            log.fine("EDGES2: contour isCW=" + isCW);
            
            if (isCW) {
                
                didReverse = true;
                
                for (int j = 0; j < result.size(); j++) {
                    
                    CurvatureScaleSpaceContour contour = result.get(j);
                                        
                    CurvatureScaleSpaceContour reversed = 
                        new CurvatureScaleSpaceContour(contour.getPeakSigma(), 
                        1.0f - contour.getPeakScaleFreeLength());
                    
                    CurvatureScaleSpaceImagePoint[] points = 
                        contour.getPeakDetails();
                    if (points.length > 1) {
                        CurvatureScaleSpaceImagePoint tmp = points[0];
                        points[0] = points[1];
                        points[1] = tmp;
                    }
                    for (int jj = 0; jj < points.length; jj++) {
                        points[jj].setScaleFreeLength(1.0f - 
                            points[jj].getScaleFreeLength());
                    }
                    reversed.setPeakDetails(points);
                    reversed.setEdgeNumber(contour.getEdgeNumber());
                                        
                    result.set(j, reversed);
                }
            }
            
            contours2.addAll(result);
        }
          
        offsetImageX2 = imgMaker.getTrimmedXOffset();
        
        offsetImageY2 = imgMaker.getTrimmedYOffset();
        
        if (contours2.isEmpty()) {
            log.info("did not find contours in image 2");
            return;
        }
        
        if ((edges2.size() > 1) || didReverse) {
            Collections.sort(contours2, new DescendingSigmaComparator());            
        }
        
        float highestPeak2 = contours2.get(0).getPeakSigma();
        
        float lowThresh2 = 0.15f * highestPeak2;
        
        for (int i = (contours2.size() - 1); i > -1; i--) {
            if (contours2.get(i).getPeakSigma() < lowThresh2) {
                contours2.remove(i);
            }
        }
      
        initialized = true;
    }
    
    public CurvatureScaleSpaceInflectionMapper(GreyscaleImage image1, 
        GreyscaleImage image2) {
        
        this.image1 = image1;
        
        this.image2 = image2;
        
    }
  
    /**
     * coordinate transformations from image 1 to image 2 are calculated from
     * matching scale space image contours.
     *
     * positive Y is down 
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |     
                 |
          180--------- 0   +X
                 |   
                 |   
                 90
                 +Y
     * </pre>
     * NOTE: this will return null if it did not find closed 
     * curves in each image for which to map between.
     * 
     * @return 
     */
    public TransformationParameters createEuclideanTransformation() {
        
        initialize();
        
        if (contours2.isEmpty() || contours1.isEmpty()) {
            return null;
        }
        
        CurvatureScaleSpaceContourMatcher matcher = 
            new CurvatureScaleSpaceContourMatcher(contours1, contours2);
        
        matcher.matchContours();
        
        List<CurvatureScaleSpaceContour> transAppliedTo1 = 
            matcher.getSolutionMatchedContours1();
        
        List<CurvatureScaleSpaceContour> transAppliedTo2 = 
            matcher.getSolutionMatchedContours2();
        
        if (transAppliedTo1 == null || transAppliedTo2 == null) {
            return null;
        }
        if (transAppliedTo1.size() != transAppliedTo2.size()) {
            throw new IllegalStateException(
                "contour matcher should have same number of contours in both lists");
        }
                
        matchedContours1.addAll(transAppliedTo1);
        
        matchedContours2.addAll(transAppliedTo2);
                
        log.info("Contour matcher solution scale=" + matcher.getSolvedScale());
        
        log.info("Contour matcher solution shift=" + matcher.getSolvedShift());
        
        // ==== adding back the image offsets removed when trimming image
        
        // ======= make a weighted sum of points to get the centroid of the edge
        // ============= the weight is sigma of total sum of sigmas
        
        
        // for each edge, we calculate the transformation from its contours
        // because the contours match, they must be from edges that match
        List<Integer> edgeNumbers1 = new ArrayList<Integer>();
        List<Integer> edgeNumbers2 = new ArrayList<Integer>();
        for (int i = 0; i < transAppliedTo1.size(); i++) {
            CurvatureScaleSpaceContour c1 = transAppliedTo1.get(i);            
            Integer e1 = Integer.valueOf(c1.getEdgeNumber());
            if (!edgeNumbers1.contains(e1)) {
                edgeNumbers1.add(e1);
                CurvatureScaleSpaceContour c2 = transAppliedTo2.get(i);
                Integer e2 = Integer.valueOf(c2.getEdgeNumber());
                edgeNumbers2.add(e2);
            }
        }
                
        PairIntArray[] xy1 = new PairIntArray[edgeNumbers1.size()];
        PairIntArray[] xy2 = new PairIntArray[edgeNumbers1.size()];
        
        List< List<Float> > weights1 = new ArrayList<List<Float> >();
        List< List<Float> > weights2 = new ArrayList<List<Float> >();
            
        TransformationParameters[] params = new 
            TransformationParameters[edgeNumbers1.size()];
        
        double[] sumS1 = new double[xy1.length];
        double[] sumS2 = new double[xy1.length];
            
        for (int i = 0; i < xy1.length; i++) {
            xy1[i] = new PairIntArray();
            xy2[i] = new PairIntArray();
            params[i] = new TransformationParameters();
            weights1.add(new ArrayList<Float>());
            weights2.add(new ArrayList<Float>());
        }
            
        for (int i = 0; i < transAppliedTo1.size(); i++) {
            
            CurvatureScaleSpaceContour c1 = transAppliedTo1.get(i);
            
            CurvatureScaleSpaceContour c2 = transAppliedTo2.get(i);
            
            Integer edge1Number = Integer.valueOf(c1.getEdgeNumber());
            int edgeIdx = edgeNumbers1.indexOf(edge1Number);
            
            float sigma1 = c1.getPeakSigma();
            float sigma2 = c2.getPeakSigma();

            StringBuilder s1 = new StringBuilder();
            StringBuilder s2 = new StringBuilder();

            if (debug) {
                s1.append(String.format("CONTOUR PEAK1: (%f, %f)", 
                    c1.getPeakSigma(), c1.getPeakScaleFreeLength()));
                s2.append(String.format("CONTOUR PEAK2: (%f, %f)", 
                    c2.getPeakSigma(), c2.getPeakScaleFreeLength()));
            }
            
            // the contours extracted from scale space images using a factor of
            // 2^(1/8) for recursive convolution tend to not have a single
            // peak, so the correction here for the single peak case is not 
            // usually needed.  for that rare case, just doubling the peak
            // for comparison to the matched contour
            if (c1.getPeakDetails().length != c2.getPeakDetails().length) {
                if (c1.getPeakDetails().length == 1) {
                    CurvatureScaleSpaceImagePoint[] p = new
                        CurvatureScaleSpaceImagePoint[]{
                        c1.getPeakDetails()[0], c1.getPeakDetails()[0]};
                    c1.setPeakDetails(p);
                    transAppliedTo1.set(i, c1);
                } else if (c2.getPeakDetails().length == 1) {
                    CurvatureScaleSpaceImagePoint[] p = new
                        CurvatureScaleSpaceImagePoint[]{
                        c2.getPeakDetails()[0], c2.getPeakDetails()[0]};                
                    c2.setPeakDetails(p);
                    transAppliedTo2.set(i, c2);
                }
            }
            
            for (int j = 0; j < c1.getPeakDetails().length; j++) {

                CurvatureScaleSpaceImagePoint spaceImagePoint = 
                    c1.getPeakDetails()[j];

                int x = spaceImagePoint.getXCoord() + offsetImageX1;
                int y = spaceImagePoint.getYCoord() + offsetImageY1;
                xy1[edgeIdx].add(x, y);
                weights1.get(edgeIdx).add(Float.valueOf(sigma1));
                sumS1[edgeIdx] += sigma1;
                   
                if (debug) {
                    s1.append(String.format(" (%d, %d)", x, y));
                }

                spaceImagePoint = c2.getPeakDetails()[j];

                x = spaceImagePoint.getXCoord() + offsetImageX2;
                y = spaceImagePoint.getYCoord() + offsetImageY2;
                xy2[edgeIdx].add(x, y);
                weights2.get(edgeIdx).add(Float.valueOf(sigma2));
                sumS2[edgeIdx] += sigma2;
                    
                if (debug) {
                    s2.append(String.format(" (%d, %d)", x, y));
                }
            }
            
            if (debug) {
                log.info(s1.toString());
                log.info(s2.toString());
            }            
        }
        
        float[][] w1 = new float[xy1.length][];
        float[][] w2 = new float[xy1.length][];
        for (int i = 0; i < xy1.length; i++) {
            List<Float> lw1 = weights1.get(i);
            List<Float> lw2 = weights1.get(i);
            w1[i] = new float[lw1.size()];
            w2[i] = new float[lw2.size()];
            for (int ii = 0; ii < lw1.size(); ii++) {
                double tmp = lw1.get(ii).floatValue() / sumS1[i];
                w1[i][ii] = Float.valueOf((float) tmp);
            }
            for (int ii = 0; ii < lw2.size(); ii++) {
                double tmp = lw2.get(ii).floatValue() / sumS2[i];
                w2[i][ii] = Float.valueOf((float) tmp);
            }
        }
              
        if (debug) {
            log.info("offsetImgX1=" + offsetImageX1 
                + " offsetImgY1=" + offsetImageY1
                + "\noffsetImgX2=" + offsetImageX2 
                + " offsetImgY2=" + offsetImageY2
            );
        }
            
        matchedXY1 = xy1;
        matchedXY2 = xy2;
        
        for (int i = 0; i < xy1.length; i++) {
            
            if (xy1[i].getN() < 3) {
                throw new IllegalStateException("need at least 3 points");
            }
            
            TransformationCalculator tc = new TransformationCalculator();
            if (debug) {
                tc.useDebugMode();
            }
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            double[] centroidsX1 = curveHelper.calculateXYCentroids(xy1[i]);
            double centroidX1 = centroidsX1[0];
            double centroidY1 = centroidsX1[1];
            double[] centroidsX2 = curveHelper.calculateXYCentroids(xy2[i]);
            double centroidX2 = centroidsX2[0];
            double centroidY2 = centroidsX2[1];
            
            params[i] = tc.calulateEuclidean(xy1[i], w1[i], xy2[i], w2[i],
                centroidX1, centroidY1, centroidX2, centroidY2);
            
            if (!doNotRefineTransformations) {
            
                params[i] = refineEuclideanSolution(params[i]);

                double mc = (float) Math.cos(params[i].getRotationInRadians());
                double ms = (float) Math.sin(params[i].getRotationInRadians());
                double scale = params[i].getScale();

                float translationX = (float) (centroidX2 - (centroidX1 * scale * mc)
                    - (centroidY1 * scale * ms));
                float translationY = (float) (centroidY2 + (centroidX1 * scale * ms)
                    - (centroidY1 * scale * mc));

                params[i].setTranslationX(translationX);
                params[i].setTranslationY(translationY);

                log.info("FINAL:\n" + params.toString());
            }
        }
        
        if (params.length == 0) {
            return null;
        } else if (params.length == 1) {
            return params[0];
        }
        
        // make a weighted average where weight is by
        // number of points in each solution or by strength of contour?
        TransformationParameters avgWeightedParams = new 
            TransformationParameters();
        
        double nTotPoints = 0;
        for (int i = 0; i < xy1.length; i++) {
            nTotPoints += xy1[i].getN();
        }
        double[] w = new double[xy1.length];
        for (int i = 0; i < xy1.length; i++) {
            w[i] = (double)xy1[i].getN()/nTotPoints;
        }
        for (int i = 0; i < xy1.length; i++) {
            float rot = params[i].getRotationInRadians();
            float scale = params[i].getScale();
            float transX = params[i].getTranslationX();
            float transY = params[i].getTranslationY();
            
            double rotA = (w[i]*avgWeightedParams.getRotationInRadians()) + rot;
            double scaleA = (w[i]*avgWeightedParams.getScale()) + scale;
            double transXA = (w[i]*avgWeightedParams.getTranslationX()) + transX;
            double transYA = (w[i]*avgWeightedParams.getTranslationY()) + transY;
            
            avgWeightedParams.setRotationInRadians((float)rotA);
            avgWeightedParams.setScale((float)scaleA);
            avgWeightedParams.setTranslationX((float)transXA);
            avgWeightedParams.setTranslationY((float)transYA);
        }
        
        return avgWeightedParams;
    }
    
    public PairIntArray[] getMatchedXY1() {
        return matchedXY1;
    }
    
    public PairIntArray[] getMatchedXY2() {
        return matchedXY2;
    }

    public List<CurvatureScaleSpaceContour> getMatchedContours1() {
        return matchedContours1;
    }
    
    public List<CurvatureScaleSpaceContour> getMatchedContours2() {
        return matchedContours2;
    }
    List<CurvatureScaleSpaceContour> getContours1() {
        return contours1;
    }
    
    List<CurvatureScaleSpaceContour> getContours2() {
        return contours2;
    }
    
    protected GreyscaleImage getImage1() {
        return image1;
    }
    protected GreyscaleImage getImage2() {
        return image2;
    }
    protected List<PairIntArray> getEdges1() {
        return edges1;
    }
    protected List<PairIntArray> getEdges2() {
        return edges2;
    }
    
    protected List<PairIntArray> getEdges1InOriginalReferenceFrame() {
        List<PairIntArray> oe = new ArrayList<PairIntArray>();
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray edge = edges1.get(i);
            PairIntArray e = new PairIntArray();
            for (int j = 0; j < edge.getN(); j++) {
                e.add(edge.getX(j) + offsetImageX1, edge.getY(j) + offsetImageY1);
            }
            oe.add(e);
        } 
        return oe;
    }
    protected List<PairIntArray> getEdges2InOriginalReferenceFrame() {
        List<PairIntArray> oe = new ArrayList<PairIntArray>();
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray edge = edges2.get(i);
            PairIntArray e = new PairIntArray();
            for (int j = 0; j < edge.getN(); j++) {
                e.add(edge.getX(j) + offsetImageX2, edge.getY(j) + offsetImageY2);
            }
            oe.add(e);
        } 
        return oe;
    }
    
    public PairInt[] getMatchedEdgesIndexes() {
         
        List<Integer> idx1 = new ArrayList<Integer>();
        List<Integer> idx2 = new ArrayList<Integer>();
        for (int i = 0; i < this.matchedContours1.size(); i++) {
            Integer edge1Idx = Integer.valueOf(
                matchedContours1.get(i).getEdgeNumber());
            Integer edge2Idx = Integer.valueOf(
                matchedContours2.get(i).getEdgeNumber());
            
            if (!idx1.contains(edge1Idx)) {
                idx1.add(edge1Idx);
                idx2.add(edge2Idx);
            }
        }
        
        PairInt[] indexes = new PairInt[idx1.size()];
        for (int i = 0; i < idx1.size(); i++) {
            indexes[i] = new PairInt(idx1.get(i), idx2.get(i));
        }
        
        return indexes;
    }

    private TransformationParameters refineEuclideanSolution(
        TransformationParameters params) {
        
        // transform edges in matchedContours1 to matchedContours2
        // and reduce the difference within max num of iterations and 
        // tolerance
       
        PairInt[] matchedEdgesIndexes = getMatchedEdgesIndexes();
        
        PairIntArray[] matchedEdges1 = new PairIntArray[matchedEdgesIndexes.length];
        PairIntArray[] matchedEdges2 = new PairIntArray[matchedEdges1.length];
        for (int i = 0; i < matchedEdgesIndexes.length; i++) {
            int edge1Idx = matchedEdgesIndexes[i].getX();
            int edge2Idx = matchedEdgesIndexes[i].getY();
            matchedEdges1[i] = edges1.get(edge1Idx);
            matchedEdges2[i] = edges2.get(edge2Idx);
        }
        
        //TODO: consider impl a non-linear conjugate gradient search to replace
        // this:
        TransformationRefiner chSqMin = new TransformationRefiner();
        
        TransformationParameters out = chSqMin.refineTransformation(matchedEdges1, 
            matchedEdges2, params);
        
        //missing translation fields are set in invoker which has necessary
        //information
        
        return out;
    }
    
}
