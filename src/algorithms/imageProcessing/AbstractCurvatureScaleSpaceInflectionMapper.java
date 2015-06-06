package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractCurvatureScaleSpaceInflectionMapper implements 
    ICurvatureScaleSpaceInflectionMapper {
    
    protected final Logger log = Logger.getLogger(this.getClass().getName());
    protected boolean debug = false;
    protected boolean useLineDrawingMode = false;
    protected boolean doRefineTransformations = false;
    protected boolean initialized = false;
    protected final ImageExt image1;
    protected final ImageExt image2;
    // for debugging, keeping a reference of originals
    protected final Image originalImage1;
    protected final Image originalImage2;
    public final int image1OriginalWidth;
    public final int image1OriginalHeight;
    protected final int image2OriginalWidth;
    protected final int image2OriginalHeight;
    protected List<PairIntArray> edges1 = null;
    protected List<PairIntArray> edges2 = null;
    protected List<CurvatureScaleSpaceContour> contours1 = new ArrayList<CurvatureScaleSpaceContour>();
    protected List<CurvatureScaleSpaceContour> contours2 = new ArrayList<CurvatureScaleSpaceContour>();
    protected Map<Integer, PairIntArray> matchedXY1ByEdgeInOrigRefFrame = null;
    protected Map<Integer, PairIntArray> matchedXY2ByEdgeInOrigRefFrame = null;
    protected Map<Integer, List<Float>> matchedXY1ByEdgeWeights = null;
    protected Map<Integer, List<Float>> matchedXY2ByEdgeWeights = null;
    /**
     * matched points from the contour lists of image 1 (matched to the same
     * in image 2) with coordinates being in the reference frames of the
     * original image 1 before any trimming.
     */
    protected PairIntArray matchedXY1 = null;
    /**
     * matched points from the contour lists of image 2 (matched to the same
     * in image 1) with coordinates being in the reference frames of the
     * original image 2 before any trimming.
     */
    protected PairIntArray matchedXY2 = null;
    /**
     * weights for points in matchedXY1 created from the peak strengths.
     */
    protected float[] matchedXY1Weights = null;
    /**
     * weights for points in matchedXY2 created from the peak strengths.
     */
    protected float[] matchedXY2Weights = null;
    /**
     * indexes for edges from edges1 which produced matching contours
     */
    protected int[] matchedEdge1Indexes = null;
    /**
     * indexes for edges from edges2 which produced matching contours
     */
    protected int[] matchedEdge2Indexes = null;
    protected int offsetImageX1 = 0;
    protected int offsetImageY1 = 0;
    protected int offsetImageX2 = 0;
    protected int offsetImageY2 = 0;
    protected boolean useOutdoorMode = false;

    public AbstractCurvatureScaleSpaceInflectionMapper(ImageExt image1, 
        ImageExt image2) {
                
        this.image1 = image1;
        this.image2 = image2;
        
        originalImage1 = (ImageExt)image1.copyImage();
        originalImage2 = (ImageExt)image2.copyImage();
        
        image1OriginalWidth = image1.getWidth();
        image1OriginalHeight = image1.getHeight();
        image2OriginalWidth = image2.getWidth();
        image2OriginalHeight = image2.getHeight();
    }

    @Override
    public void useOutdoorMode() {
        useOutdoorMode = true;
    }

    @Override
    public void useLineDrawingLineMode() {
        this.useLineDrawingMode = true;
    }

    @Override
    public void setToRefineTransformations() {
        doRefineTransformations = true;
    }

    @Override
    public void useDebugMode() {
        debug = true;
    }

    @Override
    public void initialize() {
        if (initialized) {
            return;
        }
        initialized = true;
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
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        edges1 = getEdges(imgMaker);
        
        boolean didReverse = false;
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray curve = edges1.get(i);
            Map<Float, ScaleSpaceCurve> scaleSpaceMap = imgMaker.createScaleSpaceMetricsForEdge2(curve);
            ScaleSpaceCurveImage scaleSpaceImage = 
                imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i,
                curve.getN());
            
            ContourFinder contourFinder = new ContourFinder();
            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(scaleSpaceImage, i);
            PairIntArray testContour = new PairIntArray();
            for (int j = 0; j < result.size(); j++) {
                CurvatureScaleSpaceContour c = result.get(j);
                CurvatureScaleSpaceImagePoint[] points = c.getPeakDetails();
                for (int jj = 0; jj < points.length; jj++) {
                    testContour.add(points[jj].getXCoord(), points[jj].getYCoord());
                }
            }
            boolean isCW = curveHelper.curveIsOrderedClockwise(testContour);
            log.info("EDGES1: contour isCW=" + isCW);
            if (isCW) {
                didReverse = true;
                for (int j = 0; j < result.size(); j++) {
                    CurvatureScaleSpaceContour contour = result.get(j);
                    CurvatureScaleSpaceContour reversed = 
                        new CurvatureScaleSpaceContour(contour.getPeakSigma(), 
                            1.0f - contour.getPeakScaleFreeLength());
                    CurvatureScaleSpaceImagePoint[] points = contour.getPeakDetails();
                    if (points.length > 1) {
                        CurvatureScaleSpaceImagePoint tmp = points[0];
                        points[0] = points[1];
                        points[1] = tmp;
                    }
                    for (int jj = 0; jj < points.length; jj++) {
                        points[jj].setScaleFreeLength(1.0f - points[jj].getScaleFreeLength());
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
        /*float highestPeak1 = contours1.get(0).getPeakSigma();
        float lowThresh1 = 0.15f * highestPeak1;
        for (int i = (contours1.size() - 1); i > -1; i--) {
        if (contours1.get(i).getPeakSigma() < lowThresh1) {
        contours1.remove(i);
        }
        }*/
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
        
        edges2 = getEdges(imgMaker);
        
        didReverse = false;
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray curve = edges2.get(i);
            Map<Float, ScaleSpaceCurve> scaleSpaceMap = imgMaker.createScaleSpaceMetricsForEdge2(curve);
            ScaleSpaceCurveImage scaleSpaceImage = 
                imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i,
                curve.getN());
            ContourFinder contourFinder = new ContourFinder();
            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(scaleSpaceImage, i);
            PairIntArray testContour = new PairIntArray();
            for (int j = 0; j < result.size(); j++) {
                CurvatureScaleSpaceContour c = result.get(j);
                CurvatureScaleSpaceImagePoint[] points = c.getPeakDetails();
                for (int jj = 0; jj < points.length; jj++) {
                    testContour.add(points[jj].getXCoord(), points[jj].getYCoord());
                }
            }
            boolean isCW = curveHelper.curveIsOrderedClockwise(testContour);
            log.info("EDGES2: contour isCW=" + isCW);
            if (isCW) {
                didReverse = true;
                for (int j = 0; j < result.size(); j++) {
                    CurvatureScaleSpaceContour contour = result.get(j);
                    CurvatureScaleSpaceContour reversed = 
                        new CurvatureScaleSpaceContour(contour.getPeakSigma(), 
                            1.0f - contour.getPeakScaleFreeLength());
                    CurvatureScaleSpaceImagePoint[] points = contour.getPeakDetails();
                    if (points.length > 1) {
                        CurvatureScaleSpaceImagePoint tmp = points[0];
                        points[0] = points[1];
                        points[1] = tmp;
                    }
                    for (int jj = 0; jj < points.length; jj++) {
                        points[jj].setScaleFreeLength(1.0f - points[jj].getScaleFreeLength());
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
        
debugPlot(contours1, image1, offsetImageX1, offsetImageY1);
debugPlot(contours2, image2, offsetImageX2, offsetImageY2);
       
        if (contours2.isEmpty()) {
            log.info("did not find contours in image 2");
            return;
        }
        if ((edges2.size() > 1) || didReverse) {
            Collections.sort(contours2, new DescendingSigmaComparator());
        }
        /*float highestPeak2 = contours2.get(0).getPeakSigma();
        float lowThresh2 = 0.15f * highestPeak2;
        for (int i = (contours2.size() - 1); i > -1; i--) {
        if (contours2.get(i).getPeakSigma() < lowThresh2) {
        contours2.remove(i);
        }
        }*/
    }

    protected abstract void createMatchedPointArraysFromContourPeaks();
    
    public abstract TransformationParameters createEuclideanTransformationImpl();
    
    @Override
    public TransformationParameters createEuclideanTransformation() {
        
        //TODO:  no need to check and reverse contours in the init stage
        //   so remove those
        
        initialize();
        
        if (contours2.isEmpty() || contours1.isEmpty()) {
            return null;
        }
        
        createMatchedPointArraysFromContourPeaks();
        
        return createEuclideanTransformationImpl();
    }
    
    @Override
    public PairIntArray getMatchedXY1() {
        return matchedXY1;
    }

    @Override
    public PairIntArray getMatchedXY2() {
        return matchedXY2;
    }

    @Override
    public float[] getMatchedXY1Weights() {
        return matchedXY1Weights;
    }

    @Override
    public float[] getMatchedXY2Weights() {
        return matchedXY2Weights;
    }

    @Override
    public List<CurvatureScaleSpaceContour> getContours1() {
        return contours1;
    }

    @Override
    public List<CurvatureScaleSpaceContour> getContours2() {
        return contours2;
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

    protected List<PairIntArray> getEdges1() {
        return edges1;
    }

    protected List<PairIntArray> getEdges2() {
        return edges2;
    }

    protected List<PairIntArray> getEdges1InOriginalReferenceFrame() {
        List<PairIntArray> oe = new ArrayList<PairIntArray>();
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray edge = edges1.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX1, 
                    edge.getY(j) + offsetImageY1);
            }
            oe.add(edge);
        }
        return oe;
    }

    protected PairIntArray[] getEdges1InOriginalReferenceFrameArray() {
        PairIntArray[] oe = new PairIntArray[edges1.size()];
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray edge = edges1.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX1, 
                    edge.getY(j) + offsetImageY1);
            }
            oe[i] = edge;
        }
        return oe;
    }

    protected List<PairIntArray> getEdges2InOriginalReferenceFrame() {
        List<PairIntArray> oe = new ArrayList<PairIntArray>();
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray edge = edges2.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX2, 
                    edge.getY(j) + offsetImageY2);
            }
            oe.add(edge);
        }
        return oe;
    }

    protected PairIntArray[] getEdges2InOriginalReferenceFrameArray() {
        PairIntArray[] oe = new PairIntArray[edges2.size()];
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray edge = edges2.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX2, 
                    edge.getY(j) + offsetImageY2);
            }
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getMatchedEdges1InOriginalReferenceFrameArray() {
        if (matchedEdge1Indexes == null) {
            return new PairIntArray[0];
        }
        PairIntArray[] oe = new PairIntArray[matchedEdge1Indexes.length];
        for (int i = 0; i < matchedEdge1Indexes.length; i++) {
            int eIdx = matchedEdge1Indexes[i];
            PairIntArray edge = edges1.get(eIdx).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX1, 
                    edge.getY(j) + offsetImageY1);
            }
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getMatchedEdges2InOriginalReferenceFrameArray() {
        if (matchedEdge2Indexes == null) {
            return new PairIntArray[0];
        }
        PairIntArray[] oe = new PairIntArray[matchedEdge2Indexes.length];
        for (int i = 0; i < matchedEdge2Indexes.length; i++) {
            int eIdx = matchedEdge2Indexes[i];
            PairIntArray edge = edges2.get(eIdx).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX2, 
                    edge.getY(j) + offsetImageY2);
            }
            oe[i] = edge;
        }
        return oe;
    }

    public abstract PairInt[] getMatchedEdgesIndexes();
    
    int getOffsetImageX1() {
        return offsetImageX1;
    }

    int getOffsetImageY1() {
        return offsetImageY1;
    }

    int getOffsetImageX2() {
        return offsetImageX2;
    }

    int getOffsetImageY2() {
        return offsetImageY2;
    }

    protected abstract List<PairIntArray> getEdges(
        CurvatureScaleSpaceImageMaker imgMaker);
   
    private void debugPlot(List<CurvatureScaleSpaceContour> result, ImageExt 
        img, int xOffset, int yOffset) {
        
        if (result.isEmpty()) {
            return;
        }
        
        int nExtraForDot = 1;
        int rClr = 255;
        int gClr = 0;
        int bClr = 0;
        
        for (int i = 0; i < result.size(); i++) {
            
            CurvatureScaleSpaceContour cssC = result.get(i);
            
            CurvatureScaleSpaceImagePoint[] peakDetails = cssC.getPeakDetails();
            
            for (CurvatureScaleSpaceImagePoint peakDetail : peakDetails) {
                int x = peakDetail.getXCoord() + xOffset;
                int y = peakDetail.getYCoord() + yOffset;
                for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    float xx = x + dx;
                    if ((xx > -1) && (xx < (img.getWidth() - 1))) {
                        for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                            dy++) {
                            float yy = y + dy;
                            if ((yy > -1) && (yy < (img.getHeight() - 1))) {
                                img.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                            }
                        }
                    }
                }
            }            
        }
        
        try {
            
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(dirPath + "/contours_" 
                + System.currentTimeMillis() + ".png", img);
        
        } catch (IOException e) {}
    }
}
