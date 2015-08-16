package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
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
    protected List<List<CurvatureScaleSpaceContour>> contours1 = new ArrayList<List<CurvatureScaleSpaceContour>>();
    protected List<List<CurvatureScaleSpaceContour>> contours2 = new ArrayList<List<CurvatureScaleSpaceContour>>();
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
        
        edges1 = getEdges(imgMaker);
        offsetImageX1 = imgMaker.getTrimmedXOffset();
        offsetImageY1 = imgMaker.getTrimmedYOffset();
                
        for (int i = 0; i < edges1.size(); i++) {
            
            PairIntArray curve = edges1.get(i);
            Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                imgMaker.createScaleSpaceMetricsForEdge2(curve);
            ScaleSpaceCurveImage scaleSpaceImage = 
                imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i,
                curve.getN());
/*            
try {
    MiscDebug.printScaleSpaceCurve(scaleSpaceImage,
        MiscDebug.getCurrentTimeFormatted());
} catch (IOException ex) {
    Logger.getLogger(AbstractCurvatureScaleSpaceInflectionMapper.class.getName()).log(Level.SEVERE, null, ex);
}
*/
            ContourFinder contourFinder = new ContourFinder();
            
            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(scaleSpaceImage, i);
            
            boolean reversed = contourFinder.reverseIfClockwise(result);
            
            if (reversed) {
                log.info("EDGES1: contour isCW=true");
                
                // these are extracted from contourFinder in order of decreasing
                // sigma already, so only need to be sorted if the list was
                // reversed
                Collections.sort(result, new DescendingSigmaComparator());
            }
            
MiscDebug.debugPlot(result, (ImageExt)image1.copyImage(), offsetImageX1, offsetImageY1,
    "_1_" + MiscDebug.getCurrentTimeFormatted());

            contours1.add(result);
        }
                
        if (contours1.isEmpty()) {
            log.info("no contours found in image 1");
            return;
        }
        if (edges1.size() > 1) {
            Collections.sort(contours1, new DescendingSigmaComparator2());
        }
       
        /*
        note that when modifying the contour lists in any way, one has to
        maintain decreasing order by sigma and when sigma is equal, the
        order must be by increasing scale free parameter.
        two of the search methods in the matcher depend upon those properties.
         */
        
        imgMaker = new CurvatureScaleSpaceImageMaker(image2);
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
        
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray curve = edges2.get(i);
            Map<Float, ScaleSpaceCurve> scaleSpaceMap = imgMaker.createScaleSpaceMetricsForEdge2(curve);
            ScaleSpaceCurveImage scaleSpaceImage = 
                imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i,
                curve.getN());
            
try {
    MiscDebug.printScaleSpaceCurve(scaleSpaceImage,
        MiscDebug.getCurrentTimeFormatted());
} catch (IOException ex) {
    Logger.getLogger(AbstractCurvatureScaleSpaceInflectionMapper.class.getName()).log(Level.SEVERE, null, ex);
}

            ContourFinder contourFinder = new ContourFinder();
            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(scaleSpaceImage, i);
            
            boolean reversed = contourFinder.reverseIfClockwise(result);
            
            if (reversed) {
                                
                log.info("EDGES2: contour isCW=true");
                
                // these are extracted from contourFinder in order of decreasing
                // sigma already, so only need to be sorted if the list was
                // reversed
                Collections.sort(result, new DescendingSigmaComparator());
            }
            
MiscDebug.debugPlot(result, (ImageExt)image2.copyImage(), offsetImageX2, offsetImageY2,
    "_2_" + MiscDebug.getCurrentTimeFormatted());
            
            contours2.add(result);
        }
      
        if (contours2.isEmpty()) {
            log.info("did not find contours in image 2");
            return;
        }
        if (edges2.size() > 1) {
            Collections.sort(contours2, new DescendingSigmaComparator2());
        }
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
    public List<List<CurvatureScaleSpaceContour>> getContours1() {
        return contours1;
    }

    @Override
    public List<List<CurvatureScaleSpaceContour>> getContours2() {
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
   
}
