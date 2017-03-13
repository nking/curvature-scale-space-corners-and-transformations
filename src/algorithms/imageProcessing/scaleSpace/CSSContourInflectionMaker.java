package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.SIGMA;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

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
public final class CSSContourInflectionMaker {

    protected final Logger log = Logger.getLogger(this.getClass().getName());
    protected boolean debug = false;
    protected long debugTS = 0;
    protected boolean useLineDrawingMode = false;
    protected boolean initialized = false;
    
    protected List<PairIntArray> edges = null;
    protected List<List<CurvatureScaleSpaceContour>> contourLists = new ArrayList<List<CurvatureScaleSpaceContour>>();
    
    protected final ImageExt image;
    protected final Image originalImage;

    public CSSContourInflectionMaker(ImageExt image) {

        this.image = image;

        originalImage = image.copyImage();
    }

    public void setToUseLineDrawingLineMode() {
        this.useLineDrawingMode = true;
    }

    public void setToDebug() {
        debug = true;
        debugTS = System.currentTimeMillis();
    }

    public void findContours() {

        if (initialized) {
            return;
        }

        initialized = true;

        createEdges();

        populateContours(edges, contourLists);

    }
    
    protected void createEdges() {

        CurvatureScaleSpaceImageMaker imgMaker 
            = new CurvatureScaleSpaceImageMaker(image, useLineDrawingMode);

        edges = imgMaker.getClosedCurves();
    }

    protected Image getImage() {
        return image;
    }

    Image getOriginalImage() {
        return originalImage;
    }

    private void extract(List<CurvatureScaleSpaceContour> contours,
        PairIntArray outputXY, List<Float> outputSigmaWeights) {

        float sumSigma = 0;

        for (int i = 0; i < contours.size(); i++) {

            CurvatureScaleSpaceContour c = contours.get(i);

            for (int j = 0; j < c.getPeakDetails().length; j++) {

                CurvatureScaleSpaceImagePoint spaceImagePoint =
                    c.getPeakDetails()[j];

                int x = spaceImagePoint.getXCoord();
                int y = spaceImagePoint.getYCoord();

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

    /**
     * when peak details has more than one point, this averages them and
     * replaces the details with a single point.
     * @param contours
     */
    protected void correctPeaks(List<CurvatureScaleSpaceContour> contours,
        PairIntArray edge) {

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

                int iMid = (p0.getCoordIdx() + p1.getCoordIdx())/2;

                float s = p0.getSigma();
                float tAvg = (p0.getScaleFreeLength() + p1.getScaleFreeLength())/2.f;

                int xMid = edge.getX(iMid);
                int yMid = edge.getY(iMid);

                CurvatureScaleSpaceImagePoint pMid =
                    new CurvatureScaleSpaceImagePoint(s, tAvg, xMid, yMid, iMid);
                CurvatureScaleSpaceImagePoint[] p =
                    new CurvatureScaleSpaceImagePoint[]{pMid};
                c1.setPeakDetails(p);
                contours.set(i, c1);
            }
        }
    }
    
    private void populateContours(List<PairIntArray> edges,
        List<List<CurvatureScaleSpaceContour>> contours) {

        CurvatureScaleSpaceCurvesMaker csscMaker = new CurvatureScaleSpaceCurvesMaker();

        // if use 2^(1/8) as a sigma factor should result in an error less than 10%
        // in determing the peak of a contour.  smaller factors have smaller
        // errors than that.
        float factor = (float)Math.pow(2, 1./32.);

        for (int i = 0; i < edges.size(); i++) {

            PairIntArray edge = edges.get(i);

            Map<Float, ScaleSpaceCurve> scaleSpaceMap =
                csscMaker.createScaleSpaceMetricsForEdge(edge, factor,
                SIGMA.ONE, SIGMA.TWOHUNDREDANDFIFTYSIX);
           
            // x axis is "t", that is the indexes of edge expressed as fraction
            //   of 1.
            // y axis is sigma.  
            // the x and y are the properties of inflection points.
            ScaleSpaceCurveImage scaleSpaceImage =
                csscMaker.convertScaleSpaceMapToSparseImage(
                scaleSpaceMap, i, edge.getN());

            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
                scaleSpaceImage, i);

            correctPeaks(result, edge);

            removeRedundant(result);

            boolean reversed = contourFinder.reverseIfClockwise(result, edge);

            if (reversed) {
                log.info("EDGES: contour isCW=true");

                // these are extracted from contourFinder in order of decreasing
                // sigma already, so only need to be sorted if the list was
                // reversed
                Collections.sort(result, new DescendingSigmaComparator());
            }
            
            if (debug) {
                try {
                    String fileSuffix = "edge_" + i + "_" + debugTS;
                    MiscDebug.printScaleSpaceCurve(scaleSpaceImage, fileSuffix);
                    
                    System.out.println("printing " + scaleSpaceMap.size() + 
                        " scale space contours");
                    fileSuffix = "_scalespace_" + i + "_" + debugTS;
                    MiscDebug.printScaleSpaceMap(scaleSpaceMap, fileSuffix, 
                        10);
                } catch (IOException ex) {
                }
                MiscDebug.debugPlot(result, (ImageExt) image.copyImage(), 0, 0,
                    "_" + i + "_" + debugTS);
            }

            contours.add(result);
        }

        if (contours.isEmpty()) {
            log.info("no contours found in image 1");
            return;
        }
        if (edges.size() > 1) {
            Collections.sort(contours, new DescendingSigmaComparator2());
        }

    }
 
    public List<List<CurvatureScaleSpaceContour>> getContours() {
        return contourLists;
    }

    protected List<PairIntArray> getEdges() {
        return edges;
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

    private float findMaxSigmaOfFirstPeaks(
        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatches1) {

        float maxPeakSigma = Float.MIN_VALUE;

        for (Map.Entry<Integer, List<CurvatureScaleSpaceContour>> entry : bestMatches1.entrySet()) {

            List<CurvatureScaleSpaceContour> list = entry.getValue();
            if (list.isEmpty()) {
                continue;
            }
            CurvatureScaleSpaceContour contour = list.get(0);
            float peakSigma = contour.getPeakSigma();
            if (peakSigma > maxPeakSigma) {
                maxPeakSigma = peakSigma;
            }
        }

        return maxPeakSigma;
    }
}
