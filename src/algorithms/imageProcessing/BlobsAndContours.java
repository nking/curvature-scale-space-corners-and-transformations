package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * holds the blobs, their ordered perimeters, and the scale space image
 * contours of the perimeters
 *
 * @author nichole
 */
public class BlobsAndContours {

    private final List<Set<PairInt>> blobs;

    private final List<PairIntArray> blobOrderedPerimeters;

    private final List<List<CurvatureScaleSpaceContour>> contours;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected final boolean segmentedToLineDrawing;
    
    private boolean debug = false;

    private String debugTag = "";

    /**
     * constructor, does all the work of extracting blobs, perimeter, and
     * scale space image contours.
     *
     * @param img
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @param segmentedToLineDrawing true if image contains mostly white
     * pixels and black lines for object contours (this is the result of
     * segmentation using adaptive mean thresholding, for example)
     */
    public BlobsAndContours(GreyscaleImage img, int smallestGroupLimit,
        int largestGroupLimit, boolean segmentedToLineDrawing) {

        blobs = new ArrayList<Set<PairInt>>();

        blobOrderedPerimeters = new ArrayList<PairIntArray>();

        contours = new ArrayList<List<CurvatureScaleSpaceContour>>();

        this.segmentedToLineDrawing = segmentedToLineDrawing;
        
        init(img, smallestGroupLimit, largestGroupLimit);
    }

    /**
     * constructor for using in debug mode, does all the work of extracting
     * blobs, perimeter, and scale space image contours.
     *
     * @param img
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * segmentedToLineDrawing true if image contains mostly white
     * pixels and black lines for object contours (this is the result of
     * segmentation using adaptive mean thresholding, for example)
     * @param debugTag
     */
    public BlobsAndContours(GreyscaleImage img, int smallestGroupLimit,
        int largestGroupLimit, boolean segmentedToLineDrawing, String debugTag) {

        blobs = new ArrayList<Set<PairInt>>();

        blobOrderedPerimeters = new ArrayList<PairIntArray>();

        contours = new ArrayList<List<CurvatureScaleSpaceContour>>();

        this.segmentedToLineDrawing = segmentedToLineDrawing;
        
        debug = true;

        this.debugTag = debugTag;

        init(img, smallestGroupLimit, largestGroupLimit);
    }

    protected void init(GreyscaleImage img, int smallestGroupLimit,
        int largestGroupLimit) {

        extractBlobsFromSegmentedImage(img, blobs, smallestGroupLimit,
            largestGroupLimit);

         boolean discardWhenCavityIsSmallerThanBorder = true;

        extractBoundsOfBlobs(img, blobs, blobOrderedPerimeters,
            discardWhenCavityIsSmallerThanBorder);

        populateContours(blobOrderedPerimeters, contours);
    }

    /**
     * @param segImg segmented image
     * @param outputBlobs
     * @param smallestGroupLimit
     * @param largestGroupLimit
     */
    protected void extractBlobsFromSegmentedImage(GreyscaleImage segImg,
        List<Set<PairInt>> outputBlobs, int smallestGroupLimit,
        int largestGroupLimit) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        Map<Integer, Integer> freqMap = Histogram.createAFrequencyMap(segImg);

        for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()) {

            Integer pixValue = entry.getKey();

            DFSContiguousValueFinder finder = new DFSContiguousValueFinder(segImg);
            finder.setToUse8Neighbors();
            if (smallestGroupLimit == 0) {
                finder.setMinimumNumberInCluster(1);
            } else {
                finder.setMinimumNumberInCluster(smallestGroupLimit);
            }
            
            finder.findGroups(pixValue.intValue());

            int nGroups = finder.getNumberOfGroups();

            for (int i = 0; i < nGroups; ++i) {

                PairIntArray xy = finder.getXY(i);

                if (xy.getN() < largestGroupLimit) {

                    Set<PairInt> points = Misc.convert(xy);

                    // skip blobs that are on the image boundaries because they
                    // are incomplete
                    if (!curveHelper.hasNumberOfPixelsOnImageBoundaries(3,
                        points, segImg.getWidth(), segImg.getHeight())) {

                        outputBlobs.add(points);
                    }
                }
            }
        }
        
        if (segmentedToLineDrawing) {
            // for line drawings, there may be a blob due to an objects 
            // line and to it's interior points, so we want to remove the
            // blob for the interior points and keep the exterior.  choosing
            // the exterior because later feature matching should be better
            // for points outside of the blob which is largely similar content.
            
            removeRedundantBlobs(outputBlobs);
            
        }

        if (debug) {
            Image img0 = ImageIOHelper.convertImage(segImg);
            int c = 0;
            for (int i = 0; i < outputBlobs.size(); ++i) {
                Set<PairInt> blobSet = outputBlobs.get(i);
                int clr = ImageIOHelper.getNextColorRGB(c);
                for (PairInt p : blobSet) {
                    int x = p.getX();
                    int y = p.getY();
                    ImageIOHelper.addPointToImage(x, y, img0, 0, clr);
                }
                c++;
            }
            algorithms.misc.MiscDebug.writeImageCopy(img0, "blobs_" + debugTag
                + "_" + MiscDebug.getCurrentTimeFormatted() + ".png");
        }
    }

    /**
     * given the list of blob points, extract the ordered boundaries of them
     * and remove any blobs for which the bounds were not extractable.
     * @param segImg
     * @param inOutBlobs
     * @param outputBounds
     * @param discardWhenCavityIsSmallerThanBorder
     */
    protected void extractBoundsOfBlobs(final GreyscaleImage segImg,
        final List<Set<PairInt>> inOutBlobs,
        final List<PairIntArray> outputBounds,
        boolean discardWhenCavityIsSmallerThanBorder) {

        int width = segImg.getWidth();

        int height = segImg.getHeight();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<Integer> remove = new ArrayList<Integer>();

        for (int i = 0; i < inOutBlobs.size(); ++i) {

            Set<PairInt> blob = inOutBlobs.get(i);

            EdgeExtractorForBlobBorder extractor = new EdgeExtractorForBlobBorder();

            if (debug) {
                //extractor.setToDebug();
            }

            PairIntArray closedEdge = extractor.extractAndOrderTheBorder0(
                blob, width, height,
                discardWhenCavityIsSmallerThanBorder);

            if ((closedEdge != null) &&
                (curveHelper.isAdjacent(closedEdge, 0, closedEdge.getN() - 1))) {

                /*
                int nChanged = 0;

                adjusting the edges using the unsegmented image did not work
                as well on some images, so have disable this block.
                keeping it until can review if it helps with some types of
                segmentation.

                if (blobIsDarkerThanExterior(blob, closedEdge, img)) {
                    nChanged = curveHelper.adjustEdgesTowardsBrighterPixels(
                        closedEdge, img);

                } else {
                    nChanged = curveHelper.adjustEdgesTowardsDarkerPixels(
                        closedEdge, img);
                }

                if (nChanged > 0) {

                    //TODO: this method needs to be revisited...
                    //curveHelper.removeRedundantPoints(closedEdge);

                    curveHelper.pruneAdjacentNeighborsTo2(closedEdge);

                    curveHelper.correctCheckeredSegments(closedEdge);
                }*/

                if (debug) {
                    Image img0 = ImageIOHelper.convertImage(segImg);
                    PairIntArray pa = closedEdge;
                    for (int j = 0; j < pa.getN(); ++j) {
                        int x = pa.getX(j);
                        int y = pa.getY(j);
                        ImageIOHelper.addPointToImage(x, y, img0, 0, 255, 0, 0);
                    }
                    ImageIOHelper.addPointToImage(pa.getX(0), pa.getY(0), img0, 0, 255, 255, 0);
                    ImageIOHelper.addPointToImage(pa.getX(pa.getN() - 1), pa.getY(pa.getN() - 1), img0, 0, 255, 255, 0);
                    MiscDebug.writeImageCopy(img0, "closed_edge_" + debugTag
                        + "_" + Integer.valueOf(i) +
                        "_" + MiscDebug.getCurrentTimeFormatted() + ".png");
                }

                outputBounds.add(closedEdge);

            } else {

                remove.add(Integer.valueOf(i));
            }
        }

        for (int i = (remove.size() - 1); i >  -1; --i) {
            inOutBlobs.remove(remove.get(i).intValue());
        }

        assert(inOutBlobs.size() == outputBounds.size());

        // sort by descending length
        int[] indexes = new int[inOutBlobs.size()];
        int[] lengths = new int[indexes.length];
        for (int i = 0; i < inOutBlobs.size(); ++i) {
            indexes[i] = i;
            lengths[i] = outputBounds.get(i).getN();
        }
        MultiArrayMergeSort.sortByDecr(lengths, indexes);

        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        List<PairIntArray> curves = new ArrayList<PairIntArray>();

        log.info("nBlobs before filtered to top =" + inOutBlobs.size());

        int last = (inOutBlobs.size() > 15) ? 15 : inOutBlobs.size();
        for (int i = 0; i < last; ++i) {
            int idx = indexes[i];
            blobs2.add(inOutBlobs.get(idx));
            curves.add(outputBounds.get(idx));
        }

        inOutBlobs.clear();
        outputBounds.clear();

        inOutBlobs.addAll(blobs2);
        outputBounds.addAll(curves);

        log.info("nBlobs after filtered to top =" + inOutBlobs.size());

        //TODO: put debug sections in AOP for special build after replace aspectj
        if (debug) {
            Image img0 = ImageIOHelper.convertImage(segImg);
            for (int i = 0; i < outputBounds.size(); ++i) {
                PairIntArray pa = outputBounds.get(i);
                for (int j = 0; j < pa.getN(); ++j) {
                    int x = pa.getX(j);
                    int y = pa.getY(j);
                    if (i == 0) {
                        if (j == 0 || (j == (pa.getN() - 1))) {
                            ImageIOHelper.addPointToImage(x, y, img0, 0, 200, 100, 0);
                        } else {
                            ImageIOHelper.addPointToImage(x, y, img0, 0, 255, 0, 0);
                        }
                    } else if (i == 1) {
                        ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 255, 0);
                    } else {
                        ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 0, 255);
                    }
                }
            }
            MiscDebug.writeImageCopy(img0, "blob_perimeters_" + debugTag +
                "_" + MiscDebug.getCurrentTimeFormatted() + ".png");
        }
    }

    protected void populateContours(final List<PairIntArray> closedContours,
        final List<List<CurvatureScaleSpaceContour>> outputContours) {

        outputContours.clear();

        List<ScaleSpaceCurveImage> scaleSpaceImages
            = new ArrayList<ScaleSpaceCurveImage>();

        boolean setToExtractWeakCurvesTooIfNeeded = false;

        boolean allContoursZero = true;

        for (int edgeIndex = 0; edgeIndex < closedContours.size(); ++edgeIndex) {

            PairIntArray edge = closedContours.get(edgeIndex);
            
            ScaleSpaceCurveImage sscImg =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.createScaleSpaceImage(
                    edge, edgeIndex);

            scaleSpaceImages.add(sscImg);

            List<CurvatureScaleSpaceContour> c =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);

            outputContours.add(c);

            if (!c.isEmpty()) {
                allContoursZero = false;
            }
        }

        //TODO: consider whether want to include weak contours
        if (allContoursZero) {

            setToExtractWeakCurvesTooIfNeeded = true;

            outputContours.clear();

            for (int edgeIndex = 0; edgeIndex < closedContours.size(); ++edgeIndex) {

                PairIntArray edge = closedContours.get(edgeIndex);
                
                ScaleSpaceCurveImage sscImg = scaleSpaceImages.get(edgeIndex);

                List<CurvatureScaleSpaceContour> c =
                    CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                    sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);

                outputContours.add(c);
            }
        }

        assert(closedContours.size() == outputContours.size());
    }

    /**
     * @return the blobs
     */
    public List<Set<PairInt>> getBlobs() {
        return blobs;
    }

    /**
     * @return the blobOrderedPerimeters
     */
    public List<PairIntArray> getBlobOrderedPerimeters() {
        return blobOrderedPerimeters;
    }

    /**
     * @return the contours
     */
    public List<List<CurvatureScaleSpaceContour>> getContours() {
        return contours;
    }

    /**
     * @return the debugTag
     */
    public String getDebugTag() {
        return debugTag;
    }

    /**
     * @return the debug
     */
    public boolean isDebug() {
        return debug;
    }

    // ~O(N_blob) to compare centroids.
    private void removeRedundantBlobs(List<Set<PairInt>> outputBlobs) {
        
        // for line drawings, there may be a blob due to an objects 
        // line and to it's interior points, so we want to remove the
        // blob for the interior points and keep the exterior.  choosing
        // the exterior because later feature matching should be better
        // for points outside of the blob which is largely similar content.
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        Map<PairInt, Integer> blobCenters = new HashMap<PairInt, Integer>();
        LinkedHashSet<Integer> remove = new LinkedHashSet<Integer>();
        for (int i = 0; i < outputBlobs.size(); ++i) {
            Set<PairInt> blob = outputBlobs.get(i);
            double[] xyCen = curveHelper.calculateXYCentroids(blob);
            PairInt p = new PairInt((int)Math.round(xyCen[0]), 
                (int)Math.round(xyCen[1]));
            Integer index = Integer.valueOf(i);
            
            //keep smaller blob.  it's just the line, while larger blob is all of
            //the interior points.
            if (blobCenters.containsKey(p)) {
                // compare and keep exterior perimeter
                int nCurrent = blob.size();
                int mapIdx = blobCenters.get(p).intValue();
                int nInMap = outputBlobs.get(mapIdx).size();
                if (nCurrent > nInMap) {
                    remove.add(Integer.valueOf(i));
                    continue;
                }
                remove.add(Integer.valueOf(mapIdx));
                blobCenters.put(p, index);
            } else {
                blobCenters.put(p, index);
            }
        }
        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;
        for (Entry<PairInt, Integer> entry : blobCenters.entrySet()) {
            Integer index = entry.getValue();
            PairInt p = entry.getKey();
            if (remove.contains(index)) {
                continue;
            }
            for (int nIdx = 0; nIdx < dxs8.length; ++nIdx) {
                int x2 = p.getX() + dxs8[nIdx];
                int y2 = p.getY() + dys8[nIdx];
                PairInt p2 = new PairInt(x2, y2);
                Integer index2 = blobCenters.get(p2);
                if (index2 == null) {
                    continue;
                }
                if (remove.contains(index2)) {
                    continue;
                }
                // if arrive here, there were 2 blobs with adjacent centers
                int n1 = outputBlobs.get(index.intValue()).size();
                int n2 = outputBlobs.get(index2.intValue()).size();
                //keep smaller blob.  it's just the line, while larger blob is all of
                //the interior points.
                if (n1 > n2) {
                    // keep blob 2
                    remove.add(index);
                } else {
                    // keep blob 1
                    remove.add(index2);
                }
                break;
            }
        }
        if (remove.isEmpty()) {
            return;
        }
        List<Integer> rm = new ArrayList<Integer>(remove);
        Collections.sort(rm);
        for (int i = (rm.size() - 1); i > -1; --i) {
            int idx = rm.get(i).intValue();
            outputBlobs.remove(idx);
        }
    }

}
