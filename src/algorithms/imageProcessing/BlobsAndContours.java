package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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
    
    private final int smallestGroupLimit;
    
    private final int largestGroupLimit;
    
    private final SegmentationType type;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected final boolean segmentedToLineDrawing;
    
    protected final GreyscaleImage imgGrey;
    
    protected final GreyscaleImage imgSeg;
    
    private boolean modifyBlobs = true;
    
    private boolean debug = false;

    private String debugTag = "";

    /**
     * constructor, does all the work of extracting blobs, perimeter, and
     * scale space image contours.
     *
     * @param imgGrey
     * @param imgSeg
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @param type
     * @param segmentedToLineDrawing true if image contains mostly white
     * pixels and black lines for object contours (this is the result of
     * segmentation using adaptive mean thresholding, for example)
     */
    public BlobsAndContours(GreyscaleImage imgGrey, GreyscaleImage imgSeg,
        int smallestGroupLimit, int largestGroupLimit, 
        SegmentationType type, 
        boolean segmentedToLineDrawing) {

        blobs = new ArrayList<Set<PairInt>>();

        blobOrderedPerimeters = new ArrayList<PairIntArray>();

        contours = new ArrayList<List<CurvatureScaleSpaceContour>>();

        this.segmentedToLineDrawing = segmentedToLineDrawing;
        
        this.smallestGroupLimit = smallestGroupLimit;
        
        this.largestGroupLimit = largestGroupLimit;
        
        this.type = type;
        
        this.imgGrey = imgGrey;
        
        this.imgSeg = imgSeg;
            
        init();
    }

    /**
     * constructor for using in debug mode, does all the work of extracting
     * blobs, perimeter, and scale space image contours.
     *
     * @param imgGrey
     * @param imgSeg
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * segmentedToLineDrawing true if image contains mostly white
     * pixels and black lines for object contours (this is the result of
     * segmentation using adaptive mean thresholding, for example)
     * @param type
     * @param segmentedToLineDrawing
     * @param debugTag
     */
    public BlobsAndContours(GreyscaleImage imgGrey, GreyscaleImage imgSeg,
        int smallestGroupLimit, int largestGroupLimit, 
        SegmentationType type,
        boolean segmentedToLineDrawing, String debugTag) {

        blobs = new ArrayList<Set<PairInt>>();

        blobOrderedPerimeters = new ArrayList<PairIntArray>();

        contours = new ArrayList<List<CurvatureScaleSpaceContour>>();

        this.segmentedToLineDrawing = segmentedToLineDrawing;
        
        this.type = type;
        
        debug = true;

        this.debugTag = debugTag;

        this.smallestGroupLimit = smallestGroupLimit;
        
        this.largestGroupLimit = largestGroupLimit;
        
        this.imgGrey = imgGrey;
        
        this.imgSeg = imgSeg;
            
        init();
    }

    protected void init() {

        extractBlobsFromSegmentedImage(imgSeg, blobs, smallestGroupLimit,
            largestGroupLimit);

        boolean discardWhenCavityIsSmallerThanBorder = true;

        extractBoundsOfBlobs(imgGrey, imgSeg, blobs, blobOrderedPerimeters,
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
        
        if (modifyBlobs) {
            // helps to fill in single pixel holes
            for (Set<PairInt> blob : outputBlobs) {
                growRadius(blob, segImg.getWidth(), segImg.getHeight());
            }
            /*for (Set<PairInt> blob : outputBlobs) {
                shrinkRadius(blob, segImg.getWidth(), segImg.getHeight());
            }
            
            // if shrink radius, have to redo contiguous search
            List<Integer> rmBlobIndexes = new ArrayList<Integer>();
            for (int i = 0; i < outputBlobs.size(); ++i) {
                Set<PairInt> tmp = redoContiguous(outputBlobs.get(i), 
                    segImg.getWidth(), segImg.getHeight());
                if (tmp.isEmpty()) {
                    rmBlobIndexes.add(Integer.valueOf(i));
                } else {
                    outputBlobs.set(i, tmp);
                }
            }
            for (int i = (rmBlobIndexes.size() - 1); i > -1; --i) {
                outputBlobs.remove(rmBlobIndexes.get(i).intValue());
            }
            */
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
     * @param greyImg
     * @param segImg
     * @param inOutBlobs
     * @param outputBounds
     * @param discardWhenCavityIsSmallerThanBorder
     */
    protected void extractBoundsOfBlobs(final GreyscaleImage greyImg, 
        final GreyscaleImage segImg, final List<Set<PairInt>> inOutBlobs,
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

double[] xyCen = curveHelper.calculateXYCentroids(blob);
if ((Math.abs(xyCen[0] - 124) < 50) && (Math.abs(xyCen[1] - 24) < 50)) {
    //extractor.setToDebug();
    int z = 1;
}
if ((Math.abs(xyCen[0] - 250) < 50) && (Math.abs(xyCen[1] - 50) < 50)) {
    //extractor.setToDebug();
    int z = 1;
}

            PairIntArray closedEdge = extractor.extractAndOrderTheBorder0(
                blob, width, height,
                discardWhenCavityIsSmallerThanBorder);

            if ((closedEdge != null) &&
                (curveHelper.isAdjacent(closedEdge, 0, closedEdge.getN() - 1))) {

                if (false && type.equals(SegmentationType.BINARY)) {
                
                    int nChanged = 0;

                    //adjusting the edges using the unsegmented image did not work
                    //as well on some images, so have disable this block.
                    //keeping it until can review if it helps with some types of
                    //segmentation.

                    if (blobIsDarkerThanExterior(blob, closedEdge, greyImg)) {
                        nChanged = curveHelper.adjustEdgesTowardsBrighterPixels(
                            closedEdge, greyImg);

                    } else {
                        nChanged = curveHelper.adjustEdgesTowardsDarkerPixels(
                            closedEdge, greyImg);
                    }

                    if (nChanged > 0) {

                        //TODO: this method needs to be revisited...
                        //curveHelper.removeRedundantPoints(closedEdge);

                        curveHelper.pruneAdjacentNeighborsTo2(closedEdge);

                        curveHelper.correctCheckeredSegments(closedEdge);
                    }
                }

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

        //TODO: note, results are sensitive to this limit:
        int last = (inOutBlobs.size() > 20) ? 20 : inOutBlobs.size();
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
            long ts = MiscDebug.getCurrentTimeFormatted();
            MiscDebug.writeImageCopy(segImg, "segmentation_" + debugTag +
                "_" + ts + ".png");
            MiscDebug.writeImageCopy(img0, "blob_perimeters_" + debugTag +
                "_" + ts + ".png");
        }
    }

    protected void populateContours(final List<PairIntArray> closedContours,
        final List<List<CurvatureScaleSpaceContour>> outputContours
    ) {

        outputContours.clear();

        List<ScaleSpaceCurveImage> scaleSpaceImages
            = new ArrayList<ScaleSpaceCurveImage>();

        boolean setToExtractWeakCurvesTooIfNeeded = false;
        
        if (type.equals(SegmentationType.BINARY)) {
            // might need to restrict this to the binned images
            setToExtractWeakCurvesTooIfNeeded = true;
        }

        boolean allContoursZero = true;

        for (int edgeIndex = 0; edgeIndex < closedContours.size(); ++edgeIndex) {

            PairIntArray edge = closedContours.get(edgeIndex);
    
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
double[] xycen = curveHelper.calculateXYCentroids(edge);

            ScaleSpaceCurveImage sscImg =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.createScaleSpaceImage(
                    edge, edgeIndex);

            scaleSpaceImages.add(sscImg);

            List<CurvatureScaleSpaceContour> c =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);
            
            /*if (c.isEmpty() && (somewhat large)) {
                c = CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                    sscImg, edgeIndex, true, edge);
            }*/

            //TODO: consider extracting smaller contours for each edge if needed
            // instead of only when there are no strong contours
            
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

    // worse case runtime is O(N_blob^2) to compare centroids and dimensions.
    private void removeRedundantBlobs(List<Set<PairInt>> outputBlobs) {
        
        // for line drawings, there may be a blob due to an objects 
        // line and to it's interior points, so we want to remove the
        // blob for the interior points and keep the exterior.  choosing
        // the exterior because later feature matching should be better
        // for points outside of the blob which is largely similar content.
        
        boolean preferExterior = true;
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        Map<PairInt, Integer> blobCenters = new HashMap<PairInt, Integer>();
        LinkedHashSet<Integer> remove = new LinkedHashSet<Integer>();
        for (int i = 0; i < outputBlobs.size(); ++i) {
            Set<PairInt> blob = outputBlobs.get(i);
            double[] xyCen = curveHelper.calculateXYCentroids(blob);
            PairInt p = new PairInt((int)Math.round(xyCen[0]), 
                (int)Math.round(xyCen[1]));
            Integer index = Integer.valueOf(i);
            
            //keep larger dimension blob
            if (blobCenters.containsKey(p)) {
                // compare and keep exterior perimeter
                int mapIdx = blobCenters.get(p).intValue();
                Set<PairInt> mapBlob = outputBlobs.get(mapIdx);
                 
                //xMin, xMax, yMin, yMax
                int[] minMaxXY1 = MiscMath.findMinMaxXY(blob);
                int[] minMaxXY2 = MiscMath.findMinMaxXY(mapBlob);
                int w1 = minMaxXY1[1] - minMaxXY1[0];
                int h1 = minMaxXY1[3] - minMaxXY1[2];
                int w2 = minMaxXY2[1] - minMaxXY2[0];
                int h2 = minMaxXY2[3] - minMaxXY2[2];
                if ((preferExterior && (w1 < w2) && (h1 < h2)) ||
                    (!preferExterior && (w1 > w2) && (h1 > h2))) {
                    remove.add(Integer.valueOf(i));
                    continue;
                }
                remove.add(Integer.valueOf(mapIdx));
                blobCenters.put(p, index);
            } else {
                blobCenters.put(p, index);
            }
        }
        for (Entry<PairInt, Integer> entry : blobCenters.entrySet()) {
            Integer index = entry.getValue();
            PairInt p = entry.getKey();
            if (remove.contains(index)) {
                continue;
            }
            Set<PairInt> blob1 = outputBlobs.get(index.intValue());
            int[] minMaxXY1 = MiscMath.findMinMaxXY(blob1);
            int w1 = minMaxXY1[1] - minMaxXY1[0];
            int h1 = minMaxXY1[3] - minMaxXY1[2];
            
            // look for another center within dx,dy=(3,3)
            
            for (Entry<PairInt, Integer> entry2 : blobCenters.entrySet()) {
                Integer index2 = entry2.getValue();
                if (index2.intValue() == index.intValue()) {
                    continue;
                }
                if (remove.contains(index2)) {
                    continue;
                }
                PairInt p2 = entry2.getKey();
                
                int diffX = Math.abs(p2.getX() - p.getX());
                int diffY = Math.abs(p2.getY() - p.getY());
                if ((diffX > 3) || (diffY > 3)) {
                    continue;
                }
                // if arrive here, there were 2 blobs with near centers
                Set<PairInt> blob2 = outputBlobs.get(index2.intValue());
                //keep larger dimension blob
                //xMin, xMax, yMin, yMax                
                int[] minMaxXY2 = MiscMath.findMinMaxXY(blob2);
                int w2 = minMaxXY2[1] - minMaxXY2[0];
                int h2 = minMaxXY2[3] - minMaxXY2[2];
                if ((preferExterior && (w1 < w2) && (h1 < h2)) ||
                    (!preferExterior && (w1 >w2) && (h1 > h2))) {                
                    // keep blob 2 for preferExterior
                    remove.add(index);
                } else {
                    // keep blob 1 for preferExterior
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

    private boolean blobIsDarkerThanExterior(Set<PairInt> blob, 
        PairIntArray closedEdge, GreyscaleImage greyImg) {
        
        int w = greyImg.getWidth();
        int h = greyImg.getHeight();
        
        long sumInside = 0;
        
        long sumOutside = 0;
        
        Set<PairInt> counted = new HashSet<PairInt>();
        
        // looking at a band of pixels within a distance of 2 pixels
        // on both sides of the edge.
        
        for (int i = 0; i < closedEdge.getN(); ++i) {
            int x = closedEdge.getX(i);
            int y = closedEdge.getY(i);
            
            for (int dx = -2; dx <= 2; ++dx) {
                int x2 = x + dx;
                if ((x2 < 0) || (x2 > (w - 1))) {
                    continue;
                }
                for (int dy = -2; dy <= 2; ++dy) {
                    int y2 = y + dy;
                    if ((y2 < 0) || (y2 > (h - 1))) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    if (counted.contains(p2)) {
                        continue;
                    }
                    int v = greyImg.getValue(x2, y2);
                    if (blob.contains(p2)) {
                        sumInside += v;
                    } else {
                        sumOutside += v;
                    }
                    counted.add(p2);
                }
            }
        }
        
        if (sumInside < sumOutside) {
            return true;
        }
        
        return false;
    }

    /**
     * grow blob by 1 pixel radius, but exclude growing any image boundary
     * pixels.
     * @param blob
     * @param w
     * @param h 
     */
    private void growRadius(Set<PairInt> blob, int w, int h) {
        
        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;
        
        Set<PairInt> tmp = new HashSet<PairInt>();
                
        for (PairInt p : blob) {
            for (int i = 0; i < dxs8.length; ++i) {
                int x2 = p.getX() + dxs8[i];
                int y2 = p.getY() + dys8[i];
                // do not grow image boundary pixels too
                if ((x2 < 1) || (y2 < 1) || (x2 > (w - 2)) || (y2 > (h - 2))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                tmp.add(p2);
            }
        }
        blob.addAll(tmp);
    }

    /**
     * shrink blob by 1 pixel radius
     * @param blob
     * @param w
     * @param h 
     */
    private void shrinkRadius(Set<PairInt> blob, int w, int h) {
        
        // any point in blob with a non-point neighbor gets removed
        
        Set<PairInt> tmp = new HashSet<PairInt>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (PairInt p : blob) {
            int nN = curveHelper.countNeighbors(p.getX(), p.getY(), blob, w, h);
            if (nN == 8) {
                tmp.add(p);
            }
        }
        blob.clear();
        blob.addAll(tmp);
    }

    private Set<PairInt> redoContiguous(Set<PairInt> points, int w, int h) {
        
        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        if (smallestGroupLimit == 0) {
            finder.setMinimumNumberInCluster(1);
        } else {
            finder.setMinimumNumberInCluster(smallestGroupLimit);
        }
        finder.findConnectedPointGroups(points, w, h);
            
        int nGroups = finder.getNumberOfGroups();

        int nMaxGroup = Integer.MIN_VALUE;
        int groupIdx = -1;
        
        for (int i = 0; i < nGroups; ++i) {
            Set<PairInt> group = finder.getXY(i);
            int nGroup = group.size();
            if (nGroup < largestGroupLimit) {
                if (nGroup > nMaxGroup) {
                    groupIdx = i;
                    nMaxGroup = nGroup;
                }
            }
        }
        
        if (groupIdx == -1) {
            return new HashSet<PairInt>();
        }
        
        return finder.getXY(groupIdx);
    }
}
