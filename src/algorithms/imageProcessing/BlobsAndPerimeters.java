package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class for extracting blobs and their perimeters.
 * 
 * @author nichole
 */
public class BlobsAndPerimeters {
    
    protected static boolean modifyBlobs = true;
    
    /**
     * the number of largest extracted perimeters is limited to this and
     * the number of blobs is subsequently trimmed to the same.  The variable
     * may be adjusted to a higher number for larger images.
     */
    protected static int defaultNumberLimit = 20;
          
    public static List<Set<PairInt>> extractBlobsFromSegmentedImage(
        SegmentedImageHelper imgHelper, SegmentationType type, 
        boolean useBinned) {
        
        List<Set<PairInt>> outputBlobs = new ArrayList<Set<PairInt>>();
        
        // have removed this segmentation type:
        boolean segmentedToLineDrawing = false;
        
        int smallestGroupLimit = useBinned ? 
            imgHelper.getSmallestGroupLimitBinned() : 
            imgHelper.getSmallestGroupLimit();
        
        if (smallestGroupLimit == 0) {
            throw new IllegalStateException("smallestGroupLimit must be > 0");
        }
        
        int largestGroupLimit = useBinned ? 
            imgHelper.getLargestGroupLimitBinned() : 
            imgHelper.getLargestGroupLimit();
        
        //TODO: refactor to use limits differently.  hard wiring an override here
        if (!useBinned && type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            smallestGroupLimit = 5000;
            largestGroupLimit = 100000;
        }
        
        GreyscaleImage segImg = useBinned ?
            imgHelper.getBinnedSegmentationImage(type) :
            imgHelper.getSegmentationImage(type);
            
        boolean use8Neighbors = false;
                
        Map<Integer, Integer> freqMap = Histogram.createAFrequencyMap(segImg);
        
        /*
        TODO: for low resolution, should use defaultExtractBlobs. this code
        doesn't currenly know resolution, so this section may need
        to be revised with more testing.
        */
        
        List<Set<PairInt>> outputExcludedBlobs = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> outputExcludedBoundaryBlobs = new ArrayList<Set<PairInt>>();
        
        for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()) {            
            
            Integer pixValue = entry.getKey();
            
            defaultExtractBlobs(segImg, pixValue.intValue(), 
                smallestGroupLimit, largestGroupLimit, use8Neighbors, 
                outputBlobs, outputExcludedBlobs, outputExcludedBoundaryBlobs, 
                imgHelper.debugTag);
        }
        
        int nExcl = outputExcludedBlobs.size();
        int nIncl = outputBlobs.size();
        
        boolean redo = useBinned && 
            ((nIncl == 0) || ((nExcl >= 10) && (nExcl > nIncl)));
        
        // redo with default size limit and an algorithm to separate blobs connected
        // by only 1 pixel
        if (redo) {
            smallestGroupLimit = imgHelper.getSmallestGroupLimit();
            largestGroupLimit = imgHelper.getLargestGroupLimit();
            List<Set<PairInt>> tmpOutputBlobs = new ArrayList<Set<PairInt>>();
            for (int i = 0; i < outputExcludedBlobs.size(); ++i) {
                extractBlobs2(outputExcludedBlobs.get(i), 
                    smallestGroupLimit, largestGroupLimit, tmpOutputBlobs,
                    segImg.getWidth(), segImg.getHeight());
            }
            for (int i = 0; i < outputExcludedBoundaryBlobs.size(); ++i) {
                extractBlobs2(outputExcludedBoundaryBlobs.get(i), 
                    smallestGroupLimit, largestGroupLimit, tmpOutputBlobs,
                    segImg.getWidth(), segImg.getHeight());
            }
            //TODO: may need to determine group size limits more robustly           
            outputBlobs.clear();
            outputBlobs.addAll(tmpOutputBlobs);
        }
        
        removeRedundantBlobs(outputBlobs);
        
        if (modifyBlobs && !redo) {
            // helps to fill in single pixel holes
            for (Set<PairInt> blob : outputBlobs) {
                growRadius(blob, segImg.getWidth(), segImg.getHeight());
            }
        }
        
        if (imgHelper.isInDebugMode()) {
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
            
            MiscDebug.writeImage(img0, "blobs_" + imgHelper.getDebugTag() 
                + "_" + MiscDebug.getCurrentTimeFormatted());            
        }
        
        return outputBlobs;
    }

    public static List<PairIntArray> extractBoundsOfBlobs(
        SegmentedImageHelper imgHelper, SegmentationType type, 
        final List<Set<PairInt>> inOutBlobs, boolean useBinned,
        boolean discardWhenCavityIsSmallerThanBorder) {
        
        // TODO: method arguments are vulnerable to not being consistently from
        //   same origins, so refactor this (encapsulate) at some point
        
        final List<PairIntArray> outputBounds = new ArrayList<PairIntArray>();
        
        Logger log = Logger.getLogger(BlobsAndPerimeters.class.getName());
        
        GreyscaleImage segImg = useBinned ?
            imgHelper.getBinnedSegmentationImage(type) :
            imgHelper.getSegmentationImage(type);
        
        GreyscaleImage greyImg = useBinned ?
            imgHelper.getGreyscaleImageBinned() :
            imgHelper.getGreyscaleImage();
        
        int width = segImg.getWidth();
        int height = segImg.getHeight();
        
        int numberLimit = defaultNumberLimit;
        if (useBinned) {
            int binFactor = imgHelper.getBinFactor();
            int limit = 1024/binFactor;
            if (width >= limit || height >= limit) {
                numberLimit = 30;
            }
        } else {
            if (width >= 1024 || height >= 1024) {
                numberLimit = 30;
            }
        }
        
        // sort by descending size
        int[] indexes = new int[inOutBlobs.size()];
        int[] lengths = new int[indexes.length];
        for (int i = 0; i < inOutBlobs.size(); ++i) {
            indexes[i] = i;
            lengths[i] = inOutBlobs.get(i).size();
        }
        
        MultiArrayMergeSort.sortByDecr(lengths, indexes);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        
        log.info("nBlobs before filtered =" + inOutBlobs.size());
        
        int count = 0;
        int idx;
        while ((outputBounds.size() < numberLimit) && (count < inOutBlobs.size())) {
            
            idx = indexes[count];
            
            Set<PairInt> blob = inOutBlobs.get(idx);
            
            EdgeExtractorForBlobBorder extractor = new EdgeExtractorForBlobBorder();
            
            PairIntArray closedEdge = extractor.extractAndOrderTheBorder0(blob, 
                width, height, discardWhenCavityIsSmallerThanBorder);
            
            if ((closedEdge != null) && (curveHelper.isAdjacent(closedEdge, 0, 
                closedEdge.getN() - 1))) {

                outputBounds.add(closedEdge);
                
                blobs2.add(blob);
                
            }
            
            count++;
        }
        
        inOutBlobs.clear();
        inOutBlobs.addAll(blobs2);
        
        assert (inOutBlobs.size() == outputBounds.size());
        
        log.info("nBlobs after filtered to top =" + inOutBlobs.size());
        
if (imgHelper.isInDebugMode()) {
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
    MiscDebug.writeImage(img0, "blob_perimeters_" + imgHelper.getDebugTag() + "_" + ts);
}

         return outputBounds;
    }

    // worse case runtime is O(N_blob^2) to compare centroids and dimensions.
    protected static void removeRedundantBlobs(List<Set<PairInt>> outputBlobs) {
        
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
            
            PairInt p = new PairInt((int) Math.round(xyCen[0]), 
                (int) Math.round(xyCen[1]));
            
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
        
        for (Map.Entry<PairInt, Integer> entry : blobCenters.entrySet()) {
            
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
            for (Map.Entry<PairInt, Integer> entry2 : blobCenters.entrySet()) {
                
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
                    (!preferExterior && (w1 > w2) && (h1 > h2))) {
                    
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
        
        for (int i = rm.size() - 1; i > -1; --i) {
            int idx = rm.get(i).intValue();
            outputBlobs.remove(idx);
        }
    }

    /**
     * grow blob by 1 pixel radius, but exclude growing any image boundary
     * pixels.
     * @param blob
     * @param w
     * @param h
     */
    protected static void growRadius(Set<PairInt> blob, int w, int h) {
        
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
     * grow blob by 1 pixel radius if the pixel is in the original pixel
     * set
     * @param blob set of pixels shrunk and possibly partitioned
     * @param originalSet set of pixels before shrink operation
     */
    protected static void growRadius(Set<PairInt> blob, Set<PairInt> originalSet
        ) {
        
        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;
        
        Set<PairInt> tmp = new HashSet<PairInt>();
        
        for (PairInt p : blob) {
            for (int i = 0; i < dxs8.length; ++i) {
                
                int x2 = p.getX() + dxs8[i];
                int y2 = p.getY() + dys8[i];
                
                PairInt p2 = new PairInt(x2, y2);
                
                if (originalSet.contains(p2)) {
                    tmp.add(p2);
                }
            }
        }
        blob.addAll(tmp);
    }

    /**
     * shrink blob by 1 pixel radius
     * @param blob
     */
    protected static void shrinkRadius(Set<PairInt> blob) {

        // any point in blob with a non-point neighbor gets removed
        Set<PairInt> tmp = new HashSet<PairInt>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
      
        for (PairInt p : blob) {
        
            int nN = curveHelper.countNeighbors(p.getX(), p.getY(), blob);

            if (nN == 8) {
                tmp.add(p);
            }
        }
        
        blob.clear();
        blob.addAll(tmp);        
    }
    
    /**
     * extract blobs using a dfs search to find the blobs within the image
     * (contiguous pixels with the given pixelValue) and filter those by
     * the given group size limits.
     * @param segImg
     * @param pixelValue
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @param use8Neighbors
     * @param outputBlobs contiguous pixels of value pixelValue have group
     * sizes between smallest and largestGroupLimit (exclusive)
     * @param outputExcludedBlobs the blobs excluded by the largestGroupLimit
     * @param outputExcludedBoundaryBlobs blobs excluded because they are on
     * the image bounds.  A subsequent method which shrinks the blob, may lead
     * to blobs not on the image boundaries.
     * @param debugTag
     */
    private static void defaultExtractBlobs(GreyscaleImage segImg, 
        int pixelValue, int smallestGroupLimit, int largestGroupLimit, 
        boolean use8Neighbors, List<Set<PairInt>> outputBlobs, 
        List<Set<PairInt>> outputExcludedBlobs,
        List<Set<PairInt>> outputExcludedBoundaryBlobs,
        String debugTag) {
        
        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(segImg);
            
        if (use8Neighbors) {
            finder.setToUse8Neighbors();
        }
//TODO: temporary change to plot the distribution of blob sizes for test images:        
        finder.setMinimumNumberInCluster(1);//smallestGroupLimit);
        finder.findGroups(pixelValue);
            
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        int nGroups = finder.getNumberOfGroups();
        for (int i = 0; i < nGroups; ++i) {
            PairIntArray xy = finder.getXY(i);
            Set<PairInt> points = Misc.convert(xy);
            // skip blobs that are on the image boundaries because they
            // are incomplete
            if (curveHelper.hasNumberOfPixelsOnImageBoundaries(3, 
                points, segImg.getWidth(), segImg.getHeight())) {
                
                outputExcludedBoundaryBlobs.add(points);
                
            } else {                
                if (xy.getN() < largestGroupLimit) {
                    
                    outputBlobs.add(points);
                } else {
                    outputExcludedBlobs.add(points);
                }
            }
        }
                
        if (!debugTag.equals("")) {
            List<Integer> sizes = new ArrayList<Integer>();
            for (int i = 0; i < nGroups; ++i) {
                PairIntArray xy = finder.getXY(i);
                Set<PairInt> points = Misc.convert(xy);
                if (!curveHelper.hasNumberOfPixelsOnImageBoundaries(3, 
                    points, segImg.getWidth(), segImg.getHeight())) {
                    sizes.add(Integer.valueOf(xy.getN()));
                }
            }
            HistogramHolder hist = Histogram.createSimpleHistogram(sizes);
            if (hist != null) {
                try {
                    int yMaxIdx = MiscMath.findLastZeroIndex(hist);
                    float xMax = (yMaxIdx == -1) ? 25000 : hist.getXHist()[yMaxIdx];
                    xMax = (float)Math.ceil(xMax);
                    String label = debugTag + " (0," + xMax + ") " 
                        + " def=(" + smallestGroupLimit + "," + largestGroupLimit + ")"
                        + " pixV=" + pixelValue;                    
                    hist.plotHistogram(0, xMax, label, debugTag + "_" + MiscDebug.getCurrentTimeFormatted());
                } catch (IOException ex) {
                    Logger.getLogger(BlobsAndPerimeters.class.getName()).log(
                        Level.SEVERE, null, ex);
                }
            }
        }        
    }

    /**
     * method to shrink the found blobs by '1' pixel, re-do search, then
     * grow the found blobs by '1' pixel.
     * @param segImg
     * @param pixelValue
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @param use8Neighbors
     * @param outputBlobs 
     */
    private static void extractBlobs2(Set<PairInt> inputBlobs,
        int smallestGroupLimit, int largestGroupLimit, 
        List<Set<PairInt>> outputBlobs, 
        final int imageWidth, final int imageHeight) {
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                                        
        Set<PairInt> points = inputBlobs;

        int[] minMaxXY = MiscMath.findMinMaxXY(points);

        Set<PairInt> pointsB4 = new HashSet<PairInt>(points);

        shrinkRadius(points);

        DFSConnectedGroupsFinder finder2 = new DFSConnectedGroupsFinder();

        finder2.setMinimumNumberInCluster(smallestGroupLimit - 4);

        finder2.findConnectedPointGroups(points, minMaxXY[1], minMaxXY[3]);
            
        int nGroups2 = finder2.getNumberOfGroups();

        for (int ii = 0; ii < nGroups2; ++ii) {

            Set<PairInt> blob2 = finder2.getXY(ii);

            if ((blob2.size() < smallestGroupLimit) || 
                (blob2.size() > largestGroupLimit)) {
                continue;
            }

            if (!curveHelper.hasNumberOfPixelsOnImageBoundaries(3, points, 
                imageWidth, imageHeight)) {

                // grow by '1' pixel to help maintain scale
                growRadius(blob2, pointsB4);
                
                outputBlobs.add(blob2);
            }
        }
    }

}
