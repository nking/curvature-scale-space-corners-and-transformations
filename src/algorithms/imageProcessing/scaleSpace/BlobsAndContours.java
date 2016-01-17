package algorithms.imageProcessing.scaleSpace;

import algorithms.QuickSort;
import algorithms.imageProcessing.SegmentationType;
import algorithms.imageProcessing.features.BlobPerimeterCornerHelper;
import algorithms.imageProcessing.features.BlobPerimeterRegion;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * holds the blobs, their ordered perimeters, and the scale space image
 * contours of the perimeters
 *
 * @author nichole
 */
public class BlobsAndContours  {

    protected static boolean debug = false;
    protected static String debugTag = "";

    public static List<List<CurvatureScaleSpaceContour>> populateContours(
        final BlobPerimeterCornerHelper blobPerimeterHelper,
        SegmentationType type, boolean useBinnedImage) {

        List<List<CurvatureScaleSpaceContour>> contours
            = new ArrayList<List<CurvatureScaleSpaceContour>>();

        List<ScaleSpaceCurveImage> scaleSpaceImages
            = new ArrayList<ScaleSpaceCurveImage>();

        boolean setToExtractWeakCurvesTooIfNeeded = false;

        boolean allContoursZero = true;

        List<PairIntArray> blobOrderedPerimeters =
            blobPerimeterHelper.getBlobPerimeters(type, useBinnedImage);

        for (int edgeIndex = 0; edgeIndex < blobOrderedPerimeters.size(); ++edgeIndex) {

            PairIntArray edge = blobOrderedPerimeters.get(edgeIndex);

            ScaleSpaceCurveImage sscImg =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.createScaleSpaceImage(
                    edge, edgeIndex);

            scaleSpaceImages.add(sscImg);
        }

        if (debug) {
            plot(scaleSpaceImages, "before_contours_" + debugTag);
        }

        for (int edgeIndex = 0; edgeIndex < blobOrderedPerimeters.size(); ++edgeIndex) {

            PairIntArray edge = blobOrderedPerimeters.get(edgeIndex);

            ScaleSpaceCurveImage sscImg = scaleSpaceImages.get(edgeIndex);

            List<CurvatureScaleSpaceContour> c =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);

            /*if (c.isEmpty() && (somewhat large)) {
                c = CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                    sscImg, edgeIndex, true, edge);
            }*/

            //TODO: consider extracting smaller contours for each edge if needed
            // instead of only when there are no strong contours

            contours.add(c);

            if (!c.isEmpty()) {
                allContoursZero = false;
            }
        }

        //TODO: consider whether want to include weak contours
        if (allContoursZero) {

            setToExtractWeakCurvesTooIfNeeded = true;

            contours.clear();

            for (int edgeIndex = 0; edgeIndex < blobOrderedPerimeters.size(); ++edgeIndex) {

                PairIntArray edge = blobOrderedPerimeters.get(edgeIndex);

                ScaleSpaceCurveImage sscImg = scaleSpaceImages.get(edgeIndex);

                List<CurvatureScaleSpaceContour> c =
                    CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                    sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded, edge);

                contours.add(c);
            }
        }

        if (debug) {
            plot(scaleSpaceImages, "debug_after_contours_" + debugTag);
        }

        assert(blobOrderedPerimeters.size() == contours.size());
        
        combineOverlappingPeaks(0.04f, 6.f, blobOrderedPerimeters, contours);

        if (debug) {
            plot(scaleSpaceImages, "debug_after_contours_" + debugTag);
        }

        assert(blobOrderedPerimeters.size() == contours.size());

        return contours;
    }

    private static void plot(List<ScaleSpaceCurveImage> scaleSpaceImages,
        String label) {

        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

            for (int i = 0; i < scaleSpaceImages.size(); ++i) {

                ScaleSpaceCurveImage sscImg = scaleSpaceImages.get(i);

                String str = Integer.toString(i);

                MiscDebug.plotScaleSpaceCurve(plotter, sscImg, str);
            }

            String filePath = plotter.writeFile(label);

            Logger.getLogger(BlobsAndContours.class.getName()).info(
                "write filePath=" + filePath);

            int z = 1;

        } catch (IOException ex) {
        }
    }

    /**
     * contour peaks that overlap on the scale-free axis of the scale space
     * images are combined if they are close in (x, y) space too.
     *
     * @param tolT
     * @param tolD
     * @param blobOrderedPerimeters
     * @param contoursList
     */
    protected static void combineOverlappingPeaks(final float tolT, 
        final float tolD, List<PairIntArray> blobOrderedPerimeters,
        List<List<CurvatureScaleSpaceContour>> contoursList) {

        for (int listIdx = 0; listIdx < contoursList.size(); ++listIdx) {

            PairIntArray curve = blobOrderedPerimeters.get(listIdx);
            
            List<CurvatureScaleSpaceContour> contours = contoursList.get(listIdx);

            combineOverlappingPeaks(tolT, tolD, curve, contours);
        }
    }
    
    private static double distance(float x0, float y0, float x1, float y1) {
        double dX = x0 - x1;
        double dY = y0 - y1;
        double dist = Math.sqrt(dX * dX + dY * dY);
        return dist;
    }

    private static List<Set<Integer>> findOverlappingPeaks(float[] t, 
        float tolT, int[] idxs, float[]x, float[]y, float tolD, 
        boolean isClosedCurve) {
        
        int n = t.length;
        
        List<Set<Integer>> combineIndexes = new ArrayList<Set<Integer>>();

        for (int i = 0; i < (n - 1); ++i) {
            float deltaT = t[i + 1] - t[i];
            if (deltaT > tolT) {
                continue;
            }
            Set<Integer> set = new HashSet<Integer>();
            set.add(Integer.valueOf(i));
            set.add(Integer.valueOf(i + 1));
            // check wrap around for i==0
            if ((i == 0) && isClosedCurve) {
                float deltaT2 = (1 - t[n - 1]) + t[0];
                if ((deltaT2 <= tolT) && (distance(x[idxs[n - 1]], y[idxs[n - 1]], 
                    x[idxs[0]], y[idxs[0]]) <= tolD)) {
                    
                    set.add(Integer.valueOf(n - 1));
                    for (int j = (n - 1); j > 0; --j) {
                        deltaT2 = t[j] - t[j - 1];
                        if ((deltaT2 <= tolT) 
                            && (distance(x[idxs[j]], y[idxs[j]], 
                            x[idxs[j - 1]], y[idxs[j - 1]]) <= tolD)) {
                            set.add(Integer.valueOf(j - 1));
                        } else {
                            break;
                        }
                    }
                }
            }
            // check forward until not within tolT
            for (int j = (i + 1); j < (n - 1); ++j) {
                i = j;
                float deltaT2 = t[j] - t[j + 1];
                if ((deltaT2 <= tolT) && (distance(x[idxs[j]], y[idxs[j]], 
                    x[idxs[j + 1]], y[idxs[j + 1]]) <= tolD)) {
                    set.add(Integer.valueOf(j + 1));
                } else {
                    break;
                }
            }
            combineIndexes.add(set);
        }
        
        return combineIndexes;
    }

    private static List<Set<Integer>> findOverlappingPeaks(int[] curveIndexes, 
        int tolIndexes, int[] idxs, float[]x, float[]y, float tolD, 
        boolean isClosedCurve) {
        
        int n = curveIndexes.length;
        
        List<Set<Integer>> combineIndexes = new ArrayList<Set<Integer>>();

        for (int i = 0; i < (n - 1); ++i) {
            int deltaT = curveIndexes[i + 1] - curveIndexes[i];
            if (deltaT > tolIndexes) {
                continue;
            }
            Set<Integer> set = new HashSet<Integer>();
            set.add(Integer.valueOf(i));
            set.add(Integer.valueOf(i + 1));
            // check wrap around for i==0
            if ((i == 0) && isClosedCurve) {
                int deltaT2 = (1 - curveIndexes[n - 1]) + curveIndexes[0];
                if ((deltaT2 <= tolIndexes) && (distance(x[idxs[n - 1]], y[idxs[n - 1]], 
                    x[idxs[0]], y[idxs[0]]) <= tolD)) {
                    
                    set.add(Integer.valueOf(n - 1));
                    for (int j = (n - 1); j > 0; --j) {
                        deltaT2 = curveIndexes[j] - curveIndexes[j - 1];
                        if ((deltaT2 <= tolIndexes) 
                            && (distance(x[idxs[j]], y[idxs[j]], 
                            x[idxs[j - 1]], y[idxs[j - 1]]) <= tolD)) {
                            set.add(Integer.valueOf(j - 1));
                        } else {
                            break;
                        }
                    }
                }
            }
            // check forward until not within tolIndexes
            for (int j = (i + 1); j < (n - 1); ++j) {
                i = j;
                int deltaT2 = curveIndexes[j] - curveIndexes[j + 1];
                if ((deltaT2 <= tolIndexes) && (distance(x[idxs[j]], y[idxs[j]], 
                    x[idxs[j + 1]], y[idxs[j + 1]]) <= tolD)) {
                    set.add(Integer.valueOf(j + 1));
                } else {
                    break;
                }
            }
            combineIndexes.add(set);
        }
        
        return combineIndexes;
    }

    /**
     * combine peaks which overlap in the scale free length axis and
     * are close within tolerance in (x, y) space too and return the
     * indexes that were deleted with respect to the given contours
     * list.
     * @param tolT
     * @param tolD
     * @param curve
     * @param contours
     * @return 
     */
    public static List<Integer> combineOverlappingPeaks(float tolT, float tolD, 
        PairIntArray curve, List<CurvatureScaleSpaceContour> contours) {
        
        // presumably these are all closed curves, but verify that:
        boolean isClosedCurve = false;
        if ((curve instanceof PairIntArrayWithColor) &&
            (((PairIntArrayWithColor)curve).getColor() == 1)) {
            isClosedCurve = true;
        }

        /*
        -- make an index array sorted by scale free length
        -- then traverse the sorted peak to find those withing 0.025 of
           one another on scale free axis and within distance of 5 or so in
           (x,y) space and replace them with the average.  the strongest
           sigma is retained for the average.
        */
        int n = contours.size();
        float[] x = new float[n];
        float[] y = new float[n];
        float[] t = new float[n];
        int[] idxs = new int[n];
        for (int i = 0; i < contours.size(); ++i) {

            CurvatureScaleSpaceContour contour = contours.get(i);

            CurvatureScaleSpaceImagePoint[] peakDetails = contour.getPeakDetails();

            assert (peakDetails.length == 1);

            CurvatureScaleSpaceImagePoint peakDetail = peakDetails[0];
            x[i] = peakDetail.getXCoord();
            y[i] = peakDetail.getYCoord();
            t[i] = peakDetail.getScaleFreeLength();
            idxs[i] = i;
        }
        QuickSort.sortBy1stArg(t, idxs, 0, n - 1);
            
        List<Set<Integer>> combineIndexes = findOverlappingPeaks(t, tolT, idxs, 
            x, y, tolD, isClosedCurve);

        if (combineIndexes.isEmpty()) {
            return null;
        }

        List<Integer> remove = new ArrayList<Integer>();
            
        for (Set<Integer> combineSet : combineIndexes) {
            //keep the one with largest sigma
            float sigmaMax = Float.MIN_VALUE;
            int sigmaMaxIdx = -1;
            float xAvg = 0;
            float yAvg = 0;
            float tAvg = 0;
            for (Integer index : combineSet) {
                int idx0 = idxs[index.intValue()];
                tAvg += t[index.intValue()];
                xAvg += x[idx0];
                yAvg += y[idx0];
                CurvatureScaleSpaceImagePoint[] peakDetails = contours.get(idx0).getPeakDetails();
                float sigma = peakDetails[0].getSigma();
                if (sigma > sigmaMax) {
                    sigmaMax = sigma;
                    sigmaMaxIdx = idx0;
                }
            }
            xAvg /= (float)combineSet.size();
            yAvg /= (float)combineSet.size();
            tAvg /= (float)combineSet.size();
            for (Integer index : combineSet) {
                int idx0 = idxs[index.intValue()];
                if (idx0 != sigmaMaxIdx) {
                    remove.add(Integer.valueOf(idx0));
                }
            }
            CurvatureScaleSpaceImagePoint[] peakDetails = contours.get(sigmaMaxIdx).getPeakDetails();
            peakDetails[0] = new CurvatureScaleSpaceImagePoint(sigmaMax, 
                tAvg, Math.round(xAvg), Math.round(yAvg), peakDetails[0].getCoordIdx());
        }

        Collections.sort(remove);
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i).intValue();
            contours.remove(idx);
        }
        
        return remove;
    }

    public static List<List<BlobPerimeterRegion>> populateRegions(BlobPerimeterCornerHelper 
        imgHelper, List<List<CurvatureScaleSpaceContour>> contours, 
        SegmentationType type, boolean useBinned) {
        
        List<PairIntArray> perimeterLists = imgHelper.getBlobPerimeters(type, useBinned);
        
        List<Set<PairInt>> blobs = imgHelper.getBlobs(type, useBinned);
        
        List<List<BlobPerimeterRegion>> bprs = new ArrayList<List<BlobPerimeterRegion>>();
        
        for (int i = 0; i < contours.size(); ++i) {
            
            List<BlobPerimeterRegion> bpr = extractBlobPerimeterRegions(i,
                contours.get(i), perimeterLists.get(i), blobs.get(i));
            
            bprs.add(bpr);
        }
        
        return bprs;
    }

    public static List<BlobPerimeterRegion> extractBlobPerimeterRegions(
        int theEdgeIndex, List<CurvatureScaleSpaceContour> contours,
        PairIntArray closedCurve, Set<PairInt> blob) {

        List<BlobPerimeterRegion> bprList = new ArrayList<BlobPerimeterRegion>();
                
        for (int i = 0; i < contours.size(); ++i) {
            
            CurvatureScaleSpaceContour c = contours.get(i);
            
            assert(c.getPeakDetails().length == 1);
            
            BlobPerimeterRegion bpr = extractBlobPerimeterRegion(theEdgeIndex, 
                c.getPeakDetails()[0], closedCurve, blob);
            
            bpr.setIndexWithinCurve(c.getPeakDetails()[0].getCoordIdx());
                      
            bprList.add(bpr);
        }
        
        return bprList;
    }
    
    /**
     * extract the local points surrounding (x, y) on the
     * perimeter and return an object when creating descriptors.
     * Note that the perimeter is expected to be a closed curve.
     * @param theEdgeIndex
     * @param perimeter
     * @param blob
     * @return
     */
    public static BlobPerimeterRegion extractBlobPerimeterRegion(
        int theEdgeIndex, CurvatureScaleSpaceImagePoint peakDetail, 
        PairIntArray perimeter, Set<PairInt> blob) {
        
        if (perimeter == null || (perimeter.getN() < 4)) {
            throw new IllegalArgumentException(
            "perimeter cannot be null and must have at least 4 points");
        }
        
        if (blob == null) {
            throw new IllegalArgumentException("blob cannot be null");
        }
        
        //TODO: consider asserting that this is a closed curve
        // because of averaging for some peaks, sometimes
        // (x[detailIdx, y[detailIdx]) != (peakDetail.getXCoord(), peakDetail.getYCoord())
        // so for those, have to use the preceding index and then
        // the next + 1 (instead of next)
        
        int detailIdx = peakDetail.getCoordIdx();
        float sigma = peakDetail.getSigma();

        // inflection points are found in between extremes of curvature,
        // so detailIdx must not be one of the curve endpoints.

        int xm = peakDetail.getXCoord();
        int ym = peakDetail.getYCoord();

        int xPrev, yPrev;
        int xNext, yNext;
        
        if (detailIdx == 0) {
            
            // wrap around closed curve
            xPrev = perimeter.getX(perimeter.getN() - 1);
            yPrev = perimeter.getY(perimeter.getN() - 1);
            
            if ((xm != perimeter.getX(detailIdx)) || (ym != perimeter.getY(detailIdx))) {
                xNext = perimeter.getX(detailIdx + 2);
                yNext = perimeter.getY(detailIdx + 2);
            } else {
                xNext = perimeter.getX(detailIdx + 1);
                yNext = perimeter.getY(detailIdx + 1);
            }
            
        } else if ((xm != perimeter.getX(detailIdx)) || (ym != perimeter.getY(detailIdx))) {
            
            xPrev = perimeter.getX(detailIdx - 1);
            yPrev = perimeter.getY(detailIdx - 1);
            
            if ((detailIdx + 2) < perimeter.getN()) {
                xNext = perimeter.getX(detailIdx + 2);
                yNext = perimeter.getY(detailIdx + 2);
            } else {
                // this is a closed curve, so wrap around
                xNext = perimeter.getX(0);
                yNext = perimeter.getY(0);
            }
            
        } else {
            
            xPrev = perimeter.getX(detailIdx - 1);
            yPrev = perimeter.getY(detailIdx - 1);
            
            if ((detailIdx + 1) < perimeter.getN()) {
                xNext = perimeter.getX(detailIdx + 1);
                yNext = perimeter.getY(detailIdx + 1);
            } else {
                // this is a closed curve, so wrap around
                xNext = perimeter.getX(0);
                yNext = perimeter.getY(0);
            }
            
        }
        
        BlobPerimeterRegion region = new BlobPerimeterRegion(theEdgeIndex, 
            xPrev, yPrev, xm, ym, xNext, yNext, blob, sigma);
       
        region.setIndexWithinCurve(detailIdx);

        return region;
    }
    
    /**
     * contour peaks that overlap on the scale-free axis of the scale space
     * images are combined if they are close in (x, y) space too.
     *
     * @param tolT
     * @param tolD
     * @param blobOrderedPerimeters
     * @param bprLists
     */
    protected static void combineOverlappingPeaks0(final int tolIndexes, 
        final float tolD, List<PairIntArray> blobOrderedPerimeters,
        List<List<BlobPerimeterRegion>> bprLists) {

        for (int listIdx = 0; listIdx < bprLists.size(); ++listIdx) {

            PairIntArray curve = blobOrderedPerimeters.get(listIdx);
            
            List<BlobPerimeterRegion> bprs = bprLists.get(listIdx);

            combineOverlappingPeaks0(tolIndexes, tolD, curve, bprs);
        }
    }
    
    /**
     * combine peaks which overlap in the scale free length axis and
     * are close within tolerance in (x, y) space too and return the
     * indexes that were deleted with respect to the given contours
     * list.
     * @param tolIndexes
     * @param tolD
     * @param curve
     * @param bprs
     * @return 
     */
    public static List<Integer> combineOverlappingPeaks0(int tolIndexes, float tolD, 
        PairIntArray curve, List<BlobPerimeterRegion> bprs) {
        
        // presumably these are all closed curves, but verify that:
        boolean isClosedCurve = false;
        if ((curve instanceof PairIntArrayWithColor) &&
            (((PairIntArrayWithColor)curve).getColor() == 1)) {
            isClosedCurve = true;
        }

        /*
        -- make an index array sorted by scale free length
        -- then traverse the sorted peak to find those withing 0.025 of
           one another on scale free axis and within distance of 5 or so in
           (x,y) space and replace them with the average.  the strongest
           sigma is retained for the average.
        */
        int n = bprs.size();
        float[] x = new float[n];
        float[] y = new float[n];
        int[] curveIndexes = new int[n];
        int[] idxs = new int[n];
        for (int i = 0; i < bprs.size(); ++i) {

            BlobPerimeterRegion bpr = bprs.get(i);

            x[i] = bpr.getX()[1];
            y[i] = bpr.getY()[1];
            curveIndexes[i] = bpr.getIndexWithinCurve();
            if (curveIndexes[i] == -1) {
                throw new IllegalStateException("bpr is missing curve index");
            }
            idxs[i] = i;
        }
        
        QuickSort.sortBy1stArg(curveIndexes, idxs);
                    
        List<Set<Integer>> combineIndexes = findOverlappingPeaks(curveIndexes, 
            tolIndexes, idxs, x, y, tolD, isClosedCurve);

        if (combineIndexes.isEmpty()) {
            return null;
        }

        List<Integer> remove = new ArrayList<Integer>();
            
        for (Set<Integer> combineSet : combineIndexes) {
            //keep the one with largest sigma
            float sigmaMax = Float.MIN_VALUE;
            int sigmaMaxIdx = -1;
            
            for (Integer index : combineSet) {
                int idx0 = idxs[index.intValue()];                                
                float sigma = bprs.get(idx0).getSigma();
                if (sigma > sigmaMax) {
                    sigmaMax = sigma;
                    sigmaMaxIdx = idx0;
                }
                
            }
            for (Integer index : combineSet) {
                int idx0 = idxs[index.intValue()];
                if (idx0 != sigmaMaxIdx) {
                    remove.add(Integer.valueOf(idx0));
                }
            }
        }

        Collections.sort(remove);
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i).intValue();
            bprs.remove(idx);
        }
        
        return remove;
    }

}
