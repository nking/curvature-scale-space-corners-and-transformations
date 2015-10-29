package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.misc.MiscDebug;
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
        final BlobPerimeterHelper blobPerimeterHelper,
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

}
