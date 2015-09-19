package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.SegmentedImageHelper.SegmentationType;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class BlobScaleFinder {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = false;

    public void setToDebug() {
        debug = true;
    }

    /**
     * sum the intensity of the points with an option to subtract the mean.
     * @param img
     * @param points
     * @return
     */
    protected double sumIntensity(GreyscaleImage img, Set<PairInt> points) {
        if (img == null) {
            throw new IllegalStateException("img cannot be null");
        }
        if (points == null) {
            throw new IllegalStateException("points cannot be null");
        }
        double sum = 0;
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            double v = img.getValue(x, y);
            sum += v;
        }
        return sum;
    }

    public TransformationParameters solveForScale(
        SegmentedImageHelper img1Helper, SegmentationType type1,
        boolean useBinned1,
        SegmentedImageHelper img2Helper, SegmentationType type2,
        boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {

        BlobsAndContours bc1 = img1Helper.getBlobsAndContours(type1, useBinned1);
        GreyscaleImage img1 = img1Helper.getGreyscaleImage(useBinned1);

        BlobsAndContours bc2 = img2Helper.getBlobsAndContours(type2, useBinned2);
        GreyscaleImage img2 = img2Helper.getGreyscaleImage(useBinned2);

        List<List<CurvatureScaleSpaceContour>> contours1List = bc1.getContours();
        List<List<CurvatureScaleSpaceContour>> contours2List = bc2.getContours();
        List<Set<PairInt>> blobs1 = bc1.getBlobs();
        List<Set<PairInt>> blobs2 = bc2.getBlobs();
        List<PairIntArray> perimeters1 = bc1.getBlobOrderedPerimeters();
        List<PairIntArray> perimeters2 = bc2.getBlobOrderedPerimeters();

        Map<PairInt, CSSContourMatcherWrapper> singleSolnMap =
            new HashMap<PairInt,  CSSContourMatcherWrapper>();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        /*
        Note that the matching contours are improved and filtered by using 
        feature descriptors.
        
        Finding the best match if any for each index1.
        If a subsequent match to the that index1's index2 is present, it will
        only be considered if the cost is less than the previous.
        It's a greedy best approach.
        After all of the index1 solutions are gathered, the best is found
        by adjusting the costs by the largest sigma from the first peak contours 
        from edges1 (i.e. applying the penalty from the paper to smaller peak
        sigma's costs).
        The best agreeing solutions are kept.
        
        If instead wanted a bipartite matching of best costs, would want to
        keep the top 2 or so of each index1.
        Would need to be careful to not include large peak contour false matches.
        All costs in the top would need to be adjusted by the penalty mentioned
        above.
        Then the best among more than one would be defined by cost.  
        Then more than one solutions points would be kept if similarity
        in difference of theta and scale were within a limit.
        */
       
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            index1BestMap = new HashMap<Integer, 
            FixedSizeSortedVector<IntensityFeatureComparisonStats>>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            if (contours1List.get(idx1).isEmpty()) {
                continue;
            }

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = perimeters1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            // keeping the top '2' for each index1.  comparison is by cost.
            // choosing more than one because later bipartite matching attempts
            // to match best for all index1 matchings
            FixedSizeSortedVector<IntensityFeatureComparisonStats> bestStats 
                = new FixedSizeSortedVector<IntensityFeatureComparisonStats>(2, 
                IntensityFeatureComparisonStats.class);

            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                if (contours2List.get(idx2).isEmpty()) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = perimeters2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);
                
double[] xyCen2 = curveHelper.calculateXYCentroids(curve2);
log.info("index1=" + index1.toString() + " index2=" + index2.toString()
 + " xyCen1=" + Arrays.toString(xyCen1) + " xyCen2=" + Arrays.toString(xyCen2));

                CurvatureScaleSpaceInflectionSingleEdgeMapper mapper =
                    new CurvatureScaleSpaceInflectionSingleEdgeMapper(
                    0, 0, 0, 0);

                TransformationParameters params = mapper.matchContours(
                    contours1List.get(idx1), contours2List.get(idx2));

                if ((params == null) ||
                    (mapper.getMatcher().getSolutionMatchedContours1().size() < 2)) {

                    if ((mapper.getMatcher() != null) &&
                        (mapper.getMatcher().getSolutionMatchedContours1().size() == 1)) {

                        singleSolnMap.put(new PairInt(idx1, idx2), mapper.getMatcher());
                    }

if (mapper.getMatcher() != null) {
log.info(String.format("discarding [%d] (%d,%d)  [%d] (%d,%d)  nMCs=%d",
idx1, (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
idx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
mapper.getMatcher().getSolutionMatchedContours1().size()));
} else {
log.info(String.format("discarding [%d] (%d,%d)  [%d] (%d,%d)",
idx1, (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
idx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1])));
}

                    continue;
                }

                // edit points using feature descriptors and remove outliers:
                List<FeatureComparisonStat> compStats =
                    filterContourPointsByFeatures(img1, img2, index1, index2,
                    blob1, blob2, curve1, curve2, features1, features2,
                    mapper.getMatcher());

                if (compStats.size() < 2) {
                    continue;
                }

                //TODO: consider moving this type of statistic into the
                //cost during contour matching.  wanting to avoid accepting
                //solutions which are a small number of spurious matches due
                //to one contour having many points to match to.
                int nc1 = contours1List.get(index1.intValue()).size();
                int nc2 = contours2List.get(index2.intValue()).size();
                float frac = (float)nc1/(float)nc2;
                boolean lgDiffN = ((nc1 > nc2) && frac > 2)
                    || ((nc1 < nc2) && frac < 0.5);

                if (lgDiffN) {                    
                    log.info(String.format(
                        "discarding a good match because frac of maxMatchable is low." 
                         + " nc1=%d nc2=%d (frac=%.2f) nMCs=%d",
                        nc1, nc2, frac, 
                        mapper.getMatcher().getSolutionMatchedContours1().size()));
                    continue;
                }

                //double combinedStat = calculateCombinedIntensityStat(compStats);

                IntensityFeatureComparisonStats stats = new 
                    IntensityFeatureComparisonStats(index1.intValue(), 
                    index2.intValue(), mapper.getMatcher().getSolvedCost(), 
                        mapper.getMatcher().getSolvedScale());
                stats.addAll(compStats);
                
                // bestStats keeps the top '2' smallest cost solutions added to it
                boolean added = bestStats.add(stats);
                
                if (added) {
                    log.info("  added to best for [" + index1.toString() + "] ["
                        + index2.toString() + "] cost=" + stats.getCost()
                        + " with n=" + stats.getComparisonStats().size());
                }
            }

            if (bestStats.getNumberOfItems() == 0) {
                continue;
            }

            if (debug) {
                for (int k = 0; k < bestStats.getNumberOfItems(); ++k) {
                    IntensityFeatureComparisonStats stats = bestStats.getArray()[k];
                    double[] xyCen2 = curveHelper.calculateXYCentroids(
                        blobs2.get(stats.getIndex2()));
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format(
                        "==>[%d](%d,%d) [%d](%d,%d) cost=%.2f scale=%.2f nMatched=%d(%d,%d) intSqDiff=%.1f",
                        stats.getIndex1(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
                        stats.getIndex2(), (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
                        (float)stats.getCost(), (float)stats.getScale(), 
                        stats.getComparisonStats().size(),
                        contours1List.get(stats.getIndex1()).size(), 
                        calculateCombinedIntensityStat(stats.getComparisonStats())));
                    log.info(sb.toString());
                }
            }

            index1BestMap.put(index1, bestStats);
        }
        
        List<FeatureComparisonStat> bestOverall = null;
        if (index1BestMap.isEmpty()) {
            
            // -------- process the single solution compStats ------------
            
            if (singleSolnMap.size() > 1) {
                bestOverall = processSingleSolutionsIfNoBest(img1, img2,
                    singleSolnMap, blobs1, blobs2, perimeters1, perimeters2,
                    features1, features2);
            }

            if (index1BestMap.isEmpty()) {
                return null;
            }
        } else {
            bestOverall = filterToBestConsistent(index1BestMap, contours1List,
                contours2List);
        }

        if (bestOverall == null) {
            return null;
        }

        TransformationParameters params = calculateTransformation(
            img1Helper, type1, useBinned1,
            img2Helper, type2, useBinned2,
            bestOverall, outputScaleRotTransXYStDev);

        if (params == null) {
            return null;
        }

        return params;
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
    protected BlobPerimeterRegion extractBlobPerimeterRegion(
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
            xPrev, yPrev, xm, ym, xNext, yNext, blob);

        return region;
    }

    protected List<FeatureComparisonStat> filterContourPointsByFeatures(
        GreyscaleImage img1, GreyscaleImage img2,
        Integer index1, Integer index2,
        Set<PairInt> blob1, Set<PairInt> blob2,
        PairIntArray curve1, PairIntArray curve2,
        IntensityFeatures features1, IntensityFeatures features2,
        CSSContourMatcherWrapper matcher) {

        FeatureMatcher featureMatcher = new FeatureMatcher();

        int dither = 1;

        List<FeatureComparisonStat> compStats = new ArrayList<FeatureComparisonStat>();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[] xyCen1 = curveHelper.calculateXYCentroids(curve1);
        double[] xyCen2 = curveHelper.calculateXYCentroids(curve2);

        double statSqSum = 0;

        int nMaxMatchable = matcher.getNMaxMatchable();
        int nMaxStats = 0;
        int nStats = 0;

        double cost = matcher.getSolvedCost();

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(
            "[%d](%d,%d) [%d](%d,%d) cost=%.1f scale=%.2f  nMatched=%d ",
            index1.intValue(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
            index2.intValue(), (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
            (float)cost, (float)matcher.getSolvedScale(),
            matcher.getSolutionMatchedContours1().size()));

        FeatureComparisonStat bestCompStat = null;

        for (int j = 0; j < matcher.getSolutionMatchedContours1().size(); ++j) {

            CurvatureScaleSpaceContour c1 = matcher.getSolutionMatchedContours1().get(j);
            CurvatureScaleSpaceContour c2 = matcher.getSolutionMatchedContours2().get(j);

            // the sizes of the peak details will be the same
            CurvatureScaleSpaceImagePoint[] details1 = c1.getPeakDetails();
            CurvatureScaleSpaceImagePoint[] details2 = c2.getPeakDetails();

            nMaxStats += details1.length;

            for (int jj = 0; jj < details1.length; ++jj) {

                BlobPerimeterRegion region1 = extractBlobPerimeterRegion(
                    index1.intValue(), details1[jj], curve1, blob1
                );

                BlobPerimeterRegion region2 = extractBlobPerimeterRegion(
                    index2.intValue(), details2[jj], curve2, blob2
                );

                FeatureComparisonStat compStat =
                    featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, region1, region2, dither);

sb.append(
    String.format(" (%d,%d) theta1=%d   (%d,%d) theta2=%d",
    region1.getX(), region1.getY(),
    Math.round(region1.getRelativeOrientationInDegrees()),
    region2.getX(), region2.getY(),
    Math.round(region2.getRelativeOrientationInDegrees())));

                if (compStat != null) {
                    float sumIntSqDiff = compStat.getSumIntensitySqDiff();
                    float intErrDiff = compStat.getImg2PointIntensityErr();
                    if (sumIntSqDiff < intErrDiff) {
                        if (bestCompStat == null) {
                            bestCompStat = compStat;
                        } else {
                            if (sumIntSqDiff < bestCompStat.getSumIntensitySqDiff()) {
                                bestCompStat = compStat;
                            }
                        }
                        statSqSum += (sumIntSqDiff*sumIntSqDiff);
                        sb.append(String.format("  %.1f(%.1f), ", sumIntSqDiff, intErrDiff));
                        compStats.add(compStat);
                        nStats++;
                    }
                }
            } // end details
        }// end matching contours for index1, index2

        if (bestCompStat == null) {
            return compStats;
        }

        log.info(sb.toString());

// looking at idx=2, 7 for full image and idx=0, 0 for binned
//if ((index1.intValue() == 2) && (index2.intValue() == 6)) {
if (((index1.intValue() == 0) && (index2.intValue() == 0)) ||
((index1.intValue() == 3) && (index2.intValue() == 3))) {
try {
ImageExt img1C = img1.createColorGreyscaleExt();
ImageExt img2C = img2.createColorGreyscaleExt();
for (FeatureComparisonStat compStat : compStats) {
    int x1 = compStat.getImg1Point().getX();
    int y1 = compStat.getImg1Point().getY();
    int x2 = compStat.getImg2Point().getX();
    int y2 = compStat.getImg2Point().getY();
    img1C.setRGB(x1, y1, 255, 0, 0);
    img2C.setRGB(x2, y2, 255, 0, 0);
}
String bin = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(bin + "/features1_before_redo.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/features2_before_redo.png", img2C);
int z = 1;
} catch(IOException e) {
}
}

        // if bestCompStat's difference in orientation is different than the
        // others', re-do the others to see if have an improved calculation.
        // the "re-do" should try a dither of 1 or 2
        float bestDiffTheta = AngleUtil.getAngleDifference(
            bestCompStat.getImg2PointRotInDegrees(),
            bestCompStat.getImg1PointRotInDegrees());
        if (bestDiffTheta < 0) {
            bestDiffTheta += 360;
        }

        boolean redoStats = false;
        for (FeatureComparisonStat cStat : compStats) {
            float diffTheta = AngleUtil.getAngleDifference(
                cStat.getImg2PointRotInDegrees(),
                cStat.getImg1PointRotInDegrees());
            if (diffTheta < 0) {
                diffTheta += 360;
            }
            if (Math.abs(bestDiffTheta - diffTheta) > 25) {
                redoStats = true;
                break;
            }
        }

        //TODO: may need to consider re-doing if compStats.size() is << nMaxStats too
redoStats = true;
        if (redoStats) {

            compStats = redoFilterContourPointsByFeatures(img1, img2, index1,
                index2, blob1, blob2, curve1, curve2,
                bestCompStat.getImg1PointRotInDegrees(),
                bestCompStat.getImg2PointRotInDegrees(), matcher);
            
            log.info("redone: " + printToString(compStats) + " combinedStat="
                + calculateCombinedIntensityStat(compStats));
            
            removeDiscrepantThetaDiff(compStats);

            log.info("theta diff filtered: " + printToString(compStats) + " combinedStat="
                + calculateCombinedIntensityStat(compStats));
        }

        return compStats;
    }

    List<FeatureComparisonStat> redoFilterContourPointsByFeatures(
        GreyscaleImage img1, GreyscaleImage img2,
        Integer index1, Integer index2,
        Set<PairInt> blob1, Set<PairInt> blob2,
        PairIntArray curve1, PairIntArray curve2,
        float theta1, float theta2,
        CSSContourMatcherWrapper matcher) {

        FeatureMatcher featureMatcher = new FeatureMatcher();

        // for the redo, because the orientations are set rather than found,
        // not going to re-use the invoker's instance of Features
        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        int dither = 2;

        List<FeatureComparisonStat> compStats = new ArrayList<FeatureComparisonStat>();

        int nMaxMatchable = matcher.getNMaxMatchable();
        int nMaxStats = 0;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("[%d] [%d] ", index1.intValue(), index2.intValue()));

        float theta1Radians = (float)(theta1 * Math.PI/180.);
        float theta2Radians = (float)(theta2 * Math.PI/180.);

        for (int j = 0; j < matcher.getSolutionMatchedContours1().size(); ++j) {

            CurvatureScaleSpaceContour c1 = matcher.getSolutionMatchedContours1().get(j);
            CurvatureScaleSpaceContour c2 = matcher.getSolutionMatchedContours2().get(j);

            // the sizes of the peak details will be the same
            CurvatureScaleSpaceImagePoint[] details1 = c1.getPeakDetails();
            CurvatureScaleSpaceImagePoint[] details2 = c2.getPeakDetails();

            nMaxStats += details1.length;

            for (int jj = 0; jj < details1.length; ++jj) {

                BlobPerimeterRegion region1 = extractBlobPerimeterRegion(
                    index1.intValue(), details1[jj], curve1, blob1
                );
                region1.overrideRelativeOrientation(theta1Radians);

                BlobPerimeterRegion region2 = extractBlobPerimeterRegion(
                    index2.intValue(), details2[jj], curve2, blob2
                );
                region2.overrideRelativeOrientation(theta2Radians);

                FeatureComparisonStat compStat =
                    featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, region1, region2, dither);

                if (compStat != null) {

                    float sumIntSqDiff = compStat.getSumIntensitySqDiff();
                    float intErrDiff = compStat.getImg2PointIntensityErr();

                    if (sumIntSqDiff < intErrDiff) {
                        compStats.add(compStat);
                    }

                }
            } // end details
        }// end matching contours for index1, index2

        sb.append(printToString(compStats)).append(" nMaxStats=").append(nMaxStats);

        log.info(sb.toString());

// looking at idx=2, 7 for full image and idx=0, 0 for binned
//if ((index1.intValue() == 2) && (index2.intValue() == 6)) {
if (((index1.intValue() == 0) && (index2.intValue() == 0)) ||
((index1.intValue() == 3) && (index2.intValue() == 3))) {
try {
ImageExt img1C = img1.createColorGreyscaleExt();
ImageExt img2C = img2.createColorGreyscaleExt();
MiscDebug.debugPlot(matcher.getSolutionMatchedContours1(), img1C,
    0, 0, "1_" + index1 + "_redone");
MiscDebug.debugPlot(matcher.getSolutionMatchedContours2(), img2C,
    0, 0, "2_" + index2 + "_redone");
for (FeatureComparisonStat compStat : compStats) {
    int x1 = compStat.getImg1Point().getX();
    int y1 = compStat.getImg1Point().getY();
    int x2 = compStat.getImg2Point().getX();
    int y2 = compStat.getImg2Point().getY();
    img1C.setRGB(x1, y1, 0, 255, 0);
    img2C.setRGB(x2, y2, 0, 255, 0);
}
String bin = ResourceFinder.findDirectory("bin");
ImageIOHelper.writeOutputImage(bin + "/contours_1_redone.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/contours_2_redone.png", img2C);
int z = 1;
} catch(IOException e) {
}
}

        removeOutliers(compStats);

        return compStats;
    }

    protected String printToString(List<FeatureComparisonStat> compStats) {

        StringBuilder sb = new StringBuilder();

        for (FeatureComparisonStat compStat : compStats) {

            sb.append(String.format(
                " (%d,%d) (%d,%d) theta1=%.0f theta2=%.0f intSqDiff=%.1f(%.1f)",
                compStat.getImg1Point().getX(), compStat.getImg1Point().getY(),
                compStat.getImg2Point().getX(), compStat.getImg2Point().getY(),
                compStat.getImg1PointRotInDegrees(), compStat.getImg2PointRotInDegrees(),
                compStat.getSumIntensitySqDiff(),
                compStat.getImg2PointIntensityErr()));
        }

        return sb.toString();
    }

    protected double calculateCombinedIntensityStat(
        List<FeatureComparisonStat> compStats) {

        double sum = 0;

        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }

        sum /= (double)compStats.size();

        return sum;
    }

    protected TransformationParameters calculateTransformation(
        SegmentedImageHelper img1Helper, SegmentationType type1,
        boolean useBinned1,
        SegmentedImageHelper img2Helper, SegmentationType type2,
        boolean useBinned2,
        List<FeatureComparisonStat> compStats,
        float[] outputScaleRotTransXYStDev) {

        assert(compStats.isEmpty() == false);

        int binFactor1 = img1Helper.getBinFactor(useBinned1);

        int binFactor2 = img2Helper.getBinFactor(useBinned2);


        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        int centroidX1 = 0;
        int centroidY1 = 0;

        PairIntArray matchedXY1 = new PairIntArray();
        PairIntArray matchedXY2 = new PairIntArray();

        float[] weights = new float[compStats.size()];

        double sum = 0;

        for (int i = 0; i < compStats.size(); ++i) {

            FeatureComparisonStat compStat = compStats.get(i);

            int x1 = compStat.getImg1Point().getX() * binFactor1;
            int y1 = compStat.getImg1Point().getY() * binFactor1;
            matchedXY1.add(x1, y1);

            int x2 = compStat.getImg2Point().getX() * binFactor2;
            int y2 = compStat.getImg2Point().getY() * binFactor2;
            matchedXY2.add(x2, y2);

            weights[i] = compStat.getSumIntensitySqDiff();

            sum += weights[i];
        }

        double tot = 0;
        for (int i = 0; i < compStats.size(); ++i) {
            double div = (sum - weights[i])/((compStats.size() - 1) * sum);
            weights[i] = (float)div;
            tot += div;
        }

        assert(Math.abs(tot - 1.) < 0.03);

        TransformationParameters params = tc.calulateEuclidean(
            matchedXY1, matchedXY2, weights, centroidX1, centroidY1,
            outputScaleRotTransXYStDev);

        return params;
    }

    protected void removeOutliers(List<FeatureComparisonStat> compStats) {

        if (compStats.size() < 2) {
            return;
        }

        //TODO: improve w/ a more robust outlier removal

        double[] errDivInt = new double[compStats.size()];

        float[] weights = new float[compStats.size()];
        double sum = 0;
        for (int i = 0; i < compStats.size(); ++i) {
            FeatureComparisonStat compStat = compStats.get(i);
            weights[i] = compStat.getSumIntensitySqDiff();
            sum += weights[i];
            errDivInt[i] = compStat.getImg2PointIntensityErr()/compStat.getSumIntensitySqDiff();
        }
        for (int i = 0; i < compStats.size(); ++i) {
            double div = (sum - weights[i]) / ((compStats.size() - 1) * sum);
            weights[i] = (float) div;
        }
        float[] wghtsMeanAndStDev = MiscMath.getAvgAndStDev(weights);
        float maxWeight = MiscMath.findMax(weights);

        /*
        if all stats have intensities < 5 times their errors and
        if the stdev is approx 0.15 times the mean or less, should filter here
        */
        boolean doNotFilter = true;
        if ((wghtsMeanAndStDev[1]/wghtsMeanAndStDev[0]) > 0.15) {
            doNotFilter = false;
        }
        if (doNotFilter) {
            for (int i = 0; i < errDivInt.length; ++i) {
                if (errDivInt[i] < 5) {
                    doNotFilter = false;
                    break;
                }
            }
        }
        //TODO: may need revision
        if (doNotFilter && (weights.length <= 4)) {
            return;
        }

        List<FeatureComparisonStat> filteredCompStats = new ArrayList<FeatureComparisonStat>();
        for (int i = 0; i < compStats.size(); ++i) {
            float w = weights[i];
            float diffW = Math.abs(maxWeight - w);
            if (diffW < (1.5 * wghtsMeanAndStDev[1])) {
                filteredCompStats.add(compStats.get(i));
            }
        }

        compStats.clear();

        compStats.addAll(filteredCompStats);
    }

    protected List<FeatureComparisonStat> filterToBestConsistent(
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
        index1StatsMap, List<List<CurvatureScaleSpaceContour>> contours1Lists,
        List<List<CurvatureScaleSpaceContour>> contours2Lists) {

        if (index1StatsMap == null || index1StatsMap.isEmpty()) {
            return null;
        }
        
        if (index1StatsMap.size() > 1) {
            // make corrections for cost between different edges            
            correctCostsUsingMaxSigma(index1StatsMap, contours1Lists);            
        }
        
        TreeMap<Double, List<IntensityFeatureComparisonStats>> bestMatches = 
            findBestMatchesUsingBipartite(index1StatsMap, 
            contours1Lists.size(), contours2Lists.size());

        int nMaxMatchable = 2 * Math.max(contours1Lists.size(), contours2Lists.size());
     
        Map<Integer, IntensityFeatureComparisonStats> index1Map = new
            HashMap<Integer, IntensityFeatureComparisonStats>();
            
        Map<PairInt, PairFloat> indexesDiffTheta = new HashMap<PairInt, PairFloat>();
        
        /* calculate the highest number of similar transformations and the
        lowest cost from those.  store nSimilar, indexes, cost for each iteration*/
        int[] nSimilarSummary = new int[nMaxMatchable];
        Integer[][] indexesSummary = new Integer[nMaxMatchable][];
        double[] costsSummary = new double[nMaxMatchable];
        int[] mainIndexSummary = new int[nMaxMatchable];

        int count = 0;
        
        for (Map.Entry<Double, List<IntensityFeatureComparisonStats>> entry : bestMatches.entrySet()) {

            Double adjustedCost = entry.getKey();
            
            List<IntensityFeatureComparisonStats> stats = entry.getValue();
            
            for (IntensityFeatureComparisonStats stat : stats) {
                int idx1 = stat.getIndex1();
                int idx2 = stat.getIndex2();
                assert(index1Map.get(Integer.valueOf(idx1)) == null);
                index1Map.put(Integer.valueOf(idx1), stat);
                PairInt p = new PairInt(idx1, idx2);
                PairFloat diffThetaMnStdv = indexesDiffTheta.get(p);
                if (diffThetaMnStdv == null) {
                    float[] diffThetas = new float[stat.getComparisonStats().size()];
                    for (int i = 0; i < stat.getComparisonStats().size(); ++i) {
                        FeatureComparisonStat fcs = stat.getComparisonStats().get(i);                
                        float diff = AngleUtil.getAngleDifference(
                           fcs.getImg1PointRotInDegrees(), fcs.getImg2PointRotInDegrees());
                        diffThetas[i] = diff;
                    } 
                    float[] msv = MiscMath.getAvgAndStDev(diffThetas);
                    
                    // TODO: consider whether need to exclude points further from
                    // the mean than a difference in thetaDiff of 20
                    // unless edits have changed it, this was run previously,
                    // so what is present here should already be consistent
                    // diffThetas:
                    //removeDiscrepantThetaDiff(stat.getComparisonStats());
                    
                    diffThetaMnStdv = new PairFloat(msv[0], msv[1]);
                }
                
                double scale = stat.getScale();
                
                Set<Integer> similarParamsIndexes1 = new HashSet<Integer>();
                similarParamsIndexes1.add(Integer.valueOf(idx1));
                
                // count the number of solutions similar in diffTheta and scale
                for (Map.Entry<Double, List<IntensityFeatureComparisonStats>> entry2 : bestMatches.entrySet()) {
                    Double adjustedCost2 = entry2.getKey();
                    if (adjustedCost2.equals(adjustedCost)) {
                        continue;
                    }
                    for (IntensityFeatureComparisonStats stat2 : entry2.getValue()) {
                        int idx1P = stat2.getIndex1();
                        int idx2P = stat2.getIndex2();
                        // bipartite matching should have made unique idx1 already:
                        assert(idx1 != idx1P);
                        PairInt p2 = new PairInt(idx1P, idx2P);
                        PairFloat diffThetaMnStdv2 = indexesDiffTheta.get(p2);
                        if (diffThetaMnStdv2 == null) {
                            float[] diffThetas = new float[stat2.getComparisonStats().size()];
                            for (int i = 0; i < stat2.getComparisonStats().size(); ++i) {
                                FeatureComparisonStat fcs = stat2.getComparisonStats().get(i);                
                                float diff = AngleUtil.getAngleDifference(
                                    fcs.getImg1PointRotInDegrees(), fcs.getImg2PointRotInDegrees());
                                diffThetas[i] = diff;
                            }
                            float[] msv = MiscMath.getAvgAndStDev(diffThetas);
                            diffThetaMnStdv2 = new PairFloat(msv[0], msv[1]);
                        }
                        
                        //--- compare diffTheta and scale with ---
                        if ((Math.abs(diffThetaMnStdv.getX() - diffThetaMnStdv2.getX()) < 10)
                            && (Math.abs(scale - stat2.getScale()) < 0.1)) {
                            similarParamsIndexes1.add(Integer.valueOf(idx1P));
                        }
                    }
                }
                
                // store for sort and combine
                nSimilarSummary[count] = similarParamsIndexes1.size();
                indexesSummary[count] = similarParamsIndexes1.toArray(new Integer[similarParamsIndexes1.size()]);
                costsSummary[count] = adjustedCost;
                mainIndexSummary[count] = idx1;
                
                count++;
            }
        }
        
        nSimilarSummary = Arrays.copyOf(nSimilarSummary, count);
        indexesSummary = Arrays.copyOf(indexesSummary, count);
        costsSummary = Arrays.copyOf(costsSummary, count);
        mainIndexSummary = Arrays.copyOf(mainIndexSummary, count);
        
        // sort to prefer the solution w/ largest number of similar solutions:
        //MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary, indexesSummary, mainIndexSummary);

        // OR sort to refer the solution w/ best cost and any it is similar to:
        //--- these are still sorted by costs already, so no need to resort again ---
        // MultiArrayMergeSort.sortBy1stAscThen2ndDesc(costsSummary, nSimilarSummary, indexesSummary, mainIndexSummary, 0, costsSummary.length - 1);

        Integer[] indexes = indexesSummary[0];
       
        List<FeatureComparisonStat> output = new ArrayList<FeatureComparisonStat>();
        for (Integer index1 : indexes) {
            IntensityFeatureComparisonStats stats = index1Map.get(index1);
            output.addAll(stats.getComparisonStats());
        }
        return output;
    }
    
    private List<FeatureComparisonStat> processSingleSolutionsIfNoBest(
        GreyscaleImage img1, GreyscaleImage img2,
        Map<PairInt, CSSContourMatcherWrapper> singleSolnMap,
        List<Set<PairInt>> blobs1, List<Set<PairInt>> blobs2,
        List<PairIntArray> perimeters1, List<PairIntArray> perimeters2,
        IntensityFeatures features1, IntensityFeatures features2) {

log.info("WARNING: processing single solutions... may remove these in future");

        Map<PairInt, List<FeatureComparisonStat>> compStatMap =
            new HashMap<PairInt, List<FeatureComparisonStat>>();

        TreeSet<Integer> sIndexes1Set = new TreeSet<Integer>();
        TreeSet<Integer> sIndexes2Set = new TreeSet<Integer>();
        for (Entry<PairInt, CSSContourMatcherWrapper> entry : singleSolnMap.entrySet()) {
            PairInt p = entry.getKey();
            int idx1 = p.getX();
            int idx2 = p.getY();
            sIndexes1Set.add(Integer.valueOf(idx1));
            sIndexes2Set.add(Integer.valueOf(idx2));
        }
        float[][] cost = new float[sIndexes1Set.size()][];
        Map<Integer, Integer> sIndexes1 = new HashMap<Integer, Integer>();
        Map<Integer, Integer> sIndexes2 = new HashMap<Integer, Integer>();
        int i = 0;
        for (Integer sIndex1 : sIndexes1Set) {
            cost[i] = new float[sIndexes2Set.size()];
            Arrays.fill(cost[i], Float.MAX_VALUE);
            sIndexes1.put(sIndex1, Integer.valueOf(i));
            ++i;
        }
        int j = 0;
        for (Integer sIndex2 : sIndexes2Set) {
            sIndexes2.put(sIndex2, Integer.valueOf(j));
            ++j;
        }

        int nCS = 0;

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        for (Entry<PairInt, CSSContourMatcherWrapper> entry : singleSolnMap.entrySet()) {

            PairInt p = entry.getKey();
            int idx1 = p.getX();
            int idx2 = p.getY();

            CSSContourMatcherWrapper matcher = entry.getValue();

            PairIntArray curve1 = perimeters1.get(idx1);
            Set<PairInt> blob1 = blobs1.get(idx1);
            PairIntArray curve2 = perimeters2.get(idx2);
            Set<PairInt> blob2 = blobs2.get(idx2);

            List<FeatureComparisonStat> compStats =
                filterContourPointsByFeatures(img1, img2,
                Integer.valueOf(idx1), Integer.valueOf(idx2),
                blob1, blob2, curve1, curve2, features1, features2,
                matcher);

double[] xyCen1 = curveHelper.calculateXYCentroids(curve1);
double[] xyCen2 = curveHelper.calculateXYCentroids(curve2);
log.info(String.format("single solution [%d] (%d,%d)  [%d] (%d,%d)  nCS=%d",
idx1, (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
idx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
compStats.size()));

            if (compStats.isEmpty()) {
                continue;
            }

            //TODO: combined stat needs to be recalculated so that can compare
            // different matchings

            compStatMap.put(p, compStats);

            double combStat = calculateCombinedIntensityStat(compStats);

            int cIdx1 = sIndexes1.get(Integer.valueOf(idx1)).intValue();
            int cIdx2 = sIndexes2.get(Integer.valueOf(idx2)).intValue();
            cost[cIdx1][cIdx2] = (float)combStat;

            nCS++;
        }

        if (nCS > 1) {
            List<FeatureComparisonStat> csList = new ArrayList<FeatureComparisonStat>();
            HungarianAlgorithm b = new HungarianAlgorithm();
            int[][] match = b.computeAssignments(cost);
            for (int ii = 0; ii < match.length; ii++) {
                int idx1 = match[ii][0];
                int idx2 = match[ii][1];
                if (idx1 == -1 || idx2 == -1) {
                    continue;
                }
                PairInt p = new PairInt(idx1, idx2);
                List<FeatureComparisonStat> stats = compStatMap.get(p);
                csList.addAll(stats);
            }
            removeOutliers(csList);
            if (csList.size() > 1) {
                return csList;
            }
        }

        return null;
    }

    private void removeDiscrepantThetaDiff(List<FeatureComparisonStat> compStats) {
        
        if (compStats == null || compStats.isEmpty()) {
            return;
        }
        
        float[] values = new float[compStats.size()];
        for (int i = 0; i < compStats.size(); ++i) {
            FeatureComparisonStat stat = compStats.get(i);
            float diff = AngleUtil.getAngleDifference(
                stat.getImg1PointRotInDegrees(),
                stat.getImg2PointRotInDegrees());
            values[i] = diff;
        }
        // 20 degree wide bins
        HistogramHolder hist = Histogram.createSimpleHistogram(20.f,
            values, Errors.populateYErrorsBySqrt(values));
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        float thetaDiff;
        if (yMaxIdx == -1) {
            float[] thetaDiffMeanStDev = MiscMath.getAvgAndStDev(values);
            thetaDiff = thetaDiffMeanStDev[0];
        } else {
            thetaDiff = hist.getXHist()[yMaxIdx];
        }        
        
        //TODO: consider a bin larger than 20 degrees... 25
        List<Integer> remove = new ArrayList<Integer>();
        for (int i = 0; i < values.length; ++i) {
            float diffRot = Math.abs(values[i] - thetaDiff);
            if (diffRot > 20) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i);
            compStats.remove(idx);
        }
    }

    private void correctCostsUsingMaxSigma(
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            index1StatsMap, List<List<CurvatureScaleSpaceContour>> contours1Lists) {
                
        float maxPeakSigma = Float.MIN_VALUE;
        
        Map<Integer, Float> index1PeakSigma = new HashMap<Integer, Float>();
        
        for (Entry<Integer, FixedSizeSortedVector<
            IntensityFeatureComparisonStats>> entry : index1StatsMap.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            List<CurvatureScaleSpaceContour> contours = contours1Lists.get(index1.intValue());
           
            // contours are ordered, so just need the first.
            if (contours.isEmpty()) {
                continue;
            }
            float peakSigma = contours.get(0).getPeakSigma();
            if (peakSigma > maxPeakSigma) {
                maxPeakSigma = peakSigma;
            }
            // contour 6 is not a strong blob by has maxsigma=166 which is >4*expected from 0
            index1PeakSigma.put(index1, Float.valueOf(peakSigma));
        }
        
        for (Entry<Integer, Float> entry : index1PeakSigma.entrySet()) {
            
            Integer index1 = entry.getKey();
            double penalty = maxPeakSigma - entry.getValue().floatValue();
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> index1Stats =
                index1StatsMap.get(index1);
            
            for (int i = 0; i < index1Stats.getNumberOfItems(); ++i) {
                
                IntensityFeatureComparisonStats stats = index1Stats.getArray()[i];
                                
                stats.setAdjustedCost(stats.getCost() + penalty);
            }
        }
    }

    private TreeMap<Double, List<IntensityFeatureComparisonStats>> 
    findBestMatchesUsingBipartite(Map<Integer, 
        FixedSizeSortedVector<IntensityFeatureComparisonStats>> index1StatsMap,
        int n1, int n2) {
        
        TreeMap<Double, List<IntensityFeatureComparisonStats>> matched = new
            TreeMap<Double, List<IntensityFeatureComparisonStats>>();
        
        int maxMatchable = 2 * Math.max(n1, n2);
                
        double[] costs = new double[maxMatchable];
        int[] idx1s = new int[maxMatchable];
        int[] idx2s = new int[maxMatchable];
        IntensityFeatureComparisonStats[] ics = new IntensityFeatureComparisonStats[maxMatchable];
        
        int count = 0;
        for (Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>>
            entry : index1StatsMap.entrySet()) {
            Integer index1 = entry.getKey();            
            for (IntensityFeatureComparisonStats stats : entry.getValue().getArray()) {
                if (stats == null) {
                    continue;
                }
                assert(index1.intValue() == stats.getIndex1());
                costs[count] = stats.getAdjustedCost();
                idx1s[count] = stats.getIndex1();
                idx2s[count] = stats.getIndex2();
                ics[count] = stats;
                count++;
            }
        }
        costs = Arrays.copyOf(costs, count);
        idx1s = Arrays.copyOf(idx1s, count);
        idx2s = Arrays.copyOf(idx2s, count);
        ics = Arrays.copyOf(ics, count);
        int[] lookupIndexes = new int[count];
        for (int i = 0; i < count; ++i) {
            lookupIndexes[i] = i;
        }
        
        MultiArrayMergeSort.sortByDecr(costs, lookupIndexes);
        // best answers at highest indexes
        
        Set<Integer> chosen1 = new HashSet<Integer>();
        Set<Integer> chosen2 = new HashSet<Integer>();
        
        for (int i = (count - 1); i > -1; --i) {
            
            int idx0 = lookupIndexes[i];
            
            Integer index1 = Integer.valueOf(idx1s[idx0]);
            Integer index2 = Integer.valueOf(idx2s[idx0]);
            
            if (chosen1.contains(index1)) {
                continue;
            }
            if (chosen2.contains(index2)) {
                continue;
            }
            
            Double cost = Double.valueOf(costs[i]);
            IntensityFeatureComparisonStats ifcs = ics[idx0];
                       
            assert(Math.abs(cost.doubleValue() - ifcs.getAdjustedCost()) < 0.01);
            
            /*
            TreeMap<Double, List<IntensityFeatureComparisonStats>> matched
            */
            List<IntensityFeatureComparisonStats> ifcsList = matched.get(cost);
            if (ifcsList == null) {
                ifcsList = new ArrayList<IntensityFeatureComparisonStats>();
                matched.put(cost, ifcsList);
            }
            ifcsList.add(ifcs);
                
            chosen1.add(index1);
            chosen2.add(index2);
        }
        
        return matched;
    }

}
