package algorithms.imageProcessing;

import algorithms.imageProcessing.SegmentedImageHelper.SegmentationType;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
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

        FixedSizeSortedVector<IntensityFeatureComparisonStats> topForIndex1 =
            new FixedSizeSortedVector<>(2, IntensityFeatureComparisonStats.class);

        //TODO: this section needs to be simplified

        double bestOverallStatSqSum = Double.MAX_VALUE;
        int bestOverallIdx1 = -1;
        int bestOverallC1 = 0;
        int bestOverallIdx2 = -1;
        int bestOverallC2 = 0;
        int bestOverallNMatched = -1;
        List<FeatureComparisonStat> bestOverallCompStats = null;

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            if (contours1List.get(idx1).isEmpty()) {
                continue;
            }

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = perimeters1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            double bestStatSqSum = Double.MAX_VALUE;
            double bestScale = -1;
 //double bestCost consider comparisons by cost between matches
            int bestIdx2 = -1;
            int bestC2 = 0;
            int bestNMatched = -1;
            List<FeatureComparisonStat> bestCompStats = null;

            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                if (contours2List.get(idx2).isEmpty()) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = perimeters2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);

//log.info("index1=" + index1.toString() + " index2=" + index2.toString());

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

double[] xyCen2 = curveHelper.calculateXYCentroids(curve2);
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
                    continue;
                }

                double combinedStat = calculateCombinedIntensityStat(compStats);

                //TODO: consider keeping top k instead of 1
                if (combinedStat < bestStatSqSum) {
                    bestStatSqSum = combinedStat;
                    bestIdx2 = index2.intValue();
                    bestC2 = contours2List.get(bestIdx2).size();
                    bestNMatched = mapper.getMatcher().getSolutionMatchedContours1().size();
                    bestCompStats = compStats;
                    bestScale = mapper.getMatcher().getSolvedScale();
                    log.info("  new best for [" + index1.toString() + "] ["
                        + index2.toString() + "] combinedStat=" + combinedStat 
                        + " with n=" + bestCompStats.size());
                }
            }

            if (bestCompStats == null) {
                continue;
            }

            double[] xyCen2 = curveHelper.calculateXYCentroids(blobs2.get(bestIdx2));

            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
                "==>[%d](%d,%d) [%d](%d,%d) scale=%.2f  nMatched=%d(%d,%d) intSqDiff=%.1f",
                index1.intValue(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
                bestIdx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
                bestScale, bestNMatched,
                contours1List.get(index1.intValue()).size(), bestC2,
                (float)bestStatSqSum));
            log.info(sb.toString());

            IntensityFeatureComparisonStats stats = new IntensityFeatureComparisonStats();
            stats.addAll(bestCompStats);
            topForIndex1.add(stats);

            if (bestStatSqSum < bestOverallStatSqSum) {

                bestOverallStatSqSum = bestStatSqSum;
                bestOverallIdx1 = index1.intValue();
                bestOverallIdx2 = bestIdx2;
                bestOverallC1 = contours1List.get(bestOverallIdx1).size();
                bestOverallC2 = bestC2;
                bestOverallNMatched = bestNMatched;
                bestOverallCompStats = bestCompStats;

                log.info("  best overall for [" + bestOverallIdx1 + "] [" +
                    bestOverallIdx2 + "]  combinedStat=" + bestOverallStatSqSum);
            }
        }

        addSimilarToBestOverall(bestOverallCompStats, topForIndex1);

        // -------- process the single solution compStats ------------
        if (bestOverallCompStats == null || bestOverallCompStats.isEmpty()) {

            if (singleSolnMap.size() > 1) {
                processSingleSolutionsIfNoBest(img1, img2, bestOverallCompStats,
                    singleSolnMap, blobs1, blobs2, perimeters1, perimeters2,
                    features1, features2);
            }

            if ((bestOverallCompStats == null) || bestOverallCompStats.isEmpty()) {
                return null;
            }
        }

        TransformationParameters params = calculateTransformation(
            img1Helper, type1, useBinned1,
            img2Helper, type2, useBinned2,
            bestOverallCompStats, outputScaleRotTransXYStDev);

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
    0, 0, "1_redone");
MiscDebug.debugPlot(matcher.getSolutionMatchedContours2(), img2C,
    0, 0, "2_redone");
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

    protected float calculateRotationDifferences(List<FeatureComparisonStat>
        compStats) {

        if (compStats == null || compStats.isEmpty()) {
            return Float.POSITIVE_INFINITY;
        }

        double sumDiff = 0;

        for (FeatureComparisonStat stat : compStats) {

            float diff = AngleUtil.getAngleDifference(stat.getImg1PointRotInDegrees(),
                stat.getImg2PointRotInDegrees());

            sumDiff += diff;
        }

        return (float)(sumDiff/(double)compStats.size());
    }

    private void addSimilarToBestOverall(List<FeatureComparisonStat>
        bestOverallCompStats,
        FixedSizeSortedVector<IntensityFeatureComparisonStats> topForIndex1) {

        // ---- add to bestOverallCompStats for solutions similar to best -----
        if (bestOverallCompStats != null) {

log.info("looking for solutions similar to "+bestOverallCompStats);

            float rotationBest = calculateRotationDifferences(bestOverallCompStats);

            for (int i = 0; i < topForIndex1.getNumberOfItems(); ++i) {

                IntensityFeatureComparisonStats cStats = topForIndex1.getArray()[i];

                List<FeatureComparisonStat> stats = cStats.getComparisonStats();

                if (stats.equals(bestOverallCompStats)) {
                    continue;
                }

                boolean similar = true;

                for (FeatureComparisonStat stat : stats) {

                    float rot1 = stat.getImg1PointRotInDegrees();
                    float rot2 = stat.getImg2PointRotInDegrees();

                    float rot = AngleUtil.getAngleDifference(rot1, rot2);

                    if (Math.abs(rot - rotationBest) > 30) {
                        similar = false;
                        break;
                    }
                }

                if (similar) {
log.info("  found similar:" + stats);
                    bestOverallCompStats.addAll(stats);
                }
            }
        }
    }

    private void processSingleSolutionsIfNoBest(
        GreyscaleImage img1, GreyscaleImage img2,
        List<FeatureComparisonStat> bestOverallCompStats,
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
                bestOverallCompStats = new ArrayList<FeatureComparisonStat>();
                bestOverallCompStats.addAll(csList);
            }
        }
    }
}
