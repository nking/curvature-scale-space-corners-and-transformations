package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class AbstractBlobScaleFinder {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = false;

    public void setToDebug() {
        debug = true;
    }

    /**
     * given the list of blob points, extract the ordered boundaries of them
     * and remove any blobs for which the bounds were not extractable.
     * @param inOutBlobs
     * @param outputBounds
     * @param width
     * @param height
     * @param discardWhenCavityIsSmallerThanBorder
     */
    protected void extractBoundsOfBlobs(final GreyscaleImage img,
        final List<Set<PairInt>> inOutBlobs,
        final List<PairIntArray> outputBounds, int width, int height,
        boolean discardWhenCavityIsSmallerThanBorder) {

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

                int nChanged = 0;

                if (blobIsDarkerThanExterior(blob, closedEdge, img)) {
                    nChanged = curveHelper.adjustEdgesTowardsBrighterPixels(
                        closedEdge, img);
                } else {
                    nChanged = curveHelper.adjustEdgesTowardsDarkerPixels(
                        closedEdge, img);
                }

                if (nChanged > 0) {

                    curveHelper.removeRedundantPoints(closedEdge);

                    curveHelper.pruneAdjacentNeighborsTo2(closedEdge);

                    curveHelper.correctCheckeredSegments(closedEdge);
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

        List<Set<PairInt>> blobs = new ArrayList<Set<PairInt>>();
        List<PairIntArray> curves = new ArrayList<PairIntArray>();

        int last = (inOutBlobs.size() > 15) ? 15 : inOutBlobs.size();
        for (int i = 0; i < last; ++i) {
            int idx = indexes[i];
            blobs.add(inOutBlobs.get(idx));
            curves.add(outputBounds.get(idx));
        }

        inOutBlobs.clear();
        outputBounds.clear();

        inOutBlobs.addAll(blobs);
        outputBounds.addAll(curves);

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

    protected TransformationParameters solveForScale(
        GreyscaleImage img1, GreyscaleImage img2,
        List<Set<PairInt>> blobs1, List<Set<PairInt>> blobs2,
        List<PairIntArray> bounds1, List<PairIntArray> bounds2,
        float[] outputScaleRotTransXYStDev) {
        
        List<List<CurvatureScaleSpaceContour>> contours1List = 
            populateContours(bounds1);
        
        List<List<CurvatureScaleSpaceContour>> contours2List = 
            populateContours(bounds2);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        FixedSizeSortedVector<IntensityFeatureComparisonStats> topForIndex1 =
            new FixedSizeSortedVector<>(2, IntensityFeatureComparisonStats.class);

        double bestOverallStatSqSum = Double.MAX_VALUE;
        int bestOverallIdx1 = -1;
        int bestOverallIdx2 = -1;
        int bestOverallNMatched = -1;
        List<FeatureComparisonStat> bestOverallCompStats = null;

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = bounds1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            double bestStatSqSum = Double.MAX_VALUE;
            double bestScale = -1;
            int bestIdx2 = -1;
            int bestNMatched = -1;
            List<FeatureComparisonStat> bestCompStats = null;

            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = bounds2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);

//log.info("index1=" + index1.toString() + " index2=" + index2.toString());

//TODO: need to refactor to only calculate the
// scale space curves once each!

                CurvatureScaleSpaceInflectionSingleEdgeMapper mapper =
                    new CurvatureScaleSpaceInflectionSingleEdgeMapper(
                    0, 0, 0, 0);

                TransformationParameters params = mapper.matchContours(
                    contours1List.get(idx1), contours2List.get(idx2));

                if ((params == null) ||
                    (mapper.getMatcher().getSolutionMatchedContours1().size() < 2)) {
                    
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

                double combinedStat = calculateCombinedIntensityStat(compStats);

                //TODO: consider keeping top k instead of 1
                if (combinedStat < bestStatSqSum) {
                    bestStatSqSum = combinedStat;
                    bestIdx2 = index2.intValue();
                    bestNMatched = mapper.getMatcher().getSolutionMatchedContours1().size();
                    bestCompStats = compStats;
                    bestScale = mapper.getMatcher().getSolvedScale();
                    log.info("  new best for [" + index1.toString() + "] ["
                        + index2.toString() + "] combinedStat=" + combinedStat);
                }
            }

            if (bestCompStats == null) {
                continue;
            }

            double[] xyCen2 = curveHelper.calculateXYCentroids(blobs2.get(bestIdx2));

            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
                "==>[%d](%d,%d) [%d](%d,%d) scale=%.2f  nMatched=%d  intSqDiff=%.1f",
                index1.intValue(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
                bestIdx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
                bestScale, bestNMatched, (float)bestStatSqSum));
            log.info(sb.toString());

            IntensityFeatureComparisonStats stats = new IntensityFeatureComparisonStats();
            stats.addAll(bestCompStats);
            topForIndex1.add(stats);

            if (bestStatSqSum < bestOverallStatSqSum) {
                bestOverallStatSqSum = bestStatSqSum;
                bestOverallIdx1 = index1.intValue();
                bestOverallIdx2 = bestIdx2;
                bestOverallNMatched = bestNMatched;
                bestOverallCompStats = bestCompStats;

                log.info("  best overall for [" + bestOverallIdx1 + "] [" +
                    bestOverallIdx2 + "]  combinedStat=" + bestOverallStatSqSum);
            }
        }

        if (bestOverallCompStats == null) {
            return null;
        }

        /*
        add to bestOverallCompStats for solutions similar to best
        */
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
                bestOverallCompStats.addAll(stats);
            }
        }

        TransformationParameters params = calculateTransformation(
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
     * @param perimeterIdx
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
ImageIOHelper.writeOutputImage(bin + "/contours1_before_redo.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/contours2_before_redo.png", img2C);
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

            log.info("redone: " + printToString(compStats));
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
        List<FeatureComparisonStat> compStats,
        float[] outputScaleRotTransXYStDev) {

        assert(compStats.isEmpty() == false);

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

            int x1 = compStat.getImg1Point().getX();
            int y1 = compStat.getImg1Point().getY();
            matchedXY1.add(x1, y1);

            int x2 = compStat.getImg2Point().getX();
            int y2 = compStat.getImg2Point().getY();
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

    protected TransformationParameters evaluate(TransformationParameters[]
        params, List<PairIntArray> bounds1, List<PairIntArray> bounds2,
        int img1Width, int img1Height, int img2Width, int img2Height) {

        /*
        using difference of centroids in the region of intersection defined by
        the transformation.
        */

        double tolerance = 10;

        double bestDist = Double.MAX_VALUE;
        int bestDistIdx = -1;

        for (int i = 0; i < params.length; ++i) {

            TransformationParameters p = params[i];

            List<PairIntArray> trCurves1 = applyTransformForIntersection(p,
                bounds1, img2Width, img2Height);

            List<PairIntArray> trCurves2 = reverseTransformForIntersection(p,
                bounds2, img1Width, img1Height);

            double[][] centroidXY1s = calculateCentroids(trCurves1);
            double[][] centroidXY2s = calculateCentroids(trCurves2);

            // using a greedy approach, find the SSD of the distances of the
            // centroids

            int n1 = centroidXY1s.length;
            int n2 = centroidXY2s.length;

            double distSum = 0;
            int nMatched = 0;

            Set<Integer> chosen2 = new HashSet<Integer>();

            for (int idx1 = 0; idx1 < n1; ++idx1) {

                double[] xyCen1 = centroidXY1s[idx1];

                double minDiff = Double.MAX_VALUE;
                int min2Idx = -1;

                for (int idx2 = 0; idx2 < n2; ++idx2) {

                    if (chosen2.contains(Integer.valueOf(idx2))) {
                        continue;
                    }

                    double[] xyCen2 = centroidXY2s[idx2];

                    double diffX = xyCen1[0] - xyCen2[0];
                    double diffY = xyCen1[1] - xyCen2[1];

                    double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                    if (dist > tolerance) {
                        continue;
                    }

                    if (dist < minDiff) {
                        minDiff = dist;
                        min2Idx = idx2;
                    }
                }
                if (minDiff < Double.MAX_VALUE) {
                    distSum = minDiff;
                    nMatched++;
                    chosen2.add(Integer.valueOf(min2Idx));
                }
            }

            distSum /= (double)nMatched;

            if (distSum < bestDist) {
                bestDist = distSum;
                bestDistIdx = i;
            }
        }

        return params[bestDistIdx];
    }

    protected List<PairIntArray> applyTransformForIntersection(
        TransformationParameters params, List<PairIntArray> curves,
        int imgTrWidth, int imgTrHeight) {

        Transformer transformer = new Transformer();

        List<PairIntArray> trList = new ArrayList<PairIntArray>();

        for (int ii = 0; ii < curves.size(); ++ii) {

            PairIntArray tr = transformer.applyTransformation(params,
                curves.get(ii));

            boolean allWithinBounds = true;

            for (int j = 0; j < tr.getN(); ++j) {
                int x = tr.getX(j);
                int y = tr.getY(j);
                if ((x < 0) || (y < 0) || (x > (imgTrWidth - 1))
                    || (y > (imgTrHeight - 1))) {

                    allWithinBounds = false;

                    break;
                }
            }
            if (allWithinBounds) {
                trList.add(tr);
            }
        }

        return trList;
    }

    protected List<PairIntArray> reverseTransformForIntersection(
        TransformationParameters params, List<PairIntArray> curves,
        int imgWidth, int imgHeight) {

        MatchedPointsTransformationCalculator tc =
            new MatchedPointsTransformationCalculator();

        TransformationParameters revParams = tc.swapReferenceFrames(params);

        Transformer transformer = new Transformer();

        List<Integer> remove = new ArrayList<Integer>();

        for (int ii = 0; ii < curves.size(); ++ii) {

            PairIntArray tr = transformer.applyTransformation(revParams,
                curves.get(ii));

            for (int j = 0; j < tr.getN(); ++j) {

                int x = tr.getX(j);

                int y = tr.getY(j);

                if ((x < 0) || (y < 0) || (x > (imgWidth - 1))
                    || (y > (imgHeight - 1))) {

                    remove.add(Integer.valueOf(ii));

                    break;
                }
            }
        }

        if (remove.isEmpty()) {
            return curves;
        }

        List<PairIntArray> filtered = new ArrayList<PairIntArray>(curves);

        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i).intValue();
            filtered.remove(idx);
        }

        return filtered;
    }

    protected double[][] calculateCentroids(List<PairIntArray> curves) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[][] centroids = new double[curves.size()][];

        for (int i = 0; i < curves.size(); ++i) {

            centroids[i] = curveHelper.calculateXYCentroids(curves.get(i));
        }

        return centroids;
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

    protected boolean blobIsDarkerThanExterior(Set<PairInt> blob, PairIntArray
        closedEdge, GreyscaleImage img) {

        long sumExterior = 0;

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        Set<PairInt> added = new HashSet<PairInt>();

        for (int i = 0; i < closedEdge.getN(); ++i) {

            int x = closedEdge.getX(i);
            int y = closedEdge.getY(i);

            for (int ii = 0; ii < dxs8.length; ++ii) {

                int x2 = x + dxs8[ii];
                int y2 = y + dys8[ii];

                if ((x2 < 0) || (y2 < 0) || (x2 > (img.getWidth() - 1)) ||
                    (y2 > (img.getHeight() - 1))) {
                    continue;
                }

                PairInt p2 = new PairInt(x2, y2);

                if (!blob.contains(p2) && !added.contains(p2)) {
                   sumExterior += img.getValue(x2, y2);
                   added.add(p2);
                }
            }
        }

        float avgExterior = (float)sumExterior/(float)added.size();

        float avgInterior = mean(blob, img)/(float)blob.size();

        return (avgInterior < avgExterior);
    }

    protected float mean(Set<PairInt> blob, GreyscaleImage img) {

        long sum = 0;

        for (PairInt p : blob) {
            int x = p.getX();
            int y = p.getY();
            int v = img.getValue(x, y);
            sum += v;
        }

        float mean = (float)sum/(float)blob.size();

        return mean;
    }

    protected List<List<CurvatureScaleSpaceContour>> populateContours(
        List<PairIntArray> closedContours) {
        
        List<List<CurvatureScaleSpaceContour>> list 
            = new ArrayList<List<CurvatureScaleSpaceContour>>();
        
        for (int edgeIndex = 0; edgeIndex < closedContours.size(); ++edgeIndex) {
        
            List<CurvatureScaleSpaceContour> contours = 
                CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                closedContours.get(edgeIndex), edgeIndex);
            
            list.add(contours);
        }

        return list;
    }

}
