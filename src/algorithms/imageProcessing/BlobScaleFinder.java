package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * determine scale between 2 image using blob contours.
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
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.  Note that image processing
     * such as sky masks should be applied before using this method.
     * Also note that it is expected that it will be followed by a more rigorous
     * solver for Euclidean transforms such as the FeatureMatcher and then
     * the epipolar projection solver.
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to the image to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     * @return Euclidean scale to be applied to image1 to place it in the same
     * scale reference frame as image2.  Rotation and transformation are also
     * roughly solved for.
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public TransformationParameters calculateScale(ImageExt img1,
        ImageExt img2) throws IOException, NoSuchAlgorithmException {

        /*
        (1) applies histogram equalization to greyscale images of img1, img2
        (2) bins the greyscale and color images to sizes < 300 x 300.
        */
        GreyscaleImage img1Grey = img1.copyToGreyscale();
        GreyscaleImage img2Grey = img2.copyToGreyscale();

        boolean didApplyHistEq = applyHistogramEqualizationIfNeeded(img1Grey,
            img2Grey);

        final int k = 2;
        int smallestGroupLimit = 100;
        int largestGroupLimit = 5000;

        //TransformationParameters params = calculateScale(img1Grey, img2Grey, k,
        //    smallestGroupLimit, largestGroupLimit);

        float minDimension = 300.f;//200.f
        int binFactor = (int) Math.ceil(
            Math.max((float)img1Grey.getWidth()/minDimension,
            (float)img2Grey.getHeight()/minDimension));
        int smallestGroupLimitBinned = smallestGroupLimit/(binFactor*binFactor);
        int largestGroupLimitBinned = largestGroupLimit/(binFactor*binFactor);

        log.info("binFactor=" + binFactor);

        // prevent from being smaller than needed for a convex hull
        if (smallestGroupLimitBinned < 4) {
            smallestGroupLimitBinned = 4;
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img1GreyBinned = imageProcessor.binImage(img1Grey,
            binFactor);
        GreyscaleImage img2GreyBinned = imageProcessor.binImage(img2Grey,
            binFactor);

        log.info("binned:");

        TransformationParameters paramsBinned = calculateScale(img1GreyBinned,
            img2GreyBinned, k, smallestGroupLimitBinned, largestGroupLimitBinned);

        return paramsBinned;
    }

    protected boolean applyHistogramEqualizationIfNeeded(GreyscaleImage image1,
        GreyscaleImage image2) {

        // doing this automatically for now, but can use the stats insead to
        // decide
        boolean performHistEq = true;

        boolean useStatsToDecide = false;

        if (useStatsToDecide) {
            ImageStatistics stats1 = ImageStatisticsHelper.examineImage(image1,
                true);
            ImageStatistics stats2 = ImageStatisticsHelper.examineImage(image2,
                true);
            double median1DivMedian2 = stats1.getMedian() / stats2.getMedian();
            double meanDivMedian1 = stats1.getMean() / stats1.getMedian();
            double meanDivMedian2 = stats2.getMean() / stats2.getMedian();
            if (((median1DivMedian2 > 1) && ((median1DivMedian2 - 1) > 0.2))
                || ((median1DivMedian2 < 1) && (median1DivMedian2 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian1 > 1) && ((meanDivMedian1 - 1) > 0.2))
                || ((meanDivMedian1 < 1) && (meanDivMedian1 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian2 > 1) && ((meanDivMedian2 - 1) > 0.2))
                || ((meanDivMedian2 < 1) && (meanDivMedian2 < 0.8))) {
                performHistEq = true;
            }
        }

        if (performHistEq) {
            HistogramEqualization hEq = new HistogramEqualization(image1);
            hEq.applyFilter();
            hEq = new HistogramEqualization(image2);
            hEq.applyFilter();
        }

        return performHistEq;
    }

    /**
       <pre>
        (3) extract the top 10 blobs and their contours from k=2 segmented image:
            (3a) perform segmentation for k=2
            (3b) find blobs w/ DFSContiguousValueFinder for each intensity level
            (3c) use EdgeExtractorForBlobBorder to extract closed contour for
                 each blob
            (3d) return the top 10 longest contours and the blobs
       </pre>
     * @param k
     * @param img
     * @param outputBlobs
     * @param outputBounds
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    protected void extractBlobsFromSegmentedImage(int k, GreyscaleImage img,
        List<Set<PairInt>> outputBlobs, List<PairIntArray> outputBounds,
        int smallestGroupLimit, int largestGroupLimit) throws IOException,
        NoSuchAlgorithmException {

        extractBlobsFromSegmentedImage(k, img, outputBlobs, smallestGroupLimit,
            largestGroupLimit);

        if (outputBlobs.isEmpty()) {
            return;
        }

        boolean discardWhenCavityIsSmallerThanBorder = true;

        extractBoundsOfBlobs(outputBlobs, outputBounds, img.getWidth(),
            img.getHeight(), discardWhenCavityIsSmallerThanBorder);

    }

    protected void extractBlobsFromSegmentedImage(int k, GreyscaleImage img,
        List<Set<PairInt>> outputBlobs, int smallestGroupLimit,
        int largestGroupLimit) throws IOException, NoSuchAlgorithmException {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        GreyscaleImage imgS = img.copyImage();

        ImageProcessor imageProcessor = new ImageProcessor();

        imageProcessor.applyImageSegmentation(imgS, k);

        Map<Integer, Integer> freqMap = Histogram.createAFrequencyMap(imgS);

        for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()) {

            Integer pixValue = entry.getKey();

            DFSContiguousValueFinder finder = new DFSContiguousValueFinder(imgS);
            finder.setMinimumNumberInCluster(smallestGroupLimit);
            finder.findGroups(pixValue.intValue());

            int nGroups = finder.getNumberOfGroups();

            for (int i = 0; i < nGroups; ++i) {

                PairIntArray xy = finder.getXY(i);

                if (xy.getN() < largestGroupLimit) {

                    Set<PairInt> points = Misc.convert(xy);

                    // skip blobs that are on the image boundaries because they
                    // are incomplete
                    if (!curveHelper.hasNumberOfPixelsOnImageBoundaries(3,
                        points, img.getWidth(), img.getHeight())) {

                        outputBlobs.add(points);
                    }
                }
            }
        }
    }

    /**
     * given a list of inOutBlobs, also find blobs in the image segmented into
     * k color bins and discard all but those having the same centroids as
     * inOutBlobs, then find the perimeters of those and place them in
     * outBounds.
     * @param k
     * @param img
     * @param inOutBlobs
     * @param outBounds
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    protected void extractGivenBlobsFromSegmentedImage(int k, GreyscaleImage img,
        List<Set<PairInt>> inOutBlobs, List<PairIntArray> outBounds,
        int smallestGroupLimit, int largestGroupLimit) throws IOException,
        NoSuchAlgorithmException {

        int w = img.getWidth();
        int h = img.getHeight();

        /*
        (4c) if there were results from k=2, discard any blob whose centroid
             indicates it isn't in that list.
        */
        List<Set<PairInt>> blobs = new ArrayList<Set<PairInt>>();

        extractBlobsFromSegmentedImage(k, img, blobs, smallestGroupLimit,
            largestGroupLimit);

        boolean discardWhenCavityIsSmallerThanBorder = true;

        extractBoundsOfBlobs(blobs, outBounds, w, h,
            discardWhenCavityIsSmallerThanBorder);

        inOutBlobs.clear();
        inOutBlobs.addAll(blobs);
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
    protected void extractBoundsOfBlobs(final List<Set<PairInt>> inOutBlobs,
        final List<PairIntArray> outputBounds, int width, int height,
        boolean discardWhenCavityIsSmallerThanBorder) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<Integer> remove = new ArrayList<Integer>();

        for (int i = 0; i < inOutBlobs.size(); ++i) {

            Set<PairInt> blob = inOutBlobs.get(i);

            EdgeExtractorForBlobBorder extractor = new EdgeExtractorForBlobBorder();

            if (debug) {
                extractor.setToDebug();
            }

            PairIntArray closedEdge = extractor.extractAndOrderTheBorder0(
                blob, width, height,
                discardWhenCavityIsSmallerThanBorder);

            if ((closedEdge != null) &&
                (curveHelper.isAdjacent(closedEdge, 0, closedEdge.getN() - 1))) {

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
        GreyscaleImage img1,
        GreyscaleImage img2,
        List<Set<PairInt>> blobs1, List<Set<PairInt>> blobs2,
        List<PairIntArray> bounds1, List<PairIntArray> bounds2,
        int xRelativeOffset1, int yRelativeOffset1,
        int xRelativeOffset2, int yRelativeOffset2) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        double bestOverallStatSqSum = Double.MAX_VALUE;
        int bestOverallIdx1 = -1;
        int bestOverallIdx2 = -1;
        TransformationParameters bestOverallTransformation = null;
        int bestOverallNMatched = -1;
        List<FeatureComparisonStat> bestOverallCompStats = null;

        Map<Integer, Map<Integer, TransformationParameters>> paramsMap =
            new HashMap<Integer, Map<Integer, TransformationParameters>>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = bounds1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            double bestStatSqSum = Double.MAX_VALUE;
            int bestIdx2 = -1;
            TransformationParameters bestTransformation = null;
            int bestNMatched = -1;
            List<FeatureComparisonStat> bestCompStats = null;

            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = bounds2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);

//log.info("index1=" + index1.toString() + " index2=" + index2.toString());

                CurvatureScaleSpaceInflectionSingleEdgeMapper mapper =
                    new CurvatureScaleSpaceInflectionSingleEdgeMapper(
                    curve1, index1.intValue(), curve2, index2.intValue(),
                    xRelativeOffset1, yRelativeOffset1,
                    xRelativeOffset2, yRelativeOffset2);

                TransformationParameters params = mapper.matchContours();

                if ((params == null) ||
                    (mapper.getMatcher().getSolutionMatchedContours1().size() < 3)) {

                    continue;
                }

                List<FeatureComparisonStat> compStats =
                    filterContourPointsByFeatures(img1, img2, index1, index2,
                    blob1, blob2, curve1, curve2, features1, features2,
                    mapper.getMatcher());

                if (compStats.isEmpty()) {
                    continue;
                }

                double combinedStat = calculateCombinedIntensityStat(compStats);
                
                if (combinedStat < bestStatSqSum) {
                    bestStatSqSum = combinedStat;
                    bestIdx2 = index2.intValue();
                    bestTransformation = params;
                    bestNMatched = mapper.getMatcher().getSolutionMatchedContours1().size();
                    bestCompStats = compStats;
                    
                    log.info("  new best for [" + index1.toString() + "] ["
                        + index2.toString() + "]");

                    int z = 1;
                }
            }

            if (bestTransformation == null) {
                continue;
            }

            double[] xyCen2 = curveHelper.calculateXYCentroids(blobs2.get(bestIdx2));

            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
                "==>[%d](%d,%d) [%d](%d,%d) scale=%.2f  nMatched=%d  intSqDiff=%.1f",
                index1.intValue(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
                bestIdx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
                bestTransformation.getScale(), bestNMatched, (float)bestStatSqSum));
            log.info(sb.toString());

            Map<Integer, TransformationParameters> map2 = paramsMap.get(index1);
            if (map2 == null) {
                map2 = new HashMap<Integer, TransformationParameters>();
                paramsMap.put(index1, map2);
            }
            map2.put(bestIdx2, bestTransformation);

            if (bestStatSqSum < bestOverallStatSqSum) {
                bestOverallStatSqSum = bestStatSqSum;
                bestOverallIdx1 = index1.intValue();
                bestOverallIdx2 = bestIdx2;
                bestOverallTransformation = bestTransformation;
                bestOverallNMatched = bestNMatched;
                bestOverallCompStats = bestCompStats;
                
                log.info("  best overall for [" + bestOverallIdx1 + "] [" +
                    bestOverallIdx2 + "]");
                
                int z = 1;
            }
        }

        return bestOverallTransformation;
    }

    protected TransformationParameters calculateScale(GreyscaleImage img1,
        GreyscaleImage img2, int k, int smallestGroupLimit,
        int largestGroupLimit) throws IOException, NoSuchAlgorithmException {

        /*
        extract the top 10 blobs and their contours from k=2 segmented image:
        -- perform segmentation for k=2
        -- find blobs w/ DFSContiguousValueFinder for each intensity level
        -- use EdgeExtractorForBlobBorder to extract closed contour for each
           blob
        -- return the top 10 longest contours and the blobs
        */
        List<Set<PairInt>> blobs1 = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        List<PairIntArray> bounds1 = new ArrayList<PairIntArray>();
        List<PairIntArray> bounds2 = new ArrayList<PairIntArray>();
        log.info("image1:");
        extractBlobsFromSegmentedImage(2, img1, blobs1, bounds1,
            smallestGroupLimit, largestGroupLimit);
        log.info("image2:");
        extractBlobsFromSegmentedImage(2, img2, blobs2, bounds2,
            smallestGroupLimit, largestGroupLimit);

        /*
        filter out dissimilar pairings:
        -- given blobs from image1 and blobs from image 2, use feature
           matching to rule out possible pairings, resulting in possible
           matches for each.
        */
        //ImageProcessor imageProcessor = new ImageProcessor();

if (debug) {
Image img0 = ImageIOHelper.convertImage(img1);
for (int i = 0; i < bounds1.size(); ++i) {
    PairIntArray pa = bounds1.get(i);
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
MiscDebug.writeImageCopy(img0, "blob_contours_1_" + MiscDebug.getCurrentTimeFormatted() + ".png");
img0 = ImageIOHelper.convertImage(img2);
for (int i = 0; i < bounds2.size(); ++i) {
    PairIntArray pa = bounds2.get(i);
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
MiscDebug.writeImageCopy(img0, "blob_contours_2_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}
        /*
        solve for scale:
        -- use ContourMather to get scale solutions for each pairing then
           statistical basis of combining the results and removing
           outliers.
        */
        TransformationParameters params = solveForScale(
            img1, img2,
            blobs1, blobs2, bounds1, bounds2,
            img1.getXRelativeOffset(), img1.getYRelativeOffset(),
            img2.getXRelativeOffset(), img2.getYRelativeOffset());

        return params;
    }

    /**
     * extract the local points surrounding (x, y) on the
     * perimeter and return an object when creating descriptors.
     * Note that the perimeter is expected to be a closed curve.
     * @param theEdgeIndex
     * @param x
     * @param y
     * @param perimeterIdx
     * @param perimeter
     * @param blob
     * @return
     */
    private BlobPerimeterRegion extractBlobPerimeterRegion(
        int theEdgeIndex, CurvatureScaleSpaceImagePoint peakDetail,
        PairIntArray perimeter, Set<PairInt> blob) {

        if (perimeter == null || (perimeter.getN() < 3)) {
            throw new IllegalArgumentException(
            "perimeter cannot be null and must have at least 3 points");
        }

        if (blob == null) {
            throw new IllegalArgumentException("blob cannot be null");
        }

        // because of averaging for some peaks, sometimes detailIdx and
        // (x, y) are off by 1 so using the detailIdx primarily

        int x = peakDetail.getXCoord();
        int y = peakDetail.getYCoord();
        int detailIdx = peakDetail.getCoordIdx();

        int xPrev, yPrev, xNext, yNext;
        int xm = perimeter.getX(detailIdx);
        int ym = perimeter.getY(detailIdx);

        if (detailIdx > 0) {

            xPrev = perimeter.getX(detailIdx - 1);
            yPrev = perimeter.getY(detailIdx - 1);

            if (detailIdx < (perimeter.getN() - 2)) {
                xNext = perimeter.getX(detailIdx + 1);
                yNext = perimeter.getY(detailIdx + 1);
            } else {
                // it's a closed curve, but may need to assert it here.
                xNext = perimeter.getX(perimeter.getN() - 1);
                yNext = perimeter.getY(perimeter.getN() - 1);
            }

        } else {

            assert(detailIdx == 0);

            // it's a closed curve, but may need to assert it here.
            xPrev = perimeter.getX(perimeter.getN() - 1);
            yPrev = perimeter.getY(perimeter.getN() - 1);

            xNext = perimeter.getX(detailIdx + 1);
            yNext = perimeter.getY(detailIdx + 1);
        }

        BlobPerimeterRegion region = new BlobPerimeterRegion(theEdgeIndex,
            xPrev, yPrev, xm, ym, xNext, yNext, blob);

        return region;
    }

    private List<FeatureComparisonStat> filterContourPointsByFeatures(
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

        if (redoStats) {

            compStats = redoFilterContourPointsByFeatures(img1, img2, index1,
                index2, blob1, blob2, curve1, curve2,
                bestCompStat.getImg1PointRotInDegrees(),
                bestCompStat.getImg2PointRotInDegrees(), matcher);
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
                    featureMatcher.ditherForBestLocation(
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
        
if ((index1.intValue() == 0) && (index2.intValue() == 0)) {
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
ImageIOHelper.writeOutputImage(bin + "/contours1_redone.png", img1C);
ImageIOHelper.writeOutputImage(bin + "/contours2_redone.png", img2C);
int z = 1;
} catch(IOException e) {
}
}        
        
        
        return compStats;
    }

    private String printToString(List<FeatureComparisonStat> compStats) {
        
        StringBuilder sb = new StringBuilder();
        
        for (FeatureComparisonStat compStat : compStats) {
            
            sb.append(String.format(
                " (%d,%d) (%d,%d) intSqDiff=%.1f(%.1f)",
                compStat.getImg1Point().getX(), compStat.getImg1Point().getY(),
                compStat.getImg2Point().getX(), compStat.getImg2Point().getY(),
                compStat.getSumIntensitySqDiff(), 
                compStat.getImg2PointIntensityErr()));
        }
        
        return sb.toString();
    }

    private double calculateCombinedIntensityStat(List<FeatureComparisonStat> compStats) {
        
        /*
        these are square roots of sums of squared differences.        
        */
        
        double sum = 0;
        
        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }
        
        sum /= (double)compStats.size();
        
        return sum;
    }
    
}
