package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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

        TransformationParameters params = calculateScale(img1Grey, img2Grey, k,
            smallestGroupLimit, largestGroupLimit);

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

        return params;
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

        FeatureMatcher featureMatcher = new FeatureMatcher();

        int dither = 1;

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

                double[] xyCen2 = curveHelper.calculateXYCentroids(blob2);

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

                double statSqSum = 0;

                int nMaxMatchable = mapper.getMatcher().getNMaxMatchable();
                int nMaxStats = 0;
                int nStats = 0;

                List<FeatureComparisonStat> compStats = new ArrayList<FeatureComparisonStat>();

                double cost = mapper.getMatcher().getSolvedCost();

                boolean doesMatch = true;

                StringBuilder sb = new StringBuilder();
                sb.append(String.format(
                    "[%d](%d,%d) [%d](%d,%d) cost=%.1f scale=%.2f  nMatched=%d ",
                    index1.intValue(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
                    index2.intValue(), (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
                    (float)cost, (float)mapper.getMatcher().getSolvedScale(),
                    mapper.getMatcher().getSolutionMatchedContours1().size()));

if ((index1.intValue() == 0) && (index2.intValue() == 0)) {
List<CurvatureScaleSpaceContour> list1 = new ArrayList<CurvatureScaleSpaceContour>();
for (int j = 0; j < mapper.getMatcher().getSolutionMatchedContours1().size(); ++j) {
    list1.add(mapper.getMatcher().getSolutionMatchedContours1().get(j));
}
ImageExt img1C = img1.createColorGreyscaleExt();
MiscDebug.debugPlot(list1, img1C,
img1.getXRelativeOffset(), img1.getYRelativeOffset(), "_1_ind0_ind0_");

List<CurvatureScaleSpaceContour> list2 = new ArrayList<CurvatureScaleSpaceContour>();
for (int j = 0; j < mapper.getMatcher().getSolutionMatchedContours2().size(); ++j) {
    list2.add(mapper.getMatcher().getSolutionMatchedContours2().get(j));
}
ImageExt img2C = img2.createColorGreyscaleExt();
MiscDebug.debugPlot(list2, img2C,
img2.getXRelativeOffset(), img2.getYRelativeOffset(), "_2_ind0_ind0_");
int z = 1;
}

                for (int j = 0; j < mapper.getMatcher().getSolutionMatchedContours1().size(); ++j) {

                    CurvatureScaleSpaceContour c1 =
                        mapper.getMatcher().getSolutionMatchedContours1().get(j);

                    CurvatureScaleSpaceContour c2 =
                        mapper.getMatcher().getSolutionMatchedContours2().get(j);

                    // the sizes of the peak details will be the same
                    CurvatureScaleSpaceImagePoint[] details1 = c1.getPeakDetails();
                    CurvatureScaleSpaceImagePoint[] details2 = c2.getPeakDetails();

                    nMaxStats += details1.length;

                    for (int jj = 0; jj < details1.length; ++jj) {

                        int x1 = details1[jj].getXCoord();
                        int y1 = details1[jj].getYCoord();
                        int detailIdx1 = details1[jj].getCoordIdx();

                        BlobPerimeterRegion region1 = extractBlobPerimeterRegion(
                            index1.intValue(), x1, y1, detailIdx1, curve1, blob1
                        );

                        int x2 = details2[jj].getXCoord();
                        int y2 = details2[jj].getYCoord();
                        int detailIdx2 = details2[jj].getCoordIdx();

                        BlobPerimeterRegion region2 = extractBlobPerimeterRegion(
                            index2.intValue(), x2, y2, detailIdx2, curve2, blob2
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

                            if (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr()) {

                                doesMatch = false;

                                break;

                            } else {

                                float sumIntSqDiff = compStat.getSumIntensitySqDiff();

                                statSqSum += (sumIntSqDiff*sumIntSqDiff);

                                sb.append(String.format("  %.1f(%.1f), ",
                                    sumIntSqDiff,
                                    compStat.getImg2PointIntensityErr()));

                                compStats.add(compStat);

                                nStats++;
                            }
                        }
                    } // end details

                    if (!doesMatch) {
                        break;
                    }

                }// end matching contours for index1, index2

                if (!doesMatch) {
                    continue;
                }

                log.info(sb.toString());
                sb = new StringBuilder();

                statSqSum = (nStats == 0) ? Double.MAX_VALUE : Math.sqrt(statSqSum);

                if (statSqSum < bestStatSqSum) {
                    bestStatSqSum = statSqSum;
                    bestIdx2 = index2.intValue();
                    bestTransformation = params;
                    bestNMatched = mapper.getMatcher().getSolutionMatchedContours1().size();
                    bestCompStats = compStats;

                    sb.append(String.format(
                        "   intSqDiff=%.1f nStats=(%d out of %d) ",
                        (float)bestStatSqSum, nStats, nMaxStats));
                    log.info(sb.toString());

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
        int theEdgeIndex, int x, int y, int perimeterIdx,
        PairIntArray perimeter, Set<PairInt> blob) {

        if (perimeter == null || (perimeter.getN() < 3)) {
            throw new IllegalArgumentException(
            "perimeter cannot be null and must have at least 3 points");
        }

        if (blob == null) {
            throw new IllegalArgumentException("blob cannot be null");
        }

        // because of averaging for some peaks, sometimes perimeterIdx and
        // (x, y) are off by 1 so using the perimeterIdx primarily

        int xPrev, yPrev, xNext, yNext;
        int xm = perimeter.getX(perimeterIdx);
        int ym = perimeter.getY(perimeterIdx);

        if (perimeterIdx > 0) {

            xPrev = perimeter.getX(perimeterIdx - 1);
            yPrev = perimeter.getY(perimeterIdx - 1);

            if (perimeterIdx < (perimeter.getN() - 2)) {
                xNext = perimeter.getX(perimeterIdx + 1);
                yNext = perimeter.getY(perimeterIdx + 1);
            } else {
                // it's a closed curve, but may need to assert it here.
                xNext = perimeter.getX(perimeter.getN() - 1);
                yNext = perimeter.getY(perimeter.getN() - 1);
            }

        } else {

            assert(perimeterIdx == 0);

            // it's a closed curve, but may need to assert it here.
            xPrev = perimeter.getX(perimeter.getN() - 1);
            yPrev = perimeter.getY(perimeter.getN() - 1);

            xNext = perimeter.getX(perimeterIdx + 1);
            yNext = perimeter.getY(perimeterIdx + 1);
        }

        BlobPerimeterRegion region = new BlobPerimeterRegion(theEdgeIndex,
            xPrev, yPrev, xm, ym, xNext, yNext, blob);

        return region;
    }

}
