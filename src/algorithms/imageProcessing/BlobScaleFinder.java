package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * determine scale between 2 image using blob contours.
 *
 * @author nichole
 */
public class BlobScaleFinder {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.  Note that image processing
     * such as sky masks should be applied before using this method.
     * Also not that it is expected that it will be followed by a more rigorous
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

        int smallestGroupLimit = 100;
        int largestGroupLimit = 5000;
        float minDimension = 300.f;//200.f
        int binFactor = (int) Math.ceil(
            Math.max((float)img1Grey.getWidth()/minDimension,
                (float)img2Grey.getHeight()/minDimension));
        smallestGroupLimit /= (binFactor*binFactor);
        largestGroupLimit /= (binFactor*binFactor);

        log.info("binFactor=" + binFactor);

        // prevent from being smaller than needed for a convex hull
        if (smallestGroupLimit < 4) {
            smallestGroupLimit = 4;
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img1GreyBinned = imageProcessor.binImage(img1Grey,
            binFactor);
        GreyscaleImage img2GreyBinned = imageProcessor.binImage(img2Grey,
            binFactor);

        /*
        (3) extract the top 10 blobs and their contours from k=2 segmented image:
            (3a) perform segmentation for k=2
            (3b) find blobs w/ DFSContiguousValueFinder for each intensity level
            (3c) use EdgeExtractorForBlobBorder to extract closed contour for
                 each blob
            (3d) return the top 10 longest contours and the blobs
        */
        List<Set<PairInt>> blobs1 = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        List<PairIntArray> bounds1 = new ArrayList<PairIntArray>();
        List<PairIntArray> bounds2 = new ArrayList<PairIntArray>();
        extractBlobsFromSegmentedImage(2, img1GreyBinned, blobs1, bounds1,
            smallestGroupLimit, largestGroupLimit);
        extractBlobsFromSegmentedImage(2, img2GreyBinned, blobs2, bounds2,
            smallestGroupLimit, largestGroupLimit);

        /*
        (4) extract the top 10 blobs and their contours from k=3 segmented image:
            (4a) perform segmentation for k=3
            (4b) find blobs w/ DFSContiguousValueFinder for each intensity level
            (4c) if there were results from k=2, discard any blob whose centroid
                 indicates it isn't in that list.
            (4d) use EdgeExtractorForBlobBorder to extract closed contour for
                 each blob
            (4e) if (there were no blobs returned from k=2), choose the top 10
                 longest contours
        */
        bounds1.clear();
        bounds2.clear();
        if (blobs1.size() < 10) {
            extractBlobsFromSegmentedImage(3, img1GreyBinned, blobs1, bounds1,
                smallestGroupLimit, largestGroupLimit);
        } else {
            extractGivenBlobsFromSegmentedImage(3, img1GreyBinned, blobs1,
                bounds1, smallestGroupLimit, largestGroupLimit);
        }
        if (blobs2.size() < 10) {
            extractBlobsFromSegmentedImage(3, img2GreyBinned, blobs2, bounds2,
                smallestGroupLimit, largestGroupLimit);
        } else {
            extractGivenBlobsFromSegmentedImage(3, img2GreyBinned, blobs2,
                bounds2, smallestGroupLimit, largestGroupLimit);
        }

        /*
        (5) filter out dissimilar pairings:
            (5a) given blobs from image1 and blobs from image 2, use feature
                 matching to rule out possible pairings, resulting in possible
                 matches for each.
        */

        boolean doNormalize = false;

        Map<Integer, List<Integer>> filteredBlobMatches =
            filterBlobsByFeatures(
            imageProcessor.binImage(img1Grey, binFactor),
            imageProcessor.binImage(img2Grey, binFactor), blobs1, blobs2,
            doNormalize);

        if (filteredBlobMatches.isEmpty()) {
            return null;
        }
        
        /*
        (6) solve for scale:
            (6a) use ContourMather to get scale solutions for each pairing
                 then statistical basis of combining the results and removing
                 outliers.
        */
        TransformationParameters params = solveForScale(filteredBlobMatches,
            blobs1, blobs2, bounds1, bounds2,
            img1GreyBinned.getXRelativeOffset(),
            img1GreyBinned.getYRelativeOffset(),
            img2GreyBinned.getXRelativeOffset(),
            img2GreyBinned.getYRelativeOffset());

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

        extractBoundsOfBlobs(outputBlobs, outputBounds, img.getWidth(), img.getHeight(),
            discardWhenCavityIsSmallerThanBorder);

    }

    protected void extractBlobsFromSegmentedImage(int k, GreyscaleImage img,
        List<Set<PairInt>> outputBlobs, int smallestGroupLimit,
        int largestGroupLimit) throws IOException, NoSuchAlgorithmException {

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
                    outputBlobs.add(points);
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

        int tolerance = 10;

        blobs = filterBlobsByFirstList(inOutBlobs, blobs, tolerance);

        boolean discardWhenCavityIsSmallerThanBorder = true;

        extractBoundsOfBlobs(blobs, outBounds, w, h,
            discardWhenCavityIsSmallerThanBorder);
        
        inOutBlobs.clear();
        inOutBlobs.addAll(blobs);
    }

    /**
     * given the original list of blobs origBlobs, find blobs in nextBlobs
     * having same centroids within tolerance and return those.
     *
     * @param origBlobs
     * @param nextBlobs
     * @param tolerance
     * @return
     */
    protected List<Set<PairInt>> filterBlobsByFirstList(
        List<Set<PairInt>> origBlobs, List<Set<PairInt>> nextBlobs,
        int tolerance) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<Integer> keep = new ArrayList<Integer>();

        double[][] xyCentroids0 = new double[origBlobs.size()][];
        for (int j = 0; j < origBlobs.size(); ++j) {
            xyCentroids0[j] = curveHelper.calculateXYCentroids(origBlobs.get(j));
        }

        int toleranceSq = tolerance * tolerance;

        for (int i = 0; i < nextBlobs.size(); ++i) {

            Set<PairInt> blob = nextBlobs.get(i);

            double[] xyCen = curveHelper.calculateXYCentroids(blob);

            for (int j = 0; j < origBlobs.size(); ++j) {
                double[] xyCen0 = xyCentroids0[j];
                double diffX = Math.abs(xyCen[0] - xyCen0[0]);
                double diffY = Math.abs(xyCen[1] - xyCen0[1]);
                double distSq = diffX * diffX + diffY * diffY;
                
                if (distSq < toleranceSq) {
                    keep.add(Integer.valueOf(i));
                    break;
                }
            }
        }

        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();
        for (Integer index : keep) {
            output.add(nextBlobs.get(index.intValue()));
        }

        return output;
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

        List<Integer> remove = new ArrayList<Integer>();

        for (int i = 0; i < inOutBlobs.size(); ++i) {

            Set<PairInt> blob = inOutBlobs.get(i);

            EdgeExtractorForBlobBorder extractor =
                new EdgeExtractorForBlobBorder();

            PairIntArray closedEdge = extractor.extractAndOrderTheBorder0(
                blob, width, height,
                discardWhenCavityIsSmallerThanBorder);

            if (closedEdge != null) {
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

        int last = (inOutBlobs.size() > 10) ? 10 : inOutBlobs.size();
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

    protected Map<Integer, List<Integer>> filterBlobsByFeatures(
        GreyscaleImage img1Binned, GreyscaleImage img2Binned,
        List<Set<PairInt>> blobs1, List<Set<PairInt>> blobs2, boolean doNormalize) {

MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        Map<Integer, List<Integer>> possiblePairs = new HashMap<Integer, List<Integer>>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            Set<PairInt> blob1 = blobs1.get(idx1);

            List<Integer> similar2 = new ArrayList<Integer>();

            float sumSquaredError1 = MiscMath.sumSquaredError(img1Binned, blob1);

            int sum1 = sumIntensity(img1Binned, blob1, doNormalize);

            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                Set<PairInt> blob2 = blobs2.get(idx2);

                int sum2 = sumIntensity(img2Binned, blob2, doNormalize);

                float diff = sum1 - sum2;

                float ssd = (diff * diff)/2.f;
  
double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);
double[] xyCen2 = curveHelper.calculateXYCentroids(blob2);
log.info(String.format("blob1=(%d,%d) blob2=(%d,%d) ssd=%.2f  sumSqErr=%.2f", 
(int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
(int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]), (float)ssd, (float)sumSquaredError1));

                if (ssd < sumSquaredError1) {
                    similar2.add(Integer.valueOf(idx2));
                }
            }

            if (!similar2.isEmpty()) {
                possiblePairs.put(Integer.valueOf(idx1), similar2);
            }
        }

        return possiblePairs;
    }

    /**
     * sum the intensity of the points with an option to subtract the mean.
     * @param img
     * @param points
     * @param doNormalize
     * @return
     */
    protected int sumIntensity(GreyscaleImage img, Set<PairInt> points,
        boolean doNormalize) {

        if (img == null) {
            throw new IllegalStateException("img cannot be null");
        }

        if (points == null) {
            throw new IllegalStateException("points cannot be null");
        }

        int sum = 0;

        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            int v = img.getValue(x, y);
            sum += v;
        }

        if (doNormalize) {
            double[] avgStDev = MiscMath.getAvgAndStDev(img, points);
            sum -= (points.size() * avgStDev[0]);
        }

        return sum;
    }

    protected TransformationParameters solveForScale(
        Map<Integer, List<Integer>> filteredBlobMatches,
        List<Set<PairInt>> blobs1, List<Set<PairInt>> blobs2,
        List<PairIntArray> bounds1, List<PairIntArray> bounds2,
        int xRelativeOffset1, int yRelativeOffset1,
        int xRelativeOffset2, int yRelativeOffset2) {

        Map<Integer, Map<Integer, TransformationParameters>> paramsMap =
            new HashMap<Integer, Map<Integer, TransformationParameters>>();
        
        Map<Integer, Map<Integer, Double>> costMap =
            new HashMap<Integer, Map<Integer, Double>>();

        for (Entry<Integer, List<Integer>> entry : filteredBlobMatches.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            PairIntArray curve1 = bounds1.get(index1.intValue());
                
            List<Integer> blob2Indexes = entry.getValue();
            
            double bestCost = Double.MAX_VALUE;
            int bestIdx2 = -1;
            
            for (Integer index2 : blob2Indexes) {
                
                PairIntArray curve2 = bounds2.get(index2.intValue());
                
                CurvatureScaleSpaceInflectionSingleEdgeMapper mapper = 
                    new CurvatureScaleSpaceInflectionSingleEdgeMapper(
                    curve1, index1.intValue(), curve2, index2.intValue(),
                    xRelativeOffset1, yRelativeOffset1,
                    xRelativeOffset2, yRelativeOffset2);
                
                TransformationParameters params = mapper.matchContours();
                
                if (params != null) {
                    
                    double cost = mapper.getMatcher().getSolvedCost();
                    if (cost < bestCost) {
                        bestCost = cost;
                        bestIdx2 = index2.intValue();
                    }
                
                    Map<Integer, TransformationParameters> map2 = paramsMap.get(index1);                    
                    if (map2 == null) {
                        map2 = new HashMap<Integer, TransformationParameters>();
                        paramsMap.put(index1, map2);
                    }
                    map2.put(index2, params);
                    
                    Map<Integer, Double> mapCosts2 = costMap.get(index1);
                    if (mapCosts2 == null) {
                        mapCosts2 = new HashMap<Integer, Double>();
                        costMap.put(index1, mapCosts2);
                    }
                    costMap.put(index2, mapCosts2);
                    
                }
            }
            int z = 1;
        }

        return null;
    }

}
