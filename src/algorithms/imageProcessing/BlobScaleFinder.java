package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
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

    protected boolean debug = true;
    
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

        final int k = 2;
        int smallestGroupLimit = 100;
        int largestGroupLimit = 5000;
        
        TransformationParameters params = calculateScale(img1Grey, img2Grey, k,
            smallestGroupLimit, largestGroupLimit);
        
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

        TransformationParameters paramsBinned = calculateScale(img1GreyBinned,
            img2GreyBinned, k, smallestGroupLimit, largestGroupLimit);
        
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

        Map<Integer, Map<Integer, TransformationParameters>> paramsMap =
            new HashMap<Integer, Map<Integer, TransformationParameters>>();
        
        Map<Integer, Map<Integer, Double>> costMap =
            new HashMap<Integer, Map<Integer, Double>>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {
        
            Integer index1 = Integer.valueOf(idx1);
            
            PairIntArray curve1 = bounds1.get(idx1);
            
            Set<PairInt> blob1 = blobs1.get(idx1);
            
            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {
                
                Integer index2 = Integer.valueOf(idx2);
            
                PairIntArray curve2 = bounds1.get(idx2);
            
                Set<PairInt> blob2 = blobs2.get(idx2);
            
                double bestCost = Double.MAX_VALUE;
                int bestIdx2 = -1;
                                           
log.info("index1=" + index1.toString() + " index2=" + index2.toString());

                CurvatureScaleSpaceInflectionSingleEdgeMapper mapper = 
                    new CurvatureScaleSpaceInflectionSingleEdgeMapper(
                    curve1, index1.intValue(), curve2, index2.intValue(),
                    xRelativeOffset1, yRelativeOffset1,
                    xRelativeOffset2, yRelativeOffset2);
                
                TransformationParameters params = mapper.matchContours();
                
                if ((params != null) && 
                    (mapper.getMatcher().getSolutionMatchedContours1().size() > 2)) {
                    
                    int nMaxMatchable = mapper.getMatcher().getNMaxMatchable();
                    
                    double cost = mapper.getMatcher().getSolvedCost();
                    if (cost < bestCost) {
                        bestCost = cost;
                        bestIdx2 = index2.intValue();
                    }
                
                    if ((nMaxMatchable 
                        - mapper.getMatcher().getSolutionMatchedContours1().size()) == 0) {
                        
                        String str = String.format(
                        "[%d] [%d] cost=%.1f scale=%.2f  nMatched=%d", 
                        index1.intValue(), index2.intValue(), (float)cost, 
                        (float)mapper.getMatcher().getSolvedScale(),
                        mapper.getMatcher().getSolutionMatchedContours1().size());
                        
                        log.info(str);
                    }
                    
                    /*
                    TODO: compare feature descriptors of the inflection points to help
                    remove false matches
                    */
                
if (
((index1.intValue() == 2) && (index2.intValue() == 6)) ||
((index1.intValue() == 7) && (index2.intValue() == 10))) {
    int z = 1;
    /*
    scl=0.86, cost= 13.0  nMatched=9  nMaxMatchable=9
    scl=1.07, cost=134.0          =3               =3
    */
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
        }

        // TODO: analyze map
        
        return null;
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
        extractBlobsFromSegmentedImage(2, img1, blobs1, bounds1,
            smallestGroupLimit, largestGroupLimit);
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
        TransformationParameters params = solveForScale(img1, img2,
            blobs1, blobs2, bounds1, bounds2,
            img1.getXRelativeOffset(), img1.getYRelativeOffset(),
            img2.getXRelativeOffset(), img2.getYRelativeOffset());
        
        return params;
    }
    
}
