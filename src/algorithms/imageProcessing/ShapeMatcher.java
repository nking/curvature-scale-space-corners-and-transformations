package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayDescendingComparator;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

public class ShapeMatcher {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ShapeMatcher() {
    }

    /**
    method to extract general shapes from the images and compare them in order to 
    match points.  It returns a fit to a rough Euclidean transformation.
    NOTE that the images may need pre-processing steps before using this.  For example,
    the Brown & Lowe 200? panoramic images of a mountain need to have the sky masked
    out of the image first.
     * @param image1
     * @param image2
     * @param outputMatched1
     * @param outputMatched2
    */
    public TransformationPointFit findMatchingShapes(ImageExt image1, ImageExt image2,
    PairIntArray outputMatched1, PairIntArray outputMatched2) throws 
        IOException, NoSuchAlgorithmException {
                
        GreyscaleImage img1Grey = image1.copyToGreyscale();
        GreyscaleImage img2Grey = image2.copyToGreyscale();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        final boolean performBinning = false;
        int binFactor1 = 1;
        
        int kN = 4;
        boolean performBinarySegmentation = true;
        if (performBinarySegmentation) {
            kN = 2;
        }
        
        /*
        one could start with essentially no limits here and then
        looks at the distribution of resulting contiguous group
        sizes to decide the range to keep.
        For now, choosing limits.
        */
        int smallestGroupLimit = 100;
        //TODO: consider scaling this by image size or by size and res if one
        //  day that information is passed to this method
        int largestGroupLimit = 5000;
        
        ImageExt img1Cp = (ImageExt)image1.copyImage();
        ImageExt img2Cp = (ImageExt)image2.copyImage();
        
        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(img1Grey, true); 
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(img2Grey, true); 

        log.info("stats1=" + stats1.toString());
        log.info("stats2=" + stats2.toString());
        
        boolean performHistEq = false;        
        double median1DivMedian2 = stats1.getMedian()/stats2.getMedian();
        double meanDivMedian1 = stats1.getMean()/stats1.getMedian();        
        double meanDivMedian2 = stats2.getMean()/stats2.getMedian();
        if (
            ((median1DivMedian2 > 1) && ((median1DivMedian2 - 1) > 0.2)) ||
            ((median1DivMedian2 < 1) && (median1DivMedian2 < 0.8))) {
            performHistEq = true;
        } else if (
            ((meanDivMedian1 > 1) && ((meanDivMedian1 - 1) > 0.2)) ||
            ((meanDivMedian1 < 1) && (meanDivMedian1 < 0.8))) {
            performHistEq = true;
        } else if (
            ((meanDivMedian2 > 1) && ((meanDivMedian2 - 1) > 0.2)) ||
            ((meanDivMedian2 < 1) && (meanDivMedian2 < 0.8))) {
            performHistEq = true;
        }
        if (performHistEq) {
            log.info("use histogram equalization on the greyscale images");
            HistogramEqualization hEq = new HistogramEqualization(img1Grey);
            hEq.applyFilter();
            hEq = new HistogramEqualization(img2Grey);
            hEq.applyFilter();
            /*HistogramEqualizationForColor hEqC = new HistogramEqualizationForColor(img1Cp);
            hEqC.applyFilter();
            hEqC = new HistogramEqualizationForColor(img2Cp);
            hEqC.applyFilter();*/
        }
        
        if (performBinning) {
            binFactor1 = (int) Math.ceil(
                Math.max((float)img1Grey.getWidth()/200.f,
                (float)img2Grey.getHeight()/200.));
            smallestGroupLimit /= (binFactor1*binFactor1);
            largestGroupLimit /= (binFactor1*binFactor1);
            
            log.info("binFactor1=" + binFactor1);
            
            // prevent from being smaller than needed for a convex hull
            if (smallestGroupLimit < 4) {
                smallestGroupLimit = 4;
            }
            img1Grey = imageProcessor.binImage(img1Grey, binFactor1);
            img2Grey = imageProcessor.binImage(img2Grey, binFactor1);
            img1Cp = imageProcessor.binImage(img1Cp, binFactor1);
            img2Cp = imageProcessor.binImage(img2Cp, binFactor1);
        }

        imageProcessor.applyImageSegmentation(img1Grey, kN);
        imageProcessor.applyImageSegmentation(img2Grey, kN);

        // == contiguous regions within size limits become blobs of interest,
        //    indexed by their intensity levels
        
        Map<Integer, Integer> freqMap1 = Histogram.createAFrequencyMap(img1Grey);
        Map<Integer, Integer> freqMap2 = Histogram.createAFrequencyMap(img2Grey);
        
        Map<Integer, List<PairIntArray>> contigMap1 
            = new HashMap<Integer, List<PairIntArray>>();
        Map<Integer, List<PairIntArray>> contigMap2 
            = new HashMap<Integer, List<PairIntArray>>();
        
        Map<Integer, List<GrahamScan>> hulls1 = 
            new HashMap<Integer, List<GrahamScan>>();
        Map<Integer, List<GrahamScan>> hulls2 = 
            new HashMap<Integer, List<GrahamScan>>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        Map<Integer, PairIntArray> hullCentroids1Map = 
            new HashMap<Integer, PairIntArray>();
        Map<Integer, PairIntArray> hullCentroids2Map = 
            new HashMap<Integer, PairIntArray>();
        PairIntArray allHullCentroids1 = new PairIntArray();
        PairIntArray allHullCentroids2 = new PairIntArray();
        
        for (int im = 0; im < 2; ++im) {
            
            Map<Integer, Integer> freqMap = freqMap1;
            Map<Integer, List<PairIntArray>> contigMap = contigMap1;
            Map<Integer, List<GrahamScan>> hulls = hulls1;
            GreyscaleImage imgGrey = img1Grey;
            PairIntArray hullCentroids = allHullCentroids1;
            Map<Integer, PairIntArray> hullCentroidsMap = hullCentroids1Map;
            if (im == 1) {
                freqMap = freqMap2;
                contigMap = contigMap2;
                hulls = hulls2;
                hullCentroids = allHullCentroids2;
                imgGrey = img2Grey;
                hullCentroidsMap = hullCentroids2Map;
            }
 
            for (Entry<Integer, Integer> entry : freqMap.entrySet()) {

                Integer pixValue = entry.getKey();

                DFSContiguousValueFinder finder = new DFSContiguousValueFinder(
                    imgGrey);
                finder.setMinimumNumberInCluster(smallestGroupLimit);
                finder.findGroups(pixValue.intValue());

                int nGroups = finder.getNumberOfGroups();
                List<PairIntArray> list = new ArrayList<PairIntArray>();
                for (int i = 0; i < nGroups; i++) {
                    PairIntArray xy = finder.getXY(i);
                    if (xy.getN() < largestGroupLimit) {
                        list.add(xy);
                    }
                }
                Collections.sort(list, new PairIntArrayDescendingComparator());
                
                // storing the centroids for this intensity level hulls separateley
                PairIntArray pvHullCentroids = new PairIntArray();
                
                // remove hulls with large area on image bounds
                List<Integer> rm = new ArrayList<Integer>();
                
                List<GrahamScan> listHulls = new ArrayList<GrahamScan>();
                for (int i = 0; i < list.size(); ++i) {
                    
                    PairIntArray xy = list.get(i);
                    
                    GrahamScan scan = new GrahamScan();
                    try {
                        float[] x = new float[xy.getN()];
                        float[] y = new float[x.length];
                        for (int ii = 0; ii < x.length; ++ii) {
                            x[ii] = xy.getX(ii);
                            y[ii] = xy.getY(ii);
                        }
                        
                        double[] centroidXY = 
                            curveHelper.calculateXYCentroids(x, y);

                        scan.computeHull(x, y);
                        
                        // if more than one hull point touches the bounds of the
                        // image, the hull is removed because it may be
                        // incomplete.  may need to change this rule later
                        // if it makes the solution too shallow in terms of 
                        // very close points
                        
                        int nBounds = 0;
                        for (int ii = 0; ii < scan.getXHull().length; ++ii) {
                            float xh = scan.getXHull()[ii];
                            float yh = scan.getYHull()[ii];
                            if ((xh == 0) || xh == (imgGrey.getWidth() - 1) ||
                                (yh == 0) || yh == (imgGrey.getHeight() - 1)) {
                                nBounds++;
                            }
                        }
                        if (nBounds > 3) {                            
                            rm.add(Integer.valueOf(i));
                        } else {
                            listHulls.add(scan);
                            int xh = (int)Math.round(centroidXY[0]);
                            int yh = (int)Math.round(centroidXY[1]);
                            hullCentroids.add(xh, yh);
                            pvHullCentroids.add(xh, yh);
                        }
                        
                    } catch (GrahamScanTooFewPointsException e) {
                        log.severe(e.getMessage());
                    }
                }
                
                for (int i = (rm.size() - 1); i > -1; --i) {
                    int rmIdx = rm.get(i).intValue();
                    list.remove(rmIdx);
                }

                log.info("nHulls" + (im + 1) + "=" + listHulls.size() + " for intensity=" + pixValue.toString());
                
                contigMap.put(pixValue, list);
                hulls.put(pixValue, listHulls);
                hullCentroidsMap.put(pixValue, pvHullCentroids);
            }
        }
        
        MiscDebug.writeHullImages(img1Grey, hulls1, "1_binned_hulls");
        MiscDebug.writeHullImages(img2Grey, hulls2, "2_binned_hulls");
        MiscDebug.writeImage(img1Cp, "1_binned_clr");
        MiscDebug.writeImage(img2Cp, "2_binned_clr");
       
        // make corners
        
        if (!performBinning) {
            imageProcessor.blur(img1Grey, 2);
            imageProcessor.blur(img2Grey, 2);
        }
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1Grey);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
        PairIntArray corners1 = detector.getCornersInOriginalReferenceFrame();
        
        CurvatureScaleSpaceCornerDetector detector2 = new
            CurvatureScaleSpaceCornerDetector(img2Grey);
        detector2.doNotPerformHistogramEqualization();
        detector2.findCorners();
        PairIntArray corners2 = detector2.getCornersInOriginalReferenceFrame();
        
        log.info("n1Corners=" + corners1.getN() + " n2Corners2=" 
            + corners2.getN());
        
        // experimenting with a slightly different definition for theta:
        GreyscaleImage theta1360 = imageProcessor.computeTheta360(
            detector.getGradientX(), detector.getGradientY());
        GreyscaleImage theta2360 = imageProcessor.computeTheta360(
            detector2.getGradientX(), detector2.getGradientY());
        MiscDebug.writeImage(theta1360, "1_theta360");
        MiscDebug.writeImage(theta2360, "2_theta360");

        
        //log.info("corners1=" + corners1.toString());
        //log.info("corners2=" + corners2.toString());

        MiscDebug.plotCorners(img1Grey, corners1, "1_corners");
        MiscDebug.plotCorners(img2Grey, corners2, "2_corners");
        
        /*
        If make an assumption that the histogram equalization and then color
        segmentation leaves the images in consistent state w.r.t. similar 
        colors and intensities being displayed similarly, then selection of
        blobs in this way should lead to comparable lists (which are
        subsets of the total, that is the blobs found for intensity level i0
        in image 1 are the ones to match to similar intensity level in image 2).
        
        Using the blobs, that is the hulls above, has reduced the number of
        regions to compare.
        
        Using the centroids of the blobs is appealling, but these will likely
        be in regions free of large gradients so are difficult to distinguish.
        
        Corners filtered to be only those that cross those blob hulls, however,
        should be good features to use.
        
        The problem is that the corners filtered to those lists are still too
        many to use pairwise transformation calculation from every possible
        pair (the process is approx n^4).
        
        So feature matching using patches surrounding the corners is needed and 
        this is an N^2 process but accompanied by many steps for the comparisons.
        
        Note that feature matching instead of straight pairwise transformation
        calculation unfortunately needs an assumption made about the
        scale.
        
        TODO: For sets with < a dozen or two dozen members each, could 
        consider using pairwise calculation checked by feature matches.
        
        TODO: refactor to be able to use gaussian pyramids if scale=1 is not
              known, then fit to find scale (for sets with large number of members).
             see Lindenberg 1998 and Lowe 2004
        */
        
        List<PairIntArray> filteredCornersList1 = new ArrayList<PairIntArray>();
        List<PairIntArray> filteredCornersList2 = new ArrayList<PairIntArray>();
        filterCornersAndOrderByMatchingIntensity(
            corners1, hulls1,
            corners2, hulls2, 
            filteredCornersList1, filteredCornersList2);
        
        boolean useFiltered = true;
        
        int nFiltered1Tot = 0;
        int nFiltered2Tot = 0;
        for (PairIntArray p : filteredCornersList1) {
            nFiltered1Tot += p.getN();
        }
        for (PairIntArray p : filteredCornersList2) {
            nFiltered2Tot += p.getN();
        }
        
        if ((nFiltered1Tot < 12) || (nFiltered2Tot < 12)) {
            useFiltered = false;
        }
        
        /*
        use feature description to find similar looking features within
        image2 filtered corners.
        */
        
        /*
        key = intensity level of image 1 contiguous region
        value = map with
                key = point in filteredCorners1
                value = list of similar looking points in filteredCorners2
                        for the same intensity level.
        */
        
        List<Map<PairInt, List<FeatureComparisonStat>>> similarCorners
            = new ArrayList<Map<PairInt, List<FeatureComparisonStat>>>();
        
        if (useFiltered) {
            
            for (int i = 0; i < filteredCornersList1.size(); ++i) {
                // sets with keys pixValue1 and pixValues should hold common blobs
                // and details
                PairIntArray filtered1 = filteredCornersList1.get(i);
                PairIntArray filtered2 = filteredCornersList2.get(i);

                log.info("nFiltered1=" + filtered1.getN());
                log.info("nFiltered2=" + filtered2.getN());

                log.info("filtered1=" + filtered1.toString());
                log.info("filtered2=" + filtered2.toString());

                MiscDebug.plotCorners(img1Grey, filtered1, "1_" + i + "_filtered");
                MiscDebug.plotCorners(img2Grey, filtered2, "2_" + i + "_filtered");

                Set<CornerRegion> cornerRegions1 = 
                    detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
                CornerRegion[] cr1 = findCornerRegions(filtered1, cornerRegions1);
                Set<CornerRegion> cornerRegions2 = 
                    detector2.getEdgeCornerRegionsInOriginalReferenceFrame(true);
                CornerRegion[] cr2 = findCornerRegions(filtered2, cornerRegions2);
                
                //log.info("corner region1:");
                //MiscDebug.printCornerRegion(cr1);
                //log.info("corner region2:");
                //MiscDebug.printCornerRegion(cr2);
//DEBUG before placing in tests 
log.info("corner regions:");
StringBuilder sb = new StringBuilder();
for (int iii = 0; iii < cr1.length; ++iii) {
    if ((cr1[iii] == null) || (cr2[iii] == null)) {
        continue;
    }
    /*
    try {
        MiscDebug.display(cr1[iii], cr2[iii], img1Grey, img2Grey, 
            String.valueOf(iii), 20);
    } catch (Exception e) {
        
    }*/
    sb.append("1) ").append(cr1[iii].toString()).append(" 2)").append(cr2[iii].toString());
    int z = 1;
}
log.info(sb.toString());
                
                Map<PairInt, List<FeatureComparisonStat>> similar = 
                    findSimilarFeatures(
                        img1Cp, detector.getGradientXY(), theta1360, cr1,
                        img2Cp, detector2.getGradientXY(), theta2360, cr2
                    );
                    /*findSimilarFeatures(img1Cp, 
                        img2Cp, filtered1, filtered2);*/

                similarCorners.add(similar);
            }
            
        } else {
            
            Set<CornerRegion> cornerRegions1 = 
                detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
            CornerRegion[] cr1 = findCornerRegions(corners1, cornerRegions1);
            Set<CornerRegion> cornerRegions2 = 
                detector2.getEdgeCornerRegionsInOriginalReferenceFrame(true);
            CornerRegion[] cr2 = findCornerRegions(corners2, cornerRegions2);
            
            Map<PairInt, List<FeatureComparisonStat>> similar = 
                findSimilarFeatures(
                    img1Cp, detector.getGradientXY(), theta1360, cr1,
                    img2Cp, detector2.getGradientXY(), theta2360, cr2);
                //findSimilarFeatures(img1Cp, img2Cp, corners1, corners2);

            similarCorners.add(similar);
        }
        
if (true) {
    return null;
}        
        
        /*
        if need to calculate for all hull centroids, can use the rotation 
        and translation grid based search followed by downhill simplex:
        
        boolean useGreedyMatching = true;
        boolean searchScaleToo = true;
        float scale = 1;
        TransformationPointFit[] starterPoints = pointMatcher.preSearch0Small(
            allHullCentroids1, allHullCentroids2, scale,
            0, 359, useGreedyMatching);
        TransformationPointFit[] fits2 =
            pointMatcher.refineTransformationWithDownhillSimplexWrapper(
            starterPoints, allHullCentroids1, allHullCentroids2, searchScaleToo, 
            useGreedyMatching);
        
        log.info("fit3=" + fits2[0].toString());
        */
        
        TransformationPointFit bestFit00 = null;

        float toleranceTransX = 20;//30;
        float toleranceTransY = toleranceTransX;
        boolean useGreedyMatching = true;
        boolean earlyConvergeReturn = true;
        boolean setsAreMatched = false;
        
        log.info("nAllHulls1=" + allHullCentroids1.getN() 
            + " nAllHulls2=" + allHullCentroids2.getN());
        
        PointMatcher pointMatcher = new PointMatcher();
        
        if ((bestFit00 == null) ||
            (((float)bestFit00.getNumberOfMatchedPoints()
            /(float)bestFit00.getNMaxMatchable()) < 0.5)) {
                        
            // unfortunately, the quicker solution didn't work so need to try
            // to solve it using all corners (though the internal implementation
            // does not try all permutations).
            
            long t0 = System.currentTimeMillis();
            List<TransformationPointFit> fits
                = pointMatcher.calculateEuclideanTransformationUsingPairs(
                    corners1, corners2, toleranceTransX, toleranceTransY,
                    earlyConvergeReturn, useGreedyMatching);
            long t1 = System.currentTimeMillis();
            log.info("calculateEuclideanTransformationUsingPairs seconds="
                + ((t1 - t0) * 1e-3));
            
            // the evaluation has already been tried against all corners
            // so only need to compare the fits.
            
            // this finds best that is usually a scale that is not 1, so
            // if one knows that the scale should be '1', should use the
            // method that solves mostly vertical transformations specifically
            
            for (TransformationPointFit fit : fits) {
                if (pointMatcher.fitIsBetter(bestFit00, fit)) {
                    bestFit00 = fit;
                }
            }
            log.info("best from all corners=" + bestFit00.toString());
        }
        
        // refine the solution 
        float rotHalfRangeInDegrees = 20;
        float rotDeltaInDegrees = 2;
        float transXHalfRange = 40; 
        float transXDelta = 4;
        float transYHalfRange = transXHalfRange; 
        float transYDelta = transXDelta;
        
        long t0 = System.currentTimeMillis();
        TransformationPointFit bestFit = pointMatcher.refineTheTransformation(
            bestFit00.getParameters(), corners1, corners2,
            rotHalfRangeInDegrees, rotDeltaInDegrees,
            transXHalfRange, transXDelta,
            transYHalfRange, transYDelta,
            useGreedyMatching);
        long t1 = System.currentTimeMillis();
        
        log.info(((t1-t0)*1e-3) + "seconds) refined solution=" + bestFit.toString());
        
        //TODO: use the transformation to make matching corner lists, 
        //   perhaps from larger images
        
        return bestFit;
    }

    /**
     * match the given point lists by comparing rotated and dithered blocks
     * surrounding points in the first list to blocks around points in the
     * second list.
     * The runtime complexity is O(N1 * N2), but there are many steps for
     * the rotation and dither operations.
     * 
     * @return 
     */
    protected Map<PairInt, List<FeatureComparisonStat>> findSimilarFeatures(
        ImageExt img1, GreyscaleImage gXY1, GreyscaleImage theta1,
        CornerRegion[] cornerRegions1,
        ImageExt img2, GreyscaleImage gXY2, GreyscaleImage theta2,
        CornerRegion[] cornerRegions2) {
        
        if (img1 == null) {
            throw new IllegalArgumentException("img1 cannot be null");
        }
        if (img2 == null) {
            throw new IllegalArgumentException("img2 cannot be null");
        }
        if (gXY1 == null) {
            throw new IllegalArgumentException("gXY1 cannot be null");
        }
        if (gXY2 == null) {
            throw new IllegalArgumentException("gXY2 cannot be null");
        }
        if (theta1 == null) {
            throw new IllegalArgumentException("theta1 cannot be null");
        }
        if (theta2 == null) {
            throw new IllegalArgumentException("theta2 cannot be null");
        }
        if (cornerRegions1 == null) {
            throw new IllegalArgumentException("cornerRegions1 cannot be null");
        }
        if (cornerRegions2 == null) {
            throw new IllegalArgumentException("cornerRegions2 cannot be null");
        }
        
        int blockHalfWidth = 2;
        boolean useNormalizedIntensities = false;
        
        // TODO: consider using greyscale images
        Features features1 = new Features(img1, theta1, blockHalfWidth, 
            useNormalizedIntensities);
        
        Features features2 = new Features(img2, theta2, blockHalfWidth, 
            useNormalizedIntensities);
         
        /*
        comparisons between cornerRegions1 and cornerRegions2 are made to see
        if an individual region matches in intensity and then in the gradient
        of intensity.
        
        For the comparisons, there is small dither operation around the 
        central coordinate of each cornerRegions1 region and for each dither,
        there is also a variation of the calculated orientation by +- 25.
        
        The best match among the dithers and rotations for the original 
        cornerRegions1 instance is stored if a match was found.
        */
        
        Map<PairInt, Map<PairInt, Map<Float, FeatureComparisonStat>>> 
            comparisonMap = new HashMap<PairInt, Map<PairInt, 
            Map<Float, FeatureComparisonStat>>>();
            
        float[][] offsets0 = createNeighborOffsets();
        
        int dither = 1;
        
        int nF = (2 * dither + 1) * (2 * dither + 1);
        
        for (int idx1 = 0; idx1 < cornerRegions1.length; ++idx1) {
            
            CornerRegion cornerRegion1 = cornerRegions1[idx1];
            
            if (cornerRegion1 == null) {
                continue;
            }
            
            for (int idx2 = 0; idx2 < cornerRegions1.length; ++idx2) {
                
                CornerRegion cornerRegion2 = cornerRegions2[idx2];
            
                if (cornerRegion2 == null) {
                    continue;
                }
                
                try {
                    
                    FeatureComparisonStat stat = findBestAmongDitheredRotated(
                        features1, features2, cornerRegion1, cornerRegion2,
                        dither);
                
                    if (stat != null) {
                        //storeInMap(comparisonMap, p1, p2, best.getImg1RotInDegrees(),
                        //    stat);
                    }
                    
                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                    log.log(Level.SEVERE, null, ex);
                }
            }

        }
        
        return null;
    }
   
    /**
     * create x and y offsets for the neighbor points within 2 pixel radius.
     * The result is a two-dimensional array of length 25 with the first
     * dimension being the x array and the 2nd dimension being the y array.
     * Note that the offset of (0,0) is the first given.
     * @return 
     */
    public float[][] createNeighborOffsets() {
        
        //TODO: test if this needs to be higher for higher resolution data
        int d = 8;
        
        return Misc.createNeighborOffsets(d);
    }    
    
    private void storeInMap(Map<PairInt, Map<PairInt, Map<Float, 
        FeatureComparisonStat>>> comparisonMap, PairInt p1, PairInt p2, 
        float rotInDeg, FeatureComparisonStat best) {
        
        if (best == null) {
            return;
        }
        
        if (comparisonMap == null) {
            throw new IllegalArgumentException("comparisonMap cannot be null");
        }
        
        Map<PairInt, Map<Float, FeatureComparisonStat>> p1Map = 
            comparisonMap.get(p1);
        
        if (p1Map == null) {
            p1Map = new HashMap<PairInt, Map<Float, FeatureComparisonStat>>();
            comparisonMap.put(p1, p1Map);
        }
        
        Map<Float, FeatureComparisonStat> p2Map = p1Map.get(p2);
        
        if (p2Map == null) {
            p2Map = new HashMap<Float, FeatureComparisonStat>();
            p1Map.put(p2, p2Map);
        }
        
        Float key = Float.valueOf(rotInDeg);
        
        p2Map.put(key, best);
        
    }

    private Map<Integer, PairIntArray> filterToHulls(PairIntArray corners, 
        Map<Integer, List<GrahamScan>> hulls) {
        
        PointInPolygon poly = new PointInPolygon();
        
        Map<Integer, PairIntArray> filteredCorners = new HashMap<Integer, PairIntArray>();
                
        for (Entry<Integer, List<GrahamScan>> entry : hulls.entrySet()) {
            
            Integer pixValue = entry.getKey();
            
            List<GrahamScan> hullList = entry.getValue();
           
            Set<PairInt> added = new HashSet<PairInt>();
            
            PairIntArray cornersH = new PairIntArray();
            
            for (GrahamScan scan : hullList) {
                
                for (int i = 0; i < corners.getN(); ++i) {
                    
                    int x = corners.getX(i);
                    int y = corners.getY(i);
                    
                    PairInt p = new PairInt(x, y);
                    
                    if (added.contains(p)) {
                        continue;
                    }
                    
                    if (poly.isInSimpleCurve(x, y, scan.getXHull(), 
                        scan.getYHull(), scan.getXHull().length)) {
                        
                        cornersH.add(x, y);
                        added.add(p);
                    } else {
                        for (int ii = 0; ii < scan.getXHull().length; ++ii) {
                            float diffX = Math.abs(x - scan.getXHull()[ii]);
                            float diffY = Math.abs(y - scan.getYHull()[ii]);
                            if (diffX < 2 && diffY < 2) {
                                cornersH.add(x, y);
                                added.add(p);
                                break;
                            }
                        }
                    }
                }
            }
        
            // because convex hulls are used, there may be some points that are
            // in more than one map entry
            filteredCorners.put(pixValue, cornersH);
        }
        
        return filteredCorners;
    }

    protected void filterCornersAndOrderByMatchingIntensity(
        PairIntArray corners1, Map<Integer, List<GrahamScan>> hulls1, 
        PairIntArray corners2, Map<Integer, List<GrahamScan>> hulls2, 
        List<PairIntArray> outputFilteredCornersList1, 
        List<PairIntArray> outputFilteredCornersList2) {
        
        Map<Integer, PairIntArray> filteredCorners1 = filterToHulls(
            corners1, hulls1);
        
        Map<Integer, PairIntArray> filteredCorners2 = filterToHulls(
            corners2, hulls2);
        
        Set<Integer> matchedFilter2Keys = new HashSet<Integer>();
        
        // iterate over intensity keys to find the closest matching in intensities
        // to compare those lists
        for (Entry<Integer, PairIntArray> entry : filteredCorners1.entrySet()) {
            
            Integer pixValue1 = entry.getKey();
            
            // find the closest value to it in the image2 map.  should be very
            // similar if image content is similar and pre-processing the same.
            // else, may need to use allHullCentroids1 and allHullCentroids2
            // in one calculation
            
            Integer pixValue2 = null;
            Iterator<Entry<Integer, PairIntArray>> iter = 
                filteredCorners2.entrySet().iterator();
            while (iter.hasNext()) {
                Integer pV = iter.next().getKey();
                if (matchedFilter2Keys.contains(pV)) {
                    continue;
                }
                if (pixValue2 == null) {
                    pixValue2 = pV;
                } else {
                    int diff = Math.abs(pixValue1.intValue() - pixValue2.intValue());
                    int diffV = Math.abs(pV.intValue() - pixValue2.intValue());
                    if (diffV < diff) {
                        pixValue2 = pV;
                    }
                }
            }
            
            if (pixValue2 == null) {
                continue;
            }
            matchedFilter2Keys.add(pixValue2);
            
            // sets with keys pixValue1 and pixValues should hold common blobs
            // and details
            PairIntArray filtered1 = entry.getValue();
            PairIntArray filtered2 = filteredCorners2.get(pixValue2);
            
            outputFilteredCornersList1.add(filtered1);
            outputFilteredCornersList2.add(filtered2);
        }
        
    } 

    private CornerRegion[] findCornerRegions(PairIntArray corners, 
        Set<CornerRegion> cornerRegions) {
        
        int n = corners.getN();
        
        CornerRegion[] out = new CornerRegion[n];
        
        for (int i = 0; i < corners.getN(); ++i) {
            
            int x = corners.getX(i);
            int y = corners.getY(i);
            
            CornerRegion cr = null;
            int diffSqMin = Integer.MAX_VALUE;
            for (CornerRegion cornerRegion : cornerRegions) {
                int kMaxIdx = cornerRegion.getKMaxIdx();
                int xcr = cornerRegion.getX()[kMaxIdx];
                int ycr = cornerRegion.getY()[kMaxIdx];
                int diffX = Math.abs(xcr - x);
                int diffY = Math.abs(ycr - y);
                if (diffX > 4 || diffY > 4) {
                    continue;
                }
                int diffSq = diffX * diffX + diffY * diffY;
                if (diffSq < diffSqMin) {
                    diffSqMin = diffSq;
                    cr = cornerRegion;
                }
            }
            
            out[i] = cr ;
        }
        
        return out;
    }

    protected FeatureComparisonStat findBestAmongDitheredRotated(
        Features features1, Features features2, CornerRegion cornerRegion1, 
        CornerRegion cornerRegion2, int dither) throws CornerRegion.CornerRegionDegneracyException {
        
        FeatureComparisonStat best = null;
        
        int kMaxIdx1 = cornerRegion1.getKMaxIdx();
        final int x1 = cornerRegion1.getX()[kMaxIdx1];
        final int y1 = cornerRegion1.getY()[kMaxIdx1];
        int rot1 = (int)Math.round(
            cornerRegion1.getRelativeOrientationInDegrees());
        
        int kMaxIdx2 = cornerRegion2.getKMaxIdx();
        final int x2 = cornerRegion2.getX()[kMaxIdx2];
        final int y2 = cornerRegion2.getY()[kMaxIdx2];
        int rot2 = (int)Math.round(
            cornerRegion2.getRelativeOrientationInDegrees());
        
        IntensityDescriptor desc2 = features2.extractIntensity(x2, y2, rot2);
        
        for (int x1d = (x1 - dither); x1d <= (x1 + dither); ++x1d) {
            if (!features1.isWithinXBounds(x1d)) {
                continue;
            }
            for (int y1d = (y1 - dither); y1d <= (y1 + dither); ++y1d) {
                if (!features1.isWithinYBounds(y1d)) {
                    continue;
                }
                
                for (int rotD1 = (rot1 - 25); rotD1 <= (rot1 + 25); rotD1 += 25) {
                    
                    IntensityDescriptor desc1 = features1.extractIntensity(x1d, 
                        y1d, rotD1);
                    
                    FeatureComparisonStat stat = Features.calculateIntensityStat(
                        desc1, x1d, y1d, desc2, x2, y2);
                    
                    if (stat.getSumSqDiff() > stat.getImg2PointErr()) {
                        continue;
                    }
                    
                    if (best == null) {
                        best = stat;
                    } else {
                        if (best.getSumSqDiff() > stat.getSumSqDiff()) {
                            best = stat;
                        }
                    }
                }
            }
        }
        
        return best;
    }

}
