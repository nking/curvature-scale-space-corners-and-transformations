package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
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
        
        final boolean performBinning = true;
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
                        
                        int nBounds = 0;
                        for (int ii = 0; ii < scan.getXHull().length; ++ii) {
                            float xh = scan.getXHull()[ii];
                            float yh = scan.getYHull()[ii];
                            if ((xh == 0) || xh == (imgGrey.getWidth() - 1) ||
                                (yh == 0) || yh == (imgGrey.getHeight() - 1)) {
                                nBounds++;
                            }
                        }
                        if (nBounds > 1) {
                            
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
                
                contigMap.put(pixValue, list);
                hulls.put(pixValue, listHulls);
                hullCentroidsMap.put(pixValue, pvHullCentroids);
            }
        }
        
        // ===== debug: plot the hulls =======
        Image img1W = ImageIOHelper.convertImage(img1Grey);
        Image img2W = ImageIOHelper.convertImage(img2Grey);
        
        int c = 0;
        for (Entry<Integer, List<GrahamScan>> entry : hulls1.entrySet()) {
            List<GrahamScan> hulls = entry.getValue();
            for (GrahamScan hull : hulls) {
                int[] x = new int[hull.getXHull().length];
                int[] y = new int[x.length];
                for (int i = 0; i < x.length; ++i) {
                    x[i] = Math.round(hull.getXHull()[i]);
                    y[i] = Math.round(hull.getYHull()[i]);                   
                }
                if (c == 0) {
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 255, 0, 0);
                } else if (c == 1) {
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 0, 255, 0);
                } else {
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 0, 0, 255);
                }
            }
            c++;
        }
        c = 0;
        for (Entry<Integer, List<GrahamScan>> entry : hulls2.entrySet()) {
            List<GrahamScan> hulls = entry.getValue();
            for (GrahamScan hull : hulls) {
                int[] x = new int[hull.getXHull().length];
                int[] y = new int[x.length];
                for (int i = 0; i < x.length; ++i) {
                    x[i] = Math.round(hull.getXHull()[i]);
                    y[i] = Math.round(hull.getYHull()[i]);
                }
                if (c == 0) {
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 255, 0, 0);
                } else if (c == 1) {
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 0, 255, 0);
                } else {
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 0, 0, 255);
                }
            }
            c++;
        }
       
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img1_binned.png", img1W);
            ImageIOHelper.writeOutputImage(dirPath + "/img2_binned.png", img2W);
            
            ImageIOHelper.writeOutputImage(dirPath + "/img1_binned_clr.png", img1Cp);
            ImageIOHelper.writeOutputImage(dirPath + "/img2_binned_clr.png", img2Cp);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
        // ===== end debug plot of hulls =======
        
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
        
        // ===== plot the corners ====
        img1W = ImageIOHelper.convertImage(img1Grey);
        img2W = ImageIOHelper.convertImage(img2Grey);
        ImageIOHelper.addCurveToImage(corners1, img1W, 1, 255, 0, 0);
        ImageIOHelper.addCurveToImage(corners2, img2W, 1, 255, 0, 0);
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img1_corners.png", img1W);
            ImageIOHelper.writeOutputImage(dirPath + "/img2_corners.png", img2W);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
        
        log.info("n1Corners=" + corners1.getN() + " n2Corners2=" 
            + corners2.getN());
         
        /*
        this matching using the hulls only works pretty well, but also needs
        to use information from the surrounding region to rule out fits.
        For the Brown & Lowe panoramic images, for example, the diagonal 
        features are matching to diagonal features in the 2nd image 
        (matching for ridges rather than a smaller number of adjacent features).
        
        If make an assumption that the histogram equalization and then color
        segmentation leaves the images in consistent state w.r.t. similar 
        colors and intensities being displayed similarly, then this
        should be fine to compare the surrounding intensities.  
        Can do a cross centered on the centroid, rotated for the transformation
        angle being tested and compare the intensities and gradient.
        
        This would be making an assumption about scale or re-doing surrounding 
        sum when needed for different scale.
        
        Note that if that is not the case, the images should be pre-processed
        to make that true before being given to this method.
        */
        
        float toleranceTransX = 20;//30;
        float toleranceTransY = toleranceTransX;
        boolean useGreedyMatching = true;
        boolean earlyConvergeReturn = true;
        boolean setsAreMatched = false;
        
        log.info("nAllHulls1=" + allHullCentroids1.getN() 
            + " nAllHulls2=" + allHullCentroids2.getN());
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit bestFit00 = null;
        
        for (Entry<Integer, PairIntArray> entry : hullCentroids1Map.entrySet()) {
            
            Integer pixValue1 = entry.getKey();
            
            // find the closest value to it in the image2 map.  should be very
            // similar if image content is similar and pre-processing the same.
            // else, may need to use allHullCentroids1 and allHullCentroids2
            // in one calculation
            
            Integer pixValue2 = null;
            Iterator<Entry<Integer, PairIntArray>> iter = hullCentroids2Map.entrySet().iterator();
            while (iter.hasNext()) {
                Integer pV = iter.next().getKey();
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
            
            PairIntArray hull1Centroids = entry.getValue();
            PairIntArray hull2Centroids = hullCentroids2Map.get(pixValue2);
            
            log.info("nHulls1=" + hull1Centroids.getN() + " nHulls2=" 
                + hull2Centroids.getN() + " {" + pixValue1 + "," + pixValue2 + "}");
            
            double toleranceColor = 18;
            
            long t00 = System.currentTimeMillis();
            List<TransformationPointFit> fits0
                = pointMatcher.calculateEuclideanTransformationUsingPairs(
                img1Cp, img2Cp,
                hull1Centroids, hull2Centroids,
                toleranceTransX, toleranceTransY, toleranceColor,
                earlyConvergeReturn, useGreedyMatching);
            long t10 = System.currentTimeMillis();
            log.info("calculateEuclideanTransformationForSmallSets for hull centroids seconds="
                + ((t10 - t00) * 1e-3));
            
            for (TransformationPointFit fit : fits0) {
                
                log.info("fit=" + fit.toString());
                
                TransformationPointFit fit2 = pointMatcher.transformAndEvaluateFit(
                    corners1, corners2, fit.getParameters(), 
                    toleranceTransX, toleranceTransY, useGreedyMatching,
                    setsAreMatched);

                if (pointMatcher.fitIsBetter(bestFit00, fit2)) {
                    bestFit00 = fit2;
                    log.info("bestFit0=" + bestFit00.toString());
                }
            }
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

}
