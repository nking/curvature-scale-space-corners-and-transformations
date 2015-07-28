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
        
        boolean performBinning = false;
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
        
        PairIntArray hullCentroids1 = new PairIntArray();
        PairIntArray hullCentroids2 = new PairIntArray();
        
        for (int im = 0; im < 2; ++im) {
            
            Map<Integer, Integer> freqMap = freqMap1;
            Map<Integer, List<PairIntArray>> contigMap = contigMap1;
            Map<Integer, List<GrahamScan>> hulls = hulls1;
            GreyscaleImage imgGrey = img1Grey;
            PairIntArray hullCentroids = hullCentroids1;
                        
            if (im == 1) {
                freqMap = freqMap2;
                contigMap = contigMap2;
                hulls = hulls2;
                hullCentroids = hullCentroids2;
                imgGrey = img2Grey;
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

                contigMap.put(pixValue, list);
                
                List<GrahamScan> listHulls = new ArrayList<GrahamScan>();
                for (PairIntArray xy : list) {
                    GrahamScan scan = new GrahamScan();
                    try {
                        float[] x = new float[xy.getN()];
                        float[] y = new float[x.length];
                        for (int i = 0; i < x.length; ++i) {
                            x[i] = xy.getX(i);
                            y[i] = xy.getY(i);
                        }

                        scan.computeHull(x, y);

                        listHulls.add(scan);
                        
                        double[] centroidXY = curveHelper.calculateXYCentroids(
                            scan.getXHull(), scan.getYHull()); 
                        
                        hullCentroids.add((int)Math.round(centroidXY[0]), 
                            (int)Math.round(centroidXY[1]));
                        
                    } catch (GrahamScanTooFewPointsException e) {
                        log.severe(e.getMessage());
                    }                    
                }
                hulls.put(pixValue, listHulls);                
            }
        }
        
        float toleranceTransX = 50;//30;
        float toleranceTransY = toleranceTransX;
        boolean useGreedyMatching = true;
        
        PointMatcher pointMatcher = new PointMatcher();
        List<TransformationPointFit> fits = 
            pointMatcher.calculateEuclideanTransformationForSmallSets(
                hullCentroids1, hullCentroids2, toleranceTransX, toleranceTransY, 
                useGreedyMatching);
        for (TransformationPointFit fit : fits) {
            log.info("fit=" + fit.toString());
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
                    //ImageIOHelper.addCurveToImage(x, y, img1W, 1, 255, 0, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 255, 0, 0);
                } else if (c == 1) {
                    //ImageIOHelper.addCurveToImage(x, y, img1W, 1, 0, 255, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 0, 255, 0);
                } else {
                    //ImageIOHelper.addCurveToImage(x, y, img1W, 1, 0, 0, 255);
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
                    //ImageIOHelper.addCurveToImage(x, y, img2W, 1, 255, 0, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 255, 0, 0);
                } else if (c == 1) {
                    //ImageIOHelper.addCurveToImage(x, y, img2W, 1, 0, 255, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 0, 255, 0);
                } else {
                    //ImageIOHelper.addCurveToImage(x, y, img2W, 1, 0, 0, 255);
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 0, 0, 255);
                }
            }
            c++;
        }
       
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img1_binned.png", img1W);
            ImageIOHelper.writeOutputImage(dirPath + "/img2_binned.png", img2W);
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
        
        /*
        filter the corners to keep only those within hulls to see if the numbers
        are small enough for pairwise transformation calculations.
        */
        Set<PairInt> allCorners1 = Misc.convert(corners1);
        Set<PairInt> allCorners2 = Misc.convert(corners2);
        Set<PairInt> cornersWithinHulls1 = new HashSet<PairInt>();
        Set<PairInt> cornersWithinHulls2 = new HashSet<PairInt>();
        
        log.info("n1Corners=" + allCorners1.size() + " n2Corners2=" 
            + allCorners2.size());
        
        PointInPolygon pP = new PointInPolygon();
            
        for (int im = 0; im < 2; ++im) {
            
            Map<Integer, List<GrahamScan>> hullsMap = hulls1; 
            Set<PairInt> allCorners = allCorners1;
            Set<PairInt> cornersWithinHulls = cornersWithinHulls1;
            if (im == 1) {
                hullsMap = hulls2;
                allCorners = allCorners2;
                cornersWithinHulls = cornersWithinHulls2;
            }
            
            for (Entry<Integer, List<GrahamScan>> entry : hullsMap.entrySet()) {
                List<GrahamScan> hulls = entry.getValue();
                for (GrahamScan hull : hulls) {
                    Set<PairInt> rm = new HashSet<PairInt>();
                    for (PairInt p : allCorners) {
                        // try without a tolerance
                        boolean isInHull = pP.isInSimpleCurve(p.getX(), p.getY(),
                            hull.getXHull(), hull.getYHull(), 
                            hull.getXHull().length);
                        if (isInHull) {
                            cornersWithinHulls.add(p);
                            rm.add(p);
                        }
                    }
                    for (PairInt p : rm) {
                        allCorners.remove(p);
                    }
                }
            }
        }
        
        log.info("cornersHulls1=" + cornersWithinHulls1.size() 
            + " cornersHulls2=" + cornersWithinHulls2.size());
        
        // ===== plot the corners ====
        img1W = ImageIOHelper.convertImage(img1Grey);
        img2W = ImageIOHelper.convertImage(img2Grey);
        ImageIOHelper.addToImage(cornersWithinHulls1, 0, 0, img1W, 2, 255, 0, 0);
        ImageIOHelper.addToImage(cornersWithinHulls2, 0, 0, img2W, 2, 255, 0, 0);
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img1_corners_in_hulls.png", img1W);
            ImageIOHelper.writeOutputImage(dirPath + "/img2_corners_in_hulls.png", img2W);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
                
        return null;
        //throw new UnsupportedOperationException("not yet implemented");
    
        /*
        -- segmentation by cieXY color into 3 bands
        -- dfs contiguous find of each of the 3 intensity bands
           -- convex hull of the groups
           -- compare the hulls of one image to the hulls of the other image:  
              -- by area and single pixel intensity?  (the cie xy histogram skipping
                    algorithm may assign different colors to same feature, so should
                    probably ignore intensity unless know the images are very similar).
              -- by shape
                -- furthest pairs?
                -- a description of ellipticity?
                -- points on the hull or centroid matching enough between 
                   the 2 images to be used like "corners"?
              -- by an increasing rotation and a fixed rotational range to
                 only look at hulls that are within the projected beam of the
                 hull in image1 under consideration?
                 -- need to store the best for a hull and also the fits as the
                    rotations are tried, that is, dynamic programming to avoid
                    repeating calculations.
                    -- then the best of each hull can be quickly looked up for
                       the other hulls (i.e. was rotation=20 within top fits
                       for the neighboring hulls and with consistent translation...)
                       
        */
    }

}
