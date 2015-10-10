package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayDescendingComparator;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a place holder for many steps while determining a way to create
 * correspondence lists that makes it easier to use in tests until
 * can change the design.
 *
 * @author nichole
 */
public class BinSegmentationHelper {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected final String fileName1;
    protected final String fileName2;
    protected final String filePath1;
    protected final String filePath2;

    ImageExt img1Orig = null;
    ImageExt img2Orig = null;
    GreyscaleImage img1GreyOrig = null;
    GreyscaleImage img2GreyOrig = null;
    GreyscaleImage img1Grey = null;
    GreyscaleImage img2Grey = null;
    int binFactor = 1;

    private PairIntArray corners1 = null;
    private PairIntArray corners2 = null;

    private List<PairIntArray> filteredCornersList1 = null;
    private List<PairIntArray> filteredCornersList2 = null;
    private Set<CornerRegion> cornerRegions1 = null;
    private Set<CornerRegion> cornerRegions2 = null;
    private GreyscaleImage gXY1 = null;
    private GreyscaleImage gXY2 = null;
    private GreyscaleImage theta1 = null;
    private GreyscaleImage theta2 = null;
    
    private String fileName1Root = null;
    private String fileName2Root = null;

    public BinSegmentationHelper(String fileName1, String fileName2) throws 
        IOException, Exception {
        
        this.fileName1 = fileName1;
        this.fileName2 = fileName2;

        filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        img1Orig = ImageIOHelper.readImageExt(filePath1);

        filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        img2Orig = ImageIOHelper.readImageExt(filePath2);

        int idx = fileName1.lastIndexOf(".");
        fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        fileName2Root = fileName2.substring(0, idx);
    }

    public void applySteps0(boolean performSkySubtraction) throws IOException,
        NoSuchAlgorithmException {

        ImageProcessor imageProcessor = new ImageProcessor();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        ImageHelperForTests helper = new ImageHelperForTests(img1Orig, true);
        
        if (performSkySubtraction) {
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            PairIntArray outputSkyCentroid = new PairIntArray();
            // sky are the zeros in this:
            GreyscaleImage resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img1Orig,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);

            imageProcessor.multiplyBinary(img1Orig, resultMask);

            helper = new ImageHelperForTests(img2Orig, true);
            skylineExtractor = new SkylineExtractor();
            outputSkyCentroid = new PairIntArray();
            // sky are the zeros in this:
            resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img2Orig,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);
            imageProcessor.multiplyBinary(img2Orig, resultMask);
        }
        
        TransformationParameters params90 = new TransformationParameters();
        params90.setRotationInDegrees(90);
        params90.setOriginX(0);
        params90.setOriginY(0);
        params90.setTranslationX(0);
        params90.setTranslationY(img1Orig.getWidth() - 1);

        Transformer transformer = new Transformer();

        img1Orig = (ImageExt) transformer.applyTransformation(img1Orig, params90,
            img1Orig.getHeight(), img1Orig.getWidth());

        //---------------

        img1GreyOrig = img1Orig.copyToGreyscale();
        img2GreyOrig = img2Orig.copyToGreyscale();

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

        ImageExt img1Cp = (ImageExt)img1Orig.copyImage();
        ImageExt img2Cp = (ImageExt)img2Orig.copyImage();

        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(img1GreyOrig, true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(img2GreyOrig, true);

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
            HistogramEqualization hEq = new HistogramEqualization(img1GreyOrig);
            hEq.applyFilter();
            hEq = new HistogramEqualization(img2GreyOrig);
            hEq.applyFilter();
            /*HistogramEqualizationForColor hEqC = new HistogramEqualizationForColor(img1Cp);
            hEqC.applyFilter();
            hEqC = new HistogramEqualizationForColor(img2Cp);
            hEqC.applyFilter();*/
        }

        if (performBinning) {
            float upperLimit = 300.f;
            binFactor1 = (int) Math.ceil(
                Math.max((float)img1GreyOrig.getWidth()/upperLimit,
                (float)img2GreyOrig.getHeight()/upperLimit));
            smallestGroupLimit /= (binFactor1*binFactor1);
            largestGroupLimit /= (binFactor1*binFactor1);

            log.info("binFactor1=" + binFactor1);

            // prevent from being smaller than needed for a convex hull
            if (smallestGroupLimit < 4) {
                smallestGroupLimit = 4;
            }
            img1GreyOrig = imageProcessor.binImage(img1GreyOrig, binFactor1);
            img2GreyOrig = imageProcessor.binImage(img2GreyOrig, binFactor1);
            img1Cp = imageProcessor.binImage(img1Cp, binFactor1);
            img2Cp = imageProcessor.binImage(img2Cp, binFactor1);
        }

        imageSegmentation.applyUsingKMPP(img1GreyOrig, kN);
        imageSegmentation.applyUsingKMPP(img2GreyOrig, kN);

        // == contiguous regions within size limits become blobs of interest,
        //    indexed by their intensity levels

        Map<Integer, Integer> freqMap1 = Histogram.createAFrequencyMap(img1GreyOrig);
        Map<Integer, Integer> freqMap2 = Histogram.createAFrequencyMap(img2GreyOrig);

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
            GreyscaleImage imgGrey = img1GreyOrig;
            PairIntArray hullCentroids = allHullCentroids1;
            Map<Integer, PairIntArray> hullCentroidsMap = hullCentroids1Map;
            if (im == 1) {
                freqMap = freqMap2;
                contigMap = contigMap2;
                hulls = hulls2;
                hullCentroids = allHullCentroids2;
                imgGrey = img2GreyOrig;
                hullCentroidsMap = hullCentroids2Map;
            }

            for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()) {

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

        MiscDebug.writeHullImages(img1GreyOrig, hulls1, "1_binned_hulls");
        MiscDebug.writeHullImages(img2GreyOrig, hulls2, "2_binned_hulls");
        MiscDebug.writeImage(img1Cp, "1_binned_clr");
        MiscDebug.writeImage(img2Cp, "2_binned_clr");

        // make corners

        if (!performBinning) {
            imageProcessor.blur(img1GreyOrig, 2);
            imageProcessor.blur(img2GreyOrig, 2);
        }

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1GreyOrig);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
        corners1 = detector.getCornersInOriginalReferenceFrame();
        cornerRegions1 = detector.getEdgeCornerRegions(true);
        //cornerRegions1 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY1 = detector.getEdgeFilterProducts().getGradientXY();
        img1Grey = img1GreyOrig.copyImage();
        imageProcessor.shrinkImage(img1Grey, 
            new int[]{gXY1.getXRelativeOffset(), gXY1.getYRelativeOffset(),
                gXY1.getWidth(), gXY1.getHeight()
            });
        MiscDebug.writeEdges(detector.getEdgesInOriginalReferenceFrame(), 
            img1GreyOrig, "1_edges");

        CurvatureScaleSpaceCornerDetector detector2 = new
            CurvatureScaleSpaceCornerDetector(img2GreyOrig);
        detector2.doNotPerformHistogramEqualization();
        detector2.findCorners();
        corners2 = detector2.getCornersInOriginalReferenceFrame();
        cornerRegions2 = detector2.getEdgeCornerRegions(true);
        //cornerRegions2 = detector2.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY2 = detector2.getEdgeFilterProducts().getGradientXY();
        img2Grey = img2GreyOrig.copyImage();
        imageProcessor.shrinkImage(img2Grey, 
            new int[]{gXY2.getXRelativeOffset(), gXY2.getYRelativeOffset(),
                gXY2.getWidth(), gXY2.getHeight()
            });
        MiscDebug.writeEdges(detector2.getEdgesInOriginalReferenceFrame(), 
            img2GreyOrig, "2_edges");

        log.info("n1Corners=" + corners1.getN() + " n2Corners2=" + corners2.getN());

        // experimenting with a slightly different definition for theta:
        theta1 = imageProcessor.computeTheta360(
            detector.getEdgeFilterProducts().getGradientX(), 
            detector.getEdgeFilterProducts().getGradientY());
        theta2 = imageProcessor.computeTheta360(
            detector2.getEdgeFilterProducts().getGradientX(), 
            detector2.getEdgeFilterProducts().getGradientY());
        MiscDebug.writeImage(theta1, "1_theta360");
        MiscDebug.writeImage(theta2, "2_theta360");

        //log.info("corners1=" + corners1.toString());
        //log.info("corners2=" + corners2.toString());

        MiscDebug.plotCorners(img1GreyOrig, corners1, "1_corners");
        MiscDebug.plotCorners(img2GreyOrig, corners2, "2_corners");

        ShapeMatcher shapeMatcher = new ShapeMatcher();

        filteredCornersList1 = new ArrayList<PairIntArray>();
        filteredCornersList2 = new ArrayList<PairIntArray>();
        shapeMatcher.filterCornersAndOrderByMatchingIntensity(
            corners1, hulls1,
            corners2, hulls2,
            filteredCornersList1, filteredCornersList2);

    }

    public void applySteps1(boolean subtractSky) throws IOException, NoSuchAlgorithmException {

        ImageHelperForTests helper = new ImageHelperForTests(img1Orig, true);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        if (subtractSky) {
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            PairIntArray outputSkyCentroid = new PairIntArray();
            // sky are the zeros in this:
            GreyscaleImage resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img1Orig,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);            
            imageProcessor.multiplyBinary(img1Orig, resultMask);

            helper = new ImageHelperForTests(img2Orig, true);
            skylineExtractor = new SkylineExtractor();
            outputSkyCentroid = new PairIntArray();
            // sky are the zeros in this:
            resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img2Orig,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);
            imageProcessor.multiplyBinary(img2Orig, resultMask);
        }
        
        TransformationParameters params90 = new TransformationParameters();
        params90.setRotationInDegrees(90);
        params90.setOriginX(0);
        params90.setOriginY(0);
        params90.setTranslationX(0);
        params90.setTranslationY(img1Orig.getWidth() - 1);

        Transformer transformer = new Transformer();

        img1Orig = (ImageExt) transformer.applyTransformation(img1Orig, params90,
            img1Orig.getHeight(), img1Orig.getWidth());

        //---------------

        img1GreyOrig = img1Orig.copyToGreyscale();
        img2GreyOrig = img2Orig.copyToGreyscale();

        final boolean performBinning = false;
        int binFactor1 = 1;
        
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

        ImageExt img1Cp = (ImageExt)img1Orig.copyImage();
        ImageExt img2Cp = (ImageExt)img2Orig.copyImage();

        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(img1GreyOrig, true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(img2GreyOrig, true);

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
            HistogramEqualization hEq = new HistogramEqualization(img1GreyOrig);
            hEq.applyFilter();
            hEq = new HistogramEqualization(img2GreyOrig);
            hEq.applyFilter();
            /*HistogramEqualizationForColor hEqC = new HistogramEqualizationForColor(img1Cp);
            hEqC.applyFilter();
            hEqC = new HistogramEqualizationForColor(img2Cp);
            hEqC.applyFilter();*/
        }

        if (performBinning) {
            binFactor1 = (int) Math.ceil(
                Math.max((float)img1GreyOrig.getWidth()/200.f,
                (float)img2GreyOrig.getHeight()/200.));
            smallestGroupLimit /= (binFactor1*binFactor1);
            largestGroupLimit /= (binFactor1*binFactor1);

            log.info("binFactor1=" + binFactor1);

            // prevent from being smaller than needed for a convex hull
            if (smallestGroupLimit < 4) {
                smallestGroupLimit = 4;
            }
            img1GreyOrig = imageProcessor.binImage(img1GreyOrig, binFactor1);
            img2GreyOrig = imageProcessor.binImage(img2GreyOrig, binFactor1);
            img1Cp = imageProcessor.binImage(img1Cp, binFactor1);
            img2Cp = imageProcessor.binImage(img2Cp, binFactor1);
        }

        // make corners

        if (!performBinning) {
            imageProcessor.blur(img1GreyOrig, SIGMA.ONE);
            imageProcessor.blur(img2GreyOrig, SIGMA.ONE);
        }
        
log.info("img1Grey.w=" + img1GreyOrig.getWidth() + " img1Grey.h=" + img1GreyOrig.getHeight());
log.info("img2Grey.w=" + img2GreyOrig.getWidth() + " img2Grey.h=" + img2GreyOrig.getHeight());

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1GreyOrig);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();        
        corners1 = detector.getCornersInOriginalReferenceFrame();
        cornerRegions1 = detector.getEdgeCornerRegions(true);
        //cornerRegions1 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY1 = detector.getEdgeFilterProducts().getGradientXY();
        img1Grey = img1GreyOrig.copyImage();
        imageProcessor.shrinkImage(img1Grey, 
            new int[]{gXY1.getXRelativeOffset(), gXY1.getYRelativeOffset(),
                gXY1.getWidth(), gXY1.getHeight()
            });
        MiscDebug.writeEdges(detector.getEdgesInOriginalReferenceFrame(), 
            img1GreyOrig, "1_edges");
        MiscDebug.writeImage(img1Cp, "1_clr");
        MiscDebug.plotCorners(img1GreyOrig, corners1, "1__corners");
        
        CurvatureScaleSpaceCornerDetector detector2 = new
            CurvatureScaleSpaceCornerDetector(img2GreyOrig);
        detector2.doNotPerformHistogramEqualization();
        detector2.findCorners();
        corners2 = detector2.getCornersInOriginalReferenceFrame();
        cornerRegions2 = detector2.getEdgeCornerRegions(true);
        //cornerRegions2 = detector2.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY2 = detector2.getEdgeFilterProducts().getGradientXY();
        img2Grey = img2GreyOrig.copyImage();
        imageProcessor.shrinkImage(img2Grey, 
            new int[]{gXY2.getXRelativeOffset(), gXY2.getYRelativeOffset(),
                gXY2.getWidth(), gXY2.getHeight()
            });
        MiscDebug.writeEdges(detector2.getEdgesInOriginalReferenceFrame(), 
            img2GreyOrig, "2_edges");
        MiscDebug.writeImage(img2Cp, "2_clr");
        MiscDebug.plotCorners(img2GreyOrig, corners2, "2__corners");
        
        log.info("n1Corners=" + corners1.getN() + " n2Corners2=" + corners2.getN());

        // experimenting with a slightly different definition for theta:
        theta1 = imageProcessor.computeTheta360(
            detector.getEdgeFilterProducts().getGradientX(), 
            detector.getEdgeFilterProducts().getGradientY());
        theta2 = imageProcessor.computeTheta360(
            detector2.getEdgeFilterProducts().getGradientX(), 
            detector2.getEdgeFilterProducts().getGradientY());
        
        MiscDebug.writeImage(img1Grey, "1_greyscale_trimmed");
        MiscDebug.writeImage(img2Grey, "2_greyscale_trimmed");
        
        MiscDebug.writeImage(theta1, "1_theta360_trimmed");
        MiscDebug.writeImage(theta2, "2_theta360_trimmed");
        
        MiscDebug.writeImage(gXY1, "1_gXY_trimmed");
        MiscDebug.writeImage(gXY2, "2_gXY_trimmed");
     
        /*
        //TEMP files to visualize an offset:
        GreyscaleImage theta1_tmp = theta1.copyImage();
        GreyscaleImage theta2_tmp = theta2.copyImage();
        for (int col = 0; col < theta1_tmp.getWidth(); ++col) {
            for (int row = 0; row < theta1_tmp.getHeight(); ++row) {
                int v = theta1_tmp.getValue(col, row);
                v -= 308;
                if (v < 0) {
                    v += 360;
                } else if (v > 359) {
                    v = v - 360;
                }
                theta1_tmp.setValue(col, row, v);
            }
        }
        for (int col = 0; col < theta2_tmp.getWidth(); ++col) {
            for (int row = 0; row < theta2_tmp.getHeight(); ++row) {
                int v = theta2_tmp.getValue(col, row);
                v -= 225;
                if (v < 0) {
                    v += 360;
                } else if (v > 359) {
                    v = v - 360;
                }
                theta2_tmp.setValue(col, row, v);
            }
        }
        MiscDebug.writeImage(theta1_tmp, "1_theta360_tmp");
        MiscDebug.writeImage(theta2_tmp, "2_theta360_tmp");
        */
        
            
        //log.info("corners1=" + corners1.toString());
        //log.info("corners2=" + corners2.toString());

        MiscDebug.plotCorners(img1GreyOrig, corners1, "1_corners");
        MiscDebug.plotCorners(img2GreyOrig, corners2, "2_corners");

        int z = 1;
    }

    public void applySteps2(boolean performSkySubtraction) throws IOException,
        NoSuchAlgorithmException {

        ImageProcessor imageProcessor = new ImageProcessor();

        ImageHelperForTests helper = new ImageHelperForTests(img1Orig, true);

        if (performSkySubtraction) {
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            // sky are the zeros in this:
            img1Orig = skylineExtractor.createSkyMaskedImage(
                helper.getTheta(), helper.getGradientXY(), img1Orig,
                helper.getCannyEdgeFilterSettings());

            helper = new ImageHelperForTests(img2Orig, true);
            skylineExtractor = new SkylineExtractor();
            img2Orig = skylineExtractor.createSkyMaskedImage(
                helper.getTheta(), helper.getGradientXY(), img2Orig,
                helper.getCannyEdgeFilterSettings());
        }
        
        //---------------

        img1GreyOrig = img1Orig.copyToGreyscale();
        img2GreyOrig = img2Orig.copyToGreyscale();

        final boolean performBinning = true;
        int binFactor1 = 1;

        boolean performSegmentation = true;
        int kN2 = 4;
        boolean performBinarySegmentation = true;
        if (performBinarySegmentation) {
            kN2 = 2;
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

        ImageExt img1Cp = (ImageExt)img1Orig.copyImage();
        ImageExt img2Cp = (ImageExt)img2Orig.copyImage();

        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(img1GreyOrig, true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(img2GreyOrig, true);

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
            HistogramEqualization hEq = new HistogramEqualization(img1GreyOrig);
            hEq.applyFilter();
            hEq = new HistogramEqualization(img2GreyOrig);
            hEq.applyFilter();
            /*HistogramEqualizationForColor hEqC = new HistogramEqualizationForColor(img1Cp);
            hEqC.applyFilter();
            hEqC = new HistogramEqualizationForColor(img2Cp);
            hEqC.applyFilter();*/
        }

        if (performBinning) {
            float minDimension = 300.f;//200.f
            binFactor1 = (int) Math.ceil(
                Math.max((float)img1GreyOrig.getWidth()/minDimension,
                (float)img2GreyOrig.getHeight()/minDimension));
            smallestGroupLimit /= (binFactor1*binFactor1);
            largestGroupLimit /= (binFactor1*binFactor1);

            log.info("binFactor1=" + binFactor1);

            // prevent from being smaller than needed for a convex hull
            if (smallestGroupLimit < 4) {
                smallestGroupLimit = 4;
            }
            img1GreyOrig = imageProcessor.binImage(img1GreyOrig, binFactor1);
            img2GreyOrig = imageProcessor.binImage(img2GreyOrig, binFactor1);
            img1Cp = imageProcessor.binImage(img1Cp, binFactor1);
            img2Cp = imageProcessor.binImage(img2Cp, binFactor1);
        }
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        if (performSegmentation || performBinarySegmentation) {
            imageSegmentation.applyUsingKMPP(img1GreyOrig, kN2);
            imageSegmentation.applyUsingKMPP(img2GreyOrig, kN2);
        }

        // == contiguous regions within size limits become blobs of interest,
        //    indexed by their intensity levels

        Map<Integer, Integer> freqMap1 = Histogram.createAFrequencyMap(img1GreyOrig);
        Map<Integer, Integer> freqMap2 = Histogram.createAFrequencyMap(img2GreyOrig);

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
            GreyscaleImage imgGrey = img1GreyOrig;
            PairIntArray hullCentroids = allHullCentroids1;
            Map<Integer, PairIntArray> hullCentroidsMap = hullCentroids1Map;
            if (im == 1) {
                freqMap = freqMap2;
                contigMap = contigMap2;
                hulls = hulls2;
                hullCentroids = allHullCentroids2;
                imgGrey = img2GreyOrig;
                hullCentroidsMap = hullCentroids2Map;
            }

            for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()) {

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
                            if ((xh == 0) || xh >= (imgGrey.getWidth() - 2) ||
                                (yh == 0) || yh >= (imgGrey.getHeight() - 2)) {
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

                log.info("nHulls" + (im + 1) + "=" + listHulls.size() 
                    + " for intensity=" + pixValue.toString());

                contigMap.put(pixValue, list);
                hulls.put(pixValue, listHulls);
                hullCentroidsMap.put(pixValue, pvHullCentroids);
            }
        }
        
        MiscDebug.writeHullImages(img1GreyOrig, hulls1, "1_" + fileName1Root + "_binned_hulls");
        MiscDebug.writeHullImages(img2GreyOrig, hulls2, "2_" + fileName2Root + "_binned_hulls");
        MiscDebug.writeImage(img1Cp, "1_" + fileName1Root + "_binned_clr");
        MiscDebug.writeImage(img2Cp, "2_" + fileName2Root + "_binned_clr");

        List<PairIntArray> edges1 = new ArrayList<PairIntArray>();
        List<PairIntArray> edges2 = new ArrayList<PairIntArray>();
        
        // iterate over contigMap to create scale space maps
       
        /*
        for images close in time and from the same camera, can reduce the lists
        to compare against by the keys of the contigMaps, but because the code
        may be used on other conditions, will compare all values for all keys.
        TODO: consider a flag in the final enclosing class to have streamlined
        processes when one knows the images are stereo, etc.
        */
                
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        for (int type = 0; type < 2; ++type) {
            
            Map<Integer, List<PairIntArray>> contigMap = contigMap1;
            GreyscaleImage imgGrey = img1GreyOrig;
            List<PairIntArray> edges = edges1;
            if (type == 1) {
                contigMap = contigMap2;
                imgGrey = img2GreyOrig;
                edges = edges2;
            }
                        
            for (Entry<Integer, List<PairIntArray>> entry : contigMap.entrySet()) {
               
                List<PairIntArray> contiguousPixels = entry.getValue();
                
                for (PairIntArray contiguous : contiguousPixels) {
                    
                    Set<PairInt> points = Misc.convert(contiguous);
                    
                    EdgeExtractorForBlobBorder extractor = new EdgeExtractorForBlobBorder();
/*
// perist blobs to fix the edge extractor for one as a test:
double[] xyCen = curveHelper.calculateXYCentroids(points);
String fileName = String.format("_%d_%d.dat", (int)Math.round(xyCen[0]),
(int)Math.round(xyCen[1]));
if (type == 0) {
    fileName = "blobs1" + fileName; 
} else {
    fileName = "blobs2" + fileName;
}
Misc.persistToFile(fileName, points);
*/
                    PairIntArray closedEdge = extractor.extractAndOrderTheBorder0(
                        points, imgGrey.getWidth(), imgGrey.getHeight(), 
                        discardWhenCavityIsSmallerThanBorder);
        
                    if (closedEdge != null) {
                        edges.add(closedEdge);
                    }
                }
            }
        }
        
        log.info("nEdges1=" + edges1.size() + " nEdges2=" + edges2.size());
        Collections.sort(edges1, new PairIntArrayDescendingComparator());
        Collections.sort(edges2, new PairIntArrayDescendingComparator());
        
        removeEdgesShorterThan(edges1, 3);
        removeEdgesShorterThan(edges2, 3);
        /*
        // only keep largest 510?
        int n1 = (edges1.size() > 10) ? 10 : edges1.size();
        int n2 = (edges2.size() > 10) ? 10 : edges2.size();
        edges1 = edges1.subList(0, n1);
        edges2 = edges2.subList(0, n2);
        */
        
Image img0 = ImageIOHelper.convertImage(img1GreyOrig);
for (int i = 0; i < edges1.size(); ++i) {
    PairIntArray pa = edges1.get(i);
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
img0 = ImageIOHelper.convertImage(img2GreyOrig);
for (int i = 0; i < edges2.size(); ++i) {
    PairIntArray pa = edges2.get(i);
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

        CurvatureScaleSpaceInflectionEdgeMapper mapper = 
            new CurvatureScaleSpaceInflectionEdgeMapper(edges1, edges2,
            img1GreyOrig.getXRelativeOffset(), img1GreyOrig.getYRelativeOffset(),
            img2GreyOrig.getXRelativeOffset(), img2GreyOrig.getYRelativeOffset());
mapper.debugImg1 = img1GreyOrig;
mapper.debugImg2 = img2GreyOrig;

        mapper.useDebugMode();

  //    mapper.setToRefineTransformations();

        TransformationParameters transformationParams =
            mapper.createEuclideanTransformation();

        if (transformationParams != null) {
            log.info("params from contours=" + transformationParams.toString());
        }
    }

    public void applySteps3() throws IOException, NoSuchAlgorithmException {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        TransformationParameters params90 = new TransformationParameters();
        params90.setRotationInDegrees(90);
        params90.setOriginX(0);
        params90.setOriginY(0);
        params90.setTranslationX(0);
        params90.setTranslationY(img1Orig.getWidth() - 1);

        Transformer transformer = new Transformer();

        img1Orig = (ImageExt) transformer.applyTransformation(img1Orig, params90,
            img1Orig.getHeight(), img1Orig.getWidth());

        //---------------

        img1GreyOrig = img1Orig.copyToGreyscale();
        img2GreyOrig = img2Orig.copyToGreyscale();

        final boolean performBinning = false;
        int binFactor1 = 1;
        
        int smallestGroupLimit = 100;
        //TODO: consider scaling this by image size or by size and res if one
        //  day that information is passed to this method
        int largestGroupLimit = 5000;

        ImageExt img1Cp = (ImageExt)img1Orig.copyImage();
        ImageExt img2Cp = (ImageExt)img2Orig.copyImage();

        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(img1GreyOrig, true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(img2GreyOrig, true);

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
            HistogramEqualization hEq = new HistogramEqualization(img1GreyOrig);
            hEq.applyFilter();
            hEq = new HistogramEqualization(img2GreyOrig);
            hEq.applyFilter();
            /*HistogramEqualizationForColor hEqC = new HistogramEqualizationForColor(img1Cp);
            hEqC.applyFilter();
            hEqC = new HistogramEqualizationForColor(img2Cp);
            hEqC.applyFilter();*/
        }

        if (performBinning) {
            binFactor1 = (int) Math.ceil(
                Math.max((float)img1GreyOrig.getWidth()/200.f,
                (float)img2GreyOrig.getHeight()/200.));
            smallestGroupLimit /= (binFactor1*binFactor1);
            largestGroupLimit /= (binFactor1*binFactor1);

            log.info("binFactor1=" + binFactor1);

            // prevent from being smaller than needed for a convex hull
            if (smallestGroupLimit < 4) {
                smallestGroupLimit = 4;
            }
            img1GreyOrig = imageProcessor.binImage(img1GreyOrig, binFactor1);
            img2GreyOrig = imageProcessor.binImage(img2GreyOrig, binFactor1);
            img1Cp = imageProcessor.binImage(img1Cp, binFactor1);
            img2Cp = imageProcessor.binImage(img2Cp, binFactor1);
        }

        // make corners

        if (!performBinning) {
            imageProcessor.blur(img1GreyOrig, SIGMA.ONE);
            imageProcessor.blur(img2GreyOrig, SIGMA.ONE);
        }
        
log.info("img1Grey.w=" + img1GreyOrig.getWidth() + " img1Grey.h=" + img1GreyOrig.getHeight());
log.info("img2Grey.w=" + img2GreyOrig.getWidth() + " img2Grey.h=" + img2GreyOrig.getHeight());

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1GreyOrig);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();        
        corners1 = detector.getCornersInOriginalReferenceFrame();
        cornerRegions1 = detector.getEdgeCornerRegions(true);
        //cornerRegions1 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY1 = detector.getEdgeFilterProducts().getGradientXY();
        img1Grey = img1GreyOrig.copyImage();
        imageProcessor.shrinkImage(img1Grey, 
            new int[]{gXY1.getXRelativeOffset(), gXY1.getYRelativeOffset(),
                gXY1.getWidth(), gXY1.getHeight()
            });
        MiscDebug.writeEdges(detector.getEdgesInOriginalReferenceFrame(), 
            img1GreyOrig, "1_edges");
        MiscDebug.writeImage(img1Cp, "1_clr");
        MiscDebug.plotCorners(img1GreyOrig, corners1, "1__corners");
        
        CurvatureScaleSpaceCornerDetector detector2 = new
            CurvatureScaleSpaceCornerDetector(img2GreyOrig);
        detector2.doNotPerformHistogramEqualization();
        detector2.findCorners();
        corners2 = detector2.getCornersInOriginalReferenceFrame();
        cornerRegions2 = detector2.getEdgeCornerRegions(true);
        //cornerRegions2 = detector2.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY2 = detector2.getEdgeFilterProducts().getGradientXY();
        img2Grey = img2GreyOrig.copyImage();
        imageProcessor.shrinkImage(img2Grey, 
            new int[]{gXY2.getXRelativeOffset(), gXY2.getYRelativeOffset(),
                gXY2.getWidth(), gXY2.getHeight()
            });
        MiscDebug.writeEdges(detector2.getEdgesInOriginalReferenceFrame(), 
            img2GreyOrig, "2_edges");
        MiscDebug.writeImage(img2Cp, "2_clr");
        MiscDebug.plotCorners(img2GreyOrig, corners2, "2__corners");
        
        log.info("n1Corners=" + corners1.getN() + " n2Corners2=" + corners2.getN());

        // experimenting with a slightly different definition for theta:
        theta1 = imageProcessor.computeTheta360(
            detector.getEdgeFilterProducts().getGradientX(), 
            detector.getEdgeFilterProducts().getGradientY());
        theta2 = imageProcessor.computeTheta360(
            detector2.getEdgeFilterProducts().getGradientX(), 
            detector2.getEdgeFilterProducts().getGradientY());
        
        MiscDebug.writeImage(img1Grey, "1_greyscale_trimmed");
        MiscDebug.writeImage(img2Grey, "2_greyscale_trimmed");
        
        MiscDebug.writeImage(theta1, "1_theta360_trimmed");
        MiscDebug.writeImage(theta2, "2_theta360_trimmed");
        
        MiscDebug.writeImage(gXY1, "1_gXY_trimmed");
        MiscDebug.writeImage(gXY2, "2_gXY_trimmed");
     
        /*
        //TEMP files to visualize an offset:
        GreyscaleImage theta1_tmp = theta1.copyImage();
        GreyscaleImage theta2_tmp = theta2.copyImage();
        for (int col = 0; col < theta1_tmp.getWidth(); ++col) {
            for (int row = 0; row < theta1_tmp.getHeight(); ++row) {
                int v = theta1_tmp.getValue(col, row);
                v -= 308;
                if (v < 0) {
                    v += 360;
                } else if (v > 359) {
                    v = v - 360;
                }
                theta1_tmp.setValue(col, row, v);
            }
        }
        for (int col = 0; col < theta2_tmp.getWidth(); ++col) {
            for (int row = 0; row < theta2_tmp.getHeight(); ++row) {
                int v = theta2_tmp.getValue(col, row);
                v -= 225;
                if (v < 0) {
                    v += 360;
                } else if (v > 359) {
                    v = v - 360;
                }
                theta2_tmp.setValue(col, row, v);
            }
        }
        MiscDebug.writeImage(theta1_tmp, "1_theta360_tmp");
        MiscDebug.writeImage(theta2_tmp, "2_theta360_tmp");
        */
        
            
        //log.info("corners1=" + corners1.toString());
        //log.info("corners2=" + corners2.toString());

        MiscDebug.plotCorners(img1GreyOrig, corners1, "1_corners");
        MiscDebug.plotCorners(img2GreyOrig, corners2, "2_corners");

        int z = 1;
    }
    
    protected void removeEdgesShorterThan(List<PairIntArray> output, 
        int minNumberOfPixelsInEdge) {
        
        for (int i = (output.size() - 1); i > -1; i--) {
            if (output.get(i).getN() < minNumberOfPixelsInEdge) {
                output.remove(i);
            }
        }
    }
    
    /**
     * @return the corners1
     */
    public PairIntArray getCorners1InOriginalReferenceFrame() {
        return corners1;
    }

    /**
     * @return the corners2
     */
    public PairIntArray getCorners2InOriginalReferenceFrame() {
        return corners2;
    }

    /**
     * @return the filteredCornersList1
     */
    public List<PairIntArray> getFilteredCornersList1() {
        return filteredCornersList1;
    }

    /**
     * @return the filteredCornersList2
     */
    public List<PairIntArray> getFilteredCornersList2() {
        return filteredCornersList2;
    }

    /**
     * @return the cornerRegions1
     */
    public Set<CornerRegion> getCornerRegions1() {
        return cornerRegions1;
    }

    /**
     * @return the cornerRegions2
     */
    public Set<CornerRegion> getCornerRegions2() {
        return cornerRegions2;
    }

    public void createSortedCornerRegions() {
        int n1 = cornerRegions1.size();
        int n2 = cornerRegions2.size();
        float[] k1 = new float[n1];
        float[] k2 = new float[n2];
        PairInt[] p1 = new PairInt[n1];
        PairInt[] p2 = new PairInt[n2];
        
        PairIntArray rCorners1 = new PairIntArray();
        PairIntArray rCorners2 = new PairIntArray();
    
        int count = 0;
        for (CornerRegion cr : cornerRegions1) {
            int kMaxIdx = cr.getKMaxIdx();
            k1[count] = cr.getK()[kMaxIdx];
            int x = cr.getX()[kMaxIdx];
            int y = cr.getY()[kMaxIdx];
            
            p1[count] = new PairInt(x, y);
            rCorners1.add(x, y);
            
            count++;
        }
        count = 0;
        for (CornerRegion cr : cornerRegions2) {
            int kMaxIdx = cr.getKMaxIdx();
            k2[count] = cr.getK()[kMaxIdx];
            int x = cr.getX()[kMaxIdx];
            int y = cr.getY()[kMaxIdx];
            
            p2[count] = new PairInt(x, y);
            rCorners2.add(x, y);
            
            count++;
        }

        // these are in trimmed reference frame
        
        MiscDebug.plotCorners(img1Grey, rCorners1, "1_region_corners_trimmed");
        MiscDebug.plotCorners(img2Grey, rCorners2, "2_region_corners_trimmed");
        
        QuickSort.sortBy1stArg(k1, p1);

        QuickSort.sortBy1stArg(k2, p2);

        int z = 1;
    }

    public GreyscaleImage getGreyscaleImage1() {
        return img1Grey;
    }
    public GreyscaleImage getGreyscaleImage2() {
        return img2Grey;
    }
    
    public ImageExt getImage1InOriginalReferenceFrame() {
        return img1Orig;
    }

    public GreyscaleImage getGXY1() {
        return gXY1;
    }

    public GreyscaleImage getTheta1() {
        return theta1;
    }

    public ImageExt getImage2InOriginalReferenceFrame() {
        return img2Orig;
    }

    public GreyscaleImage getGXY2() {
        return gXY2;
    }

    public GreyscaleImage getTheta2() {
        return theta2;
    }

}
