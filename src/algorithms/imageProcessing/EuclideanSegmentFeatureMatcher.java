package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 class whose goal is to find best single euclidean transformation for 
 image1 to image2.  It uses segmentation to create blobs.
 It solves for transformation of each blob in one image against the other.
 It keeps the best transformation solution of each blob in image1.
 Then evaluates all blob1 best transformations to find the best overall
 for the image corners.
 -- pros: for panorama and stereo images it does lead to a solution and it
          solves for rotation too.
 -- cons: has long runtime due to curve to curve transformation comparisons.
          also, it sometimes discards very good matching points
          because the curve needs at least 3 in order to
          estimate a euclidean transformation so any curves w/ fewer than 
          2 points are skipped.
 
 Prefer to use EuclideanSegmentFeatureMatcher2 which is faster, excepting
 for images like the checkerboard tests.  The checkerboard tests do not
 have unique features so a solution is better found by 
 EuclideanSegmentFeatureMatcher which tries many combinations from detailed
 to larger level groupings.
 
 * @author nichole
 */
public class EuclideanSegmentFeatureMatcher {

    private final FeatureMatcherSettings settings;

    // dither of 1 is fine for rotation of 0, but for larger rotations,
    // the segmentation using greyscale7 sometimes produces some differences
    // due to the order of two 1-D B3 spline operations (which are much faster
    // than a single 2-D B3 spline, so they are used even w/ the small variance
    // with rotation) so a dither=3 is necessary.
    public static final int dither = DitherDefault.dither;

    private final ImageExt img1;
    private final ImageExt img2;

    private GreyscaleImage gsImg1 = null;
    private GreyscaleImage gsImg2 = null;

    private Set<CornerRegion> cornerRegions1 = null;
    private Set<CornerRegion> cornerRegions2 = null;

    private enum State {
        DID_APPLY_HIST_EQ, COULD_NOT_DETERMINE_SCALE
    }
    private Set<State> stateSet = new HashSet<State>();

    private final boolean doDetermineScale;

    private TransformationParameters params = null;

    private float scaleTol = 0.2f;

    private float rotationInRadiansTol = (float)(20. * Math.PI/180.);

    //TODO: revise this...
    private int transXYTol = 20;

    private RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();

    private Logger log = Logger.getLogger(this.getClass().getName());

    public EuclideanSegmentFeatureMatcher(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {

        img1 = image1;
        img2 = image2;

        doDetermineScale = true;

        this.settings = settings.copy();
    }

    /**
     * constructor accepting transformation parameters.  Note, for best results,
     * the standard deviations within parameters should be populated because 
     * they are used as tolerances in matching.
     * @param image1
     * @param image2
     * @param parameters
     * @param settings
     */
    public EuclideanSegmentFeatureMatcher(ImageExt image1, ImageExt image2,
        TransformationParameters parameters, FeatureMatcherSettings settings) {

        img1 = image1;
        img2 = image2;
        doDetermineScale = false;

        this.settings = settings.copy();
    }

    public CorrespondenceList matchFeatures() throws IOException,
        NoSuchAlgorithmException {

        /*
        options:
            (1) determine scale
                (a) match remaining points derived in scale calc.
                    if resulting set spans the intersection,
                    make and return the correspondence list
                    else, follow (2)
            (2) given scale
                (b) extract corner regions from greyscale image
                (3) use feature matcher w/ scale to make the correspondence list
        */

        CorrespondenceList cl = null;

        if (doDetermineScale) {
            return solveTransformation();
        }

        applyHistEqIfNeeded();

        cl = extractAndMatch(params);

        if (settings.debug()) {
            printMatches(cl);
        }

        return cl;
    }

    private CorrespondenceList solveTransformation() throws IOException,
        NoSuchAlgorithmException {
        
        int dither2;
        if (settings.doUse2ndDerivCorners()) {
            dither2 = 1;
        } else {
            dither2 = dither;
        }

        BlobScaleFinderWrapper scaleFinder = new BlobScaleFinderWrapper(img1,
            img2, settings, rotatedOffsets, dither2);

        log.info("start calculatScale");

        params = scaleFinder.calculateScale();

        if (params == null) {
            stateSet.add(State.COULD_NOT_DETERMINE_SCALE);
            return null;
        }

        boolean didApplyHist = scaleFinder.img1Helper.didApplyHistEq();
        this.gsImg1 = scaleFinder.img1Helper.getGreyscaleImage().copyImage();
        this.gsImg2 = scaleFinder.img2Helper.getGreyscaleImage().copyImage();

        if (didApplyHist) {
            stateSet.add(State.DID_APPLY_HIST_EQ);
        }

        List<FeatureComparisonStat> stats =
            scaleFinder.getSolution().getComparisonStats();

        if (stats == null || stats.isEmpty()) {
            stateSet.add(State.COULD_NOT_DETERMINE_SCALE);
            return null;
        }

        CorrespondenceList cl = null;

        int tolXY;
        if (params.getStandardDeviations() != null) {
            tolXY = Math.round(Math.max(params.getStandardDeviations()[2],
                params.getStandardDeviations()[3]));
            if (tolXY < 3) {
                tolXY = 3;
            }
        } else {
            tolXY = 10;
        }

        int nLimit = 16;
        
        boolean covers = statsCoverIntersection(stats);
        
        boolean extractMoreCorners = (stats.size() < nLimit) || !covers;
        
        if (!extractMoreCorners) {

            List<PairInt> matched1 = new ArrayList<PairInt>();
            List<PairInt> matched2 = new ArrayList<PairInt>();
            populateLists(stats, matched1, matched2);

            cl = new CorrespondenceList(params.getScale(),
                Math.round(params.getRotationInDegrees()),
                Math.round(params.getTranslationX()),
                Math.round(params.getTranslationY()),
                Math.round(params.getStandardDeviations()[0]),
                Math.round(params.getStandardDeviations()[2]),
                Math.round(params.getStandardDeviations()[3]),
                matched1, matched2);

            return cl;
        }
            
        cl = extractAndMatch(params);
            
        if (cl != null) {
            
            addStatsToSolution(cl, stats);
            
            if (settings.debug()) {
                printMatches(cl.getPoints1(), cl.getPoints2());
                //MiscDebug.print(cl);
            }
        }
          
        return cl;
    }

    private void applyHistEqIfNeeded() {

        if (stateSet.contains(State.DID_APPLY_HIST_EQ)) {
            return;
        }

        if (gsImg1 != null) {
            // gs images were set during scale calculation
            return;
        }

        this.gsImg1 = img1.copyToGreyscale();
        this.gsImg2 = img2.copyToGreyscale();

        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(gsImg1, true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(gsImg2, true);

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
            HistogramEqualization hEq = new HistogramEqualization(gsImg1);
            hEq.applyFilter();
            hEq = new HistogramEqualization(gsImg2);
            hEq.applyFilter();
            stateSet.add(State.DID_APPLY_HIST_EQ);
        }

    }
    
    private void extractCornerRegions() {
        
        if (settings.doUse2ndDerivCorners()) {
            extract2ndDerivCornerRegions();
        } else {
            extractCannyCornerRegions();
        }
    }

    private void extractCannyCornerRegions() {

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(gsImg1, SIGMA.ONE);

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(gsImg1);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
//TODO: revisit to make sure coordinate systems are consistent:
        cornerRegions1 = detector.getEdgeCornerRegions(true);
        //cornerRegions1 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);

        if (settings.debug()) {
            List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
                Image imgCp = img1.copyImage();
                ImageIOHelper.addAlternatingColorCurvesToImage(edges, imgCp, 3);
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_1_edges_");
                imgCp = img1.copyImage();
                for (CornerRegion cr : cornerRegions1) {
                    int x = cr.getX()[cr.getKMaxIdx()];
                    int y = cr.getY()[cr.getKMaxIdx()];
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_1_cornerregions_");
                Map<Integer, Set<Integer>> junctionMap = detector.getJunctionMap();
                imgCp = img1.copyImage();
                for (Integer pixIndex : junctionMap.keySet()) {
                    int x = imgCp.getCol(pixIndex.intValue());
                    int y = imgCp.getRow(pixIndex.intValue());
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_1_junctions_");
                imgCp = img1.copyImage();
                PairIntArray corners = detector.getCornersInOriginalReferenceFrame();
                for (int ii = 0; ii < corners.getN(); ++ii) {
                    int x = corners.getX(ii);
                    int y = corners.getY(ii);
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_1_corners_");
        }

        //-------

        imageProcessor.blur(gsImg2, SIGMA.ONE);

        detector = new
            CurvatureScaleSpaceCornerDetector(gsImg2);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
        cornerRegions2 = detector.getEdgeCornerRegions(true);
        //cornerRegions2 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);

        if (settings.debug()) {
            List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
                Image imgCp = img2.copyImage();
                ImageIOHelper.addAlternatingColorCurvesToImage(edges, imgCp, 3);
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_2_edges_");
                imgCp = img2.copyImage();
                for (CornerRegion cr : cornerRegions2) {
                    int x = cr.getX()[cr.getKMaxIdx()];
                    int y = cr.getY()[cr.getKMaxIdx()];
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_2_corneregions_");
                Map<Integer, Set<Integer>> junctionMap = detector.getJunctionMap();
                imgCp = img2.copyImage();
                for (Integer pixIndex : junctionMap.keySet()) {
                    int x = imgCp.getCol(pixIndex.intValue());
                    int y = imgCp.getRow(pixIndex.intValue());
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_2_junctions_");
                imgCp = img2.copyImage();
                PairIntArray corners = detector.getCornersInOriginalReferenceFrame();
                for (int ii = 0; ii < corners.getN(); ++ii) {
                    int x = corners.getX(ii);
                    int y = corners.getY(ii);
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_2_corners_");
        }
    }
    
    private void extract2ndDerivCornerRegions() {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> points1 = imageProcessor.extract2ndDerivPoints(
            img1.copyImage().copyToGreyscale(), 200, true);
        
        Set<PairInt> points2 = imageProcessor.extract2ndDerivPoints(
            img2.copyImage().copyToGreyscale(), 200, true);

        cornerRegions1 = new HashSet<CornerRegion>();
        for (PairInt p : points1) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
            cr.setIndexWithinCurve(-1);
            cornerRegions1.add(cr);
        }
        
        cornerRegions2 = new HashSet<CornerRegion>();
        for (PairInt p : points2) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
            cr.setIndexWithinCurve(-1);
            cornerRegions2.add(cr);
        }
        
        if (settings.debug()) {
            Image imgCp = img1.copyImage();
            for (CornerRegion cr : cornerRegions1) {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
            }
            MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_1_cornerregions_");
            imgCp = img2.copyImage();
            for (CornerRegion cr : cornerRegions2) {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
            }
            MiscDebug.writeImage(imgCp, settings.getDebugTag() + "_2_cornerregions_");
        }
    }

    private CorrespondenceList findCorrespondence(TransformationParameters
        parameters) {

        FeatureMatcher featureMatcher = new FeatureMatcher();

        int tolXY;
        if (params.getStandardDeviations() != null) {
            tolXY = Math.round(Math.max(params.getStandardDeviations()[2],
                params.getStandardDeviations()[3]));
            if (tolXY < 3) {
                tolXY = 3;
            }
        } else {
            tolXY = transXYTol;
        }

        CorrespondenceList cl = featureMatcher.findSimilarFeatures(gsImg1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]),
            parameters, scaleTol, rotationInRadiansTol, tolXY,
            dither, rotatedOffsets);

        return cl;
    }

    private void populateLists(List<FeatureComparisonStat> stats,
        List<PairInt> matched1, List<PairInt> matched2) {

        for (FeatureComparisonStat stat : stats) {

            int x1 = stat.getImg1Point().getX();
            int y1 = stat.getImg1Point().getY();
            int x2 = stat.getImg2Point().getX();
            int y2 = stat.getImg2Point().getY();

            matched1.add(new PairInt(x1, y1));

            matched2.add(new PairInt(x2, y2));
        }
    }

    private CorrespondenceList extractAndMatch(
        TransformationParameters parameters) {

        extractCornerRegions();

        CorrespondenceList cl = findCorrespondence(parameters);

        return cl;
    }

    private void printMatches(CorrespondenceList cl) {
        if (cl == null) {
            return;
        }
        printMatches(cl.getPoints1(), cl.getPoints2());
    }

    private void printMatches(Collection<PairInt> m1, Collection<PairInt> m2) {

        log.info("writing " + m1.size() + " matches to images");

        int ts = MiscDebug.getCurrentTimeFormatted();
        GreyscaleImage gsImg1 = img1.copyToGreyscale();
        GreyscaleImage gsImg2 = img2.copyToGreyscale();
        String name1 = "1_" + settings.getDebugTag() + "_" + ts;
        String name2 = "2_" + settings.getDebugTag() + "_" + ts;
        name1 = name1 + "_matched";
        name2 = name2 + "_matched";
        MiscDebug.plotCorners(gsImg1, m1, name1, 2);
        MiscDebug.plotCorners(gsImg2, m2, name2, 2);
    }

    /**
     * a method to determine the intersection of transformed image 1 with
     * image 2 and then examine the distribution of stats's points in 
     * 4 quadrants of the intersection to return whether stats are present in
     * all quadrants.  A caveat of the method is that not all of the 
     * intersection necessarily has image details which could be matched, for 
     * example, clear sky does not have corners using the methods here.
     * @param stats
     * @return 
     */
    private boolean statsCoverIntersection(List<FeatureComparisonStat> stats) {
        
        /*
        calculate the intersection of the 2 images.
        divide the region into 4 parts (2 vertical and 2 horizontal) by noting
        the 4 boundary points for each and making a polygon for each.
        
        then use point in polygon tests to count the number of stats.point2's
        in each of the 4 regions.        
        */
        
        /*
       / \   ( tr )    ( tr )            (x2q2, y2q2)  d5   (x2q3, y2q3)
        |
        |                                 d2           d3             d4
        |
        0    ( tr )    ( tr )            (x2q1, y2q1)  d1   (x2q0, y2q0)
          0 -->
        */
        
        double[][] img2Intersection = MiscStats.getBoundsOfIntersectionInFrame2(
            params, img1.getWidth(), img1.getHeight(), 
            img2.getWidth(), img2.getHeight());
        
        float[] d1 = new float[]{
            (float)((img2Intersection[0][0] + img2Intersection[1][0])/2.f),
            (float)((img2Intersection[0][1] + img2Intersection[1][1])/2.f)};     
        float[] d2 = new float[]{
            (float)((img2Intersection[1][0] + img2Intersection[2][0])/2.f),
            (float)((img2Intersection[1][1] + img2Intersection[2][1])/2.f)};
        float[] d4 = new float[]{
            (float)((img2Intersection[0][0] + img2Intersection[3][0])/2.f),
            (float)((img2Intersection[0][1] + img2Intersection[3][1])/2.f)};
        float[] d5 = new float[]{
            (float)((img2Intersection[2][0] + img2Intersection[3][0])/2.f),
            (float)((img2Intersection[2][1] + img2Intersection[3][1])/2.f)};
        float[] d3 = new float[]{(d2[0] + d4[0])/2.f, (d1[1] + d5[1])/2.f};
        
        float[] xPoly0 = new float[5];
        float[] yPoly0 = new float[5];
        xPoly0[0] = (float)img2Intersection[0][0];
        yPoly0[0] = (float)img2Intersection[0][1];
        xPoly0[1] = d1[0];
        yPoly0[1] = d1[1];
        xPoly0[2] = d3[0];
        yPoly0[2] = d3[1];
        xPoly0[3] = d4[0];
        yPoly0[3] = d4[1];
        xPoly0[4] = xPoly0[0];
        yPoly0[4] = yPoly0[0];

        /*
       / \   ( tr )    ( tr )            (x2q2, y2q2)  d5   (x2q3, y2q3)
        |
        |                                 d2           d3             d4
        |
        0    ( tr )    ( tr )            (x2q1, y2q1)  d1   (x2q0, y2q0)
          0 -->
        */
        
        float[] xPoly1 = new float[5];
        float[] yPoly1 = new float[5];
        xPoly1[0] = d1[0];
        yPoly1[0] = d1[1];
        xPoly1[1] = (float)img2Intersection[1][0];
        yPoly1[1] = (float)img2Intersection[1][1];
        xPoly1[2] = d2[0];
        yPoly1[2] = d2[1];
        xPoly1[3] = d3[0];
        yPoly1[3] = d3[1];
        xPoly1[4] = xPoly1[0];
        yPoly1[4] = yPoly1[0];
        
        float[] xPoly2 = new float[5];
        float[] yPoly2 = new float[5];
        xPoly2[0] = d3[0];
        yPoly2[0] = d3[1];
        xPoly2[1] = d2[0];
        yPoly2[1] = d2[1];
        xPoly2[2] = (float)img2Intersection[2][0];
        yPoly2[2] = (float)img2Intersection[2][1];
        xPoly2[3] = d5[0];
        yPoly2[3] = d5[1];
        xPoly2[4] = xPoly2[0];
        yPoly2[4] = yPoly2[0];
        
        /*
       / \   ( tr )    ( tr )            (x2q2, y2q2)  d5   (x2q3, y2q3)
        |
        |                                 d2           d3             d4
        |
        0    ( tr )    ( tr )            (x2q1, y2q1)  d1   (x2q0, y2q0)
          0 -->
        */
        
        float[] xPoly3 = new float[5];
        float[] yPoly3 = new float[5];
        xPoly3[0] = d4[0];
        yPoly3[0] = d4[1];
        xPoly3[1] = d3[0];
        yPoly3[1] = d3[1];
        xPoly3[2] = d5[0];
        yPoly3[2] = d5[1];
        xPoly3[3] = (float)img2Intersection[3][0];
        yPoly3[3] = (float)img2Intersection[3][1];
        xPoly3[4] = xPoly3[0];
        yPoly3[4] = yPoly3[0];
        
        PointInPolygon poly = new PointInPolygon();
        
        int[] count = new int[4];
        for (FeatureComparisonStat stat : stats) {
            int x = stat.getImg2Point().getX() * stat.getBinFactor2();
            int y = stat.getImg2Point().getY() * stat.getBinFactor2();
            boolean isIn = poly.isInSimpleCurve(x, y, xPoly0, yPoly0, 5);
            if (isIn) {
                count[0]++;
            } else {
                isIn = poly.isInSimpleCurve(x, y, xPoly1, yPoly1, 5);
                if (isIn) {
                    count[1]++;
                } else {
                    isIn = poly.isInSimpleCurve(x, y, xPoly2, yPoly2, 5);
                    if (isIn) {
                        count[2]++;
                    } else {
                        isIn = poly.isInSimpleCurve(x, y, xPoly3, yPoly3, 5);
                        if (isIn) {
                            count[3]++;
                        }
                    }
                }
            }
        }
        
        int nq = 0;
        for (int c : count) {
            if (c > 0) {
                nq++;
            }
        }
        
        if (settings.debug()) {
            Image imcp = img2.copyImage();
            for (int i = 0; i < xPoly0.length; ++i) {
                ImageIOHelper.addPointToImage(xPoly0[i], yPoly0[i], imcp, 5, 0, 255, 255);
            }
            for (int i = 0; i < xPoly1.length; ++i) {
                ImageIOHelper.addPointToImage(xPoly1[i], yPoly1[i], imcp, 5, 0, 255, 0);
            }
            for (int i = 0; i < xPoly2.length; ++i) {
                ImageIOHelper.addPointToImage(xPoly2[i], yPoly2[i], imcp, 5, 0, 0, 255);
            }
            for (int i = 0; i < xPoly3.length; ++i) {
                ImageIOHelper.addPointToImage(xPoly3[i], yPoly3[i], imcp, 5, 0, 125, 125);
            }
            for (int i = 0; i < stats.size(); ++i) {
                FeatureComparisonStat stat = stats.get(i);
                PairInt p2 = stat.getImg2Point();
                ImageIOHelper.addPointToImage(p2.getX() * stat.getBinFactor2(), 
                    p2.getY() * stat.getBinFactor2(), imcp, 2, 255, 0, 0);
            }
            MiscDebug.writeImage(imcp, settings.getDebugTag() + "_scale_points");
        }
        
        return (nq == 4);
    }
    
    private void addStatsToSolution(CorrespondenceList cl, 
        List<FeatureComparisonStat> stats) {
        
        if (cl == null) {
            return;
        }
        
        List<PairInt> matched01 = cl.getPoints1();
        List<PairInt> matched02 = cl.getPoints2();
        
        Set<PairInt> added1 = new HashSet<PairInt>(matched01);
        Set<PairInt> added2 = new HashSet<PairInt>(matched02);
        
        boolean didAdd = false;
        
        for (FeatureComparisonStat stat : stats) {
            PairInt imPt1 = stat.getImg1Point();
            PairInt imPt2 = stat.getImg2Point();
            if (added1.contains(imPt1) || added2.contains(imPt2)) {
                //TODO: consider replacing the match with imPt1, imPt2 which is
                // usually better from the small first solution.
                // in that case, need to remove the added.adds below
                continue;
            }
            matched01.add(imPt1);
            matched02.add(imPt2);
            added1.add(imPt1);
            added2.add(imPt2);
            
            didAdd = true;
        }
        
        if (didAdd) {
            // TODO: redo ranges
        }
    }
}
