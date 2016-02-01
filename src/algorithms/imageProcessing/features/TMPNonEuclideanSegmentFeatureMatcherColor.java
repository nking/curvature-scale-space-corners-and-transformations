package algorithms.imageProcessing.features;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * create lists of singly matched points between 2 images.
 * It uses the criteria that matches are discarded if a point has a second
 * best match whose SSD is within 0.8*SSD of best match.
 *
 * @author nichole
 */
public class TMPNonEuclideanSegmentFeatureMatcherColor {

    private Logger log = Logger.getLogger(this.getClass().getName());

    protected final int binnedImageMaxDimension = 512;

    private final ImageExt img1;
    private final ImageExt img2;
    private final FeatureMatcherSettings settings;
    private final RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
    private final ImageExt imgBinned1;
    private final ImageExt imgBinned2;
    private final GreyscaleImage redBinnedImg1;
    private final GreyscaleImage greenBinnedImg1;
    private final GreyscaleImage blueBinnedImg1;
    private final GreyscaleImage redBinnedImg2;
    private final GreyscaleImage greenBinnedImg2;
    private final GreyscaleImage blueBinnedImg2;
    private final int binFactor1, binFactor2;

    protected final IntensityClrFeatures clrFeaturesBinned1;
    protected final IntensityClrFeatures clrFeaturesBinned2;

    private List<List<FeatureComparisonStat>> solutionStats = null;

    private List<List<PairInt>> solutionMatched1 = null;
    private List<List<PairInt>> solutionMatched2 = null;

    /**
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     * @param settings
     */
    public TMPNonEuclideanSegmentFeatureMatcherColor(ImageExt img1, ImageExt img2,
        FeatureMatcherSettings settings) {

        this.img1 = img1;
        this.img2 = img2;
        this.settings = settings;

        this.binFactor1 = (int) Math.ceil(
            Math.max((float) img1.getWidth() / (float)binnedImageMaxDimension,
            (float) img1.getHeight() / (float)binnedImageMaxDimension));

        this.binFactor2 = (int) Math.ceil(
            Math.max((float) img2.getWidth() / (float)binnedImageMaxDimension,
            (float) img2.getHeight() / (float)binnedImageMaxDimension));

        ImageProcessor imageProcessor = new ImageProcessor();

        imgBinned1 = imageProcessor.binImage(img1, binFactor1);
        imgBinned2 = imageProcessor.binImage(img2, binFactor2);

        /*
        HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(imgBinned1);
        hEq.applyFilter();
        hEq = new HistogramEqualizationForColor(imgBinned2);
        hEq.applyFilter();
        */

        redBinnedImg1 = imgBinned1.copyRedToGreyscale();
        greenBinnedImg1 = imgBinned1.copyGreenToGreyscale();
        blueBinnedImg1 = imgBinned1.copyBlueToGreyscale();

        redBinnedImg2 = imgBinned2.copyRedToGreyscale();
        greenBinnedImg2 = imgBinned2.copyGreenToGreyscale();
        blueBinnedImg2 = imgBinned2.copyBlueToGreyscale();

        GreyscaleImage gsImg1 = imgBinned1.copyToGreyscale();
        GreyscaleImage gsImg2 = imgBinned2.copyToGreyscale();

        clrFeaturesBinned1 = new IntensityClrFeatures(gsImg1, 5, rotatedOffsets);

        clrFeaturesBinned2 = new IntensityClrFeatures(gsImg2, 5, rotatedOffsets);
    }

    public boolean match() throws IOException, NoSuchAlgorithmException {

        /*boolean matchPtByPt = false;
        if (matchPtByPt) {
            return matchPtByPt();
        } else {*/
            return matchBlobByBlob();
        //}
    }

    /*
    private boolean matchPtByPt() throws IOException, NoSuchAlgorithmException {

        Set<PairInt> points1 = extractKeyPoints(redBinnedImg1,
            greenBinnedImg1, blueBinnedImg1, clrFeaturesBinned1, "1");

        Set<PairInt> points2 = extractKeyPoints(redBinnedImg2,
            greenBinnedImg2, blueBinnedImg2,  clrFeaturesBinned2, "2");
        
        List<PairInt> pointsList1 = new ArrayList<PairInt>(points1);
        List<PairInt> pointsList2 = new ArrayList<PairInt>(points2);

        int dither = 0;

        CornerMatcher matcher = new CornerMatcher(dither);

        //TODO: consider using the ransac that includes features for evaluation
                
        boolean matched = matcher.matchPoints(
            clrFeaturesBinned1, clrFeaturesBinned2, pointsList1, pointsList2,
            imgBinned1, redBinnedImg1, greenBinnedImg1, blueBinnedImg1,
            imgBinned2, redBinnedImg2, greenBinnedImg2, blueBinnedImg2,
            binFactor1, binFactor2);

        if (!matched) {
            return false;
        }

        List<FeatureComparisonStat> stats = matcher.getSolutionStats();

        log.info("nPts after SSD match (incl filter for 2nd best) = " + stats.size());

        log.info("nPts in 2nd best rejection list = " + matcher.getRejectedBy2ndBest().size());

        if (stats.isEmpty()) {
            return false;
        }

        if (settings.debug()) {
            long ts = MiscDebug.getCurrentTimeFormatted();
            MiscDebug.plotImages(stats, blueBinnedImg1.copyImage(),
                blueBinnedImg2.copyImage(), 2, "SSD_matched_" + settings.getDebugTag() + "_" + ts);
        }

        List<FeatureComparisonStat> rb2j =
            new ArrayList<FeatureComparisonStat>(matcher.getRejectedBy2ndBest());

        stats = reviseStatsForFullImages(stats, binFactor1, binFactor2, rotatedOffsets);

        rb2j = reviseStatsForFullImages(rb2j, binFactor1, binFactor2, rotatedOffsets);

        copyToInstanceVars(stats, rb2j);

        return true;
    }
    */

    private boolean matchBlobByBlob() throws IOException, NoSuchAlgorithmException {

        KeyPointsAndBounds kpab1 =
            extractKeypointsAndAssocWithBlobs(redBinnedImg1,
            greenBinnedImg1, blueBinnedImg1, imgBinned1.copyToImageExt(), "1", 
            clrFeaturesBinned1);

        KeyPointsAndBounds kpab2 =
            extractKeypointsAndAssocWithBlobs(redBinnedImg2,
            greenBinnedImg2, blueBinnedImg2, imgBinned2.copyToImageExt(), "2", 
            clrFeaturesBinned2);

        int dither = 0;
        
        /*
        TODO: need to use the homology within a point group set while matching
        to another point group set in order to get better solution.
        
        -- will compare each group in img1 to each group in img2
           -- SSD match to make a matched points list.
              -- plot and log this
           -- ransac w/ epipolar and feaures evaluation to find inliers
              -- plot and log this.
                 are the inliers true matches?
           -- possibly need another round or two of matching the points 
              that are not inliers (now that the pool of matches is smaller)
              -- then give those and the inliers back to the ransac algorithm
                 to find inliers.
                 -- is the solution improved?
        */
        
        CornerMatcher matcher = new CornerMatcher(dither);
        
        boolean useBipartiteMatching = true;
        
        List<List<FeatureComparisonStat>> fcsList = matcher.matchCornersByBlobsAndRANSAC(
            clrFeaturesBinned1, clrFeaturesBinned2, kpab1, kpab2, 
            imgBinned1, redBinnedImg1, greenBinnedImg1, blueBinnedImg1, 
            imgBinned2, redBinnedImg2, greenBinnedImg2, blueBinnedImg2, 
            binFactor1, binFactor2, useBipartiteMatching);
        
        if (fcsList == null || fcsList.isEmpty()) {
            return false;
        }
        
        fcsList = reviseStatsForFullImages(fcsList, binFactor1, binFactor2, rotatedOffsets);

        copyToInstanceVars(fcsList);

        return true;
    }

    public void copyToInstanceVars(List<List<FeatureComparisonStat>> stats) {

        this.solutionStats = new ArrayList<List<FeatureComparisonStat>>();
        this.solutionMatched1 = new ArrayList<List<PairInt>>();
        this.solutionMatched2 = new ArrayList<List<PairInt>>();

        for (List<FeatureComparisonStat> list : stats) {
        
            List<FeatureComparisonStat> list2 = new ArrayList<FeatureComparisonStat>();
            List<PairInt> m1 = new ArrayList<PairInt>();
            List<PairInt> m2 = new ArrayList<PairInt>();
            for (FeatureComparisonStat stat : list) {
                list2.add(stat.copy());
                m1.add(stat.getImg1Point().copy());
                m2.add(stat.getImg2Point().copy());
            }
            
            this.solutionStats.add(list2);
            this.solutionMatched1.add(m1);
            this.solutionMatched2.add(m2);
        }
    }

    private Set<PairInt> extract2ndDerivPoints(GreyscaleImage rImg,
        GreyscaleImage gImg, GreyscaleImage bImg) {

        int nApprox = 5000; //200

        ImageProcessor imageProcessor = new ImageProcessor();

        Set<PairInt> pR = imageProcessor.extract2ndDerivPoints(rImg, nApprox, true);
        Set<PairInt> pG = imageProcessor.extract2ndDerivPoints(gImg, nApprox, true);
        Set<PairInt> pB = imageProcessor.extract2ndDerivPoints(bImg, nApprox, true);
        Set<PairInt> pixels = new HashSet<PairInt>();
        pixels.addAll(pR);
        pixels.addAll(pG);
        pixels.addAll(pB);

        return pixels;
    }

    private Set<PairInt> extractKeyPoints(GreyscaleImage rImg,
        GreyscaleImage gImg, GreyscaleImage bImg, IntensityClrFeatures features,
        String lbl) {

        Set<PairInt> pixels = extract2ndDerivPoints(rImg, gImg, bImg);

        long ts = MiscDebug.getCurrentTimeFormatted();
        if (settings.debug()) {
            try {
                MiscDebug.writeImage(pixels, bImg.copyToColorGreyscale(),
                    "all_2ndderiv_" + settings.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
        
        // this does not remove very many points because the 2nd deriv is already
        // finding points with steep and changing intensity slopes
        boolean filterForLocalization = false;
        if (filterForLocalization) {
            log.info("before filter for localizability nPts=" + pixels.size());
            filterForLocalization0(rImg, gImg, bImg, features, pixels);
            log.info("after filter for localizability nPts=" + pixels.size());
            if (settings.debug()) {
                try {
                    MiscDebug.writeImage(pixels, rImg.copyToColorGreyscale(),
                       "all_2ndderiv_filtered_local_" + lbl + "_" + settings.getDebugTag() + "_" + ts);
                } catch (IOException ex) {
                    Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                        log(Level.SEVERE, null, ex);
                }
            }
        }

        boolean useWaveletFilter = false;
        if (useWaveletFilter) {
            log.info("before segmentation filters nPts=" + pixels.size());
            pixels = filterPointsBySegmentation0(rImg, gImg, bImg, pixels, lbl);
            log.info("after wavelet segmentation filter nPts=" + pixels.size());

            if (settings.debug()) {
                try {
                    MiscDebug.writeImage(pixels, rImg.copyToColorGreyscale(),
                        "filtered_wavelet_" + lbl + "_" + 
                        settings.getDebugTag() + "_" + ts);
                } catch (IOException ex) {
                    Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                        log(Level.SEVERE, null, ex);
                }
            }
        }

        return pixels;
    }

    private KeyPointsAndBounds extractKeypointsAndAssocWithBlobs(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg, 
        ImageExt img, String lbl, IntensityClrFeatures features) {

        Set<PairInt> pixels = extract2ndDerivPoints(rImg, gImg, bImg);

        long ts = MiscDebug.getCurrentTimeFormatted();
        if (settings.debug()) {
            try {
                MiscDebug.writeImage(pixels, bImg.copyToColorGreyscale(),
                    "all_2ndderiv_" + lbl + "_" + settings.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
        
        boolean filterForLocalization = true;
        if (filterForLocalization) {
            filterForLocalization0(rImg, gImg, bImg, features, pixels);
            if (settings.debug()) {
                try {
                    MiscDebug.writeImage(pixels, bImg.copyToColorGreyscale(),
                        "all_2ndderiv__filtered_local_" + lbl + "_" + settings.getDebugTag() + "_" + ts);
                } catch (IOException ex) {
                    Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                        log(Level.SEVERE, null, ex);
                }
            }
            log.info("number of points after localizability filter=" + 
                pixels.size());
        }
        
        // join those grouped points into larger groups where possible
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        BoundingRegions largerBoundingRegion = imageSegmentation
            .extractBlobBoundsFromLowRes(img, settings.debug(), 
            settings.getDebugTag(), features);
                
        if (settings.debug()) {
            ImageExt imgCp = (ImageExt) img.copyImage();
            ImageIOHelper.addAlternatingColorCurvesToImage(
                largerBoundingRegion.getPerimeterList(), 
                imgCp, 0, 0, 2);
            MiscDebug.writeImage(imgCp, "_larger_bounds_" + lbl + "_" + 
            settings.getDebugTag() + "_" + ts);
        }
        
        List<Set<PairInt>> pixelSetList = filterToBounds(pixels, largerBoundingRegion);
        
        KeyPointsAndBounds kpab = new KeyPointsAndBounds(pixelSetList, 
            largerBoundingRegion);
                
        if (settings.debug()) {
            ImageExt imgCp = (ImageExt) img.copyImage();
            try {
                ImageIOHelper.addAlternatingColorPointSetsToImage(pixelSetList,
                    0, 0, 2, imgCp);
                MiscDebug.writeImage(imgCp, "_points_in_larger_bounds_" + lbl + "_" + 
                settings.getDebugTag() + "_" + ts);
                
            } catch (IOException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            }
        }

        // filter by atrous based sementation only for highly textured regions (high density and large area)
        filterTexturedGroupsByWaveletSegmentation(kpab, rImg, gImg, bImg, lbl);
        
        if (settings.debug()) {
            ImageExt imgCp = (ImageExt) img.copyImage();
            try {
                ImageIOHelper.addAlternatingColorPointSetsToImage(pixelSetList,
                    0, 0, 2, imgCp);
                MiscDebug.writeImage(imgCp, "_points_filtered_textures_" + lbl + "_" + 
                settings.getDebugTag() + "_" + ts);
                
            } catch (IOException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        return kpab;
    }

    protected List<List<FeatureComparisonStat>> reviseStatsForFullImages(
        List<List<FeatureComparisonStat>> stats, int prevBinFactor1,
        int prevBinFactor2, RotatedOffsets rotatedOffsets) {

        log.info("refine stats for full image reference frames");

        if (stats.isEmpty()) {
            return stats;
        }

        //TODO: when have a work clr descriptor, put the full size feature matching backing in

        List<List<FeatureComparisonStat>> revised = new ArrayList<List<FeatureComparisonStat>>();

        for (List<FeatureComparisonStat> list : stats) {
            
            List list2 = new ArrayList<FeatureComparisonStat>();
            
            for (int i = 0; i < list.size(); ++i) {

                FeatureComparisonStat stat = list.get(i).copy();

                PairInt p1 = stat.getImg1Point();
                p1.setX(p1.getX() * prevBinFactor1);
                p1.setY(p1.getY() * prevBinFactor1);

                PairInt p2 = stat.getImg2Point();
                p2.setX(p2.getX() * prevBinFactor2);
                p2.setY(p2.getY() * prevBinFactor2);

                list2.add(stat);
            }
            
            revised.add(list2);
        }

        return revised;
    }

    protected void filterForLocalization(GreyscaleImage rImg,
        GreyscaleImage gImg, GreyscaleImage bImg,
        IntensityClrFeatures f, List<CornerRegion> corners) {

        List<Integer> remove = new ArrayList<Integer>();

        for (int i = 0; i < corners.size(); ++i) {
            CornerRegion cr = corners.get(i);

            try {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                if (f.removeDueToLocalization(rImg, gImg, bImg, x, y,
                    f.calculateOrientation(x, y))) {
                    remove.add(Integer.valueOf(i));
                }
            } catch (CornerRegion.CornerRegionDegneracyException ex) {
            }
        }

        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i);
            corners.remove(idx);
        }
    }

    protected void filterForLocalization0(GreyscaleImage rImg,
        GreyscaleImage gImg, GreyscaleImage bImg,
        IntensityClrFeatures f, Set<PairInt> keyPoints) {

        Set<PairInt> remove = new HashSet<PairInt>();

        for (PairInt kp : keyPoints) {
            
            try {
                int x = kp.getX();
                int y = kp.getY();
                if (f.removeDueToLocalization(rImg, gImg, bImg, x, y,
                    f.calculateOrientation(x, y))) {
                    remove.add(kp);
                }
            } catch (CornerRegion.CornerRegionDegneracyException ex) {
            }
        }

        keyPoints.removeAll(remove);
    }

    private Set<PairInt> filterPointsBySegmentation1(ImageExt img,
        Set<PairInt> pixels, String lbl) {

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        GreyscaleImage imgSeg1 = imageSegmentation.applyUsingCIEXYPolarThetaThenHistEq(img);
        MiscDebug.writeImage(imgSeg1, "ciexy_" +lbl + "_");

        int w = imgSeg1.getWidth();
        int h = imgSeg1.getHeight();

        List<PairInt> pixels1 = new ArrayList<PairInt>();
        List<Integer> diffs = new ArrayList<Integer>();

        int minDiffV = Integer.MAX_VALUE;

        for (PairInt p : pixels) {
            int x = p.getX();
            int y = p.getY();
            int v = -1;
            boolean skip = false;
            for (int dx = -1; dx <= 1; ++dx) {
                int x1 = x + dx;
                if ((x1 < 0) || (x1 > (w - 1))) {
                    skip = true;
                    break;
                }
                for (int dy = -1; dy <= 1; ++dy) {
                    int y1 = y + dy;
                    if (x1 == x && y1 == y) {
                        continue;
                    }
                    if ((y1 < 0) || (y1 > (h - 1))) {
                        skip = true;
                        break;
                    }
                    int v1 = imgSeg1.getValue(x1, y1);
                    if (v == -1) {
                        v = v1;
                    } else if (v1 != v) {
                        int diffV = Math.abs(v1 - v);
                        pixels1.add(new PairInt(x, y));
                        diffs.add(Integer.valueOf(diffV));
                        if (diffV < minDiffV) {
                            minDiffV = diffV;
                        }
                        skip = true;
                        break;
                    }
                }
                if (skip) {
                    break;
                }
            }
        }

        Set<PairInt> pixels2 = new HashSet<PairInt>();
        for (int i = 0; i < pixels1.size(); ++i) {
            int diffV = diffs.get(i).intValue();
            if (diffV > minDiffV) {
                pixels2.add(pixels1.get(i));
            }
        }

        return pixels2;
    }

    private List<CornerRegion> convert(Set<PairInt> points) {

        List<CornerRegion> corners = new ArrayList<CornerRegion>();
        for (PairInt p : points) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
            cr.setIndexWithinCurve(-1);
            corners.add(cr);
        }
        return corners;
    }

    private Set<PairInt> filterPointsBySegmentation0(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg,
        Set<PairInt> pixels, String lbl) {

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        GreyscaleImage segImg = imageSegmentation.createCombinedWaveletBased(rImg, gImg, bImg);

        if (settings.debug()) {
            long ts = MiscDebug.getCurrentTimeFormatted();
            MiscDebug.writeImage(segImg,
                "segmented_" +lbl + "_" + settings.getDebugTag() + "_" + ts);
        }

        Set<PairInt> pixels2 = new HashSet<PairInt>();

        for (PairInt p : pixels) {
            int v = segImg.getValue(p.getX(), p.getY());
            if (v > 0) {
                pixels2.add(p);
            }
        }

        return pixels2;
    }

    private void filterListsForLocalization(GreyscaleImage rImg,
        GreyscaleImage gImg, GreyscaleImage bImg,
        IntensityClrFeatures clrFeatures,
        List<Set<PairInt>> keyPointLists) {

        for (int i = 0; i < keyPointLists.size(); ++i) {
            Set<PairInt> points = keyPointLists.get(i);
            filterForLocalization0(rImg, gImg, bImg, clrFeatures, points);
        }
    }

    /**
     * @return the solutionStats
     */
    public List<List<FeatureComparisonStat>> getSolutionStats() {
        return solutionStats;
    }

    /**
     * @return the solutionMatched1
     */
    public List<List<PairInt>> getSolutionMatched1() {
        return solutionMatched1;
    }

    /**
     * @return the solutionMatched2
     */
    public List<List<PairInt>> getSolutionMatched2() {
        return solutionMatched2;
    }

    /**
     * 
     * @param pixelSetList list of groups of points
     * @param largerBoundingRegion a list of larger bounding regions that 
     * contains list of polynomials defining larger bounds than the
     * points within a single item of pixelSetList and medial axes for those
     * regions.  The internal contents of this structure are reordered and pruned
     * to match the alterations in pixelSetList so that the structure can
     * afterwards be used to provide complementary information for pixelSetList.
     */
    private List<Set<PairInt>> filterToBounds(Set<PairInt> pixels, 
        BoundingRegions largerBoundingRegions) {
        
        List<PairIntArray> largerBounds = largerBoundingRegions.getPerimeterList();
        
        BlobMedialAxes bma = largerBoundingRegions.getBlobMedialAxes();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
               
        // format largerBounds for PointInPolygon
        int[][] xPoly = new int[largerBounds.size()][];
        int[][] yPoly = new int[largerBounds.size()][];
        double[][] xyBoundsMinMaxXY = new double[largerBounds.size()][2];
        for (int i = 0; i < largerBounds.size(); ++i) {
            PairIntArray polynomial = largerBounds.get(i);
            xPoly[i] = new int[polynomial.getN()];
            yPoly[i] = new int[polynomial.getN()];
            for (int j = 0; j < polynomial.getN(); ++j) {
                xPoly[i][j] = polynomial.getX(j);
                yPoly[i][j] = polynomial.getY(j);
            }
            double xMin = MiscMath.findMin(xPoly[i]);
            double xMax = MiscMath.findMax(xPoly[i]);
            double yMin = MiscMath.findMin(yPoly[i]);
            double yMax = MiscMath.findMax(yPoly[i]);
            xyBoundsMinMaxXY[i] = new double[]{xMin, xMax, yMin, yMax};
        }
        
        PointInPolygon pInPoly = new PointInPolygon();
        
        // points will be removed as added to output
        Set<PairInt> points2 = new HashSet<PairInt>(pixels);
        
        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();
        
        // worst case is O(N_largerBounds * N_pixels) if no pixels are placed in output
        for (int i = 0; i < largerBounds.size(); ++i) {
                        
            Set<PairInt> add = new HashSet<PairInt>();
            
            for (PairInt p : points2) {
                if ((p.getX() < xyBoundsMinMaxXY[i][0])
                    || (p.getX() > xyBoundsMinMaxXY[i][1])
                    || (p.getY() < xyBoundsMinMaxXY[i][2])
                    || (p.getY() > xyBoundsMinMaxXY[i][3])) {
                    continue;
                }
                if (pInPoly.isInSimpleCurve(p.getX(), p.getY(),
                    xPoly[i], yPoly[i], xPoly[i].length)) {
                    add.add(p);
                }
            }
            
            output.add(add);
            
            points2.removeAll(add);
        }
        
        assert(output.size() == largerBounds.size());
        
        // anything remaining in points2 was not found within a largerBounds
        
        //TODO: could consider ways to retain those in points2 if those are
        // found to be points that are needed for better matches.
                
        List<Integer> remove = new ArrayList<Integer>();
        for (int i = 0; i < output.size(); ++i) {
            Set<PairInt> set = output.get(i);
            if (set.isEmpty()) {
                remove.add(Integer.valueOf(i));
            }
        }
        for (int i = (remove.size() - 1); i > -1; --i) {
            int rmIdx = remove.get(i);
            output.remove(rmIdx);
        }
                
        largerBoundingRegions.removeIndexes(remove);
        
        assert(output.size() == largerBounds.size());
        
        return output;        
    }

    private void filterTexturedGroupsByWaveletSegmentation(
        KeyPointsAndBounds keypointsAndBounds, 
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg,
        String lbl) {
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<Integer> highDensityLargeGroups = new ArrayList<Integer>();
    
        List<Set<PairInt>> keypoints = keypointsAndBounds.getKeyPointGroups();
            
        for (int i = 0; i < keypoints.size(); ++i) {
            
            Set<PairInt> group = keypoints.get(i);
            
            double n = group.size();
            
            double area = curveHelper.calculateArea(
                keypointsAndBounds.getBoundingRegions().getPerimeterList().get(i));
            
            double areaFraction = Math.abs(area/(double)(rImg.getWidth() * rImg.getHeight()));
            
            double[] xyCen = keypointsAndBounds.getBoundingRegions()
                .getBlobMedialAxes().getOriginalBlobXYCentroid(i);
            
            double numberDensity = Math.abs(n/area);
            
            log.info(String.format(
               "[%d] (%d,%d)  n=%d  area=%.3f fraction of image=%.3f  density=%.3f",
                i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]),
                (int)n, (float)area, (float)areaFraction, (float)numberDensity));
            
            if ((areaFraction >= 0.04) && (numberDensity >= 0.09)) {
                highDensityLargeGroups.add(Integer.valueOf(i));
            }
            
        }
        
        if (highDensityLargeGroups.isEmpty()) {
            return;
        }
        
        // for items in highDensityLargeGroups, will reduce them to only those
        // co-spatial with a segmentation point
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        GreyscaleImage segImg = imageSegmentation.createCombinedWaveletBased(rImg, gImg, bImg);

        if (settings.debug()) {
            
            long ts = System.currentTimeMillis();
            
            MiscDebug.writeImage(segImg,
                "segmented_" + lbl + "_" + settings.getDebugTag() + "_" + ts);
        }
        
        List<Integer> removeIndexes = new ArrayList<Integer>();
         
        for (Integer index : highDensityLargeGroups) {
            
            Set<PairInt> set = keypointsAndBounds.getKeyPointGroups().get(index.intValue());
            
            Set<PairInt> remove = new HashSet<PairInt>();
            
            for (PairInt p : set) {
                int v = segImg.getValue(p.getX(), p.getY());
                if (v < 1) {
                    remove.add(p);
                }
            }
            
            set.removeAll(remove);
            
            if (set.isEmpty()) {
                removeIndexes.add(index);
            }
        }
        
        for (int i = (removeIndexes.size() - 1); i > -1; --i) {
            int rmIdx = removeIndexes.get(i);
            keypointsAndBounds.getKeyPointGroups().remove(rmIdx);
        }
              
        keypointsAndBounds.getBoundingRegions().removeIndexes(removeIndexes);
        
        assert(keypointsAndBounds.getKeyPointGroups().size() == 
            keypointsAndBounds.getBoundingRegions().getPerimeterList().size());
        
    }
    
}
