package algorithms.imageProcessing.features;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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

    protected List<FeatureComparisonStat> rejectedBy2ndBest =
        new ArrayList<FeatureComparisonStat>();

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

    // trying the color space O1, O2, O3
    protected final IntensityClrFeatures clrFeaturesBinned1;
    protected final IntensityClrFeatures clrFeaturesBinned2;

    private List<FeatureComparisonStat> solutionStats = null;

    private List<PairInt> solutionMatched1 = null;
    private List<PairInt> solutionMatched2 = null;

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

        boolean matchPtByPt = false;
        if (matchPtByPt) {
            return matchPtByPt();
        } else {
            return matchBlobByBlob();
        }
    }

    private boolean matchPtByPt() throws IOException, NoSuchAlgorithmException {

        List<CornerRegion> corners1 = extractCorners(redBinnedImg1,
            greenBinnedImg1, blueBinnedImg1, "1");

        List<CornerRegion> corners2 = extractCorners(redBinnedImg2,
            greenBinnedImg2, blueBinnedImg2,  "2");

        boolean filterForLocalization = true;
        if (filterForLocalization) {
            log.info("filterForLocalization im1");
            filterForLocalization(redBinnedImg1, greenBinnedImg1, blueBinnedImg1,
                clrFeaturesBinned1, corners1);
            log.info("filterForLocalization im2");
            filterForLocalization(redBinnedImg2, greenBinnedImg2, blueBinnedImg2,
                clrFeaturesBinned2, corners2);
        }

        log.info("nPts after localization filter img1 = " + corners1.size());

        log.info("nPts after localization filter img2 = " + corners2.size());

        if (settings.debug()) {
            try {
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(corners1, imgBinned1.copyImage(),
                    "corners_filtered_" + settings.getDebugTag() + "_1_" + ts);
                MiscDebug.writeImage(corners2, imgBinned2.copyImage(),
                    "corners_filtered_" + settings.getDebugTag() + "_2_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }

        int dither = 1;

        CornerMatcher<CornerRegion> matcher = new CornerMatcher<CornerRegion>(dither);

        boolean matched = matcher.matchCorners(
            clrFeaturesBinned1, clrFeaturesBinned2, corners1, corners2,
            redBinnedImg1, greenBinnedImg1, blueBinnedImg1,
            redBinnedImg2, greenBinnedImg2, blueBinnedImg2,
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

        /*if (settings.debug()) {
            log.info("removed due to having close 2nd best:");
            for (FeatureComparisonStat stat : rb2j) {
                log.info(
                    String.format("%s %s", stat.getImg1Point().toString(),
                    stat.getImg2Point().toString()));
            }
        }*/

        copyToInstanceVars(stats, rb2j);

        return true;
    }

    private boolean matchBlobByBlob() throws IOException, NoSuchAlgorithmException {

        List<List<CornerRegion>> corners1List =
            extractCornersAndAssocWithBlobs(redBinnedImg1,
            greenBinnedImg1, blueBinnedImg1, imgBinned1.copyToImageExt(), "1");

        List<List<CornerRegion>> corners2List =
            extractCornersAndAssocWithBlobs(redBinnedImg2,
            greenBinnedImg2, blueBinnedImg2, imgBinned2.copyToImageExt(), "2");

        boolean filterForLocalization = false;
        if (filterForLocalization) {
            filterListsForLocalization(redBinnedImg1, greenBinnedImg1, blueBinnedImg1,
                clrFeaturesBinned1, corners1List);
            log.info("filterForLocalization im2");
            filterListsForLocalization(redBinnedImg2, greenBinnedImg2, blueBinnedImg2,
                clrFeaturesBinned2, corners2List);
        }

        if (settings.debug()) {
            long ts = MiscDebug.getCurrentTimeFormatted();

            ImageExt imgCp1 = (ImageExt) imgBinned1.copyImage();
            ImageIOHelper.<CornerRegion>addAlternatingColorCornerRegionListsToImage(
                corners1List, imgCp1, 0, 0, 1);
            MiscDebug.writeImage(imgCp1, "corners_1_" + settings.getDebugTag()
                + "_" + ts);

            ImageExt imgCp2 = (ImageExt) imgBinned2.copyImage();
            ImageIOHelper.<CornerRegion>addAlternatingColorCornerRegionListsToImage(
                corners2List, imgCp2, 0, 0, 1);
            MiscDebug.writeImage(imgCp2, "corners_2_" + settings.getDebugTag()
                + "_" + ts);
        }

        int dither = 1;

        CornerMatcher<CornerRegion> matcher = new CornerMatcher<CornerRegion>(dither);

        boolean matched = matcher.matchCornersByBlobs(
            clrFeaturesBinned1, clrFeaturesBinned2, corners1List, corners2List,
            redBinnedImg1, greenBinnedImg1, blueBinnedImg1,
            redBinnedImg2, greenBinnedImg2, blueBinnedImg2,
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

        /*if (settings.debug()) {
            log.info("removed due to having close 2nd best:");
            for (FeatureComparisonStat stat : rb2j) {
                log.info(
                    String.format("%s %s", stat.getImg1Point().toString(),
                    stat.getImg2Point().toString()));
            }
        }*/

        copyToInstanceVars(stats, rb2j);

        return true;
    }

    public void copyToInstanceVars(List<FeatureComparisonStat> stats) {

        this.solutionStats = new ArrayList<FeatureComparisonStat>();
        this.solutionMatched1 = new ArrayList<PairInt>();
        this.solutionMatched2 = new ArrayList<PairInt>();

        for (FeatureComparisonStat stat : stats) {
            solutionStats.add(stat.copy());
            solutionMatched1.add(stat.getImg1Point().copy());
            solutionMatched2.add(stat.getImg2Point().copy());
        }
    }

    private void copyToInstanceVars(List<FeatureComparisonStat> stats,
        List<FeatureComparisonStat> rejBy2ndBest) {

        copyToInstanceVars(stats);

        rejectedBy2ndBest.clear();

        rejectedBy2ndBest.addAll(rejBy2ndBest);
    }

    public List<FeatureComparisonStat> getRejectedBy2ndBest() {
        return rejectedBy2ndBest;
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

    private List<CornerRegion> extractCorners(GreyscaleImage rImg,
        GreyscaleImage gImg, GreyscaleImage bImg, String lbl) {

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

        log.info("before segmentation filters nPts=" + pixels.size());
        pixels = filterPointsBySegmentation0(rImg, gImg, bImg, pixels, lbl);
        log.info("after wavelet segmentation filter nPts=" + pixels.size());
        //pixels = filterPointsBySegmentation1(img, pixels, lbl);
        //log.info("after ciexy segmentation filter nPts=" + pixels.size());

        List<CornerRegion> corners = convert(pixels);

        if (settings.debug()) {
            try {
                MiscDebug.writeImage(corners, bImg.copyToColorGreyscale(),
                    "corners_" + settings.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }

        return corners;
    }

    private List<List<CornerRegion>> extractCornersAndAssocWithBlobs(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg, 
        ImageExt img, String lbl) {

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

        GreyscaleImage segImg = createCombinedSegmentation0(rImg, gImg, bImg);

        if (settings.debug()) {
            MiscDebug.writeImage(segImg,
                "segmented_" + lbl + "_" + settings.getDebugTag() + "_" + ts);
        }

        /*
        for unbinned:
        smallestGroupLimit = 100;
        largestGroupLimit = 5000;
        */
        int smallestGroupLimit = 5;
        int largestGroupLimit = Integer.MAX_VALUE;
        boolean filterOutImageBoundaryBlobs = false;
        boolean filterOutZeroPixels = true;

        List<Set<PairInt>> blobs =
            BlobsAndPerimeters.extractBlobsFromSegmentedImage(segImg,
                smallestGroupLimit, largestGroupLimit,
                filterOutImageBoundaryBlobs, filterOutZeroPixels, lbl);
        
        // look at statistics of delta E for the blobs
        //disconnectDifferentRegions(blobs, rImg, gImg, bImg);
        
        
        joinAdjacentSimilarBlobs(blobs, rImg, gImg, bImg);

        List<Set<PairInt>> pixelSetList = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < blobs.size(); ++i) {
            pixelSetList.add(new HashSet<PairInt>());
        }

        int nTot = 0;
        for (PairInt p : pixels) {
            for (int i = 0; i < blobs.size(); ++i) {
                Set<PairInt> blob = blobs.get(i);
                if (blob.contains(p)) {
                    pixelSetList.get(i).add(p);
                    nTot++;
                    break;
                }
            }
        }

        log.info("after association with blobs nPts = " + nTot);

        List<List<CornerRegion>> output = new ArrayList<List<CornerRegion>>();
        for (int i = 0; i < pixelSetList.size(); ++i) {
            Set<PairInt> pixelSet = pixelSetList.get(i);
            List<CornerRegion> corners = convert(pixelSet);
            output.add(corners);
        }

        if (settings.debug()) {
            ImageExt imgCp = (ImageExt) img.copyImage();
            try {
                ImageIOHelper.addAlternatingColorPointSetsToImage(pixelSetList,
                    0, 0, 2, imgCp);
                MiscDebug.writeImage(imgCp, "_blob_corners_" + lbl + "_" + 
                settings.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            }
            
        }

        return output;
    }

    protected List<FeatureComparisonStat> reviseStatsForFullImages(
        List<FeatureComparisonStat> stats, int prevBinFactor1,
        int prevBinFactor2, RotatedOffsets rotatedOffsets) {

        log.info("refine stats for full image reference frames");

        if (stats.isEmpty()) {
            return stats;
        }

        //TODO: when have a work clr descriptor, put the full size feature matching backing in

        List<FeatureComparisonStat> revised = new ArrayList<FeatureComparisonStat>();

        for (int i = 0; i < stats.size(); ++i) {

            FeatureComparisonStat stat = stats.get(i);

            PairInt p1 = stat.getImg1Point();
            p1.setX(p1.getX() * prevBinFactor1);
            p1.setY(p1.getY() * prevBinFactor1);

            PairInt p2 = stat.getImg2Point();
            p2.setX(p2.getX() * prevBinFactor2);
            p2.setY(p2.getY() * prevBinFactor2);

            revised.add(stat);
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

    private GreyscaleImage createCombinedSegmentation0(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg) {

        boolean use1D = true;
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        GreyscaleImage rSegImg = imageSegmentation.createGreyscale5(rImg, use1D);
        GreyscaleImage gSegImg = imageSegmentation.createGreyscale5(gImg, use1D);
        GreyscaleImage bSegImg = imageSegmentation.createGreyscale5(bImg, use1D);

        GreyscaleImage combined = rSegImg.copyImage();
        for (int i = 0; i < rSegImg.getWidth(); ++i) {
            for (int j = 0; j < rSegImg.getHeight(); ++j) {
                int g = gSegImg.getValue(i, j);
                int b = bSegImg.getValue(i, j);
                if (g > 0) {
                    combined.setValue(i, j, g);
                } else if (b > 0) {
                    combined.setValue(i, j, b);
                }
            }
        }
        return combined;
    }

    private Set<PairInt> filterPointsBySegmentation0(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg,
        Set<PairInt> pixels, String lbl) {

        GreyscaleImage segImg = createCombinedSegmentation0(rImg, gImg, bImg);

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
        List<List<CornerRegion>> cornersList) {

        for (int i = 0; i < cornersList.size(); ++i) {
            List<CornerRegion> corners = cornersList.get(i);
            filterForLocalization(rImg, gImg, bImg, clrFeatures, corners);
        }
    }

    /**
     * @return the solutionStats
     */
    public List<FeatureComparisonStat> getSolutionStats() {
        return solutionStats;
    }

    /**
     * @return the solutionMatched1
     */
    public List<PairInt> getSolutionMatched1() {
        return solutionMatched1;
    }

    /**
     * @return the solutionMatched2
     */
    public List<PairInt> getSolutionMatched2() {
        return solutionMatched2;
    }

    private void disconnectDifferentRegions(List<Set<PairInt>> blobs, 
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg) {
        
        log.info("stats for delta E 1994 of blobs:");
        CIEChromaticity cieC = new CIEChromaticity();
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < blobs.size(); ++i) {
            Set<PairInt> blob = blobs.get(i);
            List<PairInt> list = new ArrayList<PairInt>(blob);
            List<Double> deltaEs = new ArrayList<Double>();
            for (int i1 = 0; i1 < list.size(); ++i1) {
                PairInt p1 = list.get(i1);
                float[] lab1 = cieC.rgbToCIELAB(rImg.getValue(p1), 
                    gImg.getValue(p1), bImg.getValue(p1));
                for (int i2 = (i1 + 1); i2 < list.size(); ++i2) {
                    PairInt p2 = list.get(i2);
                    float[] lab2 = cieC.rgbToCIELAB(rImg.getValue(p2), 
                        gImg.getValue(p2), bImg.getValue(p2));
                    double deltaE = cieC.calcDeltaECIE94(lab1, lab2);
                    deltaEs.add(Double.valueOf(deltaE));
                }
            }
            double[] avgAndStDev = MiscMath.getAvgAndStDev(deltaEs);
            double[] xyCen = curveHelper.calculateXYCentroids(blob);
            log.info(String.format("(%d,%d)  dE avg=%.1f stdev=%.1f", 
                (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]), 
                (float)avgAndStDev[0], (float)avgAndStDev[1]));
        }
    }

    private void joinAdjacentSimilarBlobs(List<Set<PairInt>> blobs, 
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg) {
        
        // for highly textured regions like the tree leaves,
        // deltaE_1994 is a high 50 w/ st.dev. of 30
        // for uniform regions like smooth grass,
        // deltaE_1994 is below the JND of 2.3 at ~0.2 w/ st.dev. of ~3
     
        int n = 0;
        for (Set<PairInt> set : blobs) {
            n += set.size();
        }
        
        int[] blobIndexes = new int[n];
        int[] x = new int[n];
        int[] y = new int[n];
        int count = 0;
        for (int i = 0; i < blobs.size(); ++i) {
            Set<PairInt> set = blobs.get(i);
            for (PairInt p : set) {
                blobIndexes[count] = i;
                x[count] = p.getX();
                y[count] = p.getY();
                count++;
            }
        }
                
        //NearestPoints np = new NearestPoints(x, y);
        //public Set<Integer> findNeighborIndexes(int xCenter, int yCenter, float radius)
        /*
        x[0]
           x[1]  small deltaE and diff blob --> merge blobs
        */
    }

}
