package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * class encapsulating the steps from scale calculation to matching corners
 * to make correspondence lists.
 *
 * @author nichole
 */
public class FeatureMatcherWrapper {

    private final FeatureMatcherSettings settings;

    protected final int dither = 3;//4;

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

    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {

        img1 = image1;
        img2 = image2;

        doDetermineScale = true;

        this.settings = settings.copy();
    }

    /**
     * constructor accepting transformation parameters.  Note, for best results,
     * the standard deviations within parameters should be populated because they
     * are used as tolerances in matching.
     * @param image1
     * @param image2
     * @param parameters
     * @param settings
     */
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2,
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
            return solveForScale();
        }

        applyHistEqIfNeeded();

        cl = extractAndMatch(params);

        if (settings.debug()) {
            printMatches(cl);
        }

        return cl;
    }

    private CorrespondenceList solveForScale() throws IOException,
        NoSuchAlgorithmException {

        BlobScaleFinderWrapper scaleFinder = new BlobScaleFinderWrapper(img1,
            img2, settings, rotatedOffsets, dither);

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

        int binFactor1 = scaleFinder.getBinFactor1();
        int binFactor2 = scaleFinder.getBinFactor2();

        int nLimit = 16;

        if (binFactor1 != 1 || binFactor2 != 1) {

            //stats need to be revised for the location in the full size
            //image in order to be usable for correspondence
            List<FeatureComparisonStat> revisedStats = reviseStatsForFullImages(
                stats, binFactor1, binFactor2);

            stats = revisedStats;

            TransformationParameters revisedParams = null;

            if (stats.size() > 0) {

                revisedParams = MiscStats.calculateTransformation(1, 1, stats,
                    new float[4], false);
            }

            if (revisedParams != null) {
                params = revisedParams;
            } else {
                log.warning("possible ERROR in revision of stats");
            }
        }

        if (settings.debug()) {
            printMatches(stats);
        }
        
        //TODO: might add back in the check for quadrant coverage of intersection
        // and subsequent extraction of more features if coverage is low in any
        // quadrant

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

    private List<FeatureComparisonStat> reviseStatsForFullImages(
        List<FeatureComparisonStat> stats, int prevBinFactor1, int prevBinFactor2) {

        log.info("refine stats for full image reference frames");

        List<FeatureComparisonStat> revised = new ArrayList<FeatureComparisonStat>();

        FeatureMatcher featureMatcher = new FeatureMatcher();

        IntensityFeatures2 features1 = new IntensityFeatures2(5,
            settings.useNormalizedFeatures(), rotatedOffsets);
        features1.calculateGradientWithGreyscale(gsImg1);

        IntensityFeatures2 features2 = new IntensityFeatures2(5,
            settings.useNormalizedFeatures(), rotatedOffsets);
        features2.calculateGradientWithGreyscale(gsImg2);

        int dither2 = 1 * (Math.max(prevBinFactor1, prevBinFactor2));
        if (params.getTranslationX() > 5 || params.getTranslationY() > 5) {
            if (dither2 < 3) {
                dither2 = 3;
            }
        }

        int rotD = Math.round(params.getRotationInDegrees());

        final int rotationTolerance = 20;

        for (int i = 0; i < stats.size(); ++i) {

            FeatureComparisonStat stat = stats.get(i);

            int x1 = stat.getImg1Point().getX() * stat.getBinFactor1();
            int y1 = stat.getImg1Point().getY() * stat.getBinFactor1();
            int x2 = stat.getImg2Point().getX() * stat.getBinFactor2();
            int y2 = stat.getImg2Point().getY() * stat.getBinFactor2();

            // have to discard the best angles found in stat and derive new
            // for these higher resolution images
            FeatureComparisonStat compStat =
                featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, x1, y1, x2, y2, dither2, rotD,
                    rotationTolerance, gsImg1, gsImg2);

            if (compStat == null ||
                (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())) {
                continue;
            }

            revised.add(compStat);
        }
        
        return revised;
    }

    private void printMatches(List<FeatureComparisonStat> stats) {
        if (stats == null) {
            return;
        }
        List<PairInt> matched1 = new ArrayList<PairInt>();
        List<PairInt> matched2 = new ArrayList<PairInt>();
        populateLists(stats, matched1, matched2);

        printMatches(matched1, matched2);
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

}
