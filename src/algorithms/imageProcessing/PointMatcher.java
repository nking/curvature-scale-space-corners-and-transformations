package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.SubsetChooser;
import algorithms.compGeometry.PointPartitioner;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairFloat;
import algorithms.util.PairFloatArray;
import algorithms.util.PairFloatArrayUnmodifiable;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonPlotterPNG;
import algorithms.util.RangeInt;
import algorithms.util.ResourceFinder;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

/**
 * class to match the points extracted from two images.
 *
 * the transformation parameters of translation, rotation and scale are
 * found given the two sets of points.

 * The resulting transformation returns rotation in a clockwise direction.

 * @author nichole
 */
public final class PointMatcher extends AbstractPointMatcher {

    private final Logger log = Logger.getLogger(this.getClass().getName());

    // when true, there is additional logging and creation of simplex data
    // that can be viewed with resources/plot_3d_simplex_n10.html
    protected boolean debug = false;

    protected boolean usePreSearchAlt1 = true;

    protected static int minTolerance = 5;

    //NOTE: solution is sensitive to this number
    private float generalTolerance = 10;

    public static float toleranceGridFactor = 4.f;

    protected int largeSearchLimit = 20;

    protected float transDeltaPreSearch0 = 4;//10;
    protected float rotDeltaPreSearch0 = 2;

    protected float rotHalfRangeInDegreesPreSearch1Alt2 = 10;
    protected float rotDeltaInDegreesPreSearch1Alt2 = 2;
    protected float transXHalfRangePreSearch1Alt2 = 200;
    protected float transXDeltaPreSearch1Alt2 = 15;

    public void setLargeSearchLimit(int limit) {
        largeSearchLimit = limit;
    }

    /**
     * given the unmatched point sets, and the expectation that the 2 sets
     * have rotation wrt one another less than 30 degrees,
     * calculate the best fitting transformation
     * and then apply the transformation to the first set to return matched
     * points within a tolerance.
     *
     * @param skylineCorners1
     * @param skylineCorners2
     * @param corners1
     * @param corners2
     * @param outputMatchedSet1
     * @param outputMatchedSet2
     * @param useLargestToleranceForOutput use the largest tolerance for
     * applying the transformation to point sets during matching.  the largest
     * tolerance is the class variable generalTolerance.
     * If useLargestToleranceForOutput is false, the transformation's best
     * fit is used during matching (which should provide a smaller but more
     * certain matched output).  If this method is used as a precursor to
     * projection (epipolar) solvers of sets that do have projection components,
     * one might prefer to set this to true to allow more to be matched.
     * @param useGreedyMatching
     * @return best fitting transformation between unmatched points sets
     *
     * @return
     */
    public TransformationPointFit performMatchingForMostlyVerticalTranslation(
        PairIntArray skylineCorners1, PairIntArray skylineCorners2,
        PairIntArray corners1, PairIntArray corners2,
        PairIntArray outputMatchedSet1, PairIntArray outputMatchedSet2,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        float scale = 1;
        float rotationLowLimitInDegrees = 335;
        float rotationHighLimitInDegrees = 25;

        TransformationPointFit fit = calculateEuclideanTransformation(
            skylineCorners1, skylineCorners2,
            corners1, corners2,
            scale, rotationLowLimitInDegrees, rotationHighLimitInDegrees);

        if (fit == null) {
            return null;
        }

        float translationXTolerance, translationYTolerance;

        if (useLargestToleranceForOutput) {
            translationXTolerance = generalTolerance * (float)Math.sqrt(1./2);
            translationYTolerance = translationXTolerance;
        } else {
            translationXTolerance = 2 * fit.getTranslationXTolerance();
            translationYTolerance = 2 * fit.getTranslationYTolerance();
        }

        PairIntArray comb1 = skylineCorners1.copy();
        comb1.addAll(corners1);
        PairIntArray comb2 = skylineCorners2.copy();
        comb2.addAll(corners2);

        //TODO: might need tolerance to be as high as 20

        match(fit.getParameters(), comb1, comb2,
            outputMatchedSet1, outputMatchedSet2,
            translationXTolerance, translationYTolerance, useGreedyMatching);

        return fit;
    }

    /**
     * given the unmatched point sets, and the expectation that the 2 sets
     * have rotation wrt one another less than 30 degrees,
     * calculate the best fitting transformation
     * and then apply the transformation to the first set to return matched
     * points within a tolerance.
     *
     * @param set1
     * @param set2
     * @param outputMatchedSet1
     * @param outputMatchedSet2
     * @param useLargestToleranceForOutput use the largest tolerance for
     * applying the transformation to point sets during matching.  the largest
     * tolerance is the class variable generalTolerance.
     * If useLargestToleranceForOutput is false, the transformation's best
     * fit is used during matching (which should provide a smaller but more
     * certain matched output).  If this method is used as a precursor to
     * projection (epipolar) solvers of sets that do have projection components,
     * one might prefer to set this to true to allow more to be matched.
     * @param useGreedyMatching
     * @return best fitting transformation between unmatched points sets
     *
     * @return
     */
    public TransformationPointFit performMatchingForMostlyVerticalTranslation(
        PairIntArray set1, PairIntArray set2,
        PairIntArray outputMatchedSet1, PairIntArray outputMatchedSet2,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        long nPairPerm = estimateNStepsPairCalculation(set1.getN(),
            set2.getN());

        log.info("nPairPerm=" + nPairPerm);

        float scale = 1;
        float rotationLowLimitInDegrees = 335;
        float rotationHighLimitInDegrees = 25;

        TransformationPointFit fit = null;

        if (nPairPerm > 1e9) {

            log.warning("assuming scale=1 to solve." +
            "currently not able to solve for scale too in large sets in reasonable time so " +
            " you can reduce your set sizes or use vert trans methods for better solutions.");
            fit = preSearch(set1, set2, scale, rotationLowLimitInDegrees,
                rotationHighLimitInDegrees, useGreedyMatching);

        } else {
            fit = calculateEuclideanTransformation(
                set1, set2, scale, rotationLowLimitInDegrees,
                rotationHighLimitInDegrees);
        }

        if (fit == null) {
            return null;
        }

        float translationXTolerance, translationYTolerance;

        if (useLargestToleranceForOutput) {
            translationXTolerance = generalTolerance * (float)Math.sqrt(1./2);
            translationYTolerance = translationXTolerance;
        } else {
            translationXTolerance = 2 * fit.getTranslationXTolerance();
            translationYTolerance = 2 * fit.getTranslationYTolerance();
        }

        match(fit.getParameters(), set1, set2, outputMatchedSet1, outputMatchedSet2,
            translationXTolerance, translationYTolerance, useGreedyMatching);

        return fit;
    }

    /**
     * given the unmatched point sets, calculates the best fitting transformation
     * and then applies the transformation to the first set to return matched
     * points within a tolerance.
     *
     * @param unmatchedLeftXY
     * @param unmatchedRightXY
     * @param outputMatchedLeftXY
     * @param outputMatchedRightXY
     * @param useLargestToleranceForOutput use the largest tolerance for
     * applying the transformation to point sets during matching.  the largest
     * tolerance is the class variable generalTolerance.
     * If useLargestToleranceForOutput is false, the transformation's best
     * fit is used during matching (which should provide a smaller but more
     * certain matched output).  If this method is used as a precursor to
     * projection (epipolar) solvers of sets that do have projection components,
     * one might prefer to set this to true to allow more to be matched.
     * @return best fitting transformation between unmatched points sets
     *
     * @return
     */
    public TransformationPointFit performMatching(
        PairIntArray unmatchedLeftXY, PairIntArray unmatchedRightXY,
        PairIntArray outputMatchedLeftXY, PairIntArray outputMatchedRightXY,
        boolean useLargestToleranceForOutput) {

        TransformationPointFit transFit =
            calculateEuclideanTransformation(unmatchedLeftXY, unmatchedRightXY);

        if (transFit == null) {
            return null;
        }

        int nMaxMatchable = (unmatchedLeftXY.getN() < unmatchedRightXY.getN()) ?
            unmatchedLeftXY.getN() : unmatchedRightXY.getN();

        //transFit.setMaximumNumberMatchable(nMaxMatchable);

        if (nMaxMatchable == 0) {
            return transFit;
        }

        // -- transform filtered1 for matching and evaluation

        float translationXTolerance, translationYTolerance;

        if (useLargestToleranceForOutput) {
            translationXTolerance = generalTolerance * (float)Math.sqrt(1./2);
            translationYTolerance = translationXTolerance;
        } else {
            translationXTolerance = 2 * transFit.getTranslationXTolerance();
            translationYTolerance = 2 * transFit.getTranslationYTolerance();
        }

        boolean useGreedyMatching = true;

        match(transFit.getParameters(), unmatchedLeftXY, unmatchedRightXY,
            outputMatchedLeftXY, outputMatchedRightXY,
            translationXTolerance, translationYTolerance, useGreedyMatching);

        return transFit;
    }

    /**
     * calculate for Euclidean transformation from set1 to set2 for
     * the points being unmatched.
     *
     * @param set1
     * @param set2
     * @return
     */
    public TransformationPointFit calculateEuclideanTransformation(
        PairIntArray set1, PairIntArray set2) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        long nPairPerm = estimateNStepsPairCalculation(n1, n2);

        log.info("nPairPerm=" + nPairPerm + " n1=" + n1 + " n2=" + n2);

        if (nPairPerm > 10e9) {
            log.warning("assuming scale=1 to solve." +
            "currently not able to solve for scale too in large sets in reasonable time so " +
            " you can reduce your set sizes or use vert trans methods for better solutions.");
            float scale = 1;
            boolean useGreedyMatching = true;
            return preSearch(set1, set2, scale, useGreedyMatching);
        }

        boolean earlyConvergeReturn = true;

        long t0 = System.currentTimeMillis();

        List<TransformationPointFit> fits = calculateEuclideanTransformationUsingPairs(
            set1, set2, earlyConvergeReturn);

        long t1 = System.currentTimeMillis();

        log.info("calculateEuclideanTransformationUsingPairs seconds=" +
            ((t1 - t0)*1e-3));

        if (fits == null || fits.isEmpty()) {
            return null;
        }

        boolean useGreedyMatching = true;

        // find the best fit to corners1, corners2
        TransformationPointFit bestFit = null;
        for (int i = 0; i < fits.size(); ++i) {

            TransformationPointFit fit  = fits.get(i);

            TransformationPointFit fitA = evaluateForUnmatched(
                fit.getParameters(), set1, set2,
                generalTolerance, generalTolerance, useGreedyMatching);

            log.info("test fit w/ corners=" + fitA.toString());

            if (fitIsBetter(bestFit, fitA)) {
                bestFit = fitA;
            }
        }

        log.info("best fitting to all corners from best fits="
            + bestFit.toString());

        TransformationPointFit fit2 = refineAfterCalculationWithPairs(
            bestFit.getParameters().copy(), set1, set2,
            useGreedyMatching);

        if (fitIsBetter(bestFit, fit2)) {
            bestFit = fit2;
        }

        return bestFit;
    }

    /**
     * calculate for Euclidean transformation from set1 to set2 for
     * unmatched points.
     *
     * @param scale
     * @param rotationLowLimitInDegrees
     * @param rotationHighLimitInDegrees
     * @return
     */
    public TransformationPointFit calculateEuclideanTransformation(
        PairIntArray skylineCorners1, PairIntArray skylineCorners2,
        PairIntArray corners1, PairIntArray corners2,
        float scale,
        float rotationLowLimitInDegrees, float rotationHighLimitInDegrees) {

        boolean earlyConvergeReturn = false;
        boolean useLargestToleranceForOutput = true;
        boolean useGreedyMatching = true;

        long t0 = System.currentTimeMillis();

        List<TransformationPointFit> fits = new ArrayList<TransformationPointFit>();
        List<Float> tols = new ArrayList<Float>();

        //TODO: for an image with large projection effects, this may need to be a larger number.
        float[] tolerances = new float[]{10, 15};
        for (float tol : tolerances) {

            generalTolerance = tol;

            List<TransformationPointFit> fits2 =
                calculateEuclideanTransformationForSmallSets(
                skylineCorners1, skylineCorners2,
                scale, rotationLowLimitInDegrees, rotationHighLimitInDegrees,
                earlyConvergeReturn, useLargestToleranceForOutput, useGreedyMatching);

            if (fits2 != null) {
                fits.addAll(fits2);
                for (int i = 0; i < fits2.size(); ++i) {
                    tols.add(Float.valueOf(tol));
                }
            }
        }

        long t1 = System.currentTimeMillis();

        log.info("(" + ((t1 - t0)*1e-3) + " seconds for eucl w/ pairs of skyline.");

        if (fits.isEmpty()) {

            PairIntArray comb1 = skylineCorners1.copy();
            comb1.addAll(corners1);
            PairIntArray comb2 = skylineCorners2.copy();
            comb2.addAll(corners2);

            return calculateEuclideanTransformation(comb1, comb2, scale,
                rotationLowLimitInDegrees, rotationHighLimitInDegrees);
        }

        // find the best fit to corners1, corners2
        TransformationPointFit bestFit = null;
        for (int i = 0; i < fits.size(); ++i) {

            TransformationPointFit fit  = fits.get(i);

                TransformationPointFit fitA = evaluateForUnmatched(
                    fit.getParameters(), corners1, corners2,
                    generalTolerance, generalTolerance, useGreedyMatching);

                log.info("test skyline params to corners=" + fitA.toString());

                if (fitIsBetter(bestFit, fitA)) {
                    bestFit = fitA;
                    generalTolerance = tols.get(i).floatValue();
                }
        }

        log.info("best fitting to all corners from best fits="
            + bestFit.toString());

        // fit the entire datasets using skyline fit

        PairIntArray comb1 = skylineCorners1.copy();
        comb1.addAll(corners1);
        PairIntArray comb2 = skylineCorners2.copy();
        comb2.addAll(corners2);

        TransformationPointFit fit2 = refineAfterCalculationWithPairs(
            bestFit.getParameters().copy(), comb1, comb2, useGreedyMatching);

        if (fitIsBetter(bestFit, fit2)) {
            bestFit = fit2;
        }

        return bestFit;
    }

    /**
     * calculate for Euclidean transformation from set1 to set2 for
     * unmatched points.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param rotationLowLimitInDegrees
     * @param rotationHighLimitInDegrees
     * @return
     */
    public TransformationPointFit calculateEuclideanTransformation(
        PairIntArray set1, PairIntArray set2,
        float scale,
        float rotationLowLimitInDegrees, float rotationHighLimitInDegrees) {

        boolean earlyConvergeReturn = true;
        boolean useLargestToleranceForOutput = true;
        boolean useGreedyMatching = true;

        //TODO: consider using multiple tolerances if sets are small enough

        List<TransformationPointFit> fits =
            calculateEuclideanTransformationUsingPairs(set1, set2,
            scale, rotationLowLimitInDegrees, rotationHighLimitInDegrees,
            earlyConvergeReturn, useLargestToleranceForOutput, useGreedyMatching);

        if (fits.isEmpty()) {
            return null;
        }

        TransformationPointFit bestFit = null;

        for (int i = 0; i < fits.size(); ++i) {

            TransformationPointFit fit  = fits.get(i);

            log.info("compare fit=" + fit.toString());

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
        }

        log.info("best fitting to all corners from best fits="
            + bestFit.toString());

        TransformationPointFit fit2 = refineAfterCalculationWithPairs(
            bestFit.getParameters().copy(), set1, set2,
            useGreedyMatching);

        if (fitIsBetter(bestFit, fit2)) {
            bestFit = fit2;
        }

        return bestFit;
    }

    /**
     * an experimental method that creates corners to make matching points sets.
     * At a later time there will be a more formal stereoscopic matcher in the project
     * to register images and make disparity maps with.
     *
     * @param leftImg
     * @param rightImg
     * @return
     */
    public TransformationPointFit stereoscopicPointMatcher(ImageExt leftImg,
        ImageExt rightImg, PairIntArray outputMatchedLeft,
        PairIntArray outputMatchedRight) throws IOException, NoSuchAlgorithmException {

        ImageExt leftImgOrig = (ImageExt)leftImg.copyImage();
        ImageExt rightImgOrig = (ImageExt)rightImg.copyImage();

        /*
        uses color segmentation, then binary segmentation, then a gap filling
        algorithm, then a binning algorithm to make the image size < 200 x 200,
        then  creates corners, then calculates the euclidean transformation
        from the best pair's solution.
         */

        int binFactor1 = (int) Math.ceil(Math.max((float)leftImg.getWidth()/200.f,
            (float)leftImg.getHeight()/200.));
        //int binFactor2 = (int) Math.round(Math.max((float)rightImg.getWidth()/200.f,
        //    (float)rightImg.getHeight()/200.));

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        GreyscaleImage csImg1, csImg2;

        long t0 = System.currentTimeMillis();

        if (false) {
            csImg1 = imageProcessor.createGreyscaleFromColorSegmentation(leftImg, 3);
            imageSegmentation.applyImageSegmentation(csImg1, 2);
            Map<Integer, Integer> freqMap = Histogram.createAFrequencyMap(csImg1);
            int lowValue1 = Integer.MAX_VALUE;
            int highValue1 = Integer.MIN_VALUE;
            for (Integer v : freqMap.keySet()) {
                if (v.intValue() < lowValue1) {
                    lowValue1 = v.intValue();
                }
                if (v.intValue() > highValue1) {
                    highValue1 = v.intValue();
                }
            }
            imageProcessor.fillInPixels(csImg1, lowValue1, 50);
            imageProcessor.fillInPixels(csImg1, highValue1, 50);
            csImg1 = imageProcessor.binImage(csImg1, binFactor1);
            ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                + "/color_segmentation1.png", csImg1);

            csImg2 = imageProcessor.createGreyscaleFromColorSegmentation(rightImg, 3);
            imageSegmentation.applyImageSegmentation(csImg2, 2);
            freqMap = Histogram.createAFrequencyMap(csImg2);
            int lowValue2 = Integer.MAX_VALUE;
            int highValue2 = Integer.MIN_VALUE;
            for (Integer v : freqMap.keySet()) {
                if (v.intValue() < lowValue2) {
                    lowValue2 = v.intValue();
                }
                if (v.intValue() > highValue2) {
                    highValue2 = v.intValue();
                }
            }
            imageProcessor.fillInPixels(csImg2, lowValue2, 50);
            imageProcessor.fillInPixels(csImg2, highValue2, 50);
            csImg2 = imageProcessor.binImage(csImg2, binFactor1);
            ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                + "/color_segmentation2.png", csImg2);
        } else {
            ImageExt leftImgBinned = imageProcessor.binImage(leftImg, binFactor1);

            csImg1 = imageProcessor.createGreyscaleFromColorSegmentation(leftImgBinned, 3);
            imageSegmentation.applyImageSegmentation(csImg1, 2);

            ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                + "/color_segmentation1.png", csImg1);

            ImageExt rightImgBinned = imageProcessor.binImage(rightImg, binFactor1);
            csImg2 = imageProcessor.createGreyscaleFromColorSegmentation(rightImgBinned, 3);
            imageSegmentation.applyImageSegmentation(csImg2, 2);

            ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                + "/color_segmentation2.png", csImg2);
        }
        long t1 = System.currentTimeMillis();
        log.info("segmentation " + ((t1 - t0)*1e-3) + " seconds");



        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(csImg1);
        detector.doNotPerformHistogramEqualization();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        PairIntArray corners1 = detector.getCornersInOriginalReferenceFrame();

        CurvatureScaleSpaceCornerDetector detector2 = new
            CurvatureScaleSpaceCornerDetector(csImg2);
        detector2.doNotPerformHistogramEqualization();
        //detector2.findCornersIteratively(nPreferredCorners, nCrit);
        detector2.findCorners();
        PairIntArray corners2 = detector2.getCornersInOriginalReferenceFrame();

        writeImage(ImageIOHelper.convertImage(csImg1), corners1,
            "corners1_" + 0 + "_" + 0 + ".png");
        writeImage(ImageIOHelper.convertImage(csImg2), corners2,
            "corners2_" + 0 + "_" + 0 + ".png");

        t0 = System.currentTimeMillis();

        boolean useLargestToleranceForOutput = false;//true;
        boolean useGreedyMatching = true;

        log.info("n1=" + corners1.getN() + " n2=" + corners2.getN());

        TransformationPointFit fit = calculateEuclideanLeftRightTransformation(
            corners1, corners2, useLargestToleranceForOutput, useGreedyMatching);

        t1 = System.currentTimeMillis();

        log.info("(" + ((t1 - t0)*1e-3) + " seconds) fit=" + fit.toString());

        // make a matching list with small tolerance
        PairIntArray outputMatchedSet1 = new PairIntArray();
        PairIntArray outputMatchedSet2 = new PairIntArray();
        float tolTransX = 3;
        float tolTransY = tolTransX;

        match(fit.getParameters(), corners1, corners2,
            outputMatchedSet1, outputMatchedSet2, tolTransX, tolTransY,
            useGreedyMatching);

        Transformer transformer = new Transformer();

        writeTransformed(fit.getParameters(), ImageIOHelper.convertImage(csImg2),
            corners1, corners2, "transformed_binned.png");

        // then transform the solution by a scale of binFactor1.
        TransformationParameters paramsOrigScale =
            transformer.applyScaleTransformation(fit.getParameters(), binFactor1);
        paramsOrigScale.setScale(1);

        // ==== this is only the initial solution.
        // need to use it to construct
        // the rest with either stereo matching techniques such as
        // sum of absolute differences of pixel blocks from left and right
        // OR, a round of corners at full resolution of image.

        // temporarily returning the fit with modified parameters
        TransformationPointFit fit2 = new TransformationPointFit(
            paramsOrigScale, fit.getNumberOfMatchedPoints(),
            fit.getMeanDistFromModel(), fit.getStDevFromMean(),
            fit.getTranslationXTolerance(), fit.getTranslationYTolerance());

/*
        writeImage(leftImgOrig, corners1OrigScale,
                    "corners1_" + 1 + "_" + 1 + ".png");
        writeImage(rightImgOrig, corners2OrigScale,
            "corners2_" + 1 + "_" + 1 + ".png");

        writeTransformed(paramsOrigScale, rightImgOrig.copyImage(),
            corners1OrigScale, corners2OrigScale, "transformed_" + 0 + ".png");

        //writeTransformed(fit.getParameters(), rightImg.copyImage(),
        //    outputMatchedScene, outputMatchedModel,
        //    "transformedCorners_" + 0 + "_" + 0 + ".png");
 */
        return fit2;
    }

private void writeImage(Image image1, PairIntArray points1,
    String outputImageName) {
    try {
        String dirPath = ResourceFinder.findDirectory("bin");
        int nExtraForDot = 1;
        ImageIOHelper.addCurveToImage(points1, image1, nExtraForDot,
            255, 0, 0);
        ImageIOHelper.writeOutputImage(dirPath + "/" + outputImageName, image1);
    } catch (Exception e) {
         e.printStackTrace();
        log.severe("ERROR: " + e.getMessage());
    }
}
private void writeTransformed(TransformationParameters parameters,
    Image image2, PairIntArray points1,
    PairIntArray points2, String outputImageName) {
    Transformer transformer = new Transformer();
    try {
        String dirPath = ResourceFinder.findDirectory("bin");
        PairIntArray trPoints1 = transformer.applyTransformation(parameters,
            points1);
        PairIntArray transformedPoints1 = new PairIntArray();
        for (int i = 0; i < trPoints1.getN(); ++i) {
            int x = trPoints1.getX(i);
            int y = trPoints1.getY(i);
            if (x < 0 || x > (image2.getWidth() - 1)) {
                continue;
            }
            if (y < 0 || y > (image2.getHeight() - 1)) {
                continue;
            }
            transformedPoints1.add(x, y);
        }

        int nExtraForDot = 1;
        ImageIOHelper.addCurveToImage(transformedPoints1, image2, nExtraForDot,
            255, 0, 0);
        ImageIOHelper.addCurveToImage(points2, image2, nExtraForDot,
            0, 0, 255);
        ImageIOHelper.writeOutputImage(dirPath + "/" + outputImageName, image2);
    } catch (Exception e) {
         e.printStackTrace();
        log.severe("ERROR: " + e.getMessage());
    }
}

    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     *
     * It's expected that the invoker of this method is trying to solve for
     * translation for sets of points like corners in images.  This assumption
     * means that the number of point pair combinations is always far less
     * than the pixel combinations of translations over x and y.
     *
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     *
     * ALSO NOTE: if you know a better solution exists for translation
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be set to a smaller value and passed
     * to this method.
     * @param params transformation parameters to apply to matched1
     * @param matched1Transformed
     * @param matched2 set of points from image 2 that are matched to points in
     * matched1 with same indexes
     * @return
     */
    public TransformationPointFit evaluateFitForMatchedTransformed(
        TransformationParameters params, PairFloatArray matched1Transformed,
        PairIntArray matched2) {

        if (matched1Transformed == null) {
            throw new IllegalArgumentException(
                "matched1Transformed cannot be null");
        }
        if (matched2 == null) {
            throw new IllegalArgumentException("matched2 cannot be null");
        }
        if (matched1Transformed.getN() != matched2.getN()) {
            throw new IllegalArgumentException(
                "the 2 point sets must have the same length");
        }

        double[] diff = new double[matched1Transformed.getN()];

        double avg = 0;

        for (int i = 0; i < matched1Transformed.getN(); i++) {

            float transformedX = matched1Transformed.getX(i);
            float transformedY = matched1Transformed.getY(i);

            int x2 = matched2.getX(i);
            int y2 = matched2.getY(i);

            double dx = x2 - transformedX;
            double dy = y2 - transformedY;

            diff[i] = Math.sqrt(dx*dx + dy*dy);

            avg += diff[i];
        }

        avg /= (double)matched2.getN();

        double stdDev = 0;

        for (int i = 0; i < matched2.getN(); i++) {

            double d = diff[i] - avg;

            stdDev += (d * d);
        }

        stdDev = Math.sqrt(stdDev/((double)matched2.getN() - 1));

        TransformationPointFit fit = new TransformationPointFit(params.copy(),
            matched2.getN(), avg, stdDev, Float.MAX_VALUE, Float.MAX_VALUE);

        return fit;
    }

    protected TransformationPointFit evaluateFitForUnMatchedTransformedGreedy(
        TransformationParameters params, PairFloatArray unmatched1Transformed,
        PairIntArray unmatched2, float tolTransX, float tolTransY) {

        if (unmatched1Transformed == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (unmatched2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }

        PairFloatArray trFiltered1 = new PairFloatArray();
        PairIntArray filtered2 = new PairIntArray();

        reduceToIntersection(unmatched1Transformed, unmatched2, trFiltered1,
            filtered2, tolTransX, tolTransY);

        int n = trFiltered1.getN();

        Set<Integer> chosen = new HashSet<Integer>();

        double[] diffs = new double[n];
        int nMatched = 0;
        double avg = 0;

        for (int i = 0; i < n; i++) {

            float transformedX = trFiltered1.getX(i);
            float transformedY = trFiltered1.getY(i);

            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < filtered2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }

                float dx = transformedX - filtered2.getX(j);
                float dy = transformedY - filtered2.getY(j);

                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                    continue;
                }

                float diff = (float)Math.sqrt(dx*dx + dy*dy);

                if (diff < minDiff) {
                    minDiff = diff;
                    min2Idx = j;
                }
            }

            if (minDiff < Double.MAX_VALUE) {
                diffs[nMatched] = minDiff;
                nMatched++;
                chosen.add(Integer.valueOf(min2Idx));
                avg += minDiff;
            }
        }

        avg = (nMatched == 0) ? Double.MAX_VALUE :
            avg / (double)nMatched;

        double stDev = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avg;
            stDev += (d * d);
        }

        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            Math.sqrt(stDev/((double)nMatched - 1.));

        TransformationPointFit fit = new TransformationPointFit(params.copy(),
            nMatched, avg, stDev, tolTransX, tolTransY);

        int nMaxMatchable = Math.min(trFiltered1.getN(), filtered2.getN());

        fit.setMaximumNumberMatchable(nMaxMatchable);

        return fit;
    }

    protected TransformationPointFit evaluateFitForUnMatchedTransformedOptimal(
        TransformationParameters params, PairFloatArray unmatched1Transformed,
        PairIntArray unmatched2, float tolTransX, float tolTransY) {

        if (unmatched1Transformed == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (unmatched2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }

        PairFloatArray trFiltered1 = new PairFloatArray();
        PairIntArray filtered2 = new PairIntArray();

        reduceToIntersection(unmatched1Transformed, unmatched2,
            trFiltered1, filtered2, tolTransX, tolTransY);

        int n = trFiltered1.getN();

        float[][] matchedIndexesAndDiffs = calculateMatchUsingOptimal(
            trFiltered1, filtered2, tolTransX, tolTransY);

        int nMatched = matchedIndexesAndDiffs.length;
        double avg = 0;
        double stDev = 0;

        for (int i = 0; i < nMatched; i++) {
            avg += matchedIndexesAndDiffs[i][2];
        }
        avg = (nMatched == 0) ? Double.MAX_VALUE : (avg / (double)nMatched);

        for (int i = 0; i < nMatched; i++) {

            double d = matchedIndexesAndDiffs[i][2] - avg;

            stDev += (d * d);
        }
        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            (Math.sqrt(stDev/((double)nMatched - 1.)));

        TransformationPointFit fit = new TransformationPointFit(params.copy(),
            nMatched,  avg, stDev, tolTransX, tolTransY);

        int nMaxMatchable = Math.min(trFiltered1.getN(),
            filtered2.getN());

        fit.setMaximumNumberMatchable(nMaxMatchable);

        return fit;
    }

    protected TransformationPointFit evaluateFitForUnMatchedOptimal(
        PairFloatArrayUnmodifiable scaledRotatedSet1, float transX, float transY,
        float tolTransX, float tolTransY,
        PairIntArray set2, float scale, float rotationInRadians) {

        if (scaledRotatedSet1 == null) {
            throw new IllegalArgumentException(
            "unmatched1Transformed cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException(
            "unmatched2 cannot be null");
        }

        //TODO: consider filtering for the intersection

        float[][] matchedIndexesAndDiffs = calculateMatchUsingOptimal(
            scaledRotatedSet1, transX, transY, tolTransX, tolTransY,
            set2);

        int nMatched = matchedIndexesAndDiffs.length;
        double avg = 0;
        double stDev = 0;

        for (int i = 0; i < nMatched; i++) {
            avg += matchedIndexesAndDiffs[i][2];
        }
        avg = (nMatched == 0) ? Double.MAX_VALUE : (avg / (double)nMatched);

        for (int i = 0; i < nMatched; i++) {

            double d = matchedIndexesAndDiffs[i][2] - avg;

            stDev += (d * d);
        }
        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            (Math.sqrt(stDev/((double)nMatched - 1.)));

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);

        TransformationPointFit fit = new TransformationPointFit(params, nMatched,
            avg, stDev, tolTransX, tolTransY);

        fit.setMaximumNumberMatchable(Math.max(scaledRotatedSet1.getN(),
            set2.getN()));

        return fit;
    }

    /**
     * given the indexes and residuals from optimal matching, populate
     * outputMatched1 and outputMatched2;
     *
     * @param set1
     * @param set2
     * @param transTolXY
     * @param matchedIndexesAndDiffs two dimensional array holding the matched
     * indexes and the distances between the model and the point for that pair.
     * each row holds {idx1, idx2, diff}
     * @param outputMatched1 the container to hold the output matching points
     * for image 1 that are paired with outputMatched2 as a result of running
     * this method.
     * @param outputMatched2 the container to hold the output matching points
     * for image 2 that are paired with outputMatched1 as a result of running
     * this method.
     */
    public void matchPoints(
        PairIntArray set1, PairIntArray set2,
        float transTolXY,
        float[][] matchedIndexesAndDiffs,
        PairIntArray outputMatched1, PairIntArray outputMatched2) {

        if (matchedIndexesAndDiffs == null) {
            return;
        }

        for (int i = 0; i < matchedIndexesAndDiffs.length; i++) {

            int idx1 = (int)matchedIndexesAndDiffs[i][0];
            int idx2 = (int)matchedIndexesAndDiffs[i][1];
            float diff = matchedIndexesAndDiffs[i][2];

            if (diff < transTolXY) {
                outputMatched1.add(set1.getX(idx1), set1.getY(idx1));
                outputMatched2.add(set2.getX(idx2), set2.getY(idx2));
            }
        }
    }

    /**
     * given the indexes and residuals from optimal matching, populate
     * outputMatched1 and outputMatched2;
     *
     * @param set1
     * @param set2
     * @param transTolXY
     * @param matchedIndexesAndDiffs two dimensional array holding the matched
     * indexes and the distances between the model and the point for that pair.
     * each row holds {idx1, idx2, diff}
     * @param outputMatched1 the container to hold the output matching points
     * for image 1 that are paired with outputMatched2 as a result of running
     * this method.
     * @param outputMatched2 the container to hold the output matching points
     * for image 2 that are paired with outputMatched1 as a result of running
     * this method.
     */
    public void matchPoints(PairFloatArray set1, PairIntArray set2,
        float transTolXY, float[][] matchedIndexesAndDiffs,
        PairFloatArray outputMatched1, PairIntArray outputMatched2) {

        for (int i = 0; i < matchedIndexesAndDiffs.length; i++) {

            int idx1 = (int)matchedIndexesAndDiffs[i][0];
            int idx2 = (int)matchedIndexesAndDiffs[i][1];
            float diff = matchedIndexesAndDiffs[i][2];

            if (diff < transTolXY) {
                outputMatched1.add(set1.getX(idx1), set1.getY(idx1));
                outputMatched2.add(set2.getX(idx2), set2.getY(idx2));
            }
        }
    }

    protected long estimateNStepsPreSearch(PairIntArray set1,
        PairIntArray set2, float rotRange) {

        int[] xyRange = estimateTranslationXYRange(set1, set2);

        long nTotal = estimateNStepsPreSearch0(set1.getN(), set2.getN(),
            xyRange[0], xyRange[1], rotRange);

        if (usePreSearchAlt1) {

            nTotal += estimateNStepsPreSearch1Alt1(set1.getN(), set2.getN());

        } else {

            nTotal += estimateNStepsPreSearch1Alt2(set1.getN(), set2.getN());
        }

        return nTotal;
    }

    protected int[] estimateTranslationXYRange(PairIntArray set1,
        PairIntArray set2) {

        int minX2 = MiscMath.findMin(set2.getX(), set2.getN());
        int maxX2 = MiscMath.findMax(set2.getX(), set2.getN());
        int minY2 = MiscMath.findMin(set2.getY(), set2.getN());
        int maxY2 = MiscMath.findMax(set2.getY(), set2.getN());

        int minX1 = MiscMath.findMin(set1.getX(), set1.getN());
        int maxX1 = MiscMath.findMax(set1.getX(), set1.getN());
        int minY1 = MiscMath.findMin(set1.getX(), set1.getN());
        int maxY1 = MiscMath.findMax(set1.getX(), set1.getN());

        int transXStart = minX1 - maxX2;
        int transXStop = maxX1 - minX2;
        int transYStart = minY1 - maxY2;
        int transYStop = maxY1 - minY2;

        int transXRange = transXStop - transXStart;
        int transYRange = transYStop - transYStart;

        return new int[]{transXRange, transYRange};
    }

    protected long estimateNStepsPreSearch0(int n1, int n2, int xRange,
        int yRange) {

        float rotRange = 360.f;

        return estimateNStepsPreSearch0(n1, n2, xRange, yRange, rotRange);
    }

    protected long estimateNStepsPreSearch0(int n1, int n2, int xRange,
        int yRange, float rotRange) {

        int nMaxMatchable = Math.min(n1, n2);

        int nRot = (int)((rotRange/rotDeltaPreSearch0) + 1);

        long nPerFitGreedy = n1 * n2;

        int nX = (int) Math.ceil(2*xRange / transDeltaPreSearch0);

        int nY = (int) Math.ceil(2*yRange / transDeltaPreSearch0);

        int nSplit = 1;

        if (nMaxMatchable > largeSearchLimit) {

            nSplit = (int)Math.ceil((float)n1/(float)largeSearchLimit);

            nPerFitGreedy = largeSearchLimit * largeSearchLimit;
        }

        return nSplit * nRot * nX * nY * nPerFitGreedy;
    }

    protected long estimateNStepsPreSearch1Alt1(int n1, int n2) {

        int nMaxMatchable = Math.min(n1, n2);

        int nRot = (int)((360/rotDeltaPreSearch0) + 1);

        long nPerFitGreedy = n1 * n2;

        int nSplit = 1;

        if (nMaxMatchable > largeSearchLimit) {

            nSplit = (int)Math.ceil((float)n1/(float)largeSearchLimit);

            nPerFitGreedy = largeSearchLimit * largeSearchLimit;
        }

        return nSplit * nPerFitGreedy * 50;
    }

    protected long estimateNStepsPreSearch1Alt2(int n1, int n2) {

        //preSearch0 gets answer to within +- 10 degrees of rotation and +- 20 ?

        int nFits = 10;

        int nMaxMatchable = Math.min(n1, n2);

        int nRot = (int)((2*rotHalfRangeInDegreesPreSearch1Alt2/
            rotDeltaInDegreesPreSearch1Alt2) + 1);

        int nX = (int) Math.ceil(2*transXHalfRangePreSearch1Alt2
            / transXDeltaPreSearch1Alt2);
        if (nX == 0) {
            nX = 1;
        }
        int nY = nX;

        long nPerFitGreedy = n1 * n2;

        int nSplit = 1;

        if (nMaxMatchable > largeSearchLimit) {

            nSplit = (int)Math.ceil((float)n1/(float)largeSearchLimit);

            nPerFitGreedy = largeSearchLimit * largeSearchLimit;

        }

        long nTot = nSplit * nFits * nRot * nX * nY * nPerFitGreedy;

        return nTot;
    }

    /**
     * given unordered unmatched points for a transformed set1 and
     * the model set2 from image1 and image2
     * respectively, find the
     * optimal matching in set2 within tolerance and return the
     * matched information as a two dimensional array of
     * {index from set1, index from set2, diff of point in set2 from
     * model generation by point in set1}
     *
     * @param transformed1
     * @param set2
     * @param toleranceX
     * @param toleranceY
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairFloatArray transformed1, PairIntArray set2, float toleranceX,
        float toleranceY) {

        int nPoints1 = transformed1.getN();
        int nPoints2 = set2.getN();

        float[][] diffsAsCost = new float[nPoints1][nPoints2];

        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];

        // key = indexes i,j; value = diffX, diffY
        Map<PairInt, PairFloat> diffsXY = new HashMap<PairInt, PairFloat>();

        int nWithinTol = 0;

        for (int i = 0; i < transformed1.getN(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];

            float x = transformed1.getX(i);
            float y = transformed1.getY(i);

            for (int j = 0; j < set2.getN(); j++) {

                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                float diffX = x - x2;
                float diffY = y - y2;

                if ((Math.abs(diffX) > toleranceX) ||
                    (Math.abs(diffY) > toleranceY)) {

                    diffsAsCost[i][j] = Float.MAX_VALUE;
                    diffsAsCostCopy[i][j] = Float.MAX_VALUE;

                } else {

                    double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;

                    diffsXY.put(new PairInt(i, j), new PairFloat(diffX, diffY));

                    nWithinTol++;
                }
            }
        }

        if (nWithinTol == 0) {
            return new float[0][];
        }

        boolean transposed = false;
        if (nPoints1 > nPoints2) {
            diffsAsCostCopy = MatrixUtil.transpose(diffsAsCostCopy);
            transposed = true;
        }

        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(diffsAsCostCopy);

        // count the number of matches
        int count = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (idx1 == Float.MAX_VALUE || idx2 == Float.MAX_VALUE) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            count++;
        }

        float[][] output = new float[count][];

        if (count == 0) {
            return output;
        }

        count = 0;

        for (int i = 0; i < match.length; i++) {

            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }

            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            output[count] = new float[3];
            output[count][0] = idx1;
            output[count][1] = idx2;
            output[count][2] = diffsAsCost[idx1][idx2];

            count++;
        }

        return output;
    }

    /**
     * given unordered unmatched points for a transformed set1 and
     * the model set2 from image1 and image2
     * respectively, find the
     * optimal matching in set2 within tolerance and return the
     * matched information as a two dimensional array of
     * {index from set1, index from set2, diff of point in set2 from
     * model generation by point in set1}
     *
     * @return a two dimensional array holding the matched indexes and
     * the distances between the model and the point for that pair.
     * each row holds float[]{idx1, idx2, diff}
     */
    public float[][] calculateMatchUsingOptimal(
        PairFloatArrayUnmodifiable scaledRotatedSet1,
        float transX, float transY, float toleranceX, float toleranceY,
        PairIntArray set2) {

        int nPoints1 = scaledRotatedSet1.getN();
        int nPoints2 = set2.getN();

        float[][] diffsAsCost = new float[nPoints1][nPoints2];

        // the algorithm modifies diffsAsCost, so make a copy
        float[][] diffsAsCostCopy = new float[nPoints1][nPoints2];

        // key = indexes i,j; value = diffX, diffY
        Map<PairInt, PairFloat> diffsXY = new HashMap<PairInt, PairFloat>();

        int nWithinTol = 0;

        for (int i = 0; i < scaledRotatedSet1.getN(); i++) {

            diffsAsCost[i] = new float[nPoints2];
            diffsAsCostCopy[i] = new float[nPoints2];

            float x = scaledRotatedSet1.getX(i) + transX;
            float y = scaledRotatedSet1.getY(i) + transY;

            for (int j = 0; j < set2.getN(); j++) {

                int x2 = set2.getX(j);
                int y2 = set2.getY(j);

                float diffX = x - x2;
                float diffY = y - y2;

                if ((Math.abs(diffX) > toleranceX) ||
                    (Math.abs(diffY) > toleranceY)) {

                    diffsAsCost[i][j] = Float.MAX_VALUE;
                    diffsAsCostCopy[i][j] = Float.MAX_VALUE;

                } else {

                    double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                    diffsAsCost[i][j] = (float)dist;
                    diffsAsCostCopy[i][j] = (float)dist;

                    diffsXY.put(new PairInt(i, j), new PairFloat(diffX, diffY));

                    nWithinTol++;
                }
            }
        }

        if (nWithinTol == 0) {
            return new float[0][];
        }

        boolean transposed = false;
        if (nPoints1 > nPoints2) {
            diffsAsCostCopy = MatrixUtil.transpose(diffsAsCostCopy);
            transposed = true;
        }

        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(diffsAsCostCopy);

        // count the number of matches
        int count = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (idx1 == Float.MAX_VALUE || idx2 == Float.MAX_VALUE) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            count++;
        }

        float[][] output = new float[count][];

        if (count == 0) {
            return output;
        }

        count = 0;

        for (int i = 0; i < match.length; i++) {

            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }

            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }

            PairFloat diffXY = diffsXY.get(new PairInt(idx1, idx2));

            if (diffXY == null) {
                continue;
            }

            if ((diffXY.getX() > toleranceX) || (diffXY.getY() > toleranceY)) {
                continue;
            }

            output[count] = new float[3];
            output[count][0] = idx1;
            output[count][1] = idx2;
            output[count][2] = diffsAsCost[idx1][idx2];

            count++;
        }

        return output;
    }

    protected TransformationPointFit[] evaluateTranslationsOverGrid(
        PairIntArray set1, PairIntArray set2,
        final float rotationInRadians, final float scale,
        RangeInt transXStartStop, int transXDelta,
        RangeInt transYStartStop, int transYDelta,
        float tolTransX, float tolTransY,
        final boolean setsAreMatched,
        final int numberOfBestToReturn) {

        /*   _____
            |     |
            |_____| largest negative or positive translationX of set1 is the width of set2
                  _____
                 |     |
                 |_____|
        */

        if (transXDelta < 1) {
            throw new IllegalArgumentException(
            "transXDelta must be greater than 0");
        }
        if (rotationInRadians < 0 || rotationInRadians > 359) {
            throw new IllegalArgumentException(
            "rotation must be between 0 and 359, inclusive");
        }
        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");

            return null;
        }

        int nMaxMatchable = (set1.getN() < set2.getN()) ?
            set1.getN() : set2.getN();

        if (nMaxMatchable == 0) {
            return null;
        }

        int nTranslations =
            (((transXStartStop.getStop() - transXStartStop.getStart())/transXDelta) + 1) *
            (((transYStartStop.getStop() - transYStartStop.getStart())/transYDelta) + 1);

        if (nTranslations == 0) {
            return null;
        }

        Transformer transformer = new Transformer();

        int count = 0;

        FixedSizeSortedVector<TransformationPointFit> fits =
            new FixedSizeSortedVector<>(numberOfBestToReturn,
                TransformationPointFit.class);

        for (float transX = transXStartStop.getStart();
            transX <= transXStartStop.getStop(); transX += transXDelta) {

            for (float transY = transYStartStop.getStart();
                transY <= transYStartStop.getStop(); transY += transYDelta) {

                float tx0 = transX + (transXDelta/2);
                float ty0 = transY + (transYDelta/2);

                TransformationParameters params = new TransformationParameters();
                params.setRotationInRadians(rotationInRadians);
                params.setScale(scale);
                params.setTranslationX(tx0);
                params.setTranslationY(ty0);

                PairFloatArray allPoints1Tr = transformer.applyTransformation2(
                    params, set1);

                TransformationPointFit fit;

                if (setsAreMatched) {

                    fit = evaluateFitForMatchedTransformed(params,
                        allPoints1Tr, set2);

                } else {

                    // default is to use greedy matching but use optimal for small sets
                    if (nMaxMatchable < 30) {

                        fit = evaluateFitForUnMatchedTransformedOptimal(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);

                    } else {

                        fit = evaluateFitForUnMatchedTransformedGreedy(params,
                        //fit = evaluateFitForUnMatchedTransformedOptimal(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);
                    }
                }

                fits.add(fit);

                count++;
            }
        }

        return fits.getArray();
    }

    /**
     * the maximum number of iterations that the refinement of translation
     * will use in the downhill simplex.  The default value is 50.
     */
    private int dsNMaxIter = 50;
    protected void setDsNMaxIter(int n) {
        dsNMaxIter = n;
    }
    protected int getDsNMaxIter() {
        return dsNMaxIter;
    }
    float nEpsFactor = 2.0f;
    protected void setNEpsFactor(float f) {
        nEpsFactor = f;
    }
    protected float getNEpsFactor() {
        return nEpsFactor;
    }
    protected float cellFactor = 1.25f;
    protected float tolFactor = 0.67f;
    protected void setCellFactor(float f) {
        cellFactor = f;
    }
    protected void setTolFactor(float f) {
        tolFactor = f;
    }
    protected float getCellFactor() {
        return cellFactor;
    }
    protected float getTolFactor() {
        return tolFactor;
    }

    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * Note that the reference point for the rotation is the center of the
     * image 1 width and height.
     *
     * It's expected that the invoker of this method is trying to solve for
     * translation for sets of points like corners in images.  This assumption
     * means that the number of point pair combinations is always far less
     * than the pixel combinations of translations over x and y.
     *
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     *
     * ALSO NOTE: if you know a better solution exists for translation
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be set to a smaller value and passed
     * to this method.
     *
     * @param matched1 set of points from image 1 to match to image2.
     * @param matched2 set of points from image 2 to be matched with image 1
     * @param rotationInRadians given in radians with value between 0 and 2*pi, exclusive
     * @param scale
     * @param image1Width
     * @param image1Height
     * @param image2Width
     * @param image2Height
     * @return
     */
    public TransformationParameters calculateTranslationForMatched(
        PairIntArray matched1, PairIntArray matched2,
        float rotationInRadians, float scale,
        int image1Width, int image1Height, int image2Width, int image2Height) {

        if (scale < 1) {
            // numerical errors in rounding to integer can give wrong solutions
            //throw new IllegalStateException("scale cannot be smaller than 1");

            log.severe("scale cannot be smaller than 1");

            return null;
        }

        int image1CentroidX = image1Width >> 1;
        int image1CentroidY = image1Height >> 1;

        double scaleTimesCosine = scale * Math.cos(rotationInRadians);
        double scaleTimesSine = scale * Math.sin(rotationInRadians);

        double avgTransX = 0;
        double avgTransY = 0;

        for (int i = 0; i < matched1.getN(); i++) {

            int x = matched1.getX(i);
            int y = matched1.getY(i);

            double xr = image1CentroidX*scale + (
                ((x - image1CentroidX) * scaleTimesCosine) +
                ((y - image1CentroidY) * scaleTimesSine));

            double yr = image1CentroidY*scale + (
                (-(x - image1CentroidX) * scaleTimesSine) +
                ((y - image1CentroidY) * scaleTimesCosine));

            int x2 = matched2.getX(i);

            int y2 = matched2.getY(i);

            avgTransX += (int)Math.round(x2 - xr);

            avgTransY += (int)Math.round(y2 - yr);
        }

        avgTransX /= (float)matched1.getN();

        avgTransY /= (float)matched1.getN();

        TransformationParameters params =
            new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX((float)avgTransX);
        params.setTranslationY((float)avgTransY);

        return params;
    }

    /**
     * given the transformed x y that have already been scaled and rotated, add the
     * transX and transY, respectively and calculated the average residual
     * between that and set2 and the standard deviation from the average.
     * Note that set2 and (scaledRotatedX, scaledRotatedY) are NOT known to be
     * matched points so the residuals are minimized for each point in
     * the model to find the matching in set2 before computing the
     * average and standard deviation.
     *
     * @param set2
     * @param scaledRotatedSet1 the model xy points scaled and rotated
     * @param transX the x translation to apply to the model points
     * @param transY the y translation to apply to the model points
     * @return
     */
    protected TransformationPointFit evaluateFitForUnMatchedGreedy(
        PairFloatArrayUnmodifiable scaledRotatedSet1,
        float transX, float transY,
        float tolTransX, float tolTransY,
        PairIntArray set2,
        final float scale, final float rotationRadians) {

        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (scaledRotatedSet1 == null) {
            throw new IllegalArgumentException(
            "scaledRotatedSet1 cannot be null");
        }

        //TODO: consider filtering for the intersection

        int nMaxMatchable = (scaledRotatedSet1.getN() < set2.getN()) ?
            scaledRotatedSet1.getN() : set2.getN();

        Set<Integer> chosen = new HashSet<Integer>();

        double[] diffs = new double[scaledRotatedSet1.getN()];
        int nMatched = 0;
        double avg = 0;

        for (int i = 0; i < scaledRotatedSet1.getN(); i++) {

            float transformedX = scaledRotatedSet1.getX(i) + transX;
            float transformedY = scaledRotatedSet1.getY(i) + transY;

            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < set2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }

                float dx = set2.getX(j) - transformedX;
                float dy = set2.getY(j) - transformedY;

                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {

                    continue;
                }

                float diff = (float)Math.sqrt(dx*dx + dy*dy);

                if (diff < minDiff) {
                    minDiff = diff;
                    min2Idx = j;
                }
            }

            if (minDiff < Double.MAX_VALUE) {
                diffs[nMatched] = minDiff;
                nMatched++;
                chosen.add(Integer.valueOf(min2Idx));
                avg += minDiff;
            }
        }

        avg = (nMatched == 0) ? Double.MAX_VALUE :
            avg / (double)nMatched;

        double stDev = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avg;
            stDev += (d * d);
        }

        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            Math.sqrt(stDev/((double)nMatched - 1.));

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationRadians);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        params.setOriginX(0);
        params.setOriginY(0);

        TransformationPointFit fit = new TransformationPointFit(params, nMatched,
            avg, stDev, tolTransX, tolTransY);

        fit.setMaximumNumberMatchable(nMaxMatchable);

        return fit;
    }

    @SuppressWarnings("fallthrough")
    TransformationPointFit[] refineTransformationWithDownhillSimplex(
        TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2,
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        boolean searchScaleToo,
        boolean useGreedyMatching,
        boolean setsAreMatched, int nMaxIter) {

        int nMaxMatchable = (set1.getN() < set2.getN()) ?
            set1.getN() : set2.getN();

        if (nMaxMatchable == 0) {
            return null;
        }

        //TODO: revise this.  consider a sqrt(nMaxMatchable)
        double eps = Math.log(nMaxMatchable)/Math.log(10);

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.5f;
        float tau = 0.5f;

        // alternates changes in rotation, translation of X and translation of Y

        boolean go = true;

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        float rotInRadians = (fits[0] != null) ? fits[0].getRotationInRadians()
            : 0;
        float scale = (fits[0] != null) ? fits[0].getScale() : 1;

        int nIter = 0;
        int lastAlt = 0;

        int bestFitIdx = 0;
        int worstFitIdx = fits.length - 1;

        TransformationParameters[] lastParams = extractParameters(fits);

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            if (fits[bestFitIdx] == null) {
                break;
            }

            for (int i = (fits.length - 1); i > -1; --i) {
                if (fits[i] != null) {
                    worstFitIdx = i;
                    break;
                }
            }

            if (debug) {
                if ((nIter == 0) || (nIter % 5 == 0)) {
                    writeSimplex(fits, false);
                }
            }

            if (nIter > 0) {

                TransformationParameters[] currentParams = extractParameters(fits);

                boolean areTheSame = areEqual(lastParams, currentParams);

                if (areTheSame) {
                    break;
                }

                if (nIter > 50) {
                    double[] stdevs = getStandardDeviation(currentParams);
                    if (stdevs != null && (stdevs.length > 0)) {
                        boolean areNearlyZero = true;
                        for (double stdv : stdevs) {
                            if (Math.abs(stdv) > 0.01) {
                                areNearlyZero = false;
                            }
                        }
                        if (areNearlyZero) {
                            break;
                        }
                    }
                }

                lastParams = currentParams;
            }

            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            double rotSum = 0.0f;
            double scaleSum = 0.0f;
            int c = 0;
            for (int i = 0; i < (fits.length - 1); ++i) {
                if (fits[i] != null) {
                    txSum += fits[i].getTranslationX();
                    tySum += fits[i].getTranslationY();
                    rotSum += (fits[i].getRotationInRadians()/(2.*Math.PI));
                    scaleSum += (fits[i].getScale());
                    c++;
                }
            }
            transX = txSum / (float)c;
            transY = tySum / (float)c;
            rotInRadians = (float)(rotSum/(float)c)*((float)(2.*Math.PI));
            scale = (float)(scaleSum / (float)c);

            // "Reflection"

            float txReflect = 0;
            float tyReflect = 0;
            float rotReflect = 0;
            float scaleReflect = 1;

            TransformationParameters paramsReflect = new TransformationParameters();
            if (fits[worstFitIdx] != null) {
                paramsReflect.setScale(fits[worstFitIdx].getScale());
                paramsReflect.setRotationInRadians(fits[worstFitIdx].getRotationInRadians());
                paramsReflect.setTranslationX(fits[worstFitIdx].getTranslationX());
                paramsReflect.setTranslationY(fits[worstFitIdx].getTranslationY());
                paramsReflect.setOriginX(fits[worstFitIdx].getParameters().getOriginX());
                paramsReflect.setOriginY(fits[worstFitIdx].getParameters().getOriginY());
            }

            switch(lastAlt) {
                case 0: {
                    // change translation X
                    txReflect = transX + (alpha
                        * (transX - fits[worstFitIdx].getTranslationX()));
                    paramsReflect.setTranslationX(txReflect);
                    break;
                }
                case 1: {
                    // change translation Y
                    tyReflect = transY + (alpha
                        * (transY - fits[worstFitIdx].getTranslationY()));
                    paramsReflect.setTranslationY(tyReflect);
                    break;
                }
                case 2: {
                    if (searchScaleToo) {
                        scaleReflect  = scale + (alpha
                            * (scale - fits[worstFitIdx].getScale()));
                        paramsReflect.setScale(scaleReflect);
                        break;
                    }
                    // else fall through to rotation
                }
                default: {
                    // change rotation
                    rotReflect = rotInRadians + (alpha
                        * (rotInRadians - fits[worstFitIdx].getRotationInRadians()));
                    if (rotReflect < 0) {
                        rotReflect = rotReflect + (float)(2.*Math.PI);
                    } else if (rotReflect >= (2.*Math.PI)) {
                        rotReflect = rotReflect - (float)(2.*Math.PI);
                    }
                    paramsReflect.setRotationInRadians(rotReflect);
                    break;
                }
            }

            TransformationPointFit fitReflected = (fits[worstFitIdx] != null) ?
                evaluateForUnmatched(paramsReflect, set1, set2,
                tolTransX, tolTransY, useGreedyMatching) : null;

            int comp0 = compare(fits[bestFitIdx], fitReflected);
            int compLast = compare(fits[worstFitIdx], fitReflected);

            if ((comp0 < 1) && (compLast == 1)) {

                // replace last with f_refl
                fits[worstFitIdx] = fitReflected;

            } else if (comp0 == 1) {

                // reflected is better than best fit, so "expand"
                // "Expansion"

                TransformationParameters paramsExpansion = new TransformationParameters();
                if (fitReflected != null) {
                    paramsExpansion.setScale(fitReflected.getScale());
                    paramsExpansion.setRotationInRadians(fitReflected.getRotationInRadians());
                    paramsExpansion.setTranslationX(fitReflected.getTranslationX());
                    paramsExpansion.setTranslationY(fitReflected.getTranslationY());
                    paramsExpansion.setOriginX(fitReflected.getParameters().getOriginX());
                    paramsExpansion.setOriginY(fitReflected.getParameters().getOriginY());
                }

                switch (lastAlt) {
                    case 0: {
                        // change translation X
                        float txExpansion = transX + (gamma * (txReflect - transX));
                        paramsExpansion.setTranslationX(txExpansion);
                        break;
                    }
                    case 1: {
                        // change translation Y
                        float tyExpansion = transY + (gamma * (tyReflect - transY));
                        paramsExpansion.setTranslationY(tyExpansion);
                        break;
                    }
                    case 2: {
                        if (searchScaleToo) {
                            float scaleExpansion = scale + (gamma * (scaleReflect - scale));
                            paramsExpansion.setScale(scaleExpansion);
                            break;
                        }
                        // else fall through to rotation
                    }
                    default: {
                        // change rotation
                        float rotExpansion = rotInRadians + (gamma * (rotReflect - rotInRadians));
                        if (rotExpansion < 0) {
                            rotExpansion = rotExpansion + (float)(2.*Math.PI);
                        } else if (rotExpansion >= (2.*Math.PI)) {
                            rotExpansion = rotExpansion - (float)(2.*Math.PI);
                        }
                        paramsExpansion.setRotationInRadians(rotExpansion);
                        break;
                    }
                }

                TransformationPointFit fitExpansion = (fitReflected != null) ?
                    evaluateForUnmatched(paramsExpansion, set1, set2,
                    tolTransX, tolTransY, useGreedyMatching) : null;

                int compR = compare(fitReflected, fitExpansion);

                if (compR == 1) {

                    // expansion fit is better than reflected fit
                    fits[worstFitIdx] = fitExpansion;

                } else {

                    fits[worstFitIdx] = fitReflected;
                }

            } else if (compLast < 1) {

                // reflected fit is worse than the worst (last) fit, so contract
                // "Contraction"

                TransformationParameters paramsContraction = new TransformationParameters();

                paramsContraction.setScale(fits[worstFitIdx].getScale());
                paramsContraction.setRotationInRadians(fits[worstFitIdx].getRotationInRadians());
                paramsContraction.setTranslationX(fits[worstFitIdx].getTranslationX());
                paramsContraction.setTranslationY(fits[worstFitIdx].getTranslationY());
                paramsContraction.setOriginX(fits[worstFitIdx].getParameters().getOriginX());
                paramsContraction.setOriginY(fits[worstFitIdx].getParameters().getOriginY());

                switch (lastAlt) {
                    case 0: {
                        // change translation X
                        float txContraction = transX + (beta
                            * (fits[worstFitIdx].getTranslationX() - transX));
                        paramsContraction.setTranslationX(txContraction);
                        break;
                    }
                    case 1: {
                        // change translation Y
                        float tyContraction = transY + (beta
                            * (fits[worstFitIdx].getTranslationY() - transY));
                        paramsContraction.setTranslationY(tyContraction);
                        break;
                    }
                    case 2: {
                        if (searchScaleToo) {
                            float scaleContraction = scale + (beta
                                * (fits[worstFitIdx].getScale() - scale));
                            paramsContraction.setScale(scaleContraction);
                            break;
                        }
                        // else fall through to rotation
                    }
                    default: {
                        // change rotation
                        float rotContraction = rotInRadians + (beta
                            * (fits[worstFitIdx].getRotationInRadians() - rotInRadians));
                        if (rotContraction < 0) {
                            rotContraction = rotContraction + (float)(2.*Math.PI);
                        } else if (rotContraction >= (2.*Math.PI)) {
                            rotContraction = rotContraction - (float)(2.*Math.PI);
                        }
                        paramsContraction.setRotationInRadians(rotContraction);
                        break;
                    }
                }

                TransformationPointFit fitContraction = (fits[worstFitIdx] != null) ?
                    evaluateForUnmatched(paramsContraction, set1, set2,
                    tolTransX, tolTransY, useGreedyMatching) : null;

                int compC = compare(fits[worstFitIdx], fitContraction);

                if (compC > -1) {

                    fits[worstFitIdx] = fitContraction;

                } else {

                    // "Reduction"
                    for (int i = 1; i < fits.length; ++i) {

                        if (fits[i] == null || fits[bestFitIdx] == null) {
                            /*TODO: consider setting this
                            fits[i] = new TransformationPointFit(
                                new TransformationParameters(),
                                0, Double.MAX_VALUE, Double.MAX_VALUE,
                                Double.MAX_VALUE
                            );
                            */
                            continue;
                        }

                        TransformationParameters paramsI = new TransformationParameters();
                        paramsI.setScale(fits[i].getScale());
                        paramsI.setRotationInRadians(fits[i].getRotationInRadians());
                        paramsI.setTranslationX(fits[i].getTranslationX());
                        paramsI.setTranslationY(fits[i].getTranslationY());
                        paramsI.setOriginX(fits[i].getParameters().getOriginX());
                        paramsI.setOriginY(fits[i].getParameters().getOriginY());

                        switch (lastAlt) {
                            case 0: {
                                // change translation X
                                float txReduction
                                    = (fits[bestFitIdx].getTranslationX()
                                    + (tau * (fits[i].getTranslationX()
                                    - fits[bestFitIdx].getTranslationX())));
                                paramsI.setTranslationX(txReduction);
                                break;
                            }
                            case 1: {
                                // change translation Y
                                float tyReduction
                                    = (fits[bestFitIdx].getTranslationY()
                                    + (tau * (fits[i].getTranslationY()
                                    - fits[bestFitIdx].getTranslationY())));
                                paramsI.setTranslationY(tyReduction);
                                break;
                            }
                            case 2: {
                                if (searchScaleToo) {
                                    float scaleReduction
                                        = (fits[bestFitIdx].getScale()
                                        + (tau * (fits[i].getScale()
                                        - fits[bestFitIdx].getScale())));
                                    paramsI.setScale(scaleReduction);
                                    break;
                                }
                                // else fall through to rotation
                            }
                            default: {
                                // change rotation
                                float rotReduction
                                    = (fits[bestFitIdx].getRotationInRadians()
                                    + (tau * (fits[i].getRotationInRadians()
                                    - fits[bestFitIdx].getRotationInRadians())));
                                if (rotReduction < 0) {
                                    rotReduction = rotReduction + (float)(2.*Math.PI);
                                } else if (rotReduction >= (2.*Math.PI)) {
                                    rotReduction = rotReduction - (float)(2.*Math.PI);
                                }
                                paramsI.setRotationInRadians(rotReduction);
                                break;
                            }
                        }
                        fits[i] = evaluateForUnmatched(paramsI, set1, set2,
                            tolTransX, tolTransY, useGreedyMatching);
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;
            lastAlt++;
            if (searchScaleToo && (lastAlt > 3)) {
                lastAlt = 0;
            } else if (!searchScaleToo && (lastAlt > 2)) {
                lastAlt = 0;
            }

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == nMaxMatchable)
                && (fits[bestFitIdx].getMeanDistFromModel() < eps)) {
                go = false;
            } else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }
        }

        log.fine("nIter=" + Integer.toString(nIter));

        // additional step that's helpful if not enough iterations are used,
        // is to test the summed transX, transY which represent the center
        // of the simplex against the best fit
        if (fits[bestFitIdx] != null) {

            TransformationParameters p = new TransformationParameters();
            p.setScale(scale);
            p.setRotationInRadians(rotInRadians);
            p.setTranslationX(transX);
            p.setTranslationY(transY);
            p.setOriginX(fits[bestFitIdx].getParameters().getOriginX());
            p.setOriginY(fits[bestFitIdx].getParameters().getOriginY());

            TransformationPointFit fitAvg = evaluateForUnmatched(p, set1, set2,
                tolTransX, tolTransY, useGreedyMatching);

            int comp = compare(fits[bestFitIdx], fitAvg);
            if (comp == 1) {
                fits[bestFitIdx] = fitAvg;
            }
        }

        if (debug) {
            // the resulting file can be viewed with resources/plot_3d_simplex_n10.html
            writeSimplex(fits, true);
        }

        return fits;
    }

    public TransformationPointFit transformAndEvaluateFit(
        PairIntArray set1, PairIntArray set2,
        TransformationParameters params, float tolTransX, float tolTransY,
        boolean useGreedyMatching, boolean setsAreMatched) {

        if (set2 == null) {
            throw new IllegalArgumentException(
            "set2 cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException(
            "set1 cannot be null");
        }

        Transformer transformer = new Transformer();

        PairFloatArray transformedSet1 = transformer.applyTransformation2(
            params, set1);

        if (setsAreMatched) {

            return evaluateFitForMatchedTransformed(params, transformedSet1,
                set2);

        } else {

            return evaluateFitForUnMatchedTransformedGreedy(
                params, transformedSet1, set2, tolTransX, tolTransY);
        }
    }

    private TransformationParameters[] extractParameters(
        TransformationPointFit[] fits) {

        if (fits == null) {
            return new TransformationParameters[0];
        }

        TransformationParameters[] params = new TransformationParameters[fits.length];

        for (int i = 0; i < fits.length; ++i) {
            if (fits[i] != null) {
                params[i] = fits[i].getParameters();
            }
        }

        return params;
    }

    private double[] getStandardDeviation(TransformationParameters[] params) {

        if (params == null || params.length == 0) {
            return null;
        }

        double rSum = 0;
        double sSum = 0;
        double tXSum = 0;
        double tYSum = 0;

        int count = 0;
        for (int i = 0; i < params.length; ++i) {

            TransformationParameters p0 = params[i];

            if (p0 == null) {
                continue;
            }

            rSum += p0.getRotationInRadians();
            sSum += p0.getScale();
            tXSum += p0.getTranslationX();
            tYSum += p0.getTranslationY();
            count++;
        }
        rSum /= (double)count;
        sSum /= (double)count;
        tXSum /= (double)count;
        tYSum /= (double)count;

        double rStdv = 0;
        double sStdv= 0;
        double tXStdv = 0;
        double tYStdv = 0;
        for (int i = 0; i < params.length; ++i) {

            TransformationParameters p0 = params[i];

            if (p0 == null) {
                continue;
            }

            double r = p0.getRotationInRadians() - rSum;
            double s = p0.getScale() - sSum;
            double tX = p0.getTranslationX() - tXSum;
            double tY = p0.getTranslationY() - tYSum;

            rStdv += (r*r);
            sStdv += (s*s);
            tXStdv += (tX*tX);
            tYStdv += (tY*tY);
        }
        rStdv = (Math.sqrt(rStdv/(count - 1.0f)));
        sStdv = (Math.sqrt(sStdv/(count - 1.0f)));
        tXStdv = (Math.sqrt(tXStdv/(count - 1.0f)));
        tYStdv = (Math.sqrt(tYStdv/(count - 1.0f)));

        return new double[]{rStdv, sStdv, tXStdv, tYStdv};
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @return
     */
    private TransformationPointFit reevaluateForNewTolerance(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2) {

        TransformationPointFit fit2 = reevaluateForNewTolerance(
            fit, translationXTolerance, translationYTolerance,
            set1, set2, true);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @return
     */
    private TransformationPointFit reevaluateForNewToleranceOptimal(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2) {

        TransformationPointFit fit2 = reevaluateForNewTolerance(
            fit, translationXTolerance, translationYTolerance,
            set1, set2, false);

        return fit2;
    }

    /**
     * Re-evaluate the fit of the enclosed parameters, but use the new tolerance
     * for translation in x and y.
     *
     * @param fit
     * @param set1
     * @param set2
     * @return
     */
    private TransformationPointFit reevaluateForNewTolerance(
        TransformationPointFit fit, float translationXTolerance,
        float translationYTolerance, PairIntArray set1, PairIntArray set2,
        boolean useGreedyMatching) {

        if (fit == null) {
            return null;
        }

        float rotationInRadians = fit.getRotationInRadians();

        float scale = fit.getScale();

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationInRadians);
        params.setScale(scale);
        params.setTranslationX(fit.getTranslationX());
        params.setTranslationY(fit.getTranslationY());

        TransformationPointFit fit2 = evaluateForUnmatched(params,
            set1, set2, translationXTolerance, translationYTolerance,
            useGreedyMatching);

        return fit2;
    }

    private TransformationPointFit evaluateForUnmatched(
        final TransformationParameters params,
        PairIntArray set1, PairIntArray set2,
        float translationXTolerance, float translationYTolerance,
        boolean useGreedyMatching) {

        if (params == null) {
            return null;
        }

        boolean setsAreMatched = false;

        return transformAndEvaluateFit(set1, set2, params,
            translationXTolerance, translationYTolerance,
            useGreedyMatching, setsAreMatched);
    }

    /**
     * find the minima and maxima of translationX and translationY for the given
     * fits.
     * @param fits
     * @return new float[]{minTX, maxTX, minTY, maxTY}
     */
    protected float[] getTranslationMinAndMaxes(TransformationPointFit[] fits) {

        if (fits == null || fits.length == 0) {
            return null;
        }
        float minTX = Float.MAX_VALUE;
        float maxTX = Float.MIN_VALUE;
        float minTY = Float.MAX_VALUE;
        float maxTY = Float.MIN_VALUE;

        for (TransformationPointFit fit : fits) {
            if (fit == null) {
                continue;
            }
            float tx = fit.getTranslationX();
            float ty = fit.getTranslationY();
            if (tx < minTX) {
                minTX = tx;
            }
            if (ty < minTY) {
                minTY = ty;
            }
            if (tx > maxTX) {
                maxTX = tx;
            }
            if (ty > maxTY) {
                maxTY = ty;
            }
        }

        return new float[]{minTX, maxTX, minTY, maxTY};
    }

    /**
     * re-evaluate bestFit or fit, whichever has the largest translation
     * tolerance, using the smaller tolerance.  Note that there are
     * exception rules, such as when both bestFit and fit have
     * same number of points, but fit has a mean dist from model less than one
     * and bestFit has a much larger mean distance from model.  In that case,
     * even if bestFit has a smaller translation tolerance, fit will not
     * be re-evaluated because such a mean distance from model means the
     * answer has converged.
     *
     *
     * @param bestFit
     * @param fit
     * @param set1
     * @param set2
     * @param reevalFits
     * @param fitIsBetter
     */
    protected void reevaluateFitsForCommonTolerance(
        TransformationPointFit bestFit, TransformationPointFit fit,
        PairIntArray set1, PairIntArray set2,
        final TransformationPointFit[] reevalFits, final boolean[] fitIsBetter) {

        if (bestFit == null) {
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            fitIsBetter[0] = true;
            return;
        } else if (fit == null) {
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            fitIsBetter[0] = false;
            return;
        }

        if ((bestFit.getTranslationXTolerance() == fit.getTranslationXTolerance()) &&
            (bestFit.getTranslationYTolerance() == fit.getTranslationYTolerance())) {
            fitIsBetter[0] = fitIsBetter(bestFit, fit);
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            return;
        }

        /*
        check for whether fit has converged already for equal number of points
        matched.
        */
        int bestNMatches = bestFit.getNumberOfMatchedPoints();
        int compNMatches = fit.getNumberOfMatchedPoints();
        int diffEps = (int)Math.round(2.*Math.ceil(Math.max(bestNMatches,
            compNMatches)/10.));
        if (diffEps == 0) {
            diffEps = 1;
        }

        double compAvg = fit.getMeanDistFromModel();
        double bestAvg = bestFit.getMeanDistFromModel();

        double compS = fit.getStDevFromMean();
        double bestS = bestFit.getStDevFromMean();

        double aDiv = bestAvg/compAvg;

        // TODO:  this skips the equal tolerance re-eval if solutions are already
        //    alot better and have significant number of points matched
        if ((compNMatches > 10) && (bestNMatches > 10)) {
            float nDiv = (float)bestNMatches/(float)compNMatches;
            double sDiv = bestS/compS;
            double nLimit0 = 2;
            double nLimit1 = 0.5;
            if (bestNMatches > 50) {
                nLimit0 = 1.5;
            }
            if (compNMatches > 50) {
                nLimit1 = 1./1.5;
            }
            if ((nDiv > nLimit0) && (aDiv < 0.5) && (sDiv < 0.5)) {
                fitIsBetter[0] = false;
                reevalFits[0] = bestFit;
                reevalFits[1] = fit;
                return;
            } else if ((nDiv < nLimit1) && (aDiv > 2) && (sDiv > 2)) {
                fitIsBetter[0] = true;
                reevalFits[0] = bestFit;
                reevalFits[1] = fit;
                return;
            }
        }

        /*
        when tolerances are both already very small, not redoing the fit,
        just comparing as is
        */
        int limit = 7;
        if ((bestFit.getTranslationXTolerance() < limit) &&
            (bestFit.getTranslationYTolerance() < limit) &&
            (fit.getTranslationXTolerance() < limit) &&
            (fit.getTranslationYTolerance() < limit)) {

            fitIsBetter[0] = fitIsBetter(bestFit, fit);
            reevalFits[0] = bestFit;
            reevalFits[1] = fit;
            return;
        }

        /*
        -1 : both are not null and bestFit tolerances are smaller
         0 : both are not null and tolerances are same.
         1 : both are not null and fit tolerances are smaller
         2 : both are not null and the x and y fits and smaller and larger in a mix
         3 : either bestFit or fit is null
        */
        int compTol = compareTolerance(bestFit, fit);

        TransformationPointFit bestFitT = bestFit;
        TransformationPointFit fitT = fit;

        float tolX, tolY;
        if (compTol != 3) {
            if (compTol == 1) {
                tolX = fit.getTranslationXTolerance();
                tolY = fit.getTranslationYTolerance();
            } else if (compTol == -1) {
                tolX = bestFit.getTranslationXTolerance();
                tolY = bestFit.getTranslationYTolerance();
            } else {
                tolX = bestFit.getTranslationXTolerance();
                if (fit.getTranslationXTolerance() < tolX) {
                    tolX = fit.getTranslationXTolerance();
                }
                tolY = bestFit.getTranslationYTolerance();
                if (fit.getTranslationYTolerance() < tolY) {
                    tolY = fit.getTranslationYTolerance();
                }
            }
            /*
            for small tolerances, both need to be re-evaluated with optimal
            matching, and then compared
            but the results should not be saved because that would result
            in an inconsistent mix of comparisons with greedy and optimal
            matching outside of this method
            */
            int limitT = 11;
            if ((tolX < limitT) && (tolY < limitT) && (compTol != 0)) {

                bestFitT = reevaluateForNewToleranceOptimal(bestFit,
                    tolX, tolY,
                    set1, set2);

                fitT = reevaluateForNewToleranceOptimal(fit,
                    tolX, tolY,
                    set1, set2);

                float nBDiv = (float)bestFit.getNumberOfMatchedPoints()/
                    (float)bestFitT.getNumberOfMatchedPoints();

                float nFDiv = (float)fit.getNumberOfMatchedPoints()/
                    (float)fitT.getNumberOfMatchedPoints();

                if (
                    (fitT == null) || (fitT.getNumberOfMatchedPoints() < 2) ||
                    (bestFitT == null) ||
                    (bestFitT.getNumberOfMatchedPoints() < 2)
                    //|| (nBDiv > 4) || (nFDiv > 4)
                ) {
                    fitIsBetter[0] = fitIsBetter(bestFit, fit);
                    reevalFits[0] = bestFit;
                    reevalFits[1] = fit;
                } else {
                    fitIsBetter[0] = fitIsBetter(bestFitT, fitT);
                    reevalFits[0] = bestFit;
                    reevalFits[1] = fit;
                }
                return;
            }

            if (compTol == 1) {

                bestFitT = reevaluateForNewTolerance(bestFit, tolX, tolY, set1,
                    set2);

                //TODO: may need to discard only if smaller tolerances < 10
                // do not use the lower tolerance.  it may be a false fit if null here
                if (bestFitT == null) {
                    bestFitT = bestFit;
                } else if (bestFitT.getNumberOfMatchedPoints() == 0) {
                    if (bestFit.getNumberOfMatchedPoints() > 10) {
                        bestFitT = bestFit;
                    }
                }
            } else if (compTol == -1) {

                fitT = reevaluateForNewTolerance(fit, tolX, tolY, set1, set2);

                //TODO: may need to discard only if smaller tolerances < 10
                // do not use the lower tolerance if resulted in null fit
                if (fitT == null) {
                    fitT = fit;
                } else if (fitT.getNumberOfMatchedPoints() == 0) {
                    if (fit.getNumberOfMatchedPoints() > 10) {
                        fitT = fit;
                    }
                } else if (true) {
                    // need to try optimal matching for both.
                    // TODO: may be what the entire method should use for simplicity.
                    if (
                        (fit.getNumberOfMatchedPoints() > 10) &&
                        (((float)fit.getNumberOfMatchedPoints()/(float)fitT.getNumberOfMatchedPoints())
                        > 4)) {

                        bestFitT = reevaluateForNewToleranceOptimal(bestFit,
                            tolX, tolY,
                            set1, set2);

                        fitT = reevaluateForNewToleranceOptimal(fit,
                            tolX, tolY,
                            set1, set2);

                        float nBDiv = (float) bestFit.getNumberOfMatchedPoints()
                            / (float) bestFitT.getNumberOfMatchedPoints();

                        float nFDiv = (float) fit.getNumberOfMatchedPoints()
                            / (float) fitT.getNumberOfMatchedPoints();

                        if ((fitT == null) || (fitT.getNumberOfMatchedPoints() < 2)
                            || (bestFitT == null)
                            || (bestFitT.getNumberOfMatchedPoints() < 2) //|| (nBDiv > 4) || (nFDiv > 4)
                            ) {
                            fitIsBetter[0] = fitIsBetter(bestFit, fit);
                            reevalFits[0] = bestFit;
                            reevalFits[1] = fit;
                        } else {
                            fitIsBetter[0] = fitIsBetter(bestFitT, fitT);
                            reevalFits[0] = bestFit;
                            reevalFits[1] = fit;
                        }
                        return;
                    }
                }

            } else if (compTol == 2) {

                // TODO: may need to revise this
                // reduce both to smallest tolerances
                bestFitT = reevaluateForNewTolerance(bestFit,
                    tolX, tolY, set1, set2);

                fitT = reevaluateForNewTolerance(fit,
                    tolX, tolY, set1, set2);
            }
        }

        fitIsBetter[0] = fitIsBetter(bestFitT, fitT);

        reevalFits[0] = bestFit;

        if (fitIsBetter[0] && (fitT != null)) {
           reevalFits[1] = fitT;
        }

if (bestFit != null && fit != null) {
if (compTol == 1) {
    log.fine("    rot re-evaluated bestFit at lower tolerance");
} else if (compTol == 2) {
    log.fine("    rot re-evaluated bestFit and fit at common tolerance");
}
}

    }

    private PolygonPlotterPNG plotter = null;
    //TODO: put this in an aspect
    private void plotTranslationSimplex(TransformationPointFit[] fits,
        float minX, float maxX, float minY, float maxY, java.awt.Color clr) {

        try {
            if (plotter == null) {
                plotter = new PolygonPlotterPNG(minX, maxX, minY, maxY,
                    "translation simplex", "transX", "transY");
            }
            int count = 0;
            for (TransformationPointFit fit : fits) {
                if (fit == null) {
                    continue;
                }
                count++;
            }
            double[] x = new double[count];
            double[] y = new double[x.length];
            count = 0;
            for (TransformationPointFit fit : fits) {
                if (fit == null) {
                    continue;
                }
                x[count] = fit.getTranslationX();
                y[count] = fit.getTranslationY();
                count++;
            }

            if (clr == null) {
                plotter.addPolygon(x, y);
            } else {
                plotter.addPolygon(x, y, clr);
            }

        } catch (IOException e) {
            log.severe(e.getMessage());
        }
    }

    private FileWriter simplexCSVFileW = null;
    private BufferedWriter simplexCSVFileWriter = null;
    //TODO: put this in an aspect
    private void writeSimplex(TransformationPointFit[] fits, boolean closeFile) {

        //need a convention for writing for nulls.  will use max of rot, transX,
        //and transY when a fit is null
        float maxRot = Float.MIN_VALUE;
        float maxTransX = Float.MIN_VALUE;
        float maxTransY = Float.MIN_VALUE;
        for (int i = 0; i < fits.length; ++i) {
            TransformationPointFit fit = fits[i];
            if (fit != null) {
                if (fit.getParameters().getRotationInDegrees() > maxRot) {
                    maxRot = fit.getParameters().getRotationInDegrees();
                }
                if (fit.getParameters().getTranslationX() > maxTransX) {
                    maxTransX = fit.getParameters().getTranslationX();
                }
                if (fit.getParameters().getTranslationY() > maxTransY) {
                    maxTransY = fit.getParameters().getTranslationY();
                }
            }
        }

        try {
            if (simplexCSVFileWriter == null) {
                String testDirPath = ResourceFinder.findOutputTestDirectory();
                String fileName = "simplex_" + MiscDebug.getCurrentTimeFormatted() + ".csv";
                String filePath = testDirPath
                    + System.getProperty("file.separator") + fileName;
                File file = new File(filePath);
                simplexCSVFileW = new FileWriter(file);
                simplexCSVFileWriter = new BufferedWriter(simplexCSVFileW);
            }

            for (TransformationPointFit fit : fits) {

                StringBuilder sb = new StringBuilder();

                float rot, tx, ty;

                if (fit == null) {
                    rot = maxRot;
                    tx = maxTransX;
                    ty = maxTransY;
                } else {
                    rot = fit.getParameters().getRotationInDegrees();
                    tx = fit.getTranslationX();
                    ty = fit.getTranslationY();
                }
                sb.append(Float.toString(rot)).append(",")
                    .append(Float.toString(tx)).append(",")
                    .append(Float.toString(ty)).append("\n");

                String str = sb.toString();

                simplexCSVFileWriter.write(str, 0, str.length());
            }
            simplexCSVFileWriter.flush();

            if (closeFile) {
                if (simplexCSVFileW != null) {
                    simplexCSVFileW.close();
                    simplexCSVFileW = null;
                }
                if (simplexCSVFileWriter != null) {
                    simplexCSVFileWriter.close();
                    simplexCSVFileWriter = null;
                }
            }

        } catch (IOException e) {
            log.severe(e.getMessage());
        }
    }

    private void writeTranslationSimplexPlot() {

        if (plotter == null) {
            return;
        }

        try {
            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());
        } catch (IOException e) {
            log.severe(e.getMessage());
        }

    }

    /**
     * method to narrow down all parameter space to 10 fits in which at least
     * one contains a fit within +- 10 degrees of rotation and +- 20 pixels
     * of translation in x and y.  Note, still testing so these numbers
     * may change.
     * The 10 fits should be followed by a downhill simplex to further
     * constrain the transformation, or a detailed grid search of each of
     * the 10 fits with something like a delta rotation of 2 and delta
     * translation of 4 and ranges based on the error stated in the first
     * sentence.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    protected TransformationPointFit[] preSearch0(
        PairIntArray set1, PairIntArray set2, float scale,
        float startRotationInDegrees, float stopRotationInDegrees,
        boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        TransformationPointFit[] fits = null;

        // if number of points is very high, need to use alternate method
        if (false && (nMaxMatchable > largeSearchLimit)) {
            fits = preSearch0Large(set1, set2, scale, startRotationInDegrees,
                stopRotationInDegrees, useGreedyMatching);
        } else {
            fits = preSearch0Small(set1, set2, scale, startRotationInDegrees,
                stopRotationInDegrees, useGreedyMatching);
        }

        return fits;
    }

    /**
     * method to narrow down all parameter space to 10 fits in which at least
     * one contains a fit within +- 10 degrees of rotation and +- 20 pixels
     * of translation in x and y.  Note, still testing so these numbers
     * may change.
     * The 10 fits should be followed by a downhill simplex to further
     * constrain the transformation, or a detailed grid search of each of
     * the 10 fits with something like a delta rotation of 2 and delta
     * translation of 4 and ranges based on the error stated in the first
     * sentence.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    protected TransformationPointFit[] preSearch0Small(
        PairIntArray set1, PairIntArray set2, float scale,
        float rotationStartInDegrees, float rotationStopInDegrees,
        boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        /*
        The best solutions use a fine grid search of rotDelta=2 and transDelta=4
        w/ translation tolerance = transDelta.

        Such a solution takes a long time so a coarser grid is attempted here.

        Trying a fine grid in rotation and a very coarse grid in translation
        to get rotation accurate, then a finer search on translation near
        that found rotation.
        */

        /*
        editing this while testing:

        rotDelta=1 and transDelta of 100 has a smallest delta rot within the
            10 best fits which is accurate within += 20 degrees?

        rotDelta=2 and transDelta of 100 has a smallest delta rot
            within the 10 best fits which is accurate within += 25 degrees... 14, 5
            and within the best fit alone, is accurate within

        */

        float tolTransX = 0.5f*transDeltaPreSearch0 +
            (float)(Math.sin(rotDeltaPreSearch0*Math.PI/180.)*transDeltaPreSearch0);
        float tolTransY = tolTransX;

        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }

        Transformer transformer = new Transformer();

        float[] rotation = MiscMath.writeDegreeIntervals(rotationStartInDegrees,
            rotationStopInDegrees, rotDeltaPreSearch0);

        int nKeep = 10;

        long t0 = System.currentTimeMillis();
        log.info("starterPoints.length=" + nKeep);

        int minX2 = MiscMath.findMin(set2.getX(), set2.getN());
        int maxX2 = MiscMath.findMax(set2.getX(), set2.getN());
        int minY2 = MiscMath.findMin(set2.getY(), set2.getN());
        int maxY2 = MiscMath.findMax(set2.getY(), set2.getN());

        FixedSizeSortedVector<TransformationPointFit> starterPoints
            = new FixedSizeSortedVector<>(nKeep, TransformationPointFit.class);

        float[] xsr = new float[set1.getN()];
        float[] ysr = new float[set1.getN()];

        for (float rotDeg : rotation) {

            TransformationParameters rotScaleParams =
                new TransformationParameters();
            rotScaleParams.setScale(scale);
            rotScaleParams.setRotationInDegrees(rotDeg);
            rotScaleParams.setTranslationX(0);
            rotScaleParams.setTranslationY(0);
            rotScaleParams.setOriginX(0);
            rotScaleParams.setOriginY(0);

            transformer.applyTransformation(rotScaleParams, set1, xsr, ysr);

            final float minX1 = MiscMath.findMin(xsr);
            final float maxX1 = MiscMath.findMax(xsr);
            final float minY1 = MiscMath.findMin(ysr);
            final float maxY1 = MiscMath.findMax(ysr);

            float transXStart = minX1 - maxX2;
            float transXStop = maxX1 - minX2;
            float transYStart = minY1 - maxY2;
            float transYStop = maxY1 - minY2;

            float transXRange = transXStop - transXStart;
            float transYRange = transYStop - transYStart;

            int nX = (int) Math.ceil(transXRange / transDeltaPreSearch0);
            int nY = (int) Math.ceil(transYRange / transDeltaPreSearch0);

            float[] txS = new float[nX];
            txS[0] = transXStart;
            for (int i = 1; i < nX; ++i) {
                txS[i] = txS[i - 1] + transDeltaPreSearch0;
            }
            float[] tyS = new float[nY];
            tyS[0] = transYStart;
            for (int i = 1; i < nY; ++i) {
                tyS[i] = tyS[i - 1] + transDeltaPreSearch0;
            }

            for (float tx : txS) {
                for (float ty : tyS) {

                    PairFloatArray allPoints1Tr = new PairFloatArray(set1.getN());

                    for (int i = 0; i < xsr.length; ++i) {
                        allPoints1Tr.add(xsr[i] + tx, ysr[i] + ty);
                    }

                    TransformationParameters params = rotScaleParams.copy();
                    params.setTranslationX(tx);
                    params.setTranslationY(ty);

                    TransformationPointFit fit2;
                    if (useGreedyMatching) {
                        fit2 = evaluateFitForUnMatchedTransformedGreedy(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);
                    } else {
                        fit2 = evaluateFitForUnMatchedTransformedOptimal(params,
                            allPoints1Tr, set2, tolTransX, tolTransY);
                    }

                    starterPoints.add(fit2);
                }
            }
        }

        if (starterPoints.getArray()[0] == null) {
            return null;
        }

        long tm = System.currentTimeMillis() - t0;
        double ts = tm * 1e-3;
        log.info("starterPoints finished for nMaxMatchable=" + nMaxMatchable +
            " seconds=" + ts);

        if (starterPoints.getNumberOfItems() > 0) {
            log.info("best of starter points=" + starterPoints.getArray()[0].toString());
        }

        return starterPoints.getArray();
    }

    protected TransformationPointFit[]
    refineTransformationWithDownhillSimplexWrapper(TransformationPointFit[] fits,
        PairIntArray set1, PairIntArray set2, boolean searchScaleToo,
        boolean useGreedyMatching) {

        if (fits == null) {
            throw new IllegalArgumentException("fits cannot be null");
        }
        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        float centerTransX = 0;
        float centerTransY = 0;
        float minTransX = Float.MAX_VALUE;
        float minTransY = Float.MAX_VALUE;
        float maxTransX = Float.MIN_VALUE;
        float maxTransY = Float.MIN_VALUE;
        float maxTolTransX = Float.MIN_VALUE;;
        float maxTolTransY = Float.MIN_VALUE;
        int count = 0;
        for (TransformationPointFit fit : fits) {
            if (fit == null) {
                continue;
            }
            float tx = fit.getTranslationX();
            float ty = fit.getTranslationY();
            float tolX = fit.getTranslationXTolerance();
            float tolY = fit.getTranslationYTolerance();

            centerTransX += tx;
            centerTransY += ty;

            if (tx < minTransX) {
                minTransX = tx;
            }
            if (ty < minTransY) {
                minTransY = ty;
            }
            if (tx > maxTransX) {
                maxTransX = tx;
            }
            if (ty > maxTransY) {
                maxTransY = ty;
            }
            if (tolX > maxTolTransX) {
                maxTolTransX = tolX;
            }
            if (tolY > maxTolTransY) {
                maxTolTransY = tolY;
            }

            count++;
        }
        centerTransX /= (float)count;
        centerTransY /= (float)count;
        float transXHalfRange = Math.abs(maxTransX - minTransX)/2.f;
        transXHalfRange += this.transDeltaPreSearch0;
        float transYHalfRange = Math.abs(maxTransY - minTransY)/2.f;
        transYHalfRange += this.transDeltaPreSearch0;

        TransformationPointFit[] fits2 = refineTransformationWithDownhillSimplex(
            fits, set1, set2,
            centerTransX, centerTransY, maxTolTransX, maxTolTransY,
            transXHalfRange, transYHalfRange,
            searchScaleToo, useGreedyMatching, false, 50);

        return fits2;
    }

    /**
     * method to narrow down all parameter space to a fit that is within
     * +- 20 degrees in rotation and +- 50 pixels in translation in X and Y
     * of the true Euclidean transformation.
     * Note, these numbers are being learned in testing right now so may
     * change.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit preSearch(PairIntArray set1,
        PairIntArray set2, float scale, boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        float startRotationInDegrees = 0;
        float stopRotationInDegrees = 359;

        return preSearch(set1, set2, scale, startRotationInDegrees,
            stopRotationInDegrees, useGreedyMatching);
    }

    /**
     * method to narrow down all parameter space to a fit that is within
     * +- 20 degrees in rotation and +- 50 pixels in translation in X and Y
     * of the true Euclidean transformation.
     * Note, these numbers are being learned in testing right now so may
     * change.
     *
     * @param set1
     * @param set2
     * @param scale
     * @param startRotationInDegrees
     * @param useGreedyMatching
     * @param stopRotationInDegrees
     * @return
     */
    protected TransformationPointFit preSearch(PairIntArray set1,
        PairIntArray set2, float scale, float startRotationInDegrees,
        float stopRotationInDegrees, boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        TransformationPointFit[] starterPoints = preSearch0(set1, set2, scale,
            startRotationInDegrees, stopRotationInDegrees,
            useGreedyMatching);

        if (starterPoints == null) {
            return null;
        }

        boolean searchScaleToo = false;

        TransformationPointFit[] fits2 =
            refineTransformationWithDownhillSimplexWrapper(starterPoints,
            set1, set2, searchScaleToo, useGreedyMatching);

        if (fits2 == null || fits2[0] == null) {
            return null;
        }

        TransformationPointFit bestFit = fits2[0];
        /*
        TransformationPointFit fit2 = refineAfterCalculationWithPairs(
            bestFit.getParameters().copy(), set1, set2, useGreedyMatching);

        if (fitIsBetter(bestFit, fit2)) {
            bestFit = fit2;
        }
        */
        return bestFit;
    }

    protected TransformationPointFit[] boundedGridSearch(
        PairIntArray set1, PairIntArray set2,
        final TransformationParameters params,
        float scaleHalfRange, float scaleDelta,
        float rotHalfRangeInDegrees, float rotDeltaInDegrees,
        float transXHalfRange, float transXDelta,
        float transYHalfRange, float transYDelta,
        boolean useGreedyMatching) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        float tolTransX = 0.5f*transXDelta + (float)(Math.sin(rotDeltaInDegrees*Math.PI/180.)*transXDelta);
        float tolTransY = 0.5f*transYDelta + (float)(Math.sin(rotDeltaInDegrees*Math.PI/180.)*transYDelta);

        if (tolTransX < minTolerance) {
            tolTransX = generalTolerance;//minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = generalTolerance;//minTolerance;
        }

        Transformer transformer = new Transformer();

        float[] rotation = MiscMath.writeDegreeIntervals(
            params.getRotationInDegrees() - rotHalfRangeInDegrees,
            params.getRotationInDegrees() + rotHalfRangeInDegrees,
            rotDeltaInDegrees);

        log.fine("searching rotation: " + Arrays.toString(rotation));

        int nX = (int) Math.ceil(2*transXHalfRange / transXDelta);
        int nY = (int) Math.ceil(2*transYHalfRange / transYDelta);
        if (nX == 0) {
            nX = 1;
        }
        if (nY == 0) {
            nY = 1;
        }

        float[] txS = new float[nX];
        txS[0] = params.getTranslationX() - transXHalfRange;
        for (int i = 1; i < nX; ++i) {
            txS[i] = txS[i - 1] + transXDelta;
        }
        float[] tyS = new float[nY];
        tyS[0] = params.getTranslationY() - transYHalfRange;
        for (int i = 1; i < nY; ++i) {
            tyS[i] = tyS[i - 1] + transYDelta;
        }

        int nS = (int) Math.ceil(2*scaleHalfRange / scaleDelta);
        if (nS == 0) {
            nS = 1;
        }
        float[] scales = new float[nS];
        scales[0] = params.getScale() - scaleHalfRange;
        for (int i = 1; i < nS; ++i) {
            scales[i] = scales[i - 1] + scaleDelta;
        }

        int nKeep = 10;

        long t0 = System.currentTimeMillis();

        FixedSizeSortedVector<TransformationPointFit> starterPoints
            = new FixedSizeSortedVector<>(nKeep, TransformationPointFit.class);

        // add the given params so best answer is never worse than that
        if (useGreedyMatching) {
            PairFloatArray allPoints1Tr = transformer.applyTransformation2(
                params, set1);
            TransformationPointFit fit =
                evaluateFitForUnMatchedTransformedGreedy(params,
                allPoints1Tr, set2, tolTransX, tolTransY);
            starterPoints.add(fit);
        } else {
            PairFloatArray allPoints1Tr = transformer.applyTransformation2(
                params, set1);
            TransformationPointFit fit =
                evaluateFitForUnMatchedTransformedOptimal(params,
                allPoints1Tr, set2, tolTransX, tolTransY);
            starterPoints.add(fit);
        }

        float[] xsr = new float[set1.getN()];
        float[] ysr = new float[set1.getN()];

        for (float scale : scales) {
            for (float rotDeg : rotation) {

                TransformationParameters rotScaleParams =
                    new TransformationParameters();
                rotScaleParams.setScale(scale);
                rotScaleParams.setRotationInDegrees(rotDeg);
                rotScaleParams.setTranslationX(0);
                rotScaleParams.setTranslationY(0);
                rotScaleParams.setOriginX(params.getOriginX());
                rotScaleParams.setOriginY(params.getOriginY());

                transformer.applyTransformation(rotScaleParams, set1, xsr, ysr);

                for (float tx : txS) {
                    for (float ty : tyS) {

                        PairFloatArray allPoints1Tr = new PairFloatArray(set1.getN());

                        for (int i = 0; i < xsr.length; ++i) {
                            allPoints1Tr.add(xsr[i] + tx, ysr[i] + ty);
                        }

                        TransformationParameters params2 = rotScaleParams.copy();
                        params2.setTranslationX(tx);
                        params2.setTranslationY(ty);

                        TransformationPointFit fit2;
                        if (useGreedyMatching) {
                            fit2 = evaluateFitForUnMatchedTransformedGreedy(params2,
                                allPoints1Tr, set2, tolTransX, tolTransY);
                        } else {
                            fit2 = evaluateFitForUnMatchedTransformedOptimal(params2,
                                allPoints1Tr, set2, tolTransX, tolTransY);
                        }

                        starterPoints.add(fit2);
                    }
                }
            }
        }

        if (starterPoints.getArray()[0] == null) {
            return null;
        }

        long tm = System.currentTimeMillis() - t0;
        double ts = tm * 1e-3;
        log.fine("bounded grid search finished for nMaxMatchable=" + nMaxMatchable +
            " seconds=" + ts);

        if (starterPoints.getNumberOfItems() > 0) {
            log.fine("best=" + starterPoints.getArray()[0].toString());
        }

        return starterPoints.getArray();
    }

    /**
     * searches around the given params for a better fit using
     * a grid search and then a downhill simplex with the best 10 grid results
     * as starter points.
     *
     * @param set1
     * @param set2
     * @param params
     * @param rotHalfRangeInDegrees half of the range in rotation to search
     * around params.rotationInDegrees.  the value is added and subtracted
     * from params.rotationInDegrees to create the start and end of the
     * range.
     * @param rotDeltaInDegrees the interval size of rotation in degrees
     * for the search.
     * @param transXHalfRange half of the range in translation in X to search
     * around params.translationX.  the value is added and subtracted
     * from params.translationX to create the start and end of the range.
     * @param transYHalfRange half of the range in translation in Y to search
     * around params.translationY.  the value is added and subtracted
     * from params.translationY to create the start and end of the range.
     * @param transYDelta the interval size in translation of X for the search
     * @param transXDelta the interval size in translation of Y for the search
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    public TransformationPointFit refineTheTransformation(
        final TransformationParameters params, PairIntArray set1, PairIntArray set2,
        float rotHalfRangeInDegrees, float rotDeltaInDegrees,
        float transXHalfRange, float transXDelta,
        float transYHalfRange, float transYDelta,
        boolean useGreedyMatching) {

        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }

        float scaleHalfRange = 0;
        float scaleDelta = 1;

        TransformationPointFit[] fits = boundedGridSearch(
            set1, set2, params,
            scaleHalfRange, scaleDelta,
            rotHalfRangeInDegrees, rotDeltaInDegrees,
            transXHalfRange, transXDelta,
            transYHalfRange, transYDelta, useGreedyMatching);

        if (fits == null) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        if (hasConverged(fits[0], nMaxMatchable)) {
            return fits[0];
        }

        // finest grid so no need for downhill simplex search
        if ((rotDeltaInDegrees == 1) && (transXDelta == 1) && (transYDelta == 1)) {
            return fits[0];
        }

        boolean searchScaleToo = false;

        TransformationPointFit[] fits2 =
            refineTransformationWithDownhillSimplexWrapper(fits, set1, set2,
            searchScaleToo, useGreedyMatching);

        if ((fits2 == null) || (fits2.length == 0)) {
            return null;
        }

        TransformationPointFit fit2 = fits2[0];

        if (hasConverged(fit2, nMaxMatchable)) {
            return fit2;
        }

        // one last refinement
        scaleHalfRange = 0;
        scaleDelta = 1;
        rotHalfRangeInDegrees = 5;
        rotDeltaInDegrees = 1;
        transXHalfRange = 20;
        transYHalfRange = transXHalfRange;
        transXDelta = 1;
        transYDelta = transXDelta;

        TransformationPointFit[] fits3 = boundedGridSearch(
            set1, set2, fit2.getParameters(),
            scaleHalfRange, scaleDelta,
            rotHalfRangeInDegrees, rotDeltaInDegrees,
            transXHalfRange, transXDelta,
            transYHalfRange, transYDelta, useGreedyMatching);

        if ((fits3 == null) || (fits3.length == 0)) {
            return null;
        }

        return fits3[0];
    }

    /**
     * use a downhill simplex algorithm to search for a better fit to the data
     * given the starterPoints.  This method changes
     * scale, rotation, translation in X, and translation in Y.
     *
     * @param set1
     * @param set2
     * @param params
     * @param scaleHalfRange half of the range in scale to search
     * around params.scale.  the value is added and subtracted
     * from params.scale to create the start and end of the
     * range.
     * @param scaleDelta the interval size of scale for the search.
     * @param rotHalfRangeInDegrees half of the range in rotation to search
     * around params.rotationInDegrees.  the value is added and subtracted
     * from params.rotationInDegrees to create the start and end of the
     * range.
     * @param rotDeltaInDegrees the interval size of rotation in degrees
     * for the search.
     * @param transXHalfRange half of the range in translation in X to search
     * around params.translationX.  the value is added and subtracted
     * from params.translationX to create the start and end of the range.
     * @param transYHalfRange half of the range in translation in Y to search
     * around params.translationY.  the value is added and subtracted
     * from params.translationY to create the start and end of the range.
     * @param transYDelta the interval size in translation of X for the search
     * @param transXDelta the interval size in translation of Y for the search
     * @param useGreedyMatching if true, uses an N^2 algorithm to find the
     * best match made for each point in order of points in set1, else uses
     * an ~N^3 optimal bipartite matching algorithm.
     * @return
     */
    public TransformationPointFit refineTheTransformation(
        TransformationParameters params, PairIntArray set1, PairIntArray set2,
        float scaleHalfRange, float scaleDelta,
        float rotHalfRangeInDegrees, float rotDeltaInDegrees,
        float transXHalfRange, float transXDelta,
        float transYHalfRange, float transYDelta,
        boolean useGreedyMatching) {

        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (set2 == null) {
            throw new IllegalArgumentException("set2 cannot be null");
        }

        TransformationPointFit[] fits = boundedGridSearch(
            set1, set2, params,
            scaleHalfRange, scaleDelta,
            rotHalfRangeInDegrees, rotDeltaInDegrees,
            transXHalfRange, transXDelta,
            transYHalfRange, transYDelta, useGreedyMatching);

        if (fits == null) {
            return null;
        }

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        if (hasConverged(fits[0], nMaxMatchable)) {
            return fits[0];
        }

        boolean searchScaleToo = true;

        TransformationPointFit[] fits2 =
            refineTransformationWithDownhillSimplexWrapper(fits, set1, set2,
            searchScaleToo, useGreedyMatching);

        if (fits2 == null) {
            return null;
        }

        if (hasConverged(fits2[0], nMaxMatchable)) {
            return fits2[0];
        }

        // TODO: may need a small fine grid search here

        return fits2[0];

    }

    /**
     * for use when the number of data points is too high for reasonably
     * fast solution using preSearch0.
     * The method randomly reduces set1 to a smaller number of point sets
     * and finds the best solution among those.
     * @param set1
     * @param set2
     * @param scale
     * @param useGreedyMatch
     * @return
     */
    protected TransformationPointFit[] preSearch0Large(PairIntArray set1,
        PairIntArray set2, float scale,
        float startRotationInDegrees, float stopRotationInDegrees,
        boolean useGreedyMatch) {

        if ((set1 == null) || (set2 == null)) {
            return null;
        }
        if ((set1.getN() < 3) || (set2.getN() < 3)) {
            return null;
        }

        PointPartitioner partitioner = new PointPartitioner();

        PairIntArray[] set1Subsets = null;
        if (set1.getN() < 100) {
            set1Subsets = new PairIntArray[]{set1};
        } else {
            set1Subsets = partitioner.reduceByBinSampling(set1, 5, 4);
            //set1Subsets = partitioner.reduceByBinSampling(set1, 10, 2);
            //set1Subsets = partitioner.reduceByBinSampling(set1, 5, 2);

            set1Subsets = new PairIntArray[]{set1Subsets[0]};
        }

        boolean useGreedyMatching = true;

        for (PairIntArray set1Subset : set1Subsets) {

            TransformationPointFit[] fits = preSearch0Small(set1Subset, set2,
                scale, startRotationInDegrees,
                stopRotationInDegrees, useGreedyMatching);

            return fits;
        }

        return null;
    }

    /**
     Calculate the best transformation solution by matching pairs of points
     * in set1 with pairs of points in set2.
     *
     * Note that set1 and set2 should have fewer than __ members each
     * in order to use this method and get results reasonably fast.
     *
     * The number of ways to make pairs in set 1 times the number of ways to make
       pairs in set 2 =
             n_1!            n_2!
         ------------  X  ------------
         2*(n_1 - 2)!     2*(n_2 - 2)!

    After points are matched:
       Estimating scale:
           Take 2 pairs of points in both datasets and compute the distance
           between them and then take the ratio:

           scale = (distance between pair in set 1)
                   / (distance between pair in set 2)

           Note: scale is roughly determined from contour matching too in the
              inflection matcher.

       Estimating rotation:
           Take the same 2 pairs and determine the difference in their angles:
               tan(theta) = delta y / delta x

           rotation = atan((delta y between pair in set 1)
                          /(delta x between pair in set 1))
                      -
                      atan((delta y between pair in set 2)
                          /(delta x between pair in set 2))

       Estimate translation:
           Performed on one point in set 1 with its candidate match in set 2:
           From the full transformation equation, we can rewrite:
               transX = xt0 - xc*scale -
                   (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

               transY = yt0 - yc*scale -
                    ((-(x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

               where (xc, yc) is the center of the first image
               *
    * Note that if the solution converges, it will return before trying all
    * subset combinations.
    */
    protected List<TransformationPointFit>
    calculateEuclideanTransformationUsingPairs(
        PairIntArray set1, PairIntArray set2, boolean earlyConvergeReturn) {

        boolean useLargestToleranceForOutput = true;
        boolean useGreedyMatching = true;

        List<TransformationPointFit>  fits = calculateEuclideanTransformationUsingPairs(
            set1, set2, earlyConvergeReturn,
            useLargestToleranceForOutput, useGreedyMatching);

        return fits;
    }

    protected List<TransformationPointFit>
    calculateEuclideanTransformationForSmallSets(
        PairIntArray set1, PairIntArray set2,
        PairIntArray fullSet1, PairIntArray fullSet2,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        float toleranceTransX;
        float toleranceTransY;
        if (useLargestToleranceForOutput) {
            toleranceTransX = generalTolerance * (float)Math.sqrt(1./2);
            toleranceTransY = toleranceTransX;
        } else {
            toleranceTransX = 4;
            toleranceTransY = 4;
        }

        long n1P = (set1.getN() * (set1.getN() - 1))/2;

        return calculateEuclideanTransformationForSmallSets(set1, set2,
            fullSet1, fullSet2, n1P,
            toleranceTransX, toleranceTransY, useGreedyMatching);
    }

    public List<TransformationPointFit>
    calculateEuclideanTransformationForSmallSets(PairIntArray set1,
        PairIntArray set2, PairIntArray fullSet1, PairIntArray fullSet2,
        long n1Iterations,
        float toleranceTransX, float toleranceTransY,
        boolean useGreedyMatching) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        int k = 2;

        SubsetChooser s1 = new SubsetChooser(n1, k);

        int[] selected1 = new int[k];
        int[] selected2 = new int[k];

        /*
        save bestFit for scale near 1,
                   bestFit for scale > 1
                   bestFit for scale < 1
        and normalized for each.
        */
        TransformationPointFit bestFitScaleGT1 = null;
        TransformationPointFit bestFitScaleLT1 = null;
        TransformationPointFit bestFitScale1 = null;
        TransformationPointFit bestFitNormalizedScaleGT1 = null;
        TransformationPointFit bestFitNormalizedScaleLT1 = null;
        TransformationPointFit bestFitNormalizedScale1 = null;

        int maxNMatchable = Math.min(n1, n2);

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        TransformationPointFit bestFitNormalized = null;
        TransformationPointFit bestFit = null;

        List<TransformationPointFit> similarToBestFitNormalized = new ArrayList<TransformationPointFit>();
        List<TransformationPointFit> similarToBestFit = new ArrayList<TransformationPointFit>();

        long n1Count = 0;

        while (s1.getNextSubset(selected1) != -1) {

            if (n1Count > n1Iterations) {
                break;
            }

            int idx10 = selected1[0];
            int idx11 = selected1[1];

            SubsetChooser s2 = new SubsetChooser(n2, k);

            //TODO: could save the calculations for point 1 here

            while (s2.getNextSubset(selected2) != -1) {

                int idx20 = selected2[0];
                int idx21 = selected2[1];

                for (int order = 0; order < 2; ++order) {

                    if (order == 1) {
                        int swap = idx20;
                        idx20 = idx21;
                        idx21 = swap;
                    }

                    TransformationParameters params = tc.calulateEuclidean(
                        set1.getX(idx10), set1.getY(idx10),
                        set1.getX(idx11), set1.getY(idx11),
                        set2.getX(idx20), set2.getY(idx20),
                        set2.getX(idx21), set2.getY(idx21),
                        0, 0);

                    TransformationPointFit fit = evaluateForUnmatched(params,
                        fullSet1, fullSet2, toleranceTransX, toleranceTransY,
                        useGreedyMatching);

                    if ((fit == null) || (fit.getNumberOfMatchedPoints() < 2)) {
                        continue;
                    }

                    float scale = params.getScale();

                    if ((scale >= 0.95) && (scale <= 1.05)) {
                        if (fitIsBetter(bestFitScale1, fit)) {
                            bestFitScale1 = fit;
                        }
                        if (fitIsBetterNormalized(bestFitNormalizedScale1, fit)) {
                            bestFitNormalizedScale1 = fit;
                        }
                    } else if (scale > 1.05) {
                        if (fitIsBetter(bestFitScaleGT1, fit)) {
                            bestFitScaleGT1 = fit;
                        }
                        if (fitIsBetterNormalized(bestFitNormalizedScaleGT1, fit)) {
                            bestFitNormalizedScaleGT1 = fit;
                        }
                    } else if (scale < 0.95) {
                        if (fitIsBetter(bestFitScaleLT1, fit)) {
                            bestFitScaleLT1 = fit;
                        }
                        if (fitIsBetterNormalized(bestFitNormalizedScaleLT1, fit)) {
                            bestFitNormalizedScaleLT1 = fit;
                        }
                    }

                    if (fitIsBetterNormalized(bestFitNormalized, fit)) {
                        if ((bestFitNormalized != null) && (bestFitNormalized.getMeanDistFromModel() < 1)
                            && (fit.getMeanDistFromModel() < 1)) {
                            if (similarToBestFitNormalized.isEmpty()) {
                                similarToBestFitNormalized.add(bestFitNormalized);
                            }
                            similarToBestFitNormalized.add(fit);
                        } else {
                            similarToBestFitNormalized.clear();
                        }
                        bestFitNormalized = fit;
                    }
                    if (fitIsBetter(bestFit, fit)) {
                        if ((bestFit != null) && (bestFit.getMeanDistFromModel() < 1)
                            && (fit.getMeanDistFromModel() < 1)) {
                            if (similarToBestFit.isEmpty()) {
                                similarToBestFit.add(bestFit);
                            }
                            similarToBestFit.add(fit);
                        } else {
                            similarToBestFit.clear();
                        }
                        bestFit = fit;
                    }
                }
            }

            n1Count++;
        }

        log.info("similarToBestFit.size=" + similarToBestFit.size());
        if (similarToBestFit.size() > 1) {
            for (TransformationPointFit fit : similarToBestFit) {
                log.info("similarFit=" + fit.toString());
            }
        }
        log.info("similarToBestFitNormalized.size=" + similarToBestFitNormalized.size());
        if (similarToBestFitNormalized.size() > 1) {
            for (TransformationPointFit fit : similarToBestFitNormalized) {
                log.info("similarFit=" + fit.toString());
            }
        } else if (similarToBestFitNormalized.isEmpty() && (bestFitNormalized != null)) {
            log.info("bestFitNormalized=" + bestFitNormalized.toString());
            similarToBestFitNormalized.add(bestFitNormalized);
        }
        log.info("bestFit(not normalized)=" + bestFit.toString());
        similarToBestFitNormalized.add(0, bestFit);

        if (bestFitNormalized != null) {
            log.info("bestFit(normalized)=" + bestFitNormalized.toString());
        }

        if (bestFitScaleGT1 != null) {
            log.info("bestFitScaleGT1=" + bestFitScaleGT1.toString());
            similarToBestFitNormalized.add(bestFitScaleGT1);
        }

        if (bestFitScaleLT1 != null) {
            log.info("bestFitScaleLT1=" + bestFitScaleLT1.toString());
            similarToBestFitNormalized.add(bestFitScaleLT1);
        }

        if (bestFitScale1 != null) {
            log.info("bestFitScale1=" + bestFitScale1.toString());
            similarToBestFitNormalized.add(bestFitScale1);
        }

        if (bestFitNormalizedScaleGT1 != null) {
            log.info("bestFitNormalizedScaleGT1=" + bestFitNormalizedScaleGT1.toString());
            similarToBestFitNormalized.add(bestFitNormalizedScaleGT1);
        }

        if (bestFitNormalizedScaleLT1 != null) {
            log.info("bestFitNormalizedScaleLT1=" + bestFitNormalizedScaleLT1.toString());
            similarToBestFitNormalized.add(bestFitNormalizedScaleLT1);
        }

        if (bestFitNormalizedScale1 != null) {
            log.info("bestFitNormalizedScale1=" + bestFitNormalizedScale1.toString());
            similarToBestFitNormalized.add(bestFitNormalizedScale1);
        }

        return similarToBestFitNormalized;
    }

    protected List<TransformationPointFit>
    calculateEuclideanTransformationForSmallSets(
        PairIntArray set1, PairIntArray set2,
        final float scale, final float rotationLowLimitInDegrees,
        final float rotationHighLimitInDegrees,
        boolean earlyConvergeReturn,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        int k = 2;

        int[] selected1 = new int[k];
        int[] selected2 = new int[k];

        int maxNMatchable = Math.min(n1, n2);

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();
        float toleranceTransX;
        float toleranceTransY;
        if (useLargestToleranceForOutput) {
            toleranceTransX = generalTolerance * (float)Math.sqrt(1./2);
            toleranceTransY = toleranceTransX;
        } else {
            toleranceTransX = 4;
            toleranceTransY = 4;
        }

        float rotRange = AngleUtil.getAngleDifference(rotationLowLimitInDegrees,
            rotationHighLimitInDegrees);

        TransformationPointFit bestFitNormalized = null;

        List<TransformationPointFit> similarToBestFitNormalized = new ArrayList<TransformationPointFit>();

        TransformationPointFit bestFit = null;

        List<TransformationPointFit> similarToBestFit = new ArrayList<TransformationPointFit>();

        SubsetChooser s1 = new SubsetChooser(n1, k);

        while (s1.getNextSubset(selected1) != -1) {

            int idx10 = selected1[0];
            int idx11 = selected1[1];

            SubsetChooser s2 = new SubsetChooser(n2, k);

            //TODO: could save the calculations for point 1 here

            while (s2.getNextSubset(selected2) != -1) {

                int idx20 = selected2[0];
                int idx21 = selected2[1];

                for (int order = 0; order < 2; ++order) {

                    if (order == 1) {
                        int swap = idx20;
                        idx20 = idx21;
                        idx21 = swap;
                    }

                    TransformationParameters params = tc.calulateEuclidean(
                        set1.getX(idx10), set1.getY(idx10),
                        set1.getX(idx11), set1.getY(idx11),
                        set2.getX(idx20), set2.getY(idx20),
                        set2.getX(idx21), set2.getY(idx21),
                        0, 0);
                    float rotD = params.getRotationInDegrees();

                    if (Math.abs(scale - params.getScale()) > 0.05) {
                        continue;
                    }
                    if (!
                    (((rotD >= rotationLowLimitInDegrees)
                      && (Math.abs(AngleUtil.getAngleDifference(rotD, rotationLowLimitInDegrees)) < rotRange))
                    ||
                    ((rotD <= rotationHighLimitInDegrees)
                      && (Math.abs(AngleUtil.getAngleDifference(rotD, rotationHighLimitInDegrees)) < rotRange))
                    )) {
                        continue;
                    }

                    TransformationPointFit fit = evaluateForUnmatched(params,
                        set1, set2, toleranceTransX, toleranceTransY,
                        useGreedyMatching);

                    if ((fit == null) || (fit.getNumberOfMatchedPoints() < 2)) {
                        continue;
                    }

                    if (fitIsBetterNormalized(bestFitNormalized, fit)) {
                        if ((bestFitNormalized != null) && (bestFitNormalized.getMeanDistFromModel() < 1)
                            && (fit.getMeanDistFromModel() < 1)) {
                            if (similarToBestFitNormalized.isEmpty()) {
                                similarToBestFitNormalized.add(bestFitNormalized);
                            }
                            similarToBestFitNormalized.add(fit);
                        } else {
                            similarToBestFitNormalized.clear();
                        }
                        bestFitNormalized = fit;
                        if (earlyConvergeReturn
                            && hasConverged(bestFitNormalized, maxNMatchable)) {
                            if (similarToBestFitNormalized.isEmpty()) {
                                similarToBestFitNormalized.add(bestFitNormalized);
                            }
                            return similarToBestFitNormalized;
                        }
                    }
                    if (fitIsBetter(bestFit, fit)) {
                        if ((bestFit != null) && (bestFit.getMeanDistFromModel() < 1)
                            && (fit.getMeanDistFromModel() < 1)) {
                            if (similarToBestFit.isEmpty()) {
                                similarToBestFit.add(bestFit);
                            }
                            similarToBestFit.add(fit);
                        } else {
                            similarToBestFit.clear();
                        }
                        bestFit = fit;
                        if (earlyConvergeReturn
                            && hasConverged(bestFit, maxNMatchable)) {
                            if (similarToBestFit.isEmpty()) {
                                similarToBestFit.add(bestFit);
                            }
                            return similarToBestFit;
                        }
                    }
                }
            }
        }

        log.info("similarToBestFit.size=" + similarToBestFit.size());
        if (similarToBestFit.size() > 1) {
            for (TransformationPointFit fit : similarToBestFit) {
                log.info("similarFit=" + fit.toString());
            }
        }
        log.info("similarToBestFitNormalized.size=" + similarToBestFitNormalized.size());
        if (similarToBestFitNormalized.size() > 1) {
            for (TransformationPointFit fit : similarToBestFitNormalized) {
                log.info("similarFit=" + fit.toString());
            }
        } else if (similarToBestFitNormalized.isEmpty() && (bestFitNormalized != null)) {
            log.info("bestFitNormalized=" + bestFitNormalized.toString());
            similarToBestFitNormalized.add(bestFitNormalized);
        }
        if (bestFit != null) {
            log.info("bestFit(not normalized)=" + bestFit.toString());
            similarToBestFitNormalized.add(0, bestFit);
        }

        return similarToBestFitNormalized;
    }

    /**
     * this one expects a left to right ordering of images and scale near 1
     * and rotation near zero and a y translation near zero.
     *
     * runtime is O(N^2) + O(N)
     *
     * @param set1
     * @param set2
     *
     * @param useLargestToleranceForOutput
     * @param useGreedyMatching
     * @return
     */
    protected TransformationPointFit calculateEuclideanLeftRightTransformation(
        PairIntArray set1, PairIntArray set2,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        //TODO: consider that the translation in X should always be negative?

        int n1 = set1.getN();
        int n2 = set2.getN();

        float toleranceTransX;
        float toleranceTransY;
        if (useLargestToleranceForOutput) {
            toleranceTransX = generalTolerance * (float)Math.sqrt(1./2);
            toleranceTransY = toleranceTransX;
        } else {
            toleranceTransX = 4;
            toleranceTransY = 4;
        }

        List<Integer> diffXs = new ArrayList<Integer>();
        List<Integer> diffYs = new ArrayList<Integer>();

        for (int i = 0; i < n1; ++i) {
            int x1 = set1.getX(i);
            int y1 = set1.getY(i);
            for (int j = 0; j < n2; ++j) {
                int y2 = set2.getY(j);
                int diffY = y2 - y1;

                if (Math.abs(diffY) > toleranceTransY) {
                    continue;
                }

                int x2 = set2.getX(j);
                int diffX = x2 - x1;

                diffXs.add(Integer.valueOf(diffX));
                diffYs.add(Integer.valueOf(diffY));
            }
        }

        float[] values = new float[diffXs.size()];
        for (int i = 0; i < diffXs.size(); ++i) {
            values[i] = diffXs.get(i).floatValue();
        }

        HistogramHolder hist = Histogram.createSimpleHistogram(toleranceTransX,
            values, Errors.populateYErrorsBySqrt(values));

        if (hist.getXHist() == null || hist.getXHist().length == 0) {
            return null;
        }

        // sort and evaluate top 3 answers

        MultiArrayMergeSort.sortBy1stArgThen2nd(hist.getYHistFloat(), hist.getXHist());

        TransformationPointFit bestFit = null;
        TransformationPointFit bestFitNormalized = null;

        TransformationParameters params = new TransformationParameters();
        params.setScale(1);
        params.setRotationInDegrees(0);
        params.setOriginX(0);
        params.setOriginY(0);

        for (int i = 0; i < 3; i++) {

            int hIdx = hist.getYHistFloat().length - i - 1;

            float tX = hist.getXHist()[hIdx];

            double tY = 0;
            int count = 0;
            for (int ii = 0; ii < diffXs.size(); ++ii) {
                int diffX = diffXs.get(ii);
                if (Math.abs(diffX - tX) <= (toleranceTransX/2.)) {
                    tY += diffYs.get(ii);
                    count++;
                }
            }
            tY /= (double)count;

            params = params.copy();
            params.setTranslationX(tX);
            params.setTranslationY((float)tY);

            TransformationPointFit fit = evaluateForUnmatched(params,
                set1, set2, toleranceTransX, toleranceTransY,
                useGreedyMatching);

            if (fitIsBetter(bestFit, fit)) {
                bestFit = fit;
            }
            if (fitIsBetterNormalized(bestFitNormalized, fit)) {
                bestFitNormalized = fit;
            }
        }

        return bestFit;
    }

    List<TransformationPointFit> calculateEuclideanTransformationUsingPairs(
        PairIntArray set1, PairIntArray set2, boolean earlyConvergeReturn,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        float toleranceTransX;
        float toleranceTransY;
        if (useLargestToleranceForOutput) {
            toleranceTransX = generalTolerance * (float)Math.sqrt(1./2);
            toleranceTransY = toleranceTransX;
        } else {
            toleranceTransX = 4;
            toleranceTransY = 4;
        }

        return calculateEuclideanTransformationUsingPairs(
            set1, set2, toleranceTransX, toleranceTransY, earlyConvergeReturn,
            useGreedyMatching);
    }

    List<TransformationPointFit> calculateEuclideanTransformationUsingPairs(
        PairIntArray set1, PairIntArray set2,
        float toleranceTransX, float toleranceTransY,
        boolean earlyConvergeReturn, boolean useGreedyMatching) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        long np = numberOfPairPermutations(n1, n2);

        // if np is very large, need to be able to exit the algorithm
        // earlier than the full try of all permutations.
        // TODO: this may need to be set as an option to allow it to finish
        //     if user really does want all permutations tried.
        if ((np > 1e8) || ((n2 > largeSearchLimit) && earlyConvergeReturn)){
            /*
            if reduce set2 to a randomly sampled smaller set, the percentage
            of truly matchable points should remain the same.
            we then want a significant number of set1 tried against
            set 2 (but always evaluated with all points) before exiting
            the method instead of trying all , preferably, set1 would be
            randomly sampled.
            The goal would be to keep n1*(n1-1)*n2*(n2-1)*(1./4.)*n1*n2 < 1e9
            reducing n2 to 30 and randomly sampling 30 from n1 results in
            0.17 billion steps.
            */
            return calculateEuclideanTransformationUsingPairsPartitioned(
                 set1, set2, toleranceTransX, toleranceTransY,
                 useGreedyMatching);
        }

        long n1P = (set1.getN() * (set1.getN() - 1))/2;

        return calculateEuclideanTransformationForSmallSets(set1, set2,
            set1, set2, n1P,
            toleranceTransX, toleranceTransY, useGreedyMatching);
    }

    List<TransformationPointFit> calculateEuclideanTransformationUsingPairs(
        PairIntArray set1, PairIntArray set2,
        float scale,
        float rotationLowLimitInDegrees, float rotationHighLimitInDegrees,
        boolean earlyConvergeReturn,
        boolean useLargestToleranceForOutput,
        boolean useGreedyMatching) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        if ((n2 > largeSearchLimit) && earlyConvergeReturn) {
            return calculateEuclideanTransformationUsingPairsPartitioned(
                set1, set2,
                scale, rotationLowLimitInDegrees,
                rotationHighLimitInDegrees,
                useLargestToleranceForOutput, useGreedyMatching);
        }

        return calculateEuclideanTransformationForSmallSets(set1, set2,
            scale, rotationLowLimitInDegrees, rotationHighLimitInDegrees,
            earlyConvergeReturn, useLargestToleranceForOutput, useGreedyMatching);
    }

    /**
     * method to perform pairwise transformation calculations on datasets that
     * are to large to feasibly visit every permutation.
     * The name may change soon to better reflect the internal implementation
     * which is no longer partitions.
     * @param set1
     * @param set2
     * @param toleranceTransX
     * @param toleranceTransY
     * @param useGreedyMatching
     * @return
     */
    protected List<TransformationPointFit>
    calculateEuclideanTransformationUsingPairsPartitioned(
        PairIntArray set1, PairIntArray set2,
        float toleranceTransX, float toleranceTransY,
        boolean useGreedyMatching) {

        /*
        Need to be able to sample only a small part of the dataset and arrive
        at a good answer.

        Will keep all of set2 and shuffle set1 to make sure it isn't ordered.

        Will only visit the first n1Min of set1 and break from the loop.
        n1Min is such that the complexity is kept < 1e9.
            n1Min*(n1Min-1)*n2*(n2-1)*(1./4.)*n1*n2 < 1e9
            n1Min*(n1Min-1) = 1e9/(*n2*(n2-1)*(1./4.)*n1*n2)
            n1Min ~ Math.sqrt(1e9/(n2*(n2-1)*(1./4.)*n1*n2))
        */

        int n1 = set1.getN();
        int n2 = set2.getN();
        PairIntArray shuffledSet1 = null;
        try {
            shuffledSet1 = MiscMath.shuffle(set1);
        } catch (NoSuchAlgorithmException ex) {
            throw new RuntimeException(ex);
        }

        long n1Min = Math.round(Math.sqrt(1e9/(n2*(n2-1)*(1./2.)*n1*n2)));
        if (n1Min < 10) {
            //this can still take minutes if n1>150 and n2 > 150
            n1Min = 10;
        }

        log.info("n1Min=" + n1Min);

        return calculateEuclideanTransformationForSmallSets(shuffledSet1,
            set2, shuffledSet1, set2, n1Min,
            toleranceTransX, toleranceTransY, useGreedyMatching);
    }

    protected List<TransformationPointFit>
    calculateEuclideanTransformationUsingPairsPartitioned(
        PairIntArray set1, PairIntArray set2,
        final float scale,
        final float rotationLowLimitInDegrees,
        final float rotationHighLimitInDegrees,
        boolean useLargestToleranceForOutput, boolean useGreedyMatching) {

        int n1 = set1.getN();
        int n2 = set2.getN();

        PointPartitioner partitioner = new PointPartitioner();

        List<PairIntArray> set2Subsets = partitioner.randomSubsets(set2,
            largeSearchLimit);

        int k = 2;

        int[] selected1 = new int[k];
        int[] selected2 = new int[k];

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();
        float toleranceTransX;
        float toleranceTransY;
        if (useLargestToleranceForOutput) {
            toleranceTransX = generalTolerance * (float)Math.sqrt(1./2);
            toleranceTransY = toleranceTransX;
        } else {
            toleranceTransX = 4;
            toleranceTransY = 4;
        }
        double centroidX1 = 0;
        double centroidY1 = 0;

        int maxNMatchable = Math.min(n1, n2);

        float rotRange = AngleUtil.getAngleDifference(rotationLowLimitInDegrees,
            rotationHighLimitInDegrees);

        SubsetChooser s1 = new SubsetChooser(n1, k);

        TransformationPointFit bestFitNormalized = null;

        List<TransformationPointFit> similarToBestFitNormalized = new ArrayList<TransformationPointFit>();

        TransformationPointFit bestFit = null;

        List<TransformationPointFit> similarToBestFit = new ArrayList<TransformationPointFit>();

        while (s1.getNextSubset(selected1) != -1) {

            for (PairIntArray set2Subset : set2Subsets) {

                int n2S = set2Subset.getN();

                SubsetChooser s2 = new SubsetChooser(n2S, k);

                while (s2.getNextSubset(selected2) != -1) {

                    TransformationParameters params = tc.calulateEuclidean(
                        set1.getX(selected1[0]), set1.getY(selected1[0]),
                        set1.getX(selected1[1]), set1.getY(selected1[1]),
                        set2Subset.getX(selected2[0]), set2Subset.getY(selected2[0]),
                        set2Subset.getX(selected2[1]), set2Subset.getY(selected2[1]),
                        centroidX1, centroidY1);

                    float rotD = params.getRotationInDegrees();

                    if (
                        (Math.abs(scale - params.getScale()) < 0.05) &&
                        (((rotD >= rotationLowLimitInDegrees)
                          && (Math.abs(AngleUtil.getAngleDifference(rotD, rotationLowLimitInDegrees)) < rotRange))
                        ||
                        ((rotD <= rotationHighLimitInDegrees)
                          && (Math.abs(AngleUtil.getAngleDifference(rotD, rotationHighLimitInDegrees)) < rotRange))
                        )) {

                        TransformationPointFit fit = evaluateForUnmatched(params,
                            set1, set2, toleranceTransX, toleranceTransY,
                            useGreedyMatching);

                        if ((fit != null) && (fit.getNumberOfMatchedPoints() > 2)) {
                            if (fitIsBetterNormalized(bestFitNormalized, fit)) {
                                if ((bestFitNormalized != null) && (bestFitNormalized.getMeanDistFromModel() < 1)
                                    && (fit.getMeanDistFromModel() < 1)) {
                                    if (similarToBestFitNormalized.isEmpty()) {
                                        similarToBestFitNormalized.add(bestFitNormalized);
                                    }
                                    similarToBestFitNormalized.add(fit);
                                } else {
                                    similarToBestFitNormalized.clear();
                                }
                                bestFitNormalized = fit;
                            }
                            if (fitIsBetter(bestFit, fit)) {
                                if ((bestFit != null) && (bestFit.getMeanDistFromModel() < 1)
                                    && (fit.getMeanDistFromModel() < 1)) {
                                    if (similarToBestFit.isEmpty()) {
                                        similarToBestFit.add(bestFit);
                                    }
                                    similarToBestFit.add(fit);
                                } else {
                                    similarToBestFit.clear();
                                }
                                bestFit = fit;
                            }
                        }

                        params = tc.calulateEuclidean(
                            set1.getX(selected1[0]), set1.getY(selected1[0]),
                            set1.getX(selected1[1]), set1.getY(selected1[1]),
                            set2Subset.getX(selected2[1]), set2Subset.getY(selected2[1]),
                            set2Subset.getX(selected2[0]), set2Subset.getY(selected2[0]),
                            centroidX1, centroidY1);

                        rotD = params.getRotationInDegrees();

                        if (
                            (Math.abs(scale - params.getScale()) < 0.05) &&
                            (((rotD >= rotationLowLimitInDegrees)
                              && (Math.abs(AngleUtil.getAngleDifference(rotD, rotationLowLimitInDegrees)) < rotRange))
                            ||
                            ((rotD <= rotationHighLimitInDegrees)
                              && (Math.abs(AngleUtil.getAngleDifference(rotD, rotationHighLimitInDegrees)) < rotRange))
                            )) {

                            fit = evaluateForUnmatched(params, set1, set2,
                                toleranceTransX, toleranceTransY, useGreedyMatching);

                            if ((fit != null) && (fit.getNumberOfMatchedPoints() > 2)) {
                                if (fitIsBetterNormalized(bestFitNormalized, fit)) {
                                    if ((bestFitNormalized != null) && (bestFitNormalized.getMeanDistFromModel() < 1)
                                        && (fit.getMeanDistFromModel() < 1)) {
                                        if (similarToBestFitNormalized.isEmpty()) {
                                            similarToBestFitNormalized.add(bestFitNormalized);
                                        }
                                        similarToBestFitNormalized.add(fit);
                                    } else {
                                        similarToBestFitNormalized.clear();
                                    }
                                    bestFitNormalized = fit;
                                }
                                if (fitIsBetter(bestFit, fit)) {
                                    if ((bestFit != null) && (bestFit.getMeanDistFromModel() < 1)
                                        && (fit.getMeanDistFromModel() < 1)) {
                                        if (similarToBestFit.isEmpty()) {
                                            similarToBestFit.add(bestFit);
                                        }
                                        similarToBestFit.add(fit);
                                    } else {
                                        similarToBestFit.clear();
                                    }
                                    bestFit = fit;
                                }
                            }
                        }
                    }
                }
            }
        }

        log.info("similarToBestFit.size=" + similarToBestFit.size());
        if (similarToBestFit.size() > 1) {
            for (TransformationPointFit fit : similarToBestFit) {
                log.info("similarFit=" + fit.toString());
            }
        }
        log.info("similarToBestFitNormalized.size=" + similarToBestFitNormalized.size());
        if (similarToBestFitNormalized.size() > 1) {
            for (TransformationPointFit fit : similarToBestFitNormalized) {
                log.info("similarFit=" + fit.toString());
            }
        } else if (similarToBestFitNormalized.isEmpty() && (bestFitNormalized != null)) {
            log.info("bestFitNormalized=" + bestFitNormalized.toString());
            similarToBestFitNormalized.add(bestFitNormalized);
        }
        if (bestFit != null) {
            log.info("bestFit(not normalized)=" + bestFit.toString());
            similarToBestFitNormalized.add(0, bestFit);
        }

        return similarToBestFitNormalized;
    }

    protected long numberOfPairPermutations(int n1, int n2) {

        long n1p = (MiscMath.computeNDivNMinusK(n1, 2)) >> 1;

        long n2p = (MiscMath.computeNDivNMinusK(n2, 2)) >> 1;

        double lg1 = Math.log(n1p)/Math.log(2);
        double lg2 = Math.log(n2p)/Math.log(2);

        if ((lg1 + lg2) > 63) {
            return Long.MAX_VALUE;
        }

        return n1p * n2p;
    }

    protected long estimateNStepsPairCalculation(int n1, int n2) {

        long np = numberOfPairPermutations(n1, n2);

        if (np == Long.MAX_VALUE) {
            return np;
        }

        long nPerFitGreedy = n1 * n2;

        return np * 2 * nPerFitGreedy;
    }

    protected void reduceToIntersection(
        PairFloatArray set1, PairIntArray set2,
        PairFloatArray outputSet1, PairIntArray outputSet2,
        float tolTransX, float tolTransY) {

        /*
        determine minima and maxima of outputSet2

        filter set1 to place in outputSet1 only the intersection with set2
           +- tolerances

        determine minima and maxima of outputSet1

        filter set2 to place in outputSet1 only the intersection with outputSet1
           +- tolerances
        */

        int minX2 = MiscMath.findMin(set2.getX(), set2.getN());
        int maxX2 = MiscMath.findMax(set2.getX(), set2.getN());
        int minY2 = MiscMath.findMin(set2.getY(), set2.getN());
        int maxY2 = MiscMath.findMax(set2.getY(), set2.getN());

        for (int i = 0; i < set1.getN(); ++i) {
            float x = set1.getX(i);
            float y = set1.getY(i);
            if ((x < (minX2 - tolTransX)) || (x > (maxX2 + tolTransX))) {
                continue;
            }
            if ((y < (minY2 - tolTransY)) || (y > (maxY2 + tolTransY))) {
                continue;
            }
            outputSet1.add(x, y);
        }

        float minX1 = MiscMath.findMin(outputSet1.getX(), outputSet1.getN());
        float maxX1 = MiscMath.findMax(outputSet1.getX(), outputSet1.getN());
        float minY1 = MiscMath.findMin(outputSet1.getY(), outputSet1.getN());
        float maxY1 = MiscMath.findMax(outputSet1.getY(), outputSet1.getN());

        for (int i = 0; i < set2.getN(); ++i) {
            float x = set2.getX(i);
            float y = set2.getY(i);
            if ((x < (minX1 - tolTransX)) || (x > (maxX1 + tolTransX))) {
                continue;
            }
            if ((y < (minY1 - tolTransY)) || (y > (maxY1 + tolTransY))) {
                continue;
            }
            outputSet2.add(Math.round(x), Math.round(y));
        }
    }

    public void match(TransformationParameters params, PairIntArray set1,
        PairIntArray set2, PairIntArray outputMatchedSet1,
        PairIntArray outputMatchedSet2, float tolTransX, float tolTransY,
        boolean useGreedyMatching) {

        Transformer transformer = new Transformer();

        PairFloatArray transformedSet1 = transformer.applyTransformation2(
            params, set1);

        if (useGreedyMatching) {

            Set<Integer> chosen = new HashSet<Integer>();
            for (int i = 0; i < transformedSet1.getN(); ++i) {
                float transformedX = transformedSet1.getX(i);
                float transformedY = transformedSet1.getY(i);
                double minDiff = Double.MAX_VALUE;
                int min2Idx = -1;
                for (int j = 0; j < set2.getN(); ++j) {
                    if (chosen.contains(Integer.valueOf(j))) {
                        continue;
                    }
                    float dx = transformedX - set2.getX(j);
                    float dy = transformedY - set2.getY(j);
                    if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                        continue;
                    }
                    float diff = (float)Math.sqrt(dx*dx + dy*dy);
                    if (diff < minDiff) {
                        minDiff = diff;
                        min2Idx = j;
                    }
                }
                if (minDiff < Double.MAX_VALUE) {
                    chosen.add(Integer.valueOf(min2Idx));
                    outputMatchedSet1.add(set1.getX(i), set1.getY(i));
                    outputMatchedSet2.add(set2.getX(min2Idx), set2.getY(min2Idx));
                }
            }

        } else {

            float[][] matchedIndexesAndDiffs = calculateMatchUsingOptimal(
                transformedSet1, set2, tolTransX, tolTransY);

            double tolerance = Math.sqrt(tolTransX*tolTransX + tolTransY*tolTransY);

            matchPoints(set1, set2, (float)tolerance,
                matchedIndexesAndDiffs, outputMatchedSet1, outputMatchedSet2);
        }
    }

    private TransformationPointFit refineAfterCalculationWithPairs(
        TransformationParameters params, PairIntArray set1, PairIntArray set2,
        boolean useGreedyMatching) {

        float rotHalfRange = 25;
        float rotDelta = 5.f;
        float transHalfRange = 75;
        float transDelta = 15;

        int nMaxMatchable = Math.min(set1.getN(), set2.getN());

        long t0 = System.currentTimeMillis();

        TransformationPointFit fit2 = null;

        if (params.getScale() >= 1.0) {
            fit2 = refineTheTransformation(
                params, set1, set2,
                rotHalfRange, rotDelta,
                transHalfRange, transDelta, transHalfRange, transDelta,
                useGreedyMatching);
        } else {
            // scale is < 1 so reverse order of sets and fit
            MatchedPointsTransformationCalculator tc = new
                MatchedPointsTransformationCalculator();
            TransformationParameters revParams = tc.swapReferenceFrames(
                params.copy());

            TransformationPointFit revFit = refineTheTransformation(
                revParams, set2, set1,
                rotHalfRange, rotDelta,
                transHalfRange, transDelta, transHalfRange, transDelta,
                useGreedyMatching);
            // reverse the parameters.

            TransformationParameters revRevParams = tc.swapReferenceFrames(
                revFit.getParameters().copy());
            fit2 = new TransformationPointFit(revRevParams,
                revFit.getNumberOfMatchedPoints(),
                revFit.getMeanDistFromModel(),
                revFit.getStDevFromMean(),
                revFit.getTranslationXTolerance(),
                revFit.getTranslationYTolerance());
            if (fit2 != null) {
               fit2.setMaximumNumberMatchable(nMaxMatchable);
            }
        }
        long t1 = System.currentTimeMillis();
        log.info("refine seconds=" + ((t1 - t0)*1e-3));

        return fit2;
    }

    public List<TransformationPointFit> calculateEuclideanTransformationUsingPairs(
        ImageExt img1, ImageExt img2,
        PairIntArray points1, PairIntArray points2,
        float toleranceTransX, float toleranceTransY, double toleranceColor,
        boolean earlyConvergeReturn, boolean useGreedyMatching) {

        // similar to calculateEuclideanTransformationUsingPairs,
        // but additionally uses the region surrounding the
        // points to filter out some calculations

        /*
        before calculate transformation with a pair,
           -- compare the 2 points 20 neighbor summed intensity and discard
              if too different.
        after calculate transformation with a pair:
           -- compare the center pair's spokes at 0, 90, 180, 270 in the
              common reference frame and discard if too different.
              comparison should probably be total intensity.
           -- during evaluation:
                -- for each point within tolerance of greedy match,
                   -- compare the summed intensity
                   -- compare the spokes
                   only consider matched if those comparisons have similar
                   results

        for the spoke intensities, can memoize for the solid angle depending upon
        length of the spoke.
        */

        int n1 = points1.getN();
        int n2 = points2.getN();

        int k = 2;

        SubsetChooser s1 = new SubsetChooser(n1, k);

        int[] selected1 = new int[k];
        int[] selected2 = new int[k];

        TransformationPointFit bestFitNormalized = null;
        TransformationPointFit bestFit = null;

        Transformer transformer = new Transformer();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        while (s1.getNextSubset(selected1) != -1) {

            int idx10 = selected1[0];
            int idx11 = selected1[1];

            if (idx11 < idx10) {
                int swap = idx10;
                idx10 = idx11;
                idx11 = swap;
            }

            int x10 = points1.getX(idx10);
            int y10 = points1.getY(idx10);

            int x11 = points1.getX(idx11);
            int y11 = points1.getY(idx11);

            SubsetChooser s2 = new SubsetChooser(n2, k);

            while (s2.getNextSubset(selected2) != -1) {

                int idx20 = selected2[0];
                int idx21 = selected2[1];
                int x20 = points2.getX(idx20);
                int y20 = points2.getY(idx20);

                int x21 = points2.getX(idx21);
                int y21 = points2.getY(idx21);

                for (int order = 0; order < 2; ++order) {
                    if (order == 1) {
                        int swap = x20;
                        x20 = x21;
                        x21 = swap;
                        swap = y20;
                        y20 = y21;
                        y21 = swap;
                    }
                    
                    // TODO: might need a small amount of dithering to allow
                    // for finding the better centroid
                    
                    TransformationParameters params = tc.calulateEuclidean(
                        x10, y10, x11, y11, x20, y20, x21, y21, 0, 0);
                    
                    // TODO: memoize transformed patches and discarded transformations

                    boolean firstPointIsSimilar = patchesAreSimilar(params, img1, img2,
                        transformer, x10, y10);
                    
                    if (!firstPointIsSimilar) {
                        continue;
                    }

                    boolean secondPointIsSimilar = patchesAreSimilar(params, img1, img2,
                        transformer, x11, y11);

                    if (!secondPointIsSimilar) {
                        continue;
                    }

                    TransformationPointFit fit = evaluateForUnmatchedGreedy(
                        params, points1, points2, img1, img2,
                        toleranceTransX, toleranceTransY);

                    if ((fit == null) || (fit.getNumberOfMatchedPoints() < 2)) {
                        continue;
                    }

                    if (fitIsBetterNormalized(bestFitNormalized, fit)) {
                        bestFitNormalized = fit;
                    }
                    if (fitIsBetter(bestFit, fit)) {
                        bestFit = fit;
                    }
                }
            }
        }

        List<TransformationPointFit> results = new ArrayList<TransformationPointFit>();
        if (bestFitNormalized != null) {
            results.add(bestFitNormalized);
        }
        if (bestFit != null) {
            results.add(bestFit);
        }

        return results;
    }

     /**
     *
     * @param points
     * @param img
     * @return two dimensional arrays of the averaged properties of the
     * neighboring pixels to each in points.
     * output first dimension size is size of points array, 2nd dimension
     * size is size of neighbors + 1 (which is54 if use4Neighbors is chosen, else 9).
     * the first entry is for the center pixel which is the pixel in the
     * points array.
     */
    PixelColors[][] avgOfNeighborIntensities(PairIntArray points,
        ImageExt img, boolean use4Neighbors) {

        final int nNeighbors;
        int[] ncXOffsets;
        int[] ncYOffsets;
        if (use4Neighbors) {
            nNeighbors = 5;
            ncXOffsets = new int[]{0, -1, 0,  1, 0};
            ncYOffsets = new int[]{0,  0, -1, 0, 1};
        } else {
            nNeighbors = 9;
            ncXOffsets = new int[]{0, -1, -1,  0,  1, 1, 1, 0, -1};
            ncYOffsets = new int[]{0,  0, -1, -1, -1, 0, 1, 1,  1};
        }

        PixelColors[][] output = new PixelColors[points.getN()][];

        PairIntArray neighborOffsets = MiscMath.get20NeighborOffsets();

        for (int i = 0; i < points.getN(); ++i) {

            output[i] = new PixelColors[nNeighbors];

            final int x0 = points.getX(i);
            final int y0 = points.getY(i);

            for (int nIdx = 0; nIdx < ncXOffsets.length; ++nIdx) {
                int x = x0 + ncXOffsets[nIdx];
                int y = y0 + ncYOffsets[nIdx];

                PixelColors clr = avgOfNeighborIntensities(x, y, neighborOffsets,
                    img);

                output[i][nIdx] = clr;
            }
        }

        return output;
    }

    PixelColors[] avgOfNeighborIntensities(PairIntArray points, ImageExt img) {

        PixelColors[] output = new PixelColors[points.getN()];

        PairIntArray neighborOffsets = MiscMath.get20NeighborOffsets();

        for (int i = 0; i < points.getN(); ++i) {

            int x = points.getX(i);
            int y = points.getY(i);

            PixelColors clr = avgOfNeighborIntensities(x, y, neighborOffsets, img);

            output[i] = clr;
        }

        return output;
    }

    PixelColors avgOfNeighborIntensities(int x, int y, PairIntArray neighborOffsets,
        ImageExt img) {

        img.setRadiusForPopulateOnDemand(0);

        int colMax = img.getWidth() - 1;
        int rowMax = img.getHeight()- 1;

        int sumR = 0;
        int sumG = 0;
        int sumB = 0;
        float sumCIEX = 0;
        float sumCIEY = 0;
        int count = 0;
        for (int j = 0; j < neighborOffsets.getN(); ++j) {
            int xOff = neighborOffsets.getX(j);
            int x2 = x + xOff;
            if ((x2 < 0) || (x2 > colMax)) {
                continue;
            }

            int yOff = neighborOffsets.getY(j);
            int y2 = y + yOff;
            if ((y2 < 0) || (y2 > rowMax)) {
                continue;
            }
            int idx = img.getInternalIndex(x2, y2);
            sumR += img.getR(idx);
            sumG += img.getG(idx);
            sumB += img.getB(idx);
            sumCIEX += img.getCIEX(idx);
            sumCIEY += img.getCIEY(idx);

            count++;
        }

        sumR = Math.round((float)sumR/(float)count);
        sumG = Math.round((float)sumG/(float)count);
        sumB = Math.round((float)sumB/(float)count);
        sumCIEX = sumCIEX/(float)count;
        sumCIEY = sumCIEY/(float)count;

        PixelColors pixClr = new PixelColors(sumR, sumG, sumB, sumCIEX,
            sumCIEY);

        return pixClr;
    }

    private TransformationPointFit evaluateForUnmatchedGreedy(
        TransformationParameters params,
        PairIntArray points1, PairIntArray points2,
        ImageExt img1, ImageExt img2,
        float toleranceTransX, float toleranceTransY) {

        Transformer transformer = new Transformer();

        PairFloatArray transformed1 = transformer.applyTransformation2(
            params, points1);

        PairFloatArray trFiltered1 = new PairFloatArray();
        PairIntArray filtered2 = new PairIntArray();

        reduceToIntersection(transformed1, points2, trFiltered1, filtered2,
            toleranceTransX, toleranceTransY);

        int n = trFiltered1.getN();

        Set<Integer> chosen = new HashSet<Integer>();

        double[] diffs = new double[n];
        int nMatched = 0;
        double avg = 0;

        for (int i = 0; i < n; i++) {

            float transformedX = trFiltered1.getX(i);
            float transformedY = trFiltered1.getY(i);

            double minDiff = Double.MAX_VALUE;
            int min2Idx = -1;

            for (int j = 0; j < filtered2.getN(); j++) {

                if (chosen.contains(Integer.valueOf(j))) {
                    continue;
                }

                float dx = transformedX - filtered2.getX(j);
                float dy = transformedY - filtered2.getY(j);

                if ((Math.abs(dx) > toleranceTransX) ||
                    (Math.abs(dy) > toleranceTransY)) {
                    continue;
                }
                
                boolean looksSimilar = patchesAreSimilar(params, img1, img2, 
                    transformer, points1.getX(i), points1.getY(i));

                if (!looksSimilar){
                    continue;
                }

                float diff = (float)Math.sqrt(dx*dx + dy*dy);

                if (diff < minDiff) {
                    minDiff = diff;
                    min2Idx = j;
                }
            }

            if (minDiff < Double.MAX_VALUE) {
                diffs[nMatched] = minDiff;
                nMatched++;
                chosen.add(Integer.valueOf(min2Idx));
                avg += minDiff;
            }
        }

        avg = (nMatched == 0) ? Double.MAX_VALUE :
            avg / (double)nMatched;

        double stDev = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffs[i] - avg;
            stDev += (d * d);
        }

        stDev = (nMatched == 0) ? Double.MAX_VALUE :
            Math.sqrt(stDev/((double)nMatched - 1.));

        TransformationPointFit fit = new TransformationPointFit(params.copy(),
            nMatched, avg, stDev, toleranceTransX, toleranceTransY);

        int nMaxMatchable = Math.min(trFiltered1.getN(), filtered2.getN());

        fit.setMaximumNumberMatchable(nMaxMatchable);

        return fit;
    }

    private boolean containsOneSimilarWithinTolerance(PixelColors clr1,
        PixelColors[] clr2, double eps) {

        boolean oneIsSimlar = false;

        double cieXEps = 0.012;
        double cieYEps = 0.009;

        for (int i = 0; i < clr2.length; ++i) {

            int diffR = Math.abs(clr1.getRed() - clr2[i].getRed());
            int diffG = Math.abs(clr1.getGreen() - clr2[i].getGreen());
            int diffB = Math.abs(clr1.getBlue() - clr2[i].getBlue());
            float diffCIEX = Math.abs(clr1.getCIEX() - clr2[i].getCIEX());
            float diffCIEY = Math.abs(clr1.getCIEY() - clr2[i].getCIEY());

            if ((diffR > eps) || (diffG > eps) || (diffB > eps) ||
                (diffCIEX > cieXEps) || (diffCIEY > cieYEps)) {
                continue;
            }

            oneIsSimlar = true;

            break;
        }

        return oneIsSimlar;
    }

    boolean patchesAreSimilar(TransformationParameters params,
        ImageExt img1, ImageExt img2, Transformer transformer,
        final int x1, final int y1) {

        final int d = 2;

        // for each pixel (x1, y1) +-d transform the point to image 2
        // and compare the differences.
        // -- compare by division
        // -- compare by sum of block

        int count = 0;
        //int sumR = 0;
        //int sumG = 0;
        //int sumB = 0;
        
        double sumRDiv = 0;
        double sumGDiv = 0;
        double sumBDiv = 0;
        double sumCieXDiv = 0;
        double sumCieYDiv = 0;

        for (int x = (x1 - d); x <= (x1 + d); ++x) {
            if ((x < 0) || (x > (img1.getWidth() - 1))) {
                continue;
            }
            for (int y = (y1 - d); y <= (y1 + d); ++y) {
                if ((y < 0) || (y > (img1.getHeight() - 1))) {
                    continue;
                }

                int idx1 = img1.getInternalIndex(x, y);

                double[] xyT = transformer.applyTransformation(params, x, y);

                // compare to (x2, y2) given as a check
                int x2T = (int)Math.round(xyT[0]);
                int y2T = (int)Math.round(xyT[1]);
                
                if ((x2T < 0) || (x2T > (img2.getWidth() - 1)) ||
                    (y2T < 0) || (y2T > (img2.getHeight() - 1))) {
                    continue;
                }
                
                if (debug) {
                    if (y == y1 && x == x1) {
                        log.info(String.format("(%d,%d) transformed to img2 (%d,%d)",
                            x1, y1, x2T, y2T));
                    }
                }

                float cieX1 = img1.getCIEX(idx1);
                float cieY1 = img1.getCIEY(idx1);
                int r1 = img1.getR(idx1);
                int g1 = img1.getG(idx1);
                int b1 = img1.getB(idx1);
                
                int idx2 = img2.getInternalIndex(x2T, y2T);
                float cieX2 = img2.getCIEX(idx2);
                float cieY2 = img2.getCIEY(idx2);
                int r2 = img2.getR(idx2);
                int g2 = img2.getG(idx2);
                int b2 = img2.getB(idx2);

                float cieXDiv = cieX1/cieX2;
                float cieYDiv = cieY1/cieY2;
                float rDiv = (float)r1/(float)r2;
                float gDiv = (float)g1/(float)g2;
                float bDiv = (float)b1/(float)b2;

                /*
                int rDiff = Math.abs(r1 - r2);
                int gDiff = Math.abs(g1 - g2);
                int bDiff = Math.abs(b1 - b2);

                log.info("cieXDiv=+ " + cieXDiv + " cieYDiv=" + cieYDiv
                    + " rDiv=" + rDiv + " gDiv=" + gDiv + " bDiv=" + bDiv
                    + " rDiff=" + rDiff + " gDiff=" + gDiff + " bDiff=" + bDiff
                );

                sumR += rDiff;
                sumG += gDiff;
                sumB += bDiff;
                */

                sumRDiv += rDiv;
                sumGDiv += gDiv;
                sumBDiv += bDiv;
                
                sumCieXDiv += cieXDiv;
                sumCieYDiv += cieYDiv;
                
                count++;
            }
        }
        
        if (count == 0) {
            return false;
        }

        /*float rAvg = (float)sumR/(float)count;
        float gAvg = (float)sumG/(float)count;
        float bAvg = (float)sumB/(float)count;
        */
        
        double rDivAvg = sumRDiv/(double)count;
        double gDivAvg = sumGDiv/(double)count;
        double bDivAvg = sumBDiv/(double)count;
        
        double cieXDivAvg = sumCieXDiv/(double)count;
        double cieYDivAvg = sumCieYDiv/(double)count;
        
        if (debug) {
            log.info("rDivAvg=+ " + rDivAvg + " gDivAvg=" + gDivAvg + 
                " bDivAvg=" + bDivAvg 
                + " cieXDivAvg=" + cieXDivAvg + " cieYDivAvg=" + cieYDivAvg);
        }

        if ((rDivAvg < 0.84) || (gDivAvg < 0.84) || (bDivAvg < 0.84)) {
            return false;
        }
        
        if ((rDivAvg > 1.19) || (gDivAvg > 1.19) || (bDivAvg > 1.19)) {
            return false;
        }

        return true;
    }

}
