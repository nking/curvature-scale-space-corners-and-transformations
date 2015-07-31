package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairFloatArrayUnmodifiable;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.List;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import org.ejml.simple.*;
import static org.junit.Assert.fail;

/**
 *
 * @author nichole
 */
public class PointMatcher3Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public PointMatcher3Test() {
    }

    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    public void estCalculateEuclideanTransformation() throws Exception {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1437335464716L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        //image size range 250 to 10000
        //  these are altered below
        int imageWidth = 650;
        int imageHeight = 400;

        float scale = 1.0f;

        boolean spacer = true;
        
        for (int nn = 0; nn < 1; ++nn) { // repeat number of tests
        for (int rotType = 0; rotType < 5; ++rotType) {
            for (int nTest = 0; nTest < 7 /*20*/; ++nTest) { // this increases nPoints

                PointMatcher pointMatcher = new PointMatcher();
                
                PairIntArray unmatchedLeftXY = new PairIntArray();

                float rotInDegrees = (int)(sr.nextFloat() * 20);

                if (rotType == 1) {
                    rotInDegrees = sr.nextBoolean() ? (90 + rotInDegrees) :
                        (90 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(5000) + 250;
                    imageHeight = Math.abs(5000) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                } else if (rotType == 2) {
                    rotInDegrees = sr.nextBoolean() ? (180 + rotInDegrees) :
                        (180 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(250) + 250;
                    imageHeight = Math.abs(250) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                } else if (rotType == 3) {
                    rotInDegrees = sr.nextBoolean() ? (270 + rotInDegrees) :
                        (270 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(1000) + 250;
                    imageHeight = Math.abs(1000) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                } else if (rotType == 4) {
                    rotInDegrees = (360 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(500) + 250;
                    imageHeight = Math.abs(500) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                }

                int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
                int transY = (int)(0.05f * sr.nextFloat() * imageHeight);
                if (sr.nextBoolean()) {
                    transX *= -1;
                }
                if (sr.nextBoolean()) {
                    transY *= -1;
                }

                TransformationParameters params = new TransformationParameters();
                params.setRotationInDegrees(rotInDegrees);
                params.setScale(scale);
                params.setTranslationX(transX);
                params.setTranslationY(transY);
                params.setOriginX(0);
                params.setOriginY(0);

                int nPoints = (nTest + 1) * 7;

                log.info("test for nPoints=" + nPoints + " nn=" + nn 
                    + " rotType=" + rotType + " nTest=" + nTest
                    + "\nparams=" + params.toString());
                
                for (int i = 0; i < nPoints; ++i) {
                    int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                    int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                    unmatchedLeftXY.add(x, y);
                }

                // ===  transform the right points  ======

                Transformer transformer = new Transformer();

                PairIntArray unmatchedRightXY = transformer.applyTransformation(
                    params, unmatchedLeftXY);

                // consider adding small deviations to the right

                // add small number of points in left that are not in right
                //   and vice versa
                int nAdd = (int)(sr.nextFloat()*nPoints/10);
                if (nAdd == 0) {
                    nAdd = 1;
                }
                for (int i = 0; i < nAdd; ++i) {
                    int x = sr.nextInt(imageWidth);
                    int y = sr.nextInt(imageHeight);
                    unmatchedLeftXY.add(x, y);
                }
                for (int i = 0; i < nAdd; ++i) {
                    int x = sr.nextInt(imageWidth);
                    int y = sr.nextInt(imageHeight);
                    unmatchedRightXY.add(x, y);
                }

                // TODO: consider scrambling the order of the right points

                // --- TODO: in the difference between the left and right regions,
                //     need to generate points in the right

                int nMaxMatchable = nPoints;

                int nExpected = nMaxMatchable;
                int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);
                                
                long t0 = System.currentTimeMillis();
                
                boolean earlyConvergeReturn = true;
                
                List<TransformationPointFit> fits = 
                    pointMatcher.calculateEuclideanTransformationUsingPairs(
                    unmatchedLeftXY, unmatchedRightXY, earlyConvergeReturn);
                
                long t1 = System.currentTimeMillis();
                
                double timeSec = (t1 - t0) * 1e-3;
                
                assert(fits != null && !fits.isEmpty());
                
                log.info("fit for best pairwise calculation seconds=" + timeSec);
                
TransformationPointFit fit = fits.get(0);

                TransformationParameters fitParams = fit.getParameters();
                int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
                float diffRotDeg = AngleUtil.getAngleDifference(
                    fitParams.getRotationInDegrees(), rotInDegrees);
                float diffScale = Math.abs(fitParams.getScale() - scale);
                float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                String space = spacer ? "    " : "";
                log.info(space + "diff =" +
                    String.format(
                    " dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d dNPoints=%d meanDiffModel=%f",
                    diffRotDeg, diffScale, diffTransX, 
                    diffTransY, nEps, diffN, (float)fit.getMeanDistFromModel()));
            }
        }
        }
        
    }

    public void estCalculateEuclideanTransformation2()
        throws Exception {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1436749125398L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        //image size range 250 to 10000
        //  these are altered below
        int imageWidth = 650;
        int imageHeight = 400;

        // ---- random testing for stereo imaging ----

        float scale = 1.0f;

        for (int nRuns = 0; nRuns < 1; ++nRuns) { // this increases the number of tests
            for (int rotType = 0; rotType < 2; ++rotType) {
                for (int nTest = 0; nTest < 5/*20*/; ++nTest) { // this increases nPoints

                    PointMatcher pointMatcher = new PointMatcher();

                    PairIntArray unmatchedLeftXY = new PairIntArray();

                    float rotInDegrees = (int)(sr.nextFloat() * 20);

                    if (rotType == 1) {

                        rotInDegrees = 360 - rotInDegrees;
                        imageWidth = Math.abs(500) + 250;
                        imageHeight = Math.abs(500) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 2) {

                        rotInDegrees = 360 - rotInDegrees;
                        imageWidth = Math.abs(5000) + 250;
                        imageHeight = Math.abs(5000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    }

                    int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
                    int transY = (int)(0.05f * sr.nextFloat() * imageHeight);
                    if (sr.nextBoolean()) {
                        transX *= -1;
                    }
                    if (sr.nextBoolean()) {
                        transY *= -1;
                    }

                    TransformationParameters params = new TransformationParameters();
                    params.setRotationInDegrees(rotInDegrees);
                    params.setScale(scale);
                    params.setTranslationX(transX);
                    params.setTranslationY(transY);

                    int nPoints = (nTest + 1) * 7;

                    log.info("\ntest for nPoints=" + nPoints + " nTest=" + nTest
                        + " rotType=" + rotType + " nRuns=" + nRuns
                        + "\nparams=" + params.toString()
                    );

                    for (int i = 0; i < nPoints; ++i) {
                        int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                        int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                        unmatchedLeftXY.add(x, y);
                    }

                    // ===  transform the right points  ======

                    Transformer transformer = new Transformer();

                    PairIntArray unmatchedRightXY = transformer.applyTransformation(
                        params, unmatchedLeftXY);

                    // consider adding small deviations to the right

                    // add small number of points in left that are not in right
                    //   and vice versa
                    int nAdd = (int)(sr.nextFloat()*nPoints/10);
                    if (nAdd == 0) {
                        nAdd = 1;
                    }
                    for (int i = 0; i < nAdd; ++i) {
                        int x = sr.nextInt(imageWidth);
                        int y = sr.nextInt(imageHeight);
                        unmatchedLeftXY.add(x, y);
                    }
                    for (int i = 0; i < nAdd; ++i) {
                        int x = sr.nextInt(imageWidth);
                        int y = sr.nextInt(imageHeight);
                        unmatchedRightXY.add(x, y);
                    }

                    // change the order of points to make sure not biased
                    // by order
                    unmatchedRightXY.reverse();

                    float setsFractionOfImage = 1.0f;

                    int nMaxMatchable = nPoints;

                    int nExpected = nMaxMatchable;
                    int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);

                    TransformationPointFit fit =
                        pointMatcher.calculateEuclideanTransformation(
                        unmatchedLeftXY, unmatchedRightXY);

                    assert(fit != null);

                    log.info("FIT=" + fit.toString());

                    boolean converged = pointMatcher.hasConverged(fit, nMaxMatchable);

                    TransformationParameters fitParams = fit.getParameters();
                    int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
                    float diffRotDeg = AngleUtil.getAngleDifference(
                        fitParams.getRotationInDegrees(), rotInDegrees);
                    float diffScale = Math.abs(fitParams.getScale() - scale);
                    float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                    float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                    double epsTrans = 3;
                    double epsRot = 3;

                    log.info("FINAL FIT=" + fit.toString());
                    log.info("diff result and expected =" +
                        String.format(
                        " dNPoints=%d, dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d  tEps=%f",
                        diffN, diffRotDeg, diffScale, diffTransX, diffTransY, nEps, (float)epsTrans)
                    );

                    if ((diffN <= nEps) && (diffRotDeg <= epsRot) &&
                        (diffScale < 0.2) && (diffTransX <= epsTrans) &&
                        (diffTransY <= epsTrans)) {
                        converged = true;
                    }

                    assertTrue(converged);
                }
            }
        }
    }

    private void debugPlot(PairIntArray set0,
        PairIntArray set1, String label) {

        int minX0 = MiscMath.findMin(set0.getX(), set0.getN());
        int maxX0 = MiscMath.findMax(set0.getX(), set0.getN());
        int minY0 = MiscMath.findMin(set0.getY(), set0.getN());
        int maxY0 = MiscMath.findMax(set0.getY(), set0.getN());

        int minX1 = MiscMath.findMin(set1.getX(), set1.getN());
        int maxX1 = MiscMath.findMax(set1.getX(), set1.getN());
        int minY1 = MiscMath.findMin(set1.getY(), set1.getN());
        int maxY1 = MiscMath.findMax(set1.getY(), set1.getN());

        float minX = Math.min(minX0, minX1);
        float maxX = Math.max(maxX0, maxX1);

        float minY = Math.min(minY0, minY1);
        float maxY = Math.max(maxY0, maxY1);

        try {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
                minX, maxX, minY, maxY);

            plotter.addPlot(set0.getX(), set0.getY(), set1.getX(), set1.getY(), label);

            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());

        } catch (IOException e) {

        }
    }

    private void debugPlot(PairFloatArrayUnmodifiable set0,
        PairIntArray set1, String label) {

        float minX0 = MiscMath.findMin(set0.getX(), set0.getN());
        float maxX0 = MiscMath.findMax(set0.getX(), set0.getN());
        float minY0 = MiscMath.findMin(set0.getY(), set0.getN());
        float maxY0 = MiscMath.findMax(set0.getY(), set0.getN());

        int minX1 = MiscMath.findMin(set1.getX(), set1.getN());
        int maxX1 = MiscMath.findMax(set1.getX(), set1.getN());
        int minY1 = MiscMath.findMin(set1.getY(), set1.getN());
        int maxY1 = MiscMath.findMax(set1.getY(), set1.getN());

        float minX = Math.min(minX0, minX1);
        float maxX = Math.max(maxX0, maxX1);

        float minY = Math.min(minY0, minY1);
        float maxY = Math.max(maxY0, maxY1);

        try {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
                minX, maxX, minY, maxY);

            plotter.addPlot(set0.getX(), set0.getY(), set1.getX(), set1.getY(),
                label);

            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());

        } catch (IOException e) {

        }
    }

    protected class DensityTranslationResults {
        final double density;
        final int nPoints;
        final int imageWidth;
        final int imageHeight;
        final double translationDelta;
        public DensityTranslationResults(double theDensity, int numberOfPoints,
            int theImageWidth, int theImageHeight, double theTranslationDelta) {
            density = theDensity;
            nPoints = numberOfPoints;
            imageWidth = theImageWidth;
            imageHeight = theImageHeight;
            translationDelta = theTranslationDelta;
        }
    }

    /*
    https://github.com/jesolem/PCV
    adapted from code licensed under  BSD license (2-clause "Simplified BSD License").
    private SimpleMatrix calculateCameraMatrixFromFundamentalMatrix(SimpleMatrix fm) {
        SimpleMatrix u = fm.svd().getU();
        SimpleMatrix leftE = u.extractVector(true, 2);
        leftE = leftE.divide(u.get(2, 2));

        double[][] skewSymmetric = new double[3][];
        skewSymmetric[0] = new double[]{0, -leftE.get(0, 2), leftE.get(0, 1)};
        skewSymmetric[1] = new double[]{leftE.get(0, 2), 0, -leftE.get(0, 0)};
        skewSymmetric[2] = new double[]{-leftE.get(0, 1), leftE.get(0, 0), 0};

        SimpleMatrix sMatrix = new SimpleMatrix(skewSymmetric);
        double v = sMatrix.dot( fm.transpose() );


        //essential matrix from F: E = transpose(K2)*F*K1
        SimpleMatrix k = DataForTests.getVenturiCameraIntrinsics()
            .mult(DataForTests.getVenturiRotationMatrixForImage001());
        SimpleMatrix k2 = DataForTests.getVenturiCameraIntrinsics()
            .mult(DataForTests.getVenturiRotationMatrixForImage010());

        SimpleMatrix essentialMatrix = k2.transpose().mult(fm).mult(k);

        int z = 1;
    }*/

    public void estSkyline() throws Exception {

        String[] fileNames = new String[] {
            "brown_lowe_2003_image1.jpg",
            //"brown_lowe_2003_image1_rot.jpg",
            //"brown_lowe_2003_image2.jpg",
            /*"venturi_mountain_j6_0001.png",
            //"venturi_mountain_j6_0010.png",
            "seattle.jpg",
            "arches.jpg",
            "stinson_beach.jpg",
            "cloudy_san_jose.jpg",
            "stonehenge.jpg",
            "norwegian_mtn_range.jpg",
            "halfdome.jpg",
            "costa_rica.jpg",
            "new-mexico-sunrise_w725_h490.jpg",
            "arizona-sunrise-1342919937GHz.jpg",*/
            //"sky_with_rainbow.jpg",
            //"sky_with_rainbow2.jpg",
            //"patagonia_snowy_foreground.jpg",
            //"mt_rainier_snowy_field.jpg"
            //"klein_matterhorn_snowy_foreground.jpg"
            //"30.jpg",
            //"arches_sun_01.jpg",
            //"stlouis_arch.jpg",
            //"contrail.jpg"
        };

        for (String fileName : fileNames) {

            log.info("fileName=" + fileName);

            // revisit infl points.  is there a threshold removing points?
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
            int image1Width = img1.getWidth();
            int image1Height = img1.getHeight();

      /*
      ImageProcessor ip = new ImageProcessor();
      ip.testFilter(img1);
      if (true) {
          return;
      }
      */

            List<PairIntArray> edges1 = null;
            PairIntArray points1 = null;

            int nPreferredCorners = 100;
            int nCrit = 500;

            SkylineExtractor.setDebugName(fileName);

            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);

            detector.useOutdoorModeAndExtractSkyline();
            //detector.findCornersIteratively(nPreferredCorners, nCrit);
            SkylineExtractor.setDebugName(fileName);
            detector.findCorners();
            /*edges1 = detector.getEdgesInOriginalReferenceFrame();
            points1 = detector.getCornersInOriginalReferenceFrame();

            Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);

            for (PairIntArray edge : edges1) {
                ImageIOHelper.addCurveToImage(edge, image1, 2,
                    Color.YELLOW.getRed(), Color.YELLOW.getGreen(),
                    Color.YELLOW.getBlue());
            }

            ImageIOHelper.addCurveToImage(points1, image1, 1, 255, 0, 0);

            String dirPath = ResourceFinder.findDirectory("bin");
            String outFilePath = dirPath + "/tmp1_edges_and_corners_" +
                fileName + ".png";

            ImageIOHelper.writeOutputImage(outFilePath, image1);*/
        }
    }

    public void estIterativeCorners() throws Exception {

        /*
        tests skyline extractor, iterative corners, euclidean point matcher,
        then stereopojection calculation
        */

        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";

        /*
        String fileName1 = "venturi_mountain_j6_0001.png";
        String fileName2 = "venturi_mountain_j6_0010.png";
        */

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();

        //List<PairIntArray> edges1 = null;
        List<PairIntArray> edgesSkyline1 = null;
        //List<PairIntArray> edges2 = null;
        List<PairIntArray> edgesSkyline2 = null;

        int nPreferredCorners = 250;
        int nCrit = 200;

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
        //detector.useOutdoorMode();
        detector.useOutdoorModeAndExtractSkyline();
        //detector.calculateSkylineCornersOnly();
        detector.findCornersIteratively(nPreferredCorners, nCrit);
        //detector.findCorners();
        //edges1 = detector.getEdgesInOriginalReferenceFrame();
        PairIntArray skylineCorners1 = detector.getSkylineCornersInOriginalReferenceFrame();
        PairIntArray corners1 = detector.getCornersInOriginalReferenceFrame();
        edgesSkyline1 = detector.getSkylineEdgesInOriginalReferenceFrame();

        detector = new CurvatureScaleSpaceCornerDetector(img2);
        //detector.useOutdoorMode();
        detector.useOutdoorModeAndExtractSkyline();
        //detector.calculateSkylineCornersOnly();
        detector.findCornersIteratively(nPreferredCorners, nCrit);
        //detector.findCorners();
        //edges2 = detector.getEdgesInOriginalReferenceFrame();
        PairIntArray skylineCorners2 = detector.getSkylineCornersInOriginalReferenceFrame();
        PairIntArray corners2 = detector.getCornersInOriginalReferenceFrame();
        edgesSkyline2 = detector.getSkylineEdgesInOriginalReferenceFrame();

        Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        String dirPath = ResourceFinder.findDirectory("bin");
        Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        PairIntArray comb1 = skylineCorners1.copy();
        comb1.addAll(corners1);
        PairIntArray comb2 = skylineCorners2.copy();
        comb2.addAll(corners2);
        
        writeImage(image1.copyImage(), skylineCorners1,
            "skylinecorners1_" + ".png");
        writeImage(image2.copyImage(), skylineCorners2,
            "skylinecorners2_" + ".png");
        
        writeImage(image1.copyImage(), comb1,
            "corners1_" + ".png");
        writeImage(image2.copyImage(), comb2,
            "corners2_" + ".png");

        // ===== use partitioned matching =====
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        PointMatcher pointMatcher = new PointMatcher();
        //pointMatcher.setCostToNumMatchedAndDiffFromModel();
        //pointMatcher.setCostToDiffFromModel();

        boolean useLargestToleranceForOutput = true;
        boolean useGreedyMatching = true;

        TransformationPointFit trFit =
            pointMatcher.performMatchingForMostlyVerticalTranslation(
                comb1, comb2, 
                outputMatchedScene, outputMatchedModel, 
                useLargestToleranceForOutput, useGreedyMatching);

        log.info("rough euclidean from skyline alone=" + trFit.toString());
        
        writeTransformed(trFit.getParameters(), image2.copyImage(),
            outputMatchedScene, outputMatchedModel, 
            "transformedCorners_.png");
                
        if (outputMatchedScene.getN() < 7) {
            // no epipolar solution. need at least 7 points
            return;
        }

        PairFloatArray finalOutputMatchedScene = new PairFloatArray();
        PairFloatArray finalOutputMatchedModel = new PairFloatArray();

        RANSACSolver ransacSolver = new RANSACSolver();

        StereoProjectionTransformerFit sFit = ransacSolver
            .calculateEpipolarProjection(
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedScene),
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedModel),
                finalOutputMatchedScene, finalOutputMatchedModel);

        overplotEpipolarLines(sFit.getFundamentalMatrix(),
            outputMatchedScene.toPairFloatArray(), outputMatchedModel.toPairFloatArray(),
            ImageIOHelper.readImage(filePath1),
            ImageIOHelper.readImage(filePath2),
            image1Width,
            img1.getHeight(), img2.getWidth(), img2.getHeight());

        /*
        StereoProjectionTransformer st = new StereoProjectionTransformer();
        SimpleMatrix fm =
            st.calculateEpipolarProjectionForPerfectlyMatched(
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedScene),
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedModel));

        overplotEpipolarLines(fm,
            outputMatchedScene.toPairFloatArray(), outputMatchedModel.toPairFloatArray(),
            ImageIOHelper.readImage(filePath1),
            ImageIOHelper.readImage(filePath2),
            image1Width,
            img1.getHeight(), img2.getWidth(), img2.getHeight());
        */

        System.out.println("test done");
    }

    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height) throws IOException {

        String flNum = "";

        overplotEpipolarLines(fm, set1, set2, img1, img2, image1Width,
            image1Height, image2Width, image2Height, flNum);
    }
    
    public void testShapeMatching() throws Exception {

        /*
        tests skyline extractor, iterative corners, euclidean point matcher,
        then stereopojection calculation
        */
        String fileName1, fileName2;

        //-280, -14
        // Brown & Lowe 2003
        //fileName1 = "brown_lowe_2003_image1.jpg";
        //fileName2 = "brown_lowe_2003_image2.jpg";

        //-30, 0
        //https://venturi.fbk.eu/results/public-datasets/mountain-dataset/
        //fileName1 = "venturi_mountain_j6_0001.png";
        //fileName2 = "venturi_mountain_j6_0010.png";
        
        //-80,  NOTE images already rectified
        //2005 dataset from http://vision.middlebury.edu/stereo/data/ 
        fileName1 = "books_illum3_v0_695x555.png";
        fileName2 = "books_illum3_v6_695x555.png";
        
        // can find unrectified stereo images for testing here:
        //   http://www.cvlibs.net/datasets/kitti/raw_data.php
        //   http://www.robots.ox.ac.uk/NewCollegeData/index.php?n=Main.Downloads
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();

        // mask out the sky
        if (fileName1.contains("brown") || fileName1.contains("venturi")) {
            ImageHelperForTests helper = new ImageHelperForTests(img1, true);
            SkylineExtractor skylineExtractor = new SkylineExtractor();                    
            PairIntArray outputSkyCentroid = new PairIntArray();
            // sky are the zeros in this:
            GreyscaleImage resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img1,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);

            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.multiplyBinary(img1, resultMask);

            helper = new ImageHelperForTests(img2, true);
            skylineExtractor = new SkylineExtractor();                    
            outputSkyCentroid = new PairIntArray();
            // sky are the zeros in this:
            resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img2,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);
            imageProcessor.multiplyBinary(img2, resultMask);
        }
        
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        ShapeMatcher shapeMatcher = new ShapeMatcher();
        
        TransformationPointFit fit = shapeMatcher.findMatchingShapes(
            img1, img2, outputMatchedScene, outputMatchedModel);
if (true) {
    return;
}
        Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        /*
        PairIntArray comb1 = skylineCorners1.copy();
        comb1.addAll(corners1);
        PairIntArray comb2 = skylineCorners2.copy();
        comb2.addAll(corners2);
        
        writeImage(image1.copyImage(), skylineCorners1,
            "skylinecorners1_" + ".png");
        writeImage(image2.copyImage(), skylineCorners2,
            "skylinecorners2_" + ".png");
        
        GreyscaleImage img1Grey = img1.copyToGreyscale();
        GreyscaleImage img2Grey = img2.copyToGreyscale();
        
        imageProcessor.applyImageSegmentation(img1Grey, 3);
        imageProcessor.applyImageSegmentation(img2Grey, 3);
        
        writeImage(ImageIOHelper.convertImage(img1Grey), corners1,
            "segmentation1_" + ".png");
        writeImage(ImageIOHelper.convertImage(img2Grey), corners2,
            "segmentation2_" + ".png");
        
        writeImage(image1.copyImage(), comb1,
            "corners1_" + ".png");
        writeImage(image2.copyImage(), comb2,
            "corners2_" + ".png");
        */
        
        log.info("rough euclidean from shapes=" + fit.toString());
        
        writeTransformed(fit.getParameters(), image2.copyImage(),
            outputMatchedScene, outputMatchedModel, 
            "transformedPoints_.png");
                
        if (outputMatchedScene.getN() < 7) {
            // no epipolar solution. need at least 7 points
            return;
        }

        PairFloatArray finalOutputMatchedScene = new PairFloatArray();
        PairFloatArray finalOutputMatchedModel = new PairFloatArray();

        RANSACSolver ransacSolver = new RANSACSolver();

        StereoProjectionTransformerFit sFit = ransacSolver
            .calculateEpipolarProjection(
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedScene),
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedModel),
                finalOutputMatchedScene, finalOutputMatchedModel);

        overplotEpipolarLines(sFit.getFundamentalMatrix(),
            outputMatchedScene.toPairFloatArray(), outputMatchedModel.toPairFloatArray(),
            ImageIOHelper.readImage(filePath1),
            ImageIOHelper.readImage(filePath2),
            image1Width,
            img1.getHeight(), img2.getWidth(), img2.getHeight());

        /*
        StereoProjectionTransformer st = new StereoProjectionTransformer();
        SimpleMatrix fm =
            st.calculateEpipolarProjectionForPerfectlyMatched(
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedScene),
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedModel));

        overplotEpipolarLines(fm,
            outputMatchedScene.toPairFloatArray(), outputMatchedModel.toPairFloatArray(),
            ImageIOHelper.readImage(filePath1),
            ImageIOHelper.readImage(filePath2),
            image1Width,
            img1.getHeight(), img2.getWidth(), img2.getHeight());
        */

        System.out.println("test done");
    }

    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber)
        throws IOException {

        SimpleMatrix input1 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(set1);

        SimpleMatrix input2 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(set2);

        for (int ii = 0; ii < input1.numCols(); ii++) {
            double x = input1.get(0, ii);
            double y = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < input2.numCols(); ii++) {
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3,
                255, 0, 0);
        }

        StereoProjectionTransformer spTransformer = new
            StereoProjectionTransformer();

        Color clr = null;
        for (int ii = 0; ii < input2.numCols(); ii++) {
            clr = getColor(clr);
            SimpleMatrix epipolarLinesInLeft = fm.transpose().mult(input2);
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < input1.numCols(); ii++) {
            clr = getColor(clr);
            SimpleMatrix epipolarLinesInRight = fm.mult(input1);
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_1_" + outfileNumber + ".png", img1);
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_2_" + outfileNumber + ".png", img2);
    }

    public static void main(String[] args) {

        try {
            PointMatcher3Test test = new PointMatcher3Test();

            //test.testCalculateEuclideanTransformation();
            //test.testCalculateEuclideanTransformation2();
            //test.testIterativeCorners();
            test.testShapeMatching();

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    private Color getColor(Color clr) {
        if ((clr == null) || clr.equals(Color.MAGENTA)) {
            return Color.BLUE;
        }
        if (clr.equals(Color.BLUE)) {
            return Color.PINK;
        } else if (clr.equals(Color.PINK)) {
            return Color.GREEN;
        } else if (clr.equals(Color.GREEN)) {
            return Color.RED;
        } else if (clr.equals(Color.RED)) {
            return Color.CYAN;
        } else if (clr.equals(Color.CYAN)) {
            return Color.MAGENTA;
        } else if (clr.equals(Color.MAGENTA)) {
            return Color.LIGHT_GRAY;
        } else {
            return Color.ORANGE;
        }
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
}
