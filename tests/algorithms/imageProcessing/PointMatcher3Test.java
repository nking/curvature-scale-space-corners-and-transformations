package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.RangeInt;
import algorithms.util.ScatterPointPlotterPNG;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
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

    public void estPerformVerticalPartitionedMatching() throws Exception {

        PointMatcher pointMatcher = new PointMatcher();

        int imageWidth = 2048;
        int imageHeight = 1256;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1435259647733L;
        //seed = 1435616626667L;
        //seed = 1435634400691L;
        //seed = 1435697643135L;
        //seed = 1435799127020L;
        //seed = 1435802308803L;
        //seed = 1435857338349L;
        seed = 1435898464646L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        // ---- random testing for stereo imaging ----

        int numberOfPartitions = sr.nextInt(3) + 1;

        float scale = 1.0f;

        for (int rotType = 0; rotType < 4; ++rotType) {

            for (int nTest = 0; nTest < 15; ++nTest) {

                PairIntArray unmatchedLeftXY = new PairIntArray();

                int transX = (int)(0.25f * sr.nextFloat() * imageWidth);
                int transY = (int)(0.05f * sr.nextFloat() * imageHeight);
                if (sr.nextBoolean()) {
                    transX *= -1;
                }
                if (sr.nextBoolean()) {
                    transY *= -1;
                }

                float rotInDegrees = (int)(sr.nextFloat() * 20);

                if (rotType == 1) {
                    rotInDegrees = sr.nextBoolean() ? (270 + rotInDegrees) :
                        (270 - rotInDegrees);
                } else if (rotType == 2) {
                    rotInDegrees = sr.nextBoolean() ? (180 + rotInDegrees) :
                        (180 - rotInDegrees);
                } else if (rotType == 3) {
                    rotInDegrees = sr.nextBoolean() ? (90 + rotInDegrees) :
                        (90 - rotInDegrees);
                }

                TransformationParameters params = new TransformationParameters();
                params.setRotationInDegrees(rotInDegrees);
                params.setScale(scale);
                params.setTranslationX(transX);
                params.setTranslationY(transY);

                int nPoints = (nTest + 1) * 7;

                log.info("test for nPoints=" + nPoints + "\nparams=" + params.toString()
                    + "\nnumber of vertical partitions=" + numberOfPartitions);

                for (int i = 0; i < nPoints; ++i) {
                    int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                    int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                    unmatchedLeftXY.add(x, y);
                }

                // ===  transform the right points  ======

                Transformer transformer = new Transformer();

                PairIntArray unmatchedRightXY = transformer.applyTransformation(
                    params, unmatchedLeftXY, imageWidth >> 1, imageHeight >> 1);

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

/*if (nPoints != 49) {
    continue;
}*/
                PairIntArray outputMatchedLeftXY = new PairIntArray();
                PairIntArray outputMatchedRightXY = new PairIntArray();

                boolean useLargestToleranceForOutput = false;

                TransformationPointFit fit =
                    pointMatcher.performVerticalPartitionedMatching(
                    numberOfPartitions, unmatchedLeftXY, unmatchedRightXY,
                    imageWidth, imageHeight, imageWidth, imageHeight,
                    outputMatchedLeftXY, outputMatchedRightXY,
                    useLargestToleranceForOutput);

                assert(fit != null);
                TransformationParameters fitParams = fit.getParameters();
                int diffN = Math.abs(nPoints - fit.getNumberOfMatchedPoints());
                float diffRotDeg = getAngleDifference(
                    fitParams.getRotationInDegrees(), rotInDegrees);
                float diffScale = Math.abs(fitParams.getScale() - scale);
                float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                log.info("FINAL FIT=" + fit.toString());
                log.info("diff result and expected =" +
                    String.format(
                    " dNPoints=%d, dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f",
                    diffN, diffRotDeg, diffScale, diffTransX, diffTransY)
                );

                double epsTrans = 3;
                if (nPoints < 10) {
                    epsTrans = 10;
                } else if (nPoints < 40) {
                    epsTrans = 7;
                }

                assertTrue(diffN < 0.5*nPoints);
                assertTrue(diffRotDeg <= 10.0);
                assertTrue(diffScale < 0.2);
                assertTrue(diffTransX <= epsTrans);
                assertTrue(diffTransY <= epsTrans);

                numberOfPartitions++;
                if (numberOfPartitions > 3) {
                    numberOfPartitions = 1;
                }
            }
        }
    }

    public void testCalculateTranslationFromGridThenDownhillSimplex()
        throws Exception {

        PointMatcher pointMatcher = new PointMatcher();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1436162255215L;
        //seed = 1436200575218L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        //image size range 250 to 10000
        //  these are altered below
        int imageWidth = 650;
        int imageHeight = 400;

        // ---- random testing for stereo imaging ----

        int numberOfPartitions = sr.nextInt(3) + 1;

        float scale = 1.0f;

        int transXDelta = 1;
        int transYDelta = 1;
        int halfRange = 4;

        List<DensityTranslationResults> converged0 =
            new ArrayList<DensityTranslationResults>();
        List<DensityTranslationResults> converged1 =
            new ArrayList<DensityTranslationResults>();

        for (int nRuns = 0; nRuns < 10;/*50;*/ ++nRuns) {
            for (int rotType = 0; rotType < 4; ++rotType) {
                for (int nTest = 0; nTest < 15; ++nTest) {

                    PairIntArray unmatchedLeftXY = new PairIntArray();

                    float rotInDegrees = (int)(sr.nextFloat() * 20);

                    if (rotType == 1) {
                        rotInDegrees = sr.nextBoolean() ? (270 + rotInDegrees) :
                            (270 - rotInDegrees);
                        transXDelta = 2;
                        transYDelta = 2;
                        halfRange = 10;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(250) + 250;
                        imageHeight = Math.abs(250) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 2) {
                        rotInDegrees = sr.nextBoolean() ? (180 + rotInDegrees) :
                            (180 - rotInDegrees);
                        transXDelta = 3;
                        transYDelta = 3;
                        halfRange = 20;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(1000) + 250;
                        imageHeight = Math.abs(1000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 3) {
                        rotInDegrees = sr.nextBoolean() ? (90 + rotInDegrees) :
                            (90 - rotInDegrees);
                        transXDelta = 4;
                        transYDelta = 4;
                        halfRange = 100;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(2000) + 250;
                        imageHeight = Math.abs(2000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 3) {
                        rotInDegrees = sr.nextBoolean() ? (45 + rotInDegrees) :
                            (45 - rotInDegrees);
                        transXDelta = 5;
                        transYDelta = 5;
                        halfRange = 115;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(5000) + 250;
                        imageHeight = Math.abs(5000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    }

                    int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(2048)));
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

                    log.info("\ntest for nPoints=" + nPoints + "\nparams=" + params.toString()
                        + "\ntransXDelta=" + transXDelta + " transYDelta=" + transYDelta
                        + " translation range=" + (2*halfRange));

                    for (int i = 0; i < nPoints; ++i) {
                        int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                        int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                        unmatchedLeftXY.add(x, y);
                    }

                    // ===  transform the right points  ======

                    Transformer transformer = new Transformer();

                    PairIntArray unmatchedRightXY = transformer.applyTransformation(
                        params, unmatchedLeftXY, imageWidth >> 1, imageHeight >> 1);

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

    /*if (nPoints != 49) {
        continue;
    }*/
                    /*
                    private int dsNMaxIter = 50;
                    pointMatcher.setDsNMaxIter(int n);
                    float nEpsFactor = 2.0f;
                    pointMatcher.setNEpsFactor(float f);
                    protected float cellFactor = 1.25f;
                    protected float tolFactor = 0.75f;
                    pointMatcher.setCellFactor(float f);
                    pointMatcher.setTolFactor(float f);
                    */
                    float tolTransX = pointMatcher.getTolFactor() * transXDelta;
                    float tolTransY = pointMatcher.getTolFactor() * transYDelta;
                    if (tolTransX < 1) {
                        tolTransX = 1;
                    }
                    if (tolTransY < 1) {
                        tolTransY = 1;
                    }

                    boolean setsAreMatched = false;
                    float setsFractionOfImage = 1.0f;

                    RangeInt transXStartStop = new RangeInt(
                        (int)(params.getTranslationX() - halfRange),
                        (int)(params.getTranslationX() + halfRange));
                    RangeInt transYStartStop = new RangeInt(
                        (int)(params.getTranslationY() - halfRange),
                        (int)(params.getTranslationY() + halfRange));

                    PairFloatArray scaledRotatedLeft = pointMatcher.scaleAndRotate(
                        unmatchedLeftXY, params.getRotationInRadians(), params.getScale(),
                        imageWidth, imageHeight);

                    int nIter = 0;
                    int nMaxIter = 10;

                    int dsMaxIter = pointMatcher.getDsNMaxIter();
                    double dens = (double)nPoints/(double)imageWidth;

                    boolean converged = false;

                    while ((nIter == 0) ||
                        (!converged && (nIter < nMaxIter) && (transXDelta > 1))) {

                        if (nIter > 0) {
                            transXDelta--;
                            transYDelta--;
                        }

                        TransformationPointFit fit =
                            pointMatcher.calculateTranslationFromGridThenDownhillSimplex(
                                scaledRotatedLeft,  unmatchedLeftXY, unmatchedRightXY,
                                imageWidth, imageHeight, imageWidth, imageHeight,
                                params.getRotationInRadians(), params.getScale(),
                                transXStartStop, transXDelta,
                                transYStartStop, transYDelta,
                                tolTransX, tolTransY, setsAreMatched,
                                setsFractionOfImage);

                        assert(fit != null);
                        TransformationParameters fitParams = fit.getParameters();
                        int diffN = Math.abs(nPoints - fit.getNumberOfMatchedPoints());
                        float diffRotDeg = getAngleDifference(
                            fitParams.getRotationInDegrees(), rotInDegrees);
                        float diffScale = Math.abs(fitParams.getScale() - scale);
                        float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                        float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                        log.info("FINAL FIT=" + fit.toString());
                        log.info("diff result and expected =" +
                            String.format(
                            " dNPoints=%d, dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f",
                            diffN, diffRotDeg, diffScale, diffTransX, diffTransY)
                        );

                        double epsTrans = 10;
                        double epsRot = 5;

                        if ((diffN < 0.5*nPoints) && (diffRotDeg <= epsRot) &&
                            (diffScale < 0.2) && (diffTransX <= epsTrans) &&
                            (diffTransY <= epsTrans)) {
                            converged = true;
                        }

                        nIter++;
                    }
                    assertTrue(converged);
                    if (nIter > 1) {
                        log.info("needed change for translation delta for dens=" + dens
                            + " nPoints=" + nPoints + " w=" + imageWidth + " h=" + imageHeight
                            + " transXDelta=" + transXDelta
                            + " transYDelta=" + transYDelta
                        );
                        converged1.add(
                            new DensityTranslationResults(dens, nPoints,
                                imageWidth, imageHeight, transXDelta));
                    } else {
                        log.info("density=" + dens + " nPoints=" + nPoints
                            + " w=" + imageWidth + " h=" + imageHeight
                            + " transXDelta=" + transXDelta
                            + " transYDelta=" + transYDelta
                        );
                        converged0.add(
                            new DensityTranslationResults(dens, nPoints,
                                imageWidth, imageHeight, transXDelta));
                    }
                    /* For having rotation and scale correct already,
                    a transXDelta and transYDelta of values 1 are needed
                    when the point density is >= 0.04 nPoints/dimension in 
                    pixels in order to assure convergence (diff from expected <= 3).

                    That's 41 points or more in a 1024 x 1024 image that would lead
                    to requirement of translation delta = 1, which means trying every
                    pixel combination of translation in the starter points, and then
                    the downhill simplex portion isn't needed, so that
                    would be equiv to brute force imageDimension^2 for just translation,
                    (and then factors for the rotation and scale iterations).
                    
Would like to see if I can change any parameters in the downhill simplex
such that starter points created w/ delta translation = 2 or 3 at smallest,
will still produce a converging result in the downhill simplex.
                    
                    */
                }
            }
        }

        //need to plot density vs transDelta for converged and those
        //that took more than one iteration to converge.
        float[] x = new float[converged0.size()];
        float[] y = new float[x.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = (float)converged0.get(i).density;
            y[i] = (float)converged0.get(i).translationDelta;
        }
        ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
        plotter.plot(x, y, "converged", "density", "trDelta");
        plotter.writeFile(MiscDebug.getCurrentTimeFormatted());

        x = new float[converged1.size()];
        y = new float[x.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = (float)converged1.get(i).density;
            y[i] = (float)converged1.get(i).translationDelta;
        }
        plotter = new ScatterPointPlotterPNG();
        plotter.plot(x, y, "converged after decr trDelta", "density", "trDelta");
        plotter.writeFile(MiscDebug.getCurrentTimeFormatted());

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

    private void examineIterativeCorners() throws Exception {


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
        PairIntArray points1 = null;
        PairIntArray points2 = null;

        int nPreferredCorners = 100;
        int nCrit = 500;

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
        //detector.useOutdoorMode();
        //detector.useOutdoorModeAndExtractSkyline();
        detector.calculateSkylineCornersOnly();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        //edges1 = detector.getEdgesInOriginalReferenceFrame();
        points1 = detector.getSkylineCornersInOriginalReferenceFrame();
        edgesSkyline1 = detector.getSkylineEdgesInOriginalReferenceFrame();

        detector = new CurvatureScaleSpaceCornerDetector(img2);
        //detector.useOutdoorMode();
        //detector.useOutdoorModeAndExtractSkyline();
        detector.calculateSkylineCornersOnly();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        //edges2 = detector.getEdgesInOriginalReferenceFrame();
        points2 = detector.getSkylineCornersInOriginalReferenceFrame();
        edgesSkyline2 = detector.getSkylineEdgesInOriginalReferenceFrame();

        Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        String dirPath = ResourceFinder.findDirectory("bin");
        Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        // ===== use partitioned matching =====
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        PointMatcher pointMatcher = new PointMatcher();
        //pointMatcher.setCostToNumMatchedAndDiffFromModel();
        //pointMatcher.setCostToDiffFromModel();

        boolean useLargestToleranceForOutput = true;

        TransformationPointFit trFit =
            pointMatcher.performPartitionedMatching(points1, points2,
                image1Width, image1Height, image2Width, image2Height,
                outputMatchedScene, outputMatchedModel,
                useLargestToleranceForOutput);

        log.info("rough euclidean from skyline alone=" + trFit.toString());

/*
the skyline is hopefully useful in getting the transformation solution
into the local search region

now should try rough euclidean match using the transformation from
skyline as a start of solution with whole image corners and prefer the
later if better.

then make a larger output matched scene and model set of matched points to
determine the epipolar projection

(for brown & lee, the skyline-only region of overlap is small so
difficult to make a precise epipolar projection solution.  it's
better to have points of correspondence spread over as much of the
images as possible)
*/

        /*
        pointMatcher.performMatching(points1, points2,
        image1Width >> 1, image1Height >> 1,
            image2Width >> 1, image2Height >> 1,
            outputMatchedScene, outputMatchedModel, 1.0f);
        */

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

    private void adjustPointsOfInterest() throws Exception {

        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        // revisit infl points.  is there a threshold removing points?
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();

        List<PairIntArray> edges1 = null;
        List<PairIntArray> edges2 = null;
        PairIntArray points1 = null;
        PairIntArray points2 = null;

        boolean makeInflectionPoints = false;

        if (makeInflectionPoints) {

            CurvatureScaleSpaceInflectionMapperForOpenCurves inflMapper = new
                CurvatureScaleSpaceInflectionMapperForOpenCurves(img1, img2);

            PairIntArray[] xyPeaks = inflMapper.createUnmatchedXYFromContourPeaks();
            points1 = xyPeaks[0];
            points2 = xyPeaks[1];

            // there may be a couple redundant points per set that should be removed
            /*
            DataForTests.writePointsToTestResources(points1,
                "brown_lowe_2003_image1_infl_pts.tsv");
            DataForTests.writePointsToTestResources(points2,
                "brown_lowe_2003_image2_infl_pts.tsv");
            */


            edges1 = inflMapper.getEdges1InOriginalReferenceFrame();

            edges2 = inflMapper.getEdges2InOriginalReferenceFrame();

        } else {

            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);

            detector.useOutdoorMode();

            detector.findCorners();

            edges1 = detector.getEdgesInOriginalReferenceFrame();

            points1 = detector.getCornersInOriginalReferenceFrame();

            detector = new
                CurvatureScaleSpaceCornerDetector(img2);

            detector.useOutdoorMode();

            detector.findCorners();

            edges2 = detector.getEdgesInOriginalReferenceFrame();

            points2 = detector.getCornersInOriginalReferenceFrame();
        }

        int nPointsEdges1 = 0;
        for (PairIntArray edge1 : edges1) {
            nPointsEdges1 += edge1.getN();
        }
        int nPointsEdges2 = 0;
        for (PairIntArray edge2 : edges2) {
            nPointsEdges2 += edge2.getN();
        }
        log.info("total number of points in edges=" + nPointsEdges1 + " , "
            + nPointsEdges2);

        Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);

        for (PairIntArray edge : edges1) {
            ImageIOHelper.addCurveToImage(edge, image1, 2,
                Color.YELLOW.getRed(), Color.YELLOW.getGreen(),
                Color.YELLOW.getBlue());
        }

        ImageIOHelper.addCurveToImage(points1, image1, 1, 255, 0, 0);

        String dirPath = ResourceFinder.findDirectory("bin");
        String outFilePath = dirPath + "/tmp1_edges_infl.png";

        ImageIOHelper.writeOutputImage(outFilePath, image1);

        Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        for (PairIntArray edge : edges2) {
            ImageIOHelper.addCurveToImage(edge, image2, 2,
                Color.YELLOW.getRed(), Color.YELLOW.getGreen(),
                Color.YELLOW.getBlue());
        }

        ImageIOHelper.addCurveToImage(points2, image2, 1, 255, 0, 0);

        outFilePath = dirPath + "/tmp2_edges_infl.png";

        ImageIOHelper.writeOutputImage(outFilePath, image2);

        log.info("POINTS1: ");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < points1.getN(); i++) {
            String str = String.format("%d %d\n", points1.getX(i), points1.getY(i));
            sb.append(str);
        }
        log.info(sb.toString());
        log.info("POINTS2: ");
        sb = new StringBuilder();
        for (int i = 0; i < points2.getN(); i++) {
            String str = String.format("%d %d\n", points2.getX(i), points2.getY(i));
            sb.append(str);
        }
        log.info(sb.toString());

    }

    //@Test
    public void est155() throws Exception {

        PairIntArray scene = new PairIntArray();
        PairIntArray model = new PairIntArray();
        DataForTests.readBrownAndLoweMatches(scene, model);

        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int xSceneCentroid = img1.getWidth() >> 1;
        int ySceneCentroid = img1.getHeight() >> 1;
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int xModelCentroid = img2.getWidth() >> 1;
        int yModelCentroid = img2.getHeight() >> 1;

        StereoProjectionTransformer st = new StereoProjectionTransformer();

        SimpleMatrix left =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(scene);
        SimpleMatrix right =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(model);

        SimpleMatrix fm =
            st.calculateEpipolarProjectionForPerfectlyMatched(left, right);

        overplotEpipolarLines(fm, scene.toPairFloatArray(),
            model.toPairFloatArray(), img1, img2, img1.getWidth(),
            img1.getHeight(), img2.getWidth(), img2.getHeight());

    }

    //@Test
    public void est156() throws Exception {

        PairIntArray scene = new PairIntArray();
        PairIntArray model = new PairIntArray();
        DataForTests.readBrownAndLoweCorners(scene, model);

        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int xSceneCentroid = img1.getWidth() >> 1;
        int ySceneCentroid = img1.getHeight() >> 1;
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int xModelCentroid = img2.getWidth() >> 1;
        int yModelCentroid = img2.getHeight() >> 1;

        StereoProjectionTransformer st = new StereoProjectionTransformer();

        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        StereoProjectionTransformerFit fit =
            st.calculateEpipolarProjectionForUnmatched(scene, model,
            xSceneCentroid, ySceneCentroid,
            xModelCentroid, yModelCentroid,
            outputMatchedScene, outputMatchedModel);

        overplotEpipolarLines(fit.getFundamentalMatrix(),
            scene.toPairFloatArray(), model.toPairFloatArray(),
            img1, img2, img1.getWidth(),
            img1.getHeight(), img2.getWidth(), img2.getHeight());
    }

    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height) throws IOException {

        String flNum = "";

        overplotEpipolarLines(fm, set1, set2, img1, img2, image1Width,
            image1Height, image2Width, image2Height, flNum);
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

            //test.testPerformVerticalPartitionedMatching();
            test.testCalculateTranslationFromGridThenDownhillSimplex();

            /*
            tests for :
            -- for same set w/ projection
            -- for same set w/ projection and noise
            tests for scales which are close to 1 and less than 2
            */

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

    private TransformationPointFit[] simulateNearestGridSearchResults(
        int nFits, int imageWidth2, int imageHeight2,
        TransformationParameters params,
        final PairFloatArray unmatched1ScaledRotated,
        final PairIntArray set2, float tolTransX, float tolTransY) {

        /*
        simulating the results of a grid search that divides
        the space of
            -imageWidth2 to +imageWidth2 by intervals of 10
        and
            -imageHeight2 to +imageHeight2 by intervals of 10
        and finds that the top 10 of those are the 10 closest to the
        real transX, transY given in params.
        */

        int transXStart = -1*imageWidth2 + 1;
        int transXStop = imageWidth2 - 1;
        int transYStart = -1*imageHeight2 + 1;
        int transYStop = imageHeight2 - 1;

        //2047*2/10=409
        /*
        -2047  -1638  -1229  -820  -411  -2  407  816  1225  1634  2043
                                       -49
        */
        int dx = (transXStop - transXStart)/10;
        int dy = (transYStop - transYStart)/10;

        // these are usually a small value such as 8, but need the evaluation
        // to still consider points and calculate diffs for points
        // within the size of the bin

        int transXBinNumber = (int)((params.getTranslationX() - transXStart)/dx);
        int transYBinNumber = (int)((params.getTranslationY() - transYStart)/dy);

        PointMatcher pointMatcher = new PointMatcher();

        TransformationPointFit[] fits = new TransformationPointFit[nFits];

        int count = 0;
        for (int i = (transXBinNumber - 1); i <= (transXBinNumber + 1); ++i) {
            if (i < 0 || i > 10) { continue;}
            for (int j = (transYBinNumber - 1); j <= (transYBinNumber + 1); ++j) {
                if (j < 0 || j > 10) { continue;}

                int tx = transXStart + (i * dx) + (dx/2);
                int ty = transYStart + (j * dy) + (dy/2);

                fits[count] =
                    pointMatcher.evaluateFitForUnMatchedGreedy(
                        unmatched1ScaledRotated, tx, ty, tolTransX, tolTransY,
                        set2, params.getScale(), params.getRotationInRadians());

                count++;
            }
        }

        int[] xs = new int[]{transXBinNumber - 2, transXBinNumber + 2};
        int[] ys = new int[]{transYBinNumber - 2, transYBinNumber + 2};
        for (int i : xs) {
            if (i < 0 || i > 9) { continue;}
            if (count > 9) {break;}

            for (int j : ys) {
                if (j < 0 || j > 9) { continue;}
                if (count > 9) {break;}

                int tx = transXStart + (i * dx) + (dx/2);
                int ty = transYStart + (j * dy) + (dy/2);

                TransformationParameters params1 = new TransformationParameters();
                params1.setRotationInRadians(params.getRotationInRadians());
                params1.setScale(params.getScale());
                params1.setTranslationX(tx);
                params1.setTranslationY(ty);

                fits[count] =
                    pointMatcher.evaluateFitForUnMatchedTransformedGreedy(
                        params1, unmatched1ScaledRotated, set2,
                        tolTransX, tolTransY);

                count++;
            }
        }

        pointMatcher.sortByDescendingMatches(fits, 0, (fits.length - 1));

        return fits;
    }

    private float getAngleDifference(float rotDegrees0, float rotDegrees1) {
         /*
         I  |  0
        ---------
         II | III
        */
        int q0 = 0;
        if (rotDegrees0 >= 270) {
            q0 = 3;
        } else if (rotDegrees0 >= 180) {
            q0 = 2;
        } else if (rotDegrees0 >= 90) {
            q0 = 1;
        }
        int q1 = 0;
        if (rotDegrees1 >= 270) {
            q0 = 3;
        } else if (rotDegrees1 >= 180) {
            q0 = 2;
        } else if (rotDegrees1 >= 90) {
            q0 = 1;
        }

        /*
         I  |  0
        ---------
         II | III
        */
        float angleDiff = -1;
        if (q0 == 0){
            if (q1 == 0) {
                angleDiff = Math.abs(rotDegrees1 - rotDegrees0);
            } else if (q1 == 1) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 2) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else {
                angleDiff = Math.abs(360 - rotDegrees1 + rotDegrees0);
            }
        } else if (q0 == 1) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 1) {
                angleDiff = Math.abs(rotDegrees0 - rotDegrees1);
            } else if (q1 == 2) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            }
        } else if (q0 == 2) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 1) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else if (q1 == 2) {
                angleDiff = Math.abs(rotDegrees0 - rotDegrees1);
            } else {
                angleDiff = (rotDegrees1 - rotDegrees0);
            }
        } else if (q0 == 3) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (360 - rotDegrees0 + rotDegrees1);
            } else if (q1 == 1) {
                float diff = (rotDegrees0 - rotDegrees1);
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 2) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else {
                angleDiff = Math.abs(rotDegrees0 - rotDegrees1);
            }
        }

        if (angleDiff >= 360) {
            angleDiff = 360 - angleDiff;
        }

        return angleDiff;
    }

}
