package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.IOException;
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
public class PointMatcher2Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public PointMatcher2Test() {
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

    public void estPerformMatchingForMostlyVerticalTranslation() throws Exception {
        
        String fileName1, fileName2;
        
        for (int nMethod = 0; nMethod < 2; ++nMethod) {
            for (int nData = 0; nData < 2; ++nData) {

                if (nData == 0) {
                    // transX ~ -280 and transY ~ -20
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                } else {
                    // transX ~ -34 and transY ~ 0
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                }

                String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
                ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
                int image1Width = img1.getWidth();
                int image1Height = img1.getHeight();

                String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
                ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
                int image2Width = img2.getWidth();
                int image2Height = img2.getHeight();

                int nPreferredCorners = 200;
                int nCrit = 100;

                /*
                // quick look at simplified image w/ color segmentation to explore contour matching visually
                ImageProcessor imageProcessor = new ImageProcessor();
                ImageExt clrImg1 = ImageIOHelper.readImageExt(filePath1);
                GreyscaleImage csImg1 = imageProcessor.createGreyscaleFromColorSegmentation(clrImg1, 3);
                ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                    + "/color_seg1.png", csImg1);
                ImageExt clrImg2 = ImageIOHelper.readImageExt(filePath2);
                GreyscaleImage csImg2 = imageProcessor.createGreyscaleFromColorSegmentation(clrImg2, 3);
                ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                    + "/color_seg2.png", csImg2);
                */

                CurvatureScaleSpaceCornerDetector detector = new
                    CurvatureScaleSpaceCornerDetector(img1);
                detector.useOutdoorModeAndExtractSkyline();
                detector.findCornersIteratively(nPreferredCorners, nCrit);
                //detector.findCorners();
                PairIntArray corners1 = detector.getCornersInOriginalReferenceFrame();
                PairIntArray skylineCorners1 = detector.getSkylineCornersInOriginalReferenceFrame();
                List<PairIntArray> skylineEdges1 = detector.getSkylineEdgesInOriginalReferenceFrame();

                CurvatureScaleSpaceCornerDetector detector2 = new
                    CurvatureScaleSpaceCornerDetector(img2);
                detector2.useOutdoorModeAndExtractSkyline();
                //detector2.findCorners();
                detector2.findCornersIteratively(nPreferredCorners, nCrit);
                PairIntArray corners2 = detector2.getCornersInOriginalReferenceFrame();
                PairIntArray skylineCorners2 = detector2.getSkylineCornersInOriginalReferenceFrame();
                List<PairIntArray> skylineEdges2 = detector2.getSkylineEdgesInOriginalReferenceFrame();

                Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);
                String dirPath = ResourceFinder.findDirectory("bin");
                Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

                PairIntArray outputMatchedScene = new PairIntArray();
                PairIntArray outputMatchedModel = new PairIntArray();

                PointMatcher pointMatcher = new PointMatcher();

                log.info("nSkylineCorners1=" + skylineCorners1.getN() + " nSkylineCorners2="
                    + skylineCorners2.getN());

                float scale = 1;

                PairIntArray comb1 = skylineCorners1.copy();
                comb1.addAll(corners1);
                PairIntArray comb2 = skylineCorners2.copy();
                comb2.addAll(corners2);

                boolean useLargestToleranceForOutput= true;
                boolean useGreedyMatching = true;

                TransformationPointFit fit = null;
                
                long t0 = System.currentTimeMillis();
                
                if (nMethod == 0) {
                    fit = pointMatcher.performMatchingForMostlyVerticalTranslation(
                        skylineCorners1, skylineCorners2,
                        corners1, corners2,
                        outputMatchedScene, outputMatchedModel,
                        useLargestToleranceForOutput, useGreedyMatching);
                } else {
                    fit = pointMatcher.performMatching(comb1, comb2, 
                        outputMatchedScene, outputMatchedModel, 
                        useLargestToleranceForOutput);
                }

                long t1 = System.currentTimeMillis();
                double t0Sec = (t1 - t0) * 1e-3;

                log.info("(" + t0Sec +  " seconds)"
                    + " euclidean tr=" + fit.toString());

                float rotInDegrees, transX, transY;
                if (nData == 0) {
                    rotInDegrees = 0;
                    transX = -280;
                    transY = -20;
                } else {
                    rotInDegrees = 0;
                    transX = -34;
                    transY = 0;
                }

                TransformationParameters fitParams = fit.getParameters();
                float diffRotDeg = AngleUtil.getAngleDifference(
                fitParams.getRotationInDegrees(), rotInDegrees);
                float diffScale = Math.abs(fitParams.getScale() - scale);
                float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                assertTrue(diffScale < 0.1);
                assertTrue(diffRotDeg <= 20);
                assertTrue(diffTransX <= 50);
                assertTrue(diffTransY <= 50);

                writeTransformed(fit.getParameters(), image2.copyImage(),
                    comb1, comb2, "transformedSkyline_" + nData + ".png");

                writeImage(image1.copyImage(), comb1, 
                    "corners1_" + nData + "_" + nMethod + ".png");
                writeImage(image2.copyImage(), comb2, 
                    "corners2_" + nData + "_" + nMethod + ".png");

                writeTransformed(fit.getParameters(), image2.copyImage(),
                    outputMatchedScene, outputMatchedModel, 
                    "transformedCorners_" + nData + "_" + nMethod + ".png");

                assertTrue(outputMatchedScene.getN() >= 7);

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
            }
        }

        System.out.println("test done");
    }

    public void testPerformMatchingForMostlyVerticalTranslation2() throws Exception {
        
        String fileName1, fileName2;
        
        for (int nMethod = 0; nMethod < 1; ++nMethod) {
            for (int nData = 0; nData < 1; ++nData) {

                //if (nData == 0) {
                    // transX ~ -280 and transY ~ -20
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                //}

                String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
                ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
                int image1Width = img1.getWidth();
                int image1Height = img1.getHeight();

                String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
                ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
                int image2Width = img2.getWidth();
                int image2Height = img2.getHeight();

                int nPreferredCorners = 200;
                int nCrit = 100;

                
                // quick look at simplified image w/ color segmentation to explore contour matching visually
                ImageProcessor imageProcessor = new ImageProcessor();
                ImageExt clrImg1 = ImageIOHelper.readImageExt(filePath1);
                GreyscaleImage csImg1 = imageProcessor.createGreyscaleFromColorSegmentation(clrImg1, 3);
                imageProcessor.applyImageSegmentation(csImg1, 2);
                imageProcessor.fillInGaps(csImg1, 50);
                ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                    + "/color_seg1.png", csImg1);
                ImageExt clrImg2 = ImageIOHelper.readImageExt(filePath2);
                GreyscaleImage csImg2 = imageProcessor.createGreyscaleFromColorSegmentation(clrImg2, 3);
                imageProcessor.applyImageSegmentation(csImg2, 2);
                imageProcessor.finInGaps(csImg2, 50);
                ImageIOHelper.writeOutputImage(ResourceFinder.findDirectory("bin")
                    + "/color_seg2.png", csImg2);
if (true){return;}
                CurvatureScaleSpaceCornerDetector detector = new
                    CurvatureScaleSpaceCornerDetector(csImg1);
                detector.doNotPerformHistogramEqualization();
                detector.findCornersIteratively(nPreferredCorners, nCrit);
                //detector.findCorners();
                PairIntArray corners1 = detector.getCornersInOriginalReferenceFrame();

                CurvatureScaleSpaceCornerDetector detector2 = new
                    CurvatureScaleSpaceCornerDetector(csImg2);
                detector2.doNotPerformHistogramEqualization();
                detector2.findCornersIteratively(nPreferredCorners, nCrit);
                PairIntArray corners2 = detector2.getCornersInOriginalReferenceFrame();

                Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);
                String dirPath = ResourceFinder.findDirectory("bin");
                Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

                PairIntArray outputMatchedScene = new PairIntArray();
                PairIntArray outputMatchedModel = new PairIntArray();

                PointMatcher pointMatcher = new PointMatcher();

                log.info("corners1=" + corners1.getN() + " corners2="
                    + corners2.getN());
                
                writeImage(image1.copyImage(), corners1, 
                    "corners1_" + nData + "_" + nMethod + ".png");
                writeImage(image2.copyImage(), corners2, 
                    "corners2_" + nData + "_" + nMethod + ".png");

                float scale = 1;

                boolean useLargestToleranceForOutput= true;
                boolean useGreedyMatching = true;

                TransformationPointFit fit = null;
                
                long t0 = System.currentTimeMillis();
                
                if (nMethod == 0) {
                    fit = pointMatcher.performMatchingForMostlyVerticalTranslation(
                        corners1, corners2,
                        outputMatchedScene, outputMatchedModel,
                        useLargestToleranceForOutput, useGreedyMatching);
                } else {
                    fit = pointMatcher.performMatching(corners1, corners2, 
                        outputMatchedScene, outputMatchedModel, 
                        useLargestToleranceForOutput);
                }

                long t1 = System.currentTimeMillis();
                double t0Sec = (t1 - t0) * 1e-3;

                log.info("(" + t0Sec +  " seconds)"
                    + " euclidean tr=" + fit.toString());

                //190,218  268,217
                float rotInDegrees, transX, transY;
                //if (nData == 0) {
                    rotInDegrees = 0;
                    transX = +78;
                    transY = 1;
                //}

                TransformationParameters fitParams = fit.getParameters();
                float diffRotDeg = AngleUtil.getAngleDifference(
                fitParams.getRotationInDegrees(), rotInDegrees);
                float diffScale = Math.abs(fitParams.getScale() - scale);
                float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                assertTrue(diffScale < 0.1);
                assertTrue(diffRotDeg <= 20);
                assertTrue(diffTransX <= 50);
                assertTrue(diffTransY <= 50);

                writeTransformed(fit.getParameters(), image2.copyImage(),
                    corners1, corners2, "transformedSkyline_" + nData + ".png");

                writeTransformed(fit.getParameters(), image2.copyImage(),
                    outputMatchedScene, outputMatchedModel, 
                    "transformedCorners_" + nData + "_" + nMethod + ".png");

                assertTrue(outputMatchedScene.getN() >= 7);

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
            }
        }

        System.out.println("test done");
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

    public static void main(String[] args) {

        try {
            PointMatcher2Test test = new PointMatcher2Test();

            //test.testPerformMatchingForMostlyVerticalTranslation();
            test.testPerformMatchingForMostlyVerticalTranslation2();

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

}
