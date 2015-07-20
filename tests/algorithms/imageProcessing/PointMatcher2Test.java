package algorithms.imageProcessing;

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
   
    public void testPerformMatchingForMostlyVerticalTranslation() throws Exception {

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

        int nPreferredCorners = 200;
        int nCrit = 100;

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

        // ===== use partitioned matching =====
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        PointMatcher pointMatcher = new PointMatcher();

        log.info("nSkylineCorners1=" + skylineCorners1.getN() + " nSkylineCorners2=" 
            + skylineCorners2.getN());

        float scale = 1;
        float rotationLowLimitInDegrees = 335;
        float rotationHighLimitInDegrees = 25;
        
        TransformationPointFit skyLineFit =
            pointMatcher.calculateEuclideanTransformation(
                skylineCorners1, skylineCorners2,
                scale, rotationLowLimitInDegrees,
                rotationHighLimitInDegrees); 

        log.info("euclidean tr from skyline alone=" + skyLineFit.toString());
        
        writeTransformed(skyLineFit.getParameters(), image2.copyImage(),
            skylineCorners1, skylineCorners2, "transformedSkyline.png");
        
        writeImage(image1.copyImage(), corners1, "corners1.png");
        writeImage(image2.copyImage(), corners2, "corners2.png");
        
        log.info("nCorners1=" + corners1.getN() + " nCorners2=" + corners2.getN());

        TransformationPointFit allFit =
            pointMatcher.calculateEuclideanTransformation(corners1, corners2,
                scale, rotationLowLimitInDegrees, rotationHighLimitInDegrees); 
        
        if (allFit != null) {
            log.info("all fit from all corners=" + allFit.toString());
        }
        
        float rotHalfRange = 20;
        float rotDelta = 2f;
        float transHalfRange = 200;
        float transDelta = 4;
        
        boolean useGreedyMatching = true;
                
        TransformationPointFit trFit2 = pointMatcher.refineTheTransformation(
            skyLineFit.getParameters(), corners1, corners2,
            rotHalfRange, rotDelta,
            transHalfRange, transDelta, transHalfRange, transDelta,
            useGreedyMatching);

        if (trFit2 != null) {
            log.info("fit from skyline and all corners=" + trFit2.toString());
        }
        
        writeTransformed(trFit2.getParameters(), image2.copyImage(),
            corners1, corners2, "transformedCorners.png");
        
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

            test.testPerformMatchingForMostlyVerticalTranslation();

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
