package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.CannyEdgeFilter;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.HistogramEqualization;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.imageProcessing.SegmentedCellMerger;
import algorithms.imageProcessing.WaterShed;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MedianSmooth;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntPair;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.ejml.data.Complex64F;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class AndroidStatuesTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public AndroidStatuesTest() {
    }

    public void test0() throws Exception {

        String fileName1, fileName2;

        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
//trees and grass contributing too many 2nd deriv pts.  does assoc w/ blobs retain enough remaining pts?
        for (int i = 0; i < 2; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_02.jpg";
                    fileName2 = "android_statues_04.jpg";
                    //fileName1 = "android_statues_02_gingerbreadman.jpg";
                    //fileName2 = "android_statues_04_gingerbreadman.jpg";
                    settings.setUseNormalizedFeatures(true);
                    settings.setToUse2ndDerivCorners();
                    break;
                }
                default: {
                    fileName1 = "android_statues_01.jpg";
                    fileName2 = "android_statues_03.jpg";
                    settings.setUseNormalizedFeatures(true);
                    settings.setToUse2ndDerivCorners();
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, settings, false);
        }
    }

    public void estRot90() throws Exception {

        String fileName1, fileName2;

        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);

        for (int i = 0; i < 5; ++i) {
            fileName1 = null;
            fileName2 = null;
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_02.jpg";
                    fileName2 = "android_statues_04.jpg";
                    //fileName1 = "campus_010.jpg";
                    //fileName2 = "campus_011.jpg";
                    settings.setUseNormalizedFeatures(true);
                    settings.setToUse2ndDerivCorners();
                    break;
                }
                default: {
                    fileName1 = "android_statues_01.jpg";
                    fileName2 = "android_statues_03.jpg";
                    settings.setUseNormalizedFeatures(true);
                    settings.setToUse2ndDerivCorners();
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, settings, true);
        }
    }

    private void runCorrespondenceList(String fileName1, String fileName2,
        FeatureMatcherSettings settings, boolean rotateBy90) throws Exception {

        if (fileName1 == null) {
            return;
        }

        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);

        settings.setDebugTag(fileName1Root);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();
        
        int maxDimension = 350;
        int binFactor1 = (int) Math.ceil(Math.max((float)w1/maxDimension,
            (float)h1/ maxDimension));
        int binFactor2 = (int) Math.ceil(Math.max((float)w2/maxDimension,
            (float)h2/ maxDimension));

        ImageExt img1Binned = imageProcessor.binImage(img1, binFactor1);
        ImageExt img2Binned = imageProcessor.binImage(img2, binFactor2);
                
        Set<PairIntPair> similarClass = new HashSet<PairIntPair>();
        Set<PairIntPair> differentClass = new HashSet<PairIntPair>();
        
        GreyscaleImage img1Ws = imageSegmentation.createGreyscaleO1Watershed(img1Binned, 
            "_" + fileName1Root + "_1_binned", img1.getWidth(), img1.getHeight());
        SegmentedCellMerger scm = new SegmentedCellMerger(img1Binned, img1Ws,
            0, "_" + fileName1Root + "_1_binned");
        populateClasses(similarClass, differentClass, fileName1Root);
        scm.setClassPairs(similarClass, differentClass, fileName1Root);
        scm.merge();
        
        /*
        GreyscaleImage img2Ws = imageSegmentation.createGreyscaleO1Watershed(img2Binned, 
            "_" + fileName2Root + "_2_binned",
            img2.getWidth(), img2.getHeight());
        
        SegmentedCellMerger scm = new SegmentedCellMerger(img2Binned, img2Ws,
            0, "_" + fileName2Root + "_2_binned");
        populateClasses(similarClass, differentClass, fileName2Root);
        scm.setClassPairs(similarClass, differentClass, fileName2Root);
        scm.merge();
        */
        
if (true) {
    return;
}
        if (rotateBy90) {
            TransformationParameters params90 = new TransformationParameters();
            params90.setRotationInDegrees(90);
            params90.setOriginX(0);
            params90.setOriginY(0);
            params90.setTranslationX(0);
            params90.setTranslationY(img1.getWidth() - 1);
            Transformer transformer = new Transformer();
            img1 = (ImageExt) transformer.applyTransformation(img1,
                params90, img1.getHeight(), img1.getWidth());

            /*
            venturi:
                tx += -50
                ty += -1
                rot += 0
            books:
                tx += -74 ---> total is 620 when rot=270
                ty += -0
                rot += 0
            */

            /*
            MatchedPointsTransformationCalculator tc =
                new MatchedPointsTransformationCalculator();

            TransformationParameters revParams = tc.swapReferenceFrames(params90);
            transformer.transformToOrigin(0, 0, revParams);
            revParams.setTranslationX(revParams.getTranslationX() + -74);
            revParams.setTranslationY(revParams.getTranslationY() + -0);
            revParams.setRotationInDegrees(revParams.getRotationInDegrees() - 0);
            log.info("revParams: " + revParams.toString());

            ImageExt img1RevTr = img1.copyToImageExt();
            img1RevTr = (ImageExt) transformer.applyTransformation(img1RevTr,
                revParams, img1RevTr.getHeight(), img1RevTr.getWidth());
            MiscDebug.writeImage(img1RevTr, "rot90_rev_trans");
            */
        }

        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();

        log.info("fileName1Root=" + fileName1Root);

        EpipolarColorSegmentedSolver solver = new EpipolarColorSegmentedSolver(img1, img2, settings);

        boolean solved = solver.solve();

        assertTrue(solved);

        //MiscDebug.writeImagesInAlternatingColor(img1, img2, stats,
        //    fileName1Root + "_matched_non_euclid", 2);

    }

     public static void main(String[] args) {

        try {
            AndroidStatuesTest test = new AndroidStatuesTest();
            //test.test0();
            //test.testRot90();

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    public void estColor() throws Exception {

        //String fileName1 = "android_statues_02.jpg";
        //String fileName2 = "android_statues_04.jpg";
        String fileName1 = "android_statues_02_gingerbreadman.jpg";
        String fileName2 = "android_statues_04_gingerbreadman.jpg";
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        CIEChromaticity cieC = new CIEChromaticity();
        int binnedImageMaxDimension = 512;
        int binFactor1 = (int) Math.ceil(
            Math.max((float) img1.getWidth() / (float)binnedImageMaxDimension,
            (float) img1.getHeight() / (float)binnedImageMaxDimension));
        int binFactor2 = (int) Math.ceil(
            Math.max((float) img2.getWidth() / (float)binnedImageMaxDimension,
            (float) img2.getHeight() / (float)binnedImageMaxDimension));
        ImageProcessor imageProcessor = new ImageProcessor();
        ImageExt imgBinned1 = imageProcessor.binImage(img1, binFactor1);
        ImageExt imgBinned2 = imageProcessor.binImage(img2, binFactor2);

        /*
        HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(imgBinned1);
        hEq.applyFilter();
        hEq = new HistogramEqualizationForColor(imgBinned2);
        hEq.applyFilter();
        */
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        GreyscaleImage redBinnedImg1 = imgBinned1.copyRedToGreyscale();
        GreyscaleImage greenBinnedImg1 = imgBinned1.copyGreenToGreyscale();
        GreyscaleImage blueBinnedImg1 = imgBinned1.copyBlueToGreyscale();
        GreyscaleImage redBinnedImg2 = imgBinned2.copyRedToGreyscale();
        GreyscaleImage greenBinnedImg2 = imgBinned2.copyGreenToGreyscale();
        GreyscaleImage blueBinnedImg2 = imgBinned2.copyBlueToGreyscale();
        GreyscaleImage gsImg1 = imgBinned1.copyToGreyscale();
        GreyscaleImage gsImg2 = imgBinned2.copyToGreyscale();
        IntensityClrFeatures clrFeaturesBinned1 = new IntensityClrFeatures(gsImg1.copyImage(),
            5, rotatedOffsets);
        IntensityClrFeatures clrFeaturesBinned2 = new IntensityClrFeatures(gsImg2.copyImage(),
            5, rotatedOffsets);
        IntensityFeatures features1 = new IntensityFeatures(5, true, rotatedOffsets);
        IntensityFeatures features2 = new IntensityFeatures(5, true, rotatedOffsets);

        /*
        looking at trace/determinant of autocorrelation
        and eigenvalues of greyscale, autocorrelation, and lab colors for
        selected points in both images

        statues subsets:

        0 (64, 100) (96, 109)
        1 (67, 103) (103, 111)
        2 (68, 78)  (113, 86)
        3 (66, 49)  (106, 50)
        4 (92, 108) (157, 118)
        5 (92, 111) (160, 122)   delta e = 28.9
        6 (69, 129) (108, 142)   delta e = 26.4

        is the edelta for the gingerbread man's white stripes and shadows
        the same for shadow and higher illumination?
        */
        // single values of edelta
        List<PairInt> points1 = new ArrayList<PairInt>();
        List<PairInt> points2 = new ArrayList<PairInt>();
        points1.add(new PairInt(64, 100)); points2.add(new PairInt(96, 109));
        points1.add(new PairInt(67, 103)); points2.add(new PairInt(103, 111));
        points1.add(new PairInt(68, 78)); points2.add(new PairInt(113, 86));
        points1.add(new PairInt(66, 49)); points2.add(new PairInt(106, 50));
        points1.add(new PairInt(92, 108)); points2.add(new PairInt(157, 118));
        points1.add(new PairInt(92, 111)); points2.add(new PairInt(160, 122));
        points1.add(new PairInt(69, 129)); points2.add(new PairInt(108, 142));
        // this one is too look at localizability:
        points1.add(new PairInt(46, 65)); points2.add(new PairInt(6, 4));

        int n = points1.size();
        for (int i = 0; i < n; ++i) {

            StringBuilder sb = new StringBuilder();

            PairInt p1 = points1.get(i);
            PairInt p2 = points2.get(i);
            sb.append(p1.toString()).append(p2.toString());

            int r1 = redBinnedImg1.getValue(p1.getX(), p1.getY());
            int g1 = greenBinnedImg1.getValue(p1.getX(), p1.getY());
            int b1 = blueBinnedImg1.getValue(p1.getX(), p1.getY());
            int r2 = redBinnedImg2.getValue(p2.getX(), p2.getY());
            int g2 = greenBinnedImg2.getValue(p2.getX(), p2.getY());
            int b2 = blueBinnedImg2.getValue(p2.getX(), p2.getY());

            float[] lab1 = cieC.rgbToCIELAB(r1, g1, b1);
            float[] lab2 = cieC.rgbToCIELAB(r2, g2, b2);
            float[] cieXY1 = cieC.rgbToCIEXYZ(r1, g1, b1);
            float[] cieXY2 = cieC.rgbToCIEXYZ(r2, g2, b2);
            double deltaE = cieC.calcDeltaECIE94(lab1, lab2);

            sb.append(String.format("  dE=%.1f", (float)deltaE));

            int rot1 = clrFeaturesBinned1.calculateOrientation(p1.getX(), p1.getY());
            int rot2 = clrFeaturesBinned2.calculateOrientation(p2.getX(), p2.getY());

            IntensityDescriptor desc_l1 = clrFeaturesBinned1.extractIntensityLOfCIELAB(
                redBinnedImg1, greenBinnedImg1, blueBinnedImg1, p1.getX(), p1.getY(),
                rot1);

            IntensityDescriptor desc_a1 = clrFeaturesBinned1.extractIntensityAOfCIELAB(
                redBinnedImg1, greenBinnedImg1, blueBinnedImg1, p1.getX(), p1.getY(),
                rot1);

            IntensityDescriptor desc_b1 = clrFeaturesBinned1.extractIntensityBOfCIELAB(
                redBinnedImg1, greenBinnedImg1, blueBinnedImg1, p1.getX(), p1.getY(),
                rot1);

            IntensityDescriptor desc_l2 = clrFeaturesBinned2.extractIntensityLOfCIELAB(
                redBinnedImg2, greenBinnedImg2, blueBinnedImg2, p2.getX(), p2.getY(),
                rot2);

            IntensityDescriptor desc_a2 = clrFeaturesBinned2.extractIntensityAOfCIELAB(
                redBinnedImg2, greenBinnedImg2, blueBinnedImg2, p2.getX(), p2.getY(),
                rot2);

            IntensityDescriptor desc_b2 = clrFeaturesBinned2.extractIntensityBOfCIELAB(
                redBinnedImg2, greenBinnedImg2, blueBinnedImg2, p2.getX(), p2.getY(),
                rot2);

            IntensityDescriptor desc1 = features1.extractIntensity(gsImg1,
                p1.getX(), p1.getY(), rot1);

            IntensityDescriptor desc2 = features2.extractIntensity(gsImg2,
                p2.getX(), p2.getY(), rot2);

            double det, trace;
            SimpleMatrix a_l1 = clrFeaturesBinned1.createAutoCorrelationMatrix(desc_l1);
            det = a_l1.determinant();
            trace = a_l1.trace();
            sb.append(String.format("\n  L1_det(A)/trace=%.1f", (float)(det/trace)));
            SimpleMatrix a_l2 = clrFeaturesBinned2.createAutoCorrelationMatrix(desc_l2);
            det = a_l2.determinant();
            trace = a_l2.trace();
            sb.append(String.format("  L2_det(A)/trace=%.1f", (float)(det/trace)));

            try {
                sb.append("\n  eigen values:\n");
                SimpleEVD eigen1 = a_l1.eig();
                for (int j = 0; j < eigen1.getNumberOfEigenvalues(); ++j) {
                    Complex64F eigen = eigen1.getEigenvalue(j);
                    sb.append(String.format("    [1] %d %.1f %.1f\n", j,
                        (float)eigen.getReal(), (float)eigen.getMagnitude()));
                }
                sb.append("\n");
                SimpleEVD eigen2 = a_l2.eig();
                for (int j = 0; j < eigen2.getNumberOfEigenvalues(); ++j) {
                    Complex64F eigen = eigen2.getEigenvalue(j);
                    sb.append(String.format("    [1] %d %.1f %.1f\n", j,
                        (float)eigen.getReal(), (float)eigen.getMagnitude()));
                }
            } catch (Throwable t) {
            }

            if (desc1 != null && desc2 != null) {
                SimpleMatrix gs_l1 = features1.createAutoCorrelationMatrix(desc1);
                det = gs_l1.determinant();
                trace = gs_l1.trace();
                sb.append(String.format("\n  Grey_det(A)/trace=%.1f", (float)(det/trace)));
                SimpleMatrix gs_l2 = features2.createAutoCorrelationMatrix(desc2);
                det = gs_l2.determinant();
                trace = gs_l2.trace();
                sb.append(String.format("  Grey_det(A)/trace=%.1f", (float)(det/trace)));

                try {
                    sb.append("\n  eigen values:\n");
                    SimpleEVD eigen1 = gs_l1.eig();
                    for (int j = 0; j < eigen1.getNumberOfEigenvalues(); ++j) {
                        Complex64F eigen = eigen1.getEigenvalue(j);
                        sb.append(String.format("    [1] %d %.1f %.1f\n", j,
                            (float) eigen.getReal(), (float) eigen.getMagnitude()));
                    }
                    sb.append("\n");
                    SimpleEVD eigen2 = gs_l2.eig();
                    for (int j = 0; j < eigen2.getNumberOfEigenvalues(); ++j) {
                        Complex64F eigen = eigen2.getEigenvalue(j);
                        sb.append(String.format("    [2] %d %.1f %.1f\n", j,
                            (float) eigen.getReal(), (float) eigen.getMagnitude()));
                    }
                } catch (Throwable t) {
                }

            }

            log.info(sb.toString());
        }
    }

    private void populateClasses(Set<PairIntPair> similarClass, 
        Set<PairIntPair> differentClass, String fileNameRoot) {
        
        if (fileNameRoot.contains("android_statues_02")) {
            
            similarClass.add(new PairIntPair(278,133,285,135));
            similarClass.add(new PairIntPair(285,135,290,131));
            similarClass.add(new PairIntPair(290,131,289,125));
            similarClass.add(new PairIntPair(289,125,282,120));
            similarClass.add(new PairIntPair(282,120,276,123));
            similarClass.add(new PairIntPair(278,133,276,123));
            similarClass.add(new PairIntPair(278,112,282,120));
            similarClass.add(new PairIntPair(267,110,278,112));
            similarClass.add(new PairIntPair(267,110,269,119));
            similarClass.add(new PairIntPair(267,110,273,103));
            similarClass.add(new PairIntPair(267,110,257,112));
            similarClass.add(new PairIntPair(273,103,273,96));
            similarClass.add(new PairIntPair(267,110,262,120));
            similarClass.add(new PairIntPair(262,96,273,96));
            similarClass.add(new PairIntPair(269,89,273,96));
            similarClass.add(new PairIntPair(269,89,271,80));
            similarClass.add(new PairIntPair(269,89,275,86));
            similarClass.add(new PairIntPair(275,86,278,91));
            similarClass.add(new PairIntPair(278,91,282,87));
            similarClass.add(new PairIntPair(275,86,282,87));
            similarClass.add(new PairIntPair(269,89,271,80));
            similarClass.add(new PairIntPair(271,80,275,86));
            similarClass.add(new PairIntPair(262,79,271,80));
            similarClass.add(new PairIntPair(262,79,269,71));
            similarClass.add(new PairIntPair(282,87,290,75));
            similarClass.add(new PairIntPair(290,75,288,66));
            similarClass.add(new PairIntPair(288,66,280,67));
            similarClass.add(new PairIntPair(265,63,269,71));
            similarClass.add(new PairIntPair(265,63,260,65));
            similarClass.add(new PairIntPair(269,71,272,56));
            similarClass.add(new PairIntPair(260,54,272,56));
            similarClass.add(new PairIntPair(272,56,276,46));
            similarClass.add(new PairIntPair(276,46,276,40));
            similarClass.add(new PairIntPair(261,41,264,35));
            similarClass.add(new PairIntPair(70,82,84,93));
            similarClass.add(new PairIntPair(59,169,86,171));
            similarClass.add(new PairIntPair(84,140,72,137));
            similarClass.add(new PairIntPair(72,137,61,132));
            similarClass.add(new PairIntPair(44,123,61,132));
            similarClass.add(new PairIntPair(51,39,62,28));
            similarClass.add(new PairIntPair(31,71,51,39));
            
            differentClass.add(new PairIntPair(19,168,34,166));
            differentClass.add(new PairIntPair(19,168,39,173));
            differentClass.add(new PairIntPair(163,79,148,76));
            differentClass.add(new PairIntPair(163,79,176,67));
            differentClass.add(new PairIntPair(181,86,163,79));
            differentClass.add(new PairIntPair(233,8,0,0));
             
        } else if (fileNameRoot.contains("android_statues_04")) {
            
            similarClass.add(new PairIntPair(105,143,137,145));
            similarClass.add(new PairIntPair(125,130,105,143));
            similarClass.add(new PairIntPair(125,130,137,145));
            similarClass.add(new PairIntPair(117,125,105,143));
            similarClass.add(new PairIntPair(117,125,125,130));
            similarClass.add(new PairIntPair(117,125,114,115));
            similarClass.add(new PairIntPair(114,115,125,130));
            similarClass.add(new PairIntPair(128,114,131,119));
            similarClass.add(new PairIntPair(128,114,114,115));
            similarClass.add(new PairIntPair(131,119,125,130));
            similarClass.add(new PairIntPair(98,128,103,130));
            similarClass.add(new PairIntPair(98,128,93,127));
            similarClass.add(new PairIntPair(92,120,83,121));
            similarClass.add(new PairIntPair(92,120,93,127));
            similarClass.add(new PairIntPair(93,127,83,121));
            similarClass.add(new PairIntPair(144,116,144,128));
            similarClass.add(new PairIntPair(144,116,136,125));
            similarClass.add(new PairIntPair(144,128,137,145));
            similarClass.add(new PairIntPair(190,105,206,119));
            similarClass.add(new PairIntPair(240,107,223,103));
            similarClass.add(new PairIntPair(240,107,225,134));
            similarClass.add(new PairIntPair(240,107,259,105));
            similarClass.add(new PairIntPair(240,107,253,133));
            similarClass.add(new PairIntPair(11,52,15,60));
            
            differentClass.add(new PairIntPair(138,29,245,54));
            differentClass.add(new PairIntPair(240,107,245,54));
            differentClass.add(new PairIntPair(190,105,245,54));
            differentClass.add(new PairIntPair(54,73,46,75));
            
        } else if (fileNameRoot.contains("android_statues_01")) {
            
            similarClass.add(new PairIntPair(121,68,117,63));
            similarClass.add(new PairIntPair(121,68,113,67));
            similarClass.add(new PairIntPair(117,63,113,67));
            similarClass.add(new PairIntPair(114,82,120,82));
            similarClass.add(new PairIntPair(120,82,123,78));
            similarClass.add(new PairIntPair(121,68,128,76));
            similarClass.add(new PairIntPair(166,68,170,65));
            similarClass.add(new PairIntPair(166,68,160,81));
            similarClass.add(new PairIntPair(179,78,160,81));
            similarClass.add(new PairIntPair(178,87,160,81));
            similarClass.add(new PairIntPair(302,111,302,120));
            similarClass.add(new PairIntPair(302,111,310,120));
            similarClass.add(new PairIntPair(45,50,49,45));
            similarClass.add(new PairIntPair(93,27,101,26));
            similarClass.add(new PairIntPair(202,69,212,69));
                
            differentClass.add(new PairIntPair(121,68,133,68));
            differentClass.add(new PairIntPair(111,56,117,63));
            differentClass.add(new PairIntPair(265,130,281,122));
            differentClass.add(new PairIntPair(147,36,152,22));
            differentClass.add(new PairIntPair(147,36,153,47));
        
        } else if (fileNameRoot.contains("android_statues_03")) {   
            
            similarClass.add(new PairIntPair(48,89,21,89));
            similarClass.add(new PairIntPair(48,89,28,127));
            similarClass.add(new PairIntPair(48,89,68,125));
            similarClass.add(new PairIntPair(48,89,70,90));
            similarClass.add(new PairIntPair(145,106,149,141));
            similarClass.add(new PairIntPair(145,106,164,139));
            similarClass.add(new PairIntPair(145,106,172,105));
            similarClass.add(new PairIntPair(145,106,154,72));
            similarClass.add(new PairIntPair(86,58,120,58));
        
            differentClass.add(new PairIntPair(48,89,65,106));
            differentClass.add(new PairIntPair(48,89,94,84));
            differentClass.add(new PairIntPair(172,105,182,90));
        }
        
    }

}
