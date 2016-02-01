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
import algorithms.imageProcessing.WaterShed;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MedianSmooth;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
                    //fileName1 = "campus_010.jpg";
                    //fileName2 = "campus_011.jpg";
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

        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();
        
        //createTmpImages(img1, "_1");
        //createTmpImages(img2, "_2");
        
        int maxDimension = 350;
        int binFactor1 = (int) Math.ceil(Math.max((float)w1/maxDimension,
            (float)h1/ maxDimension));
        int binFactor2 = (int) Math.ceil(Math.max((float)w2/maxDimension,
            (float)h2/ maxDimension));

        ImageExt img1Binned = imageProcessor.binImage(img1, binFactor1);
        ImageExt img2Binned = imageProcessor.binImage(img2, binFactor2);
        //createTmpImages(img1Binned, "_1_binned");
        //createTmpImages(img2Binned, "_2_binned");
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        imageSegmentation.extractObjectEdges(img1Binned, "_" + fileName1Root + "_1_binned", 
            img1.getWidth(), img1.getHeight());
        imageSegmentation.extractObjectEdges(img2Binned, "_" + fileName2Root + "_2_binned",
            img2.getWidth(), img2.getHeight());
        
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

    private void createTmpImages(ImageExt img1, String lbl) {
        
        // O1 is (R-G)/sqrt(2)
        // O2 (R+G-2B)/sqrt(6)
        // O3 (R+G+B)/sqrt(2)
        
        CIEChromaticity cieC = new CIEChromaticity();
        GreyscaleImage rImg1 = img1.copyRedToGreyscale();
        GreyscaleImage gImg1 = img1.copyGreenToGreyscale();
        GreyscaleImage bImg1 = img1.copyBlueToGreyscale();
        GreyscaleImage o1 = rImg1.createFullRangeIntWithDimensions();
        GreyscaleImage o2 = rImg1.createFullRangeIntWithDimensions();
        GreyscaleImage o3 = rImg1.createFullRangeIntWithDimensions();
        GreyscaleImage lImg = rImg1.createFullRangeIntWithDimensions();
        GreyscaleImage aImg = rImg1.createFullRangeIntWithDimensions();
        GreyscaleImage bImg = rImg1.createFullRangeIntWithDimensions();
        // finding mode of l, a, and b in order to make a deltaE image
        float[] ls = new float[rImg1.getNPixels()];
        float[] as = new float[rImg1.getNPixels()];
        float[] bs = new float[rImg1.getNPixels()];
        float[] hueAngle = new float[rImg1.getNPixels()];
        float[] cieXYAngle = new float[rImg1.getNPixels()];
        for (int i = 0; i < rImg1.getNPixels(); ++i) {
            int r = rImg1.getValue(i);
            int g = gImg1.getValue(i);
            int b = bImg1.getValue(i);
            o1.setValue(i, (int)Math.round((double)(r - g)/Math.sqrt(2)));
            o2.setValue(i, (int)Math.round((double)(r + g - 2*b)/Math.sqrt(6)));
            o3.setValue(i, (int)Math.round((double)(r + g + b)/Math.sqrt(2)));
            float[] lab = cieC.rgbToCIELAB(r, g, b);
            lImg.setValue(i, (int)Math.round(lab[0]));
            aImg.setValue(i, (int)Math.round(lab[1]));
            bImg.setValue(i, (int)Math.round(lab[2]));
            ls[i] = lab[0];
            as[i] = lab[1];
            bs[i] = lab[2];
            cieXYAngle[i] = (float)(Math.atan(img1.getCIEY(i)/img1.getCIEX(i)) * 180. / Math.PI);
            if (cieXYAngle[i] < 0) {
                cieXYAngle[i] += 360.;
            }
            hueAngle[i] = (float)(Math.atan(bs[i]/as[i]) * 180. / Math.PI);
            if (hueAngle[i] < 0) {
                hueAngle[i] += 360.;
            }
        }
        HistogramEqualization hEq = new HistogramEqualization(o1);
        hEq.applyFilter();
        hEq = new HistogramEqualization(o2);
        hEq.applyFilter();
        hEq = new HistogramEqualization(o3);
        hEq.applyFilter();
        MiscDebug.writeImage(o1, "_o1_" + lbl);
        MiscDebug.writeImage(o2, "_o2_" + lbl);
        MiscDebug.writeImage(o3, "_o3_" + lbl);
        
        GreyscaleImage o3MinusO2 = rImg1.createFullRangeIntWithDimensions();
        for (int i = 0; i < rImg1.getNPixels(); ++i) {
            o3MinusO2.setValue(i, o3.getValue(i) - o2.getValue(i));
        }
        hEq = new HistogramEqualization(o3MinusO2);
        hEq.applyFilter();
        MiscDebug.writeImage(o3MinusO2, "_o3MinusO2" + lbl);
        
        Arrays.sort(ls);
        Arrays.sort(as);
        Arrays.sort(bs);
        float medianL = ls[ls.length/2];
        float medianA = as[ls.length/2];
        float medianB = bs[ls.length/2];
        GreyscaleImage deltaEImg = rImg1.createFullRangeIntWithDimensions();
        for (int i = 0; i < rImg1.getNPixels(); ++i) {
            float l1 = lImg.getValue(i);
            float a1 = aImg.getValue(i);
            float b1 = bImg.getValue(i);
            double deltaE = cieC.calcDeltaECIE94(l1, a1, b1, medianL, medianA, 
                medianB);
            deltaEImg.setValue(i, (int)Math.round(deltaE));
        }
        
        hEq = new HistogramEqualization(lImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(aImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(bImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(deltaEImg);
        hEq.applyFilter();
        
        MiscDebug.writeImage(lImg, "_l_" + lbl);
        MiscDebug.writeImage(aImg, "_a_" + lbl);
        MiscDebug.writeImage(bImg, "_b_" + lbl);
        MiscDebug.writeImage(bImg, "_deltaE_" + lbl);    
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        GreyscaleImage segImg = imageSegmentation.createGreyscale5(bImg, true);
        MiscDebug.writeImage(segImg, "_b_seg_" + lbl);
                
        WaterShed ws = new WaterShed();
        int[][] labelled1 = ws.createLabelledImage(aImg.copyImage());
        GreyscaleImage wsImg1 = new GreyscaleImage(aImg.getWidth(), aImg.getHeight(),
            GreyscaleImage.Type.Bits32FullRangeInt);
        for (int j = 0; j < aImg.getHeight(); ++j) {
            for (int i = 0; i < aImg.getWidth(); ++i) {
                int v = labelled1[i][j];
                wsImg1.setValue(i, j, v);
            }
        }
        MiscDebug.writeImage(wsImg1, "_a_watershed_" + lbl);

        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage aImgAM = aImg.copyImage();
        imageProcessor.applyAdaptiveMeanThresholding(aImgAM, 1);
        MiscDebug.writeImage(aImgAM, "_a_adaptive_median_" + lbl);
        
        GreyscaleImage bImgAM = bImg.copyImage();
        imageProcessor.applyAdaptiveMeanThresholding(bImgAM, 1);
        MiscDebug.writeImage(bImgAM, "_b_adaptive_median_" + lbl);
        

        ws = new WaterShed();
        labelled1 = ws.createLabelledImage(o1.copyImage());
        wsImg1 = new GreyscaleImage(o1.getWidth(), o1.getHeight(),
            GreyscaleImage.Type.Bits32FullRangeInt);
        for (int j = 0; j < o1.getHeight(); ++j) {
            for (int i = 0; i < o1.getWidth(); ++i) {
                int v = labelled1[i][j];
                wsImg1.setValue(i, j, v);
            }
        }
        MiscDebug.writeImage(wsImg1, "_o1_watershed_" + lbl);

        GreyscaleImage o1AM = o1.copyImage();
        imageProcessor.applyAdaptiveMeanThresholding(o1AM, 1);
        MiscDebug.writeImage(o1AM, "_o1_adaptive_median_" + lbl);
        
        
        ws = new WaterShed();
        labelled1 = ws.createLabelledImage(o2.copyImage());
        wsImg1 = new GreyscaleImage(o2.getWidth(), o2.getHeight(),
            GreyscaleImage.Type.Bits32FullRangeInt);
        for (int j = 0; j < o2.getHeight(); ++j) {
            for (int i = 0; i < o2.getWidth(); ++i) {
                int v = labelled1[i][j];
                wsImg1.setValue(i, j, v);
            }
        }
        MiscDebug.writeImage(wsImg1, "_o2_watershed_" + lbl);

        GreyscaleImage o2AM = o2.copyImage();
        imageProcessor.applyAdaptiveMeanThresholding(o2AM, 1);
        MiscDebug.writeImage(o2AM, "_o2_adaptive_median_" + lbl);
        
         
        ws = new WaterShed();
        labelled1 = ws.createLabelledImage(o3.copyImage());
        wsImg1 = new GreyscaleImage(o3.getWidth(), o3.getHeight(),
            GreyscaleImage.Type.Bits32FullRangeInt);
        for (int j = 0; j < o3.getHeight(); ++j) {
            for (int i = 0; i < o3.getWidth(); ++i) {
                int v = labelled1[i][j];
                wsImg1.setValue(i, j, v);
            }
        }
        MiscDebug.writeImage(wsImg1, "_o3_watershed_" + lbl);

        GreyscaleImage o3AM = o3.copyImage();
        imageProcessor.applyAdaptiveMeanThresholding(o3AM, 1);
        MiscDebug.writeImage(o3AM, "_o3_adaptive_median_" + lbl);
    }
}
