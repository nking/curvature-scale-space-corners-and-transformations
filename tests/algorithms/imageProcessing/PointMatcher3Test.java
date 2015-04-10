package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import static algorithms.imageProcessing.StereoProjectionTransformer.rewriteInto3ColumnMatrix;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.ejml.simple.*;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

/**
 *
 * @author nichole
 */
public class PointMatcher3Test {

    private void smallestSubsets() throws Exception {
        
        //String fileName1 = "brown_lowe_2003_image1.jpg";
        //String fileName2 = "brown_lowe_2003_image2.jpg";
       
        String fileName1 = "venturi_mountain_j6_0001.png";
        String fileName2 = "venturi_mountain_j6_0010.png";
        
        // revisit infl points.  is there a threshold removing points?
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleB(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
       
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleB(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        /*
        PairIntArray bl2003points1 = new PairIntArray();
        PairIntArray bl2003points2 = new PairIntArray();
        bl2003points1.add(471, 156); bl2003points2.add(194, 156);
        bl2003points1.add(411, 185); bl2003points2.add(138, 175);
        bl2003points1.add(331, 185); bl2003points2.add(54, 159);
        bl2003points1.add(353, 265); bl2003points2.add(63, 246);
        bl2003points1.add(352, 306); bl2003points2.add(52, 288);
        bl2003points1.add(502, 360); bl2003points2.add(188, 348);
        bl2003points1.add(384, 357); bl2003points2.add(77, 341);
        PairIntArray points1 = bl2003points1;
        PairIntArray points2 = bl2003points2;
        */
        
        PairIntArray venturipoints1 = new PairIntArray();
        PairIntArray venturipoints2 = new PairIntArray();
        venturipoints1.add(142, 240);  venturipoints2.add(106, 243);
        venturipoints1.add(285, 254);  venturipoints2.add(252, 257);
        venturipoints1.add(221, 350);  venturipoints2.add(184, 355);
        //venturipoints1.add(277, 356);  venturipoints2.add(244, 358);
        venturipoints1.add(589, 329);  venturipoints2.add(550, 326);
        //venturipoints1.add(589, 188);  venturipoints2.add(550, 188);
        //venturipoints1.add(588, 201);  venturipoints2.add(551, 201);
        venturipoints1.add(264, 201);  venturipoints2.add(230, 202);
        
        //venturipoints1.add(545, 325);  venturipoints2.add(509, 321);
        venturipoints1.add(83, 411); venturipoints2.add(42, 419);
        
        venturipoints1.add(456, 206); venturipoints2.add(422, 205);
        
        PairIntArray points1 = venturipoints1;
        PairIntArray points2 = venturipoints2;
        
        Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);

        ImageIOHelper.addCurveToImage(points1, image1, 1, 255, 0, 0);

        String dirPath = ResourceFinder.findDirectory("bin");
        String outFilePath = dirPath + "/tmp1_edges_infl.png";

        ImageIOHelper.writeOutputImage(outFilePath, image1);

        Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

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
       
        StereoProjectionTransformer st = new StereoProjectionTransformer();
        
        SimpleMatrix[] fms =
            st.calculateEpipolarProjectionFor7Points(
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(points1),
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(points2));
        
        for (SimpleMatrix fm : fms) {
                        
            overplotEpipolarLines(fm,
                points1.toPairFloatArray(), points2.toPairFloatArray(), 
                ImageIOHelper.readImage(filePath1),
                ImageIOHelper.readImage(filePath2), 
                image1Width, 
                img1.getHeight(), img2.getWidth(), img2.getHeight());
        }
        
        System.out.println("test done");
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
    
    private void examineInvPointLists() throws Exception {
        
        //TODO: implement the code for this, including inverting the image.
        
        // not cheking in the images for the temporary change
        /*
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName1Inv = "brown_lowe_2003_image1_inv.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String fileName2Inv = "brown_lowe_2003_image2_inv.jpg";
        */
        String fileName1 = "venturi_mountain_j6_0001.png";
        String fileName1Inv = "venturi_mountain_j6_0001_inv.png";
        String fileName2 = "venturi_mountain_j6_0010.png";
        String fileName2Inv = "venturi_mountain_j6_0010_inv.png";
        
        // revisit infl points.  is there a threshold removing points?
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath1Inv = ResourceFinder.findFileInTestResources(fileName1Inv);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        Image img1Inv = ImageIOHelper.readImage(filePath1Inv);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
       
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        String filePath2Inv = ResourceFinder.findFileInTestResources(fileName2Inv);
        Image img2Inv = ImageIOHelper.readImage(filePath2Inv);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        
        List<PairIntArray> edges1 = null;
        List<PairIntArray> edges2 = null;
        PairIntArray points1 = null;
        PairIntArray points2 = null;
        
        boolean makeInflectionPoints = false;
        
        int dist = 8;//10;
        
        PairIntArray tmp1 = null;
        PairIntArray tmp2 = null;
        
        if (makeInflectionPoints) {
            
            CurvatureScaleSpaceInflectionMapperForOpenCurves inflMapper = new
                CurvatureScaleSpaceInflectionMapperForOpenCurves(img1, img2);

            PairIntArray[] xyPeaks = inflMapper.createUnmatchedXYFromContourPeaks();        
            points1 = xyPeaks[0];
            points2 = xyPeaks[1];

            edges1 = inflMapper.getEdges1InOriginalReferenceFrame();
            
            edges2 = inflMapper.getEdges2InOriginalReferenceFrame();

            tmp1 = points1.copy();
            tmp2 = points2.copy();
            
            /*
            inflMapper = new
                CurvatureScaleSpaceInflectionMapperForOpenCurves(img1Inv, img2Inv);

            xyPeaks = inflMapper.createUnmatchedXYFromContourPeaks();        
            PairIntArray tmp1Inv = xyPeaks[0];
            PairIntArray tmp2Inv = xyPeaks[1];
            
            // points similar within 2 pixels in both sets:
            PairIntArray common1 = new PairIntArray();
            for (int i = 0; i < tmp1.getN(); i++) {
                int x = tmp1.getX(i);
                int y = tmp1.getY(i);
                for (int ii = 0; ii < tmp1Inv.getN(); ii++) {
                    int xInv = tmp1Inv.getX(ii);
                    int yInv = tmp1Inv.getY(ii);
                    float dx = Math.abs(xInv - x);
                    float dy = Math.abs(yInv - y);
                    if ((dx <= dist) && (dy <= dist)) {
                        common1.add(x, y);
                        break;
                    }
                }
            }
            points1 = common1;

            // points similar within 2 pixels in both sets:
            PairIntArray common2 = new PairIntArray();
            for (int i = 0; i < tmp2.getN(); i++) {
                int x = tmp2.getX(i);
                int y = tmp2.getY(i);
                for (int ii = 0; ii < tmp2Inv.getN(); ii++) {
                    int xInv = tmp2Inv.getX(ii);
                    int yInv = tmp2Inv.getY(ii);
                    float dx = Math.abs(xInv - x);
                    float dy = Math.abs(yInv - y);
                    if ((dx <= dist) && (dy <= dist)) {
                        common2.add(x, y);
                        break;
                    }
                }
            }
            points2 = common2;
            */
        } else {
            
            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);
            detector.useOutdoorMode();            
            detector.findCorners();
            edges1 = detector.getEdgesInOriginalReferenceFrame();
            points1 = detector.getCornersInOriginalReferenceFrame();
            
            tmp1 = points1.copy();
            
            /*
            detector = new CurvatureScaleSpaceCornerDetector(img1Inv);
            detector.useOutdoorMode();            
            detector.findCorners();
            PairIntArray tmp1Inv = detector.getCornersInOriginalReferenceFrame();
            
            // points similar within 2 pixels in both sets:
            PairIntArray common1 = new PairIntArray();
            for (int i = 0; i < tmp1.getN(); i++) {
                int x = tmp1.getX(i);
                int y = tmp1.getY(i);
                for (int ii = 0; ii < tmp1Inv.getN(); ii++) {
                    int xInv = tmp1Inv.getX(ii);
                    int yInv = tmp1Inv.getY(ii);
                    float dx = Math.abs(xInv - x);
                    float dy = Math.abs(yInv - y);
                    if ((dx <= dist) && (dy <= dist)) {
                        common1.add(x, y);
                        break;
                    }
                }
            }
            points1 = common1;
            */
            
            detector = new
                CurvatureScaleSpaceCornerDetector(img2);
            detector.useOutdoorMode();            
            detector.findCorners();
            edges2 = detector.getEdgesInOriginalReferenceFrame();
            points2 = detector.getCornersInOriginalReferenceFrame();
            
            tmp2 = points2.copy();
            
            /*
            detector = new CurvatureScaleSpaceCornerDetector(img2Inv);
            detector.useOutdoorMode();            
            detector.findCorners();
            PairIntArray tmp2Inv = detector.getCornersInOriginalReferenceFrame();
            
            // points similar within 2 pixels in both sets:
            PairIntArray common2 = new PairIntArray();
            for (int i = 0; i < tmp2.getN(); i++) {
                int x = tmp2.getX(i);
                int y = tmp2.getY(i);
                for (int ii = 0; ii < tmp2Inv.getN(); ii++) {
                    int xInv = tmp2Inv.getX(ii);
                    int yInv = tmp2Inv.getY(ii);
                    float dx = Math.abs(xInv - x);
                    float dy = Math.abs(yInv - y);
                    if ((dx <= dist) && (dy <= dist)) {
                        common2.add(x, y);
                        break;
                    }
                }
            }
            points2 = common2;
            */
        }
        
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
                
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();
            
        PointMatcher pointMatcher = new PointMatcher();
                
        pointMatcher.performPartitionedMatching(tmp1, tmp2,
        image1Width >> 1, image1Height >> 1,
            image2Width >> 1, image2Height >> 1,
            outputMatchedScene, outputMatchedModel);
        
        if (outputMatchedScene.getN() < 7) {
            // no solution
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
    
    @Test
    public void testSkyline() throws Exception {
        
        String[] fileNames = new String[] {
            "brown_lowe_2003_image1.jpg",
            //"brown_lowe_2003_image1_rot.jpg",
            //"brown_lowe_2003_image2.jpg",
            "venturi_mountain_j6_0001.png",
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
            "arizona-sunrise-1342919937GHz.jpg",
            "sky_with_rainbow.jpg",
            "sky_with_rainbow2.jpg",
            //"30.jpg",
            //"arches_sun_01.jpg",
            //"stlouis_arch.jpg", 
            //"contrail.jpg"
        };
        /*
         "brown_lowe_2003_image1.jpg",         fwhm hue=0.14, fwhm saturation=0.16       nsadded=
         "venturi_mountain_j6_0001.png",       fwhm hue=0.12, fwhm saturation=0.31       nsadded=528 to 118228 points <==== ? contig would revert
         "seattle.jpg",                        fwhm hue=0.11, fwhm saturation=0.12       nsadded=
         "arches.jpg",                         fwhm hue=0.11, fwhm saturation=0.29       nsadded=3739 to 217049 points;  contig passes, OK
         "stinson_beach.jpg",                  fwhm hue=0.11, fwhm saturation=0.37       nsadded=19112 to 63105 points;  contig fixes, OK
         "cloudy_san_jose.jpg",                fwhm hue=0.11, fwhm saturation=0.15       nsadded=31017 to 163605 points; contig fixes, OK
         "30.jpg",                             fwhm hue=0.  , fwhm saturation=0.         nsadded=
         "sky_with_rainbow.jpg",               fwhm hue=0.15, fwhm saturation=0.06,0.46  nsadded=42887 to 117339 points;  contig passes, OK;  other problems
         "sky_with_rainbow2.jpg",              fwhm hue=0.12, fwhm saturation=0.07       nsadded=547 to 63561 points; contig passes, OK;  other problems
         "stonehenge.jpg",                     fwhm hue=0.06, fwhm saturation=0.22,0.74  nsadded=178 to 162882 points; contig fixes, OK
         "norwegian_mtn_range.jpg",            fwhm hue=0.13, fwhm saturation=0.21       nsadded=902 to 124404 points; contig passes, OK;
         "halfdome.jpg",                       fwhm hue=0.12, fwhm saturation=0.15       nsadded=
         "costa_rica.jpg",                     fwhm hue=0.06, fwhm saturation=0.33       nsadded=
         "new-mexico-sunrise_w725_h490.jpg",   fwhm hue=0.16, fwhm saturation=0.15       nsadded=1588 to 295376 points;contig passes, OK;
         "arizona-sunrise-1342919937GHz.jpg"   fwhm hue=0.11, fwhm saturation=0.21       nsadded=49902 to 231634;contig passes, OK;
        */
        
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

            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);
            detector.useOutdoorModeAndExtractSkyline();
            //detector.findCornersIteratively(nPreferredCorners, nCrit);
            detector.findCorners();
            edges1 = detector.getEdgesInOriginalReferenceFrame();
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

            ImageIOHelper.writeOutputImage(outFilePath, image1);
        }
    }
    
    private void examineIterativeCorners() throws Exception {
        
        
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        
        /*
        String fileName1 = "venturi_mountain_j6_0001.png";
        String fileName2 = "venturi_mountain_j6_0010.png";
        */
        
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
        
        int nPreferredCorners = 100;
        int nCrit = 500;
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
        detector.useOutdoorMode();
        //detector.useOutdoorModeAndExtractSkyline();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        edges1 = detector.getEdgesInOriginalReferenceFrame();
        points1 = detector.getCornersInOriginalReferenceFrame();
        
        detector = new CurvatureScaleSpaceCornerDetector(img2);
        //detector.useOutdoorMode();
        detector.useOutdoorModeAndExtractSkyline();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        edges2 = detector.getEdgesInOriginalReferenceFrame();
        points2 = detector.getCornersInOriginalReferenceFrame();
        
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
            
        log.info("POINTS1: " + points1.getN() + " points");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < points1.getN(); i++) {
            String str = String.format("%d %d\n", points1.getX(i), points1.getY(i));
            sb.append(str);
        }
        log.info(sb.toString());
        log.info("POINTS2: " + points2.getN() + " points");
        sb = new StringBuilder();
        for (int i = 0; i < points2.getN(); i++) {
            String str = String.format("%d %d\n", points2.getX(i), points2.getY(i));
            sb.append(str);
        }
        log.info(sb.toString());
                
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();
            
        PointMatcher pointMatcher = new PointMatcher();
        //pointMatcher.setCostToNumMatchedAndDiffFromModel();
        //pointMatcher.setCostToDiffFromModel();
        
        pointMatcher.performPartitionedMatching(points1, points2,
        image1Width >> 1, image1Height >> 1,
            image2Width >> 1, image2Height >> 1,
            outputMatchedScene, outputMatchedModel);
        /*
        pointMatcher.performMatching(points1, points2,
        image1Width >> 1, image1Height >> 1,
            image2Width >> 1, image2Height >> 1,
            outputMatchedScene, outputMatchedModel, 1.0f);
        */
        if (outputMatchedScene.getN() < 7) {
            // no solution
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
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public PointMatcher3Test() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
   
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
  
    //@Test
    public void est1() throws Exception {

        // test for dataset which already matches exactly
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
        
        double scale = 1.;
        double rotation = 0.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 1);
    }
    
    //@Test
    public void est2() throws Exception {

        // test for exact match plus noise
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);

        double scale = 1.;
        double rotation = 0.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 2);
    }
    
    //@Test
    public void est3() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 3);
        
    }
    
    //@Test
    public void est4() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
       
        int nModelPoints = (int)(1.3f * nScenePoints);
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 0;
                
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 4);
    }
    
    //@Test
    public void est5() throws Exception {

        // test for exact match translated in X
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1419991550083L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 110;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 5);
        
    }
    
    //@Test
    public void est6() throws Exception {

        // test for exact match translated in X plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
                
        double scale = 1.;
        double rotation = 0.;
        int translateX = 200;
        int translateY = 110;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 6);
    }
    
    //@Test
    public void est7() throws Exception {

        // test for exact match rotated by 30 degrees
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1420149550374L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = nScenePoints;
                
        double rotation = 30. * Math.PI/180.;
        double scale = 1.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 7);
        
    }
  
    //@Test
    public void est8() throws Exception {

        // test for rotation plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
                
        double rotation = 30. * Math.PI/180.;
        double scale = 1;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 8);
    }
    
    //@Test
    public void est9() throws Exception {

        // test for scale plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 0.;
        double scale = 4.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 9);
    }
    
    //@Test
    public void est10() throws Exception {

        // test for scale smaller than 1, plus random points
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420179838620L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 0.;
        double scale = 1./4.;
        int translateX = 0;
        int translateY = 0;
        
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 10);
    }
    
    //@Test
    public void est11() throws Exception {

        // test for rotation and translation, plus random points
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420159107635L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 30.*Math.PI/180.;
        double scale = 1.;
        int translateX = 200;
        int translateY = -10;
     
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 11);
    }

    //@Test
    public void est12() throws Exception {

        // test for scale, rotation and translation, plus random points
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420159107635L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 30.*Math.PI/180.;
        double scale = 2.;
        int translateX = 200;
        int translateY = -10;
     
        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 12);
    }
    
    //@Test
    public void est13() throws Exception {

        // test for scale, rotation and translation, plus random points
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1420159107635L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 14.*Math.PI/180.;
        double scale = 1.;
        int translateX = 280;
        int translateY = -14;

        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 13);
    }
    
    //@Test
    public void est14() throws Exception {

        // test for scale close to 1
       
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1420187783647L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        int nScenePoints = 70;
        
        int xRange = 400;
        int yRange = 300;
        
        int nModelPoints = (int)(1.3f * nScenePoints);
        
        double rotation = 0.;
        double scale = 1.2;
        int translateX = 100;
        int translateY = -14;

        runTest(sr, nScenePoints, nModelPoints, xRange, yRange,
            scale, rotation, translateX, translateY, 13);
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
    
    private void runTest(SecureRandom sr, int nScenePoints, int nModelPoints, 
        int xRange, int yRange,
        double scale, double rotation, int translationX, int translationY,
        int testNumber) throws Exception {       
        
        int xModelMin = 0;
        int yModelMin = 0;
        
        PairIntArray model = createRandomPoints(sr, nScenePoints,
            xModelMin, yModelMin, xRange, yRange);
        
        PairIntArray scene = model.copy();
                        
        populateWithRandomPoints(sr, model, nScenePoints, nModelPoints, 
            xModelMin, yModelMin, xRange, yRange);
        
        int xModelCentroid = (xModelMin + xRange) >> 1;
        int yModelCentroid = (yModelMin + yRange) >> 1;
        
        scaleAndRotate(scene, 1./scale, -1*rotation, xModelCentroid, 
            yModelCentroid);
        int tx = (translationX == 0) ? 0 : (int)(-1*translationX/scale);
        int ty = (translationY == 0) ? 0 : (int)(-1*translationY/scale);
        translateX(scene, tx);
        translateY(scene, ty);
        
        // transform the centroid point from model to use for calculations
        TransformationParameters revParams = new TransformationParameters();
        revParams.setScale((float)(1./scale));
        revParams.setRotationInRadians(-1.f*(float)rotation);
        revParams.setTranslationX(tx);
        revParams.setTranslationY(ty);
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        double[] xySceneCen = tc.applyTransformation(revParams, 
            xModelCentroid, yModelCentroid, xModelCentroid, yModelCentroid);
           
        int xSceneCentroid = (int)xySceneCen[0];
        int ySceneCentroid = (int)xySceneCen[1];
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateProjectiveTransformationWrapper(
            scene, model, xSceneCentroid, ySceneCentroid, 
            xModelCentroid, yModelCentroid, 1.0f);

        System.out.println("=> " + fit.toString());
        
        TransformationParameters params = fit.getParameters();

        Transformer transformer = new Transformer();

        PairIntArray transformed = transformer.applyTransformation(
            params, scene, xSceneCentroid, ySceneCentroid);

        overplotTransformed(transformed, model, xRange, yRange, testNumber);

        int count = 0;
        for (int i = 0; i < transformed.getN(); i++) {
            int x = transformed.getX(i);
            // tolerance?
            if ((x < 0) || (x > 2*xModelCentroid)) {
                continue;
            }
            int y = transformed.getY(i);
            if ((y < 0) || (y > 2*yModelCentroid)) {
                continue;
            }
            count++;
        }
        
        System.out.println("=> Number of transformed scene points within bounds = " 
            + count);
        
        int nExpected = (nScenePoints > count) ? count : nScenePoints; 
        
        assertTrue(Math.abs(nExpected - fit.getNumberOfMatchedPoints()) 
            < 0.1*nScenePoints);
        
        assertTrue(Math.abs(params.getRotationInRadians() - rotation) <= 10.0);
        assertTrue(Math.abs(params.getScale() - scale) < 1.0);
        assertTrue(Math.abs(params.getTranslationX() - translationX) <= 1.0);
        assertTrue(Math.abs(params.getTranslationY() - translationY) <= 1.0);
        
        //solveForProjective(fit, scene, model, xSceneCentroid, ySceneCentroid);
    }
    
    private void solveForProjective(TransformationPointFit fit,
        PairIntArray scene, PairIntArray model, int xSceneCentroid, int 
            ySceneCentroid) {
        
        // or use params mean dist?
        double tolerance = 2;
        tolerance = fit.getMeanDistFromModel();
                
        Transformer transformer = new Transformer();
        
        PairFloatArray transformed = transformer.applyTransformation2(
            fit.getParameters(), scene, xSceneCentroid, ySceneCentroid);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        float[][] matched = pointMatcher.calculateMatchUsingOptimal(
            transformed, model, tolerance);
        
        PairIntArray matched1 = new PairIntArray();
        PairIntArray matched2 = new PairIntArray();
        
        pointMatcher.matchPoints(scene, model, tolerance, matched,
            matched1, matched2);
        
        /*
        from here, one could use RANSAC epipolar projection solver to
        calculate the epipolar projection from the true matched points,
        but the projection may have been large enough that there are few
        matched points still.
        
        So, wanting to look at transformed and model points overplotted in
        x vs y
        and the differences as x vs dx and y vs dy to see if a projective trend
        is visible.
        presumably, the x vs dx should be compared to the same for a larger
        tolerance to see the trend.
        
        the brown & lowe 2003 manually matched points should show this.
        */
        
    }
    
    private void overplotTransformed(double[][] transformed, double[][] model,
        int width, int height, int testNumber) throws IOException {
        
        Image image = new Image(width, height);
        
        for (int ii = 0; ii < transformed.length; ii++) {
            double x2 = transformed[ii][0];
            double y2 = transformed[ii][1];
            ImageIOHelper.addPointToImage((float) x2, (float) y2, image, 3, 
                0, 0, 255);
        }
        for (int ii = 0; ii < model.length; ii++) {
            double x = model[ii][0];
            double y = model[ii][1];
            ImageIOHelper.addPointToImage((float) x, (float) y, image, 2, 
                255, 0, 0);
        }
        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_t2_" + testNumber + ".png", image);
        
    }
    
    private void populateWithRandomPoints( SecureRandom sr, double[][] m, 
        int nPoints1, int nPoints2, int xMin, int yMin,
        double xRange, double yRange) {
        
        for (int i = nPoints1; i < nPoints2; i++) {
            m[i] = new double[2];
            double x = xMin + (sr.nextDouble() * xRange);
            double y = yMin + (sr.nextDouble() * yRange);
            m[i][0] = x;
            m[i][1] = y;
        }
    }
    
    private void populateWithRandomPoints(SecureRandom sr, PairIntArray m, 
        int nPoints1, int nPoints2, int xMin, int yMin,
        int xRange, int yRange) {
        
        for (int i = nPoints1; i < nPoints2; i++) {
            int x = xMin + sr.nextInt(xRange);
            int y = yMin + sr.nextInt(yRange);
            m.add(x, y);
        }
    }
    
    private void translateX(double[][] m, double translateX) {
        
        for (int i = 0; i < m.length; i++) {
            double x = m[i][0];
            m[i][0] = x + translateX;
        }
    }
    private void translateY(double[][] m, double translateY) {
        
        for (int i = 0; i < m.length; i++) {
            double y = m[i][1];
            m[i][1] = y + translateY;
        }
    }
    
    private void translateX(PairIntArray m, int translateX) {
        for (int i = 0; i < m.getN(); i++) {
            int x = m.getX(i);
            int y = m.getY(i);
            m.set(i, (x + translateX), y);
        }
    }
    
    private void translateY(PairIntArray m, int translateY) {
        for (int i = 0; i < m.getN(); i++) {
            int x = m.getX(i);
            int y = m.getY(i);
            m.set(i, x, (y + translateY));
        }
    }
    
    private void scaleAndRotate(double[][] m, 
        double scale, double rotationInRadians, 
        double centroidX, double centroidY) {
        /*
        xr[i] = centroidX1*s + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
        yr[i] = centroidY1*s + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
        */
        double scaleTimesCosine = scale * Math.cos(rotationInRadians);
        double scaleTimesSine = scale * Math.sin(rotationInRadians);
        
        for (int i = 0; i < m.length; i++) {
            double x = m[i][0];
            double y = m[i][1];
            double rx = centroidX*scale + ( 
                ((x - centroidX) * scaleTimesCosine) +
                ((y - centroidY) * scaleTimesSine));
            double ry = centroidY*scale + ( 
                (-(x - centroidX) * scaleTimesSine) +
                ((y - centroidY) * scaleTimesCosine));
            m[i][0] = rx;
            m[i][1] = ry;
        }
    }
    
    private void scaleAndRotate(PairIntArray m, 
        double scale, double rotationInRadians, 
        double centroidX, double centroidY) {
        
        //rotationInRadians *= -1;
        
        /*
        xr[i] = centroidX1*s + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
        yr[i] = centroidY1*s + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
        */
        double scaleTimesCosine = scale * Math.cos(rotationInRadians);
        double scaleTimesSine = scale * Math.sin(rotationInRadians);
        
        for (int i = 0; i < m.getN(); i++) {
            double x = m.getX(i);
            double y = m.getY(i);
            double rx = centroidX*scale + ( 
                ((x - centroidX) * scaleTimesCosine) +
                ((y - centroidY) * scaleTimesSine));
            double ry = centroidY*scale + ( 
                (-(x - centroidX) * scaleTimesSine) +
                ((y - centroidY) * scaleTimesCosine));
            
            m.set(i, (int)Math.round(rx), (int)Math.round(ry));
        }
    }
      
    private double[][] createRandomPoints(SecureRandom sr, int nPoints,
        int xMin, int yMin, double xRange, double yRange) {
     
        double[][] sMatrix = new double[nPoints][2];
        for (int i = 0; i < nPoints; i++) {
            sMatrix[i] = new double[2];
            double x = xMin + (sr.nextDouble()*xRange);
            double y = yMin + (sr.nextDouble()*yRange);
            sMatrix[i][0] = x;
            sMatrix[i][1] = y;
        }
        
        return sMatrix;
    }
    
    private PairIntArray createRandomPoints(SecureRandom sr, int nPoints,
        int xMin, int yMin, int xRange, int yRange) {
     
        PairIntArray output = new PairIntArray(nPoints);
        for (int i = 0; i < nPoints; i++) {
            int x = xMin + sr.nextInt(xRange);
            int y = yMin + sr.nextInt(yRange);
            output.add(x, y);
        }
        
        return output;
    }
    
    private double[][] createCopyOfSize(double[][] m, int sizeToCreate) {
         
        int end = m.length;
        if (sizeToCreate < end) {
            end = sizeToCreate;
        }
        
        double[][] copy = new double[sizeToCreate][2];
        for (int i = 0; i < end; i++) {
            copy[i] = new double[2];
            copy[i][0] = m[i][0];
            copy[i][1] = m[i][1];
        }
        
        return copy;
    }
    
    private double[] transform(SimpleMatrix fm, double x, double y) {
        
        double d = (fm.get(2, 0) * x) + (fm.get(2, 1) * y) 
            + fm.get(2, 2);
        
        double x2 = (fm.get(0, 0) * x) + (fm.get(0, 1) * y) 
            + fm.get(0, 2);
        
        x2 /= d;
        
        double y2 = (fm.get(1, 0) * x) + (fm.get(1, 1) * y) 
            + fm.get(1, 2);
        
        y2 /= d;
        
        return new double[]{x2, y2};
    }
    
    public static void main(String[] args) {
        
        try {
            PointMatcher3Test test = new PointMatcher3Test();

            /*
            test.test1();
            test.test2();
            test.test3();
            test.test4();
            test.test5();
            test.test6();
            test.test7();      
            test.test8();
            test.test9();
            test.test10();          
            test.test11();
            test.test12();
            test.test13();
            test.test14();
            */
            //test.test15();
            //test.test155();
            
            //test.test156();
            
            //test.adjustPointsOfInterest();
            //test.examineInvPointLists();
            //test.smallestSubsets();
            //test.examineIterativeCorners();
            
            test.testSkyline();
            
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

    private void overplotTransformed(PairIntArray transformed, 
        PairIntArray model, int width, int height, int testNumber) 
        throws IOException {
        
        Image image = new Image(width, height);
        
        for (int ii = 0; ii < transformed.getN(); ii++) {
            double x2 = transformed.getX(ii);
            double y2 = transformed.getY(ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, image, 3, 
                0, 0, 255);
        }
        for (int ii = 0; ii < model.getN(); ii++) {
            double x = model.getX(ii);
            double y = model.getY(ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, image, 2, 
                255, 0, 0);
        }
        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_t2_" + testNumber + ".png", image);
    }
    
    private void plot(PairIntArray points, int width, int height, int plotNumber) 
        throws IOException {
        
        Image image = new Image(width, height);
        
        for (int ii = 0; ii < points.getN(); ii++) {
            double x2 = points.getX(ii);
            double y2 = points.getY(ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, image, 2, 
                0, 0, 255);
        }
        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_" + plotNumber + ".png", image);
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
