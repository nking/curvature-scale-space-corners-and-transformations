package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.*;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FeatureMatcher2Test extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public FeatureMatcher2Test() {
    }

    public void testCorrespondenceFromFeatures() throws Exception {
        PairIntArray points1, points2;
        String fileName1, fileName2;
        
        for (int i = 0; i < 3; ++i) {
            points1 = new PairIntArray();
            points2 = new PairIntArray();
            switch(i) {
                case 0: {
                    getBrownAndLoweFeatureCenters90(points1, points2);
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 1: {
                    getVenturiFeatureCenters90(points1, points2);
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    break;
                }
                default: {
                    getBooksFeatureCenters90(points1, points2);
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
            }
            
            runCorrespondenceFromFeatures(fileName1, fileName2, points1, points2);
        }
    }
    
    protected void runCorrespondenceFromFeatures(
        String fileName1, String fileName2,
        PairIntArray points1, PairIntArray points2) throws Exception {
        
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        //testing that a known list of matches are found as "matching" by the algorithms.

        BinSegmentationHelper helper = new BinSegmentationHelper(fileName1, fileName2);
        helper.applySteps3();

        helper.createSortedCornerRegions();

        final GreyscaleImage gsImg1 = helper.getGreyscaleImage1();
        final GreyscaleImage gXY1 = helper.getGXY1();
        final GreyscaleImage theta1 = helper.getTheta1();
        final GreyscaleImage gsImg2 = helper.getGreyscaleImage2();
        final GreyscaleImage gXY2 = helper.getGXY2();
        final GreyscaleImage theta2 = helper.getTheta2();

        assertTrue(gsImg1.getWidth() == gXY1.getWidth());
        assertTrue(gsImg1.getHeight() == gXY1.getHeight());
        assertTrue(gsImg2.getWidth() == gXY2.getWidth());
        assertTrue(gsImg2.getHeight() == gXY2.getHeight());
        assertTrue(gsImg1.getWidth() == theta1.getWidth());
        assertTrue(gsImg1.getHeight() == theta1.getHeight());
        assertTrue(gsImg2.getWidth() == theta2.getWidth());
        assertTrue(gsImg2.getHeight() == theta2.getHeight());

        TransformationParameters params = tc.calulateEuclidean(
            points1.getX(0), points1.getY(0), points1.getX(1), points1.getY(1),
            points2.getX(0), points2.getY(0), points2.getX(1), points2.getY(1),
            0, 0);

        log.info("expected params=" + params);
        
        BlobScaleFinderWrapper scaleFinder = new BlobScaleFinderWrapper(
            helper.img1Orig, helper.img2Orig);
            
        scaleFinder.setToDebug();

        TransformationParameters params0 = scaleFinder.calculateScale();
        log.info("params from scale finder=" + params0);
        
        Set<CornerRegion> cornerRegions1 = helper.getCornerRegions1();
        Set<CornerRegion> cornerRegions2 = helper.getCornerRegions2();
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        float scale = params0.getScale();
        
        CorrespondenceList cl = matcher.findSimilarFeatures(gsImg1, gXY1, theta1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2, gXY2, theta2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]), 
            scale);
        
        Collection<PairInt> m1 = cl.getPoints1();
        Collection<PairInt> m2 = cl.getPoints2();
        MiscDebug.plotCorners(gsImg1.copyImage(), m1, "1_" + fileName1Root  + "_matched", 2);
        MiscDebug.plotCorners(gsImg2.copyImage(), m2, "2_" + fileName2Root + "_matched", 2);
        
        log.info(" " + m1.size() + " points for first correspondence");
        
        assertTrue(Math.abs(params.getRotationInDegrees() 
            - cl.getRotationInDegrees()) < cl.getRangeRotation());
        
        assertTrue(Math.abs(params.getTranslationX()
            - cl.getTranslationX()) < cl.getRangeTranslationX());
        
        assertTrue(Math.abs(params.getTranslationY()
            - cl.getTranslationY()) < cl.getRangeTranslationY());
        
        MiscDebug.writeImage(cornerRegions1, ImageIOHelper.convertImage(gsImg1),
            "1_" + fileName1Root + "_cr");
        MiscDebug.writeImage(cornerRegions2, ImageIOHelper.convertImage(gsImg2), 
            "2_" + fileName1Root + "_cr");
    }
    
    protected static void getBrownAndLoweFeatureCentersBinned(PairIntArray out1,
        PairIntArray out2) {
        out1.add(144, 67);
        out2.add(52, 64);
        out1.add(150, 68);
        out2.add(56, 66);
        out1.add(148, 79);
        out2.add(53, 76);
        out1.add(165, 85);
        out2.add(66, 83);
        out1.add(128, 87);
        out2.add(32, 82);
    }

    protected static void getBrownAndLoweFeatureCenters90(PairIntArray out1,
        PairIntArray out2) {

        out1.add(206, 65);
        out2.add(170, 200);

        //difficult test point to try after method improved:
        //out1.add(165, 186);
        //out2.add(55, 139);
        out1.add(161, 182);//best match is 162,183
        out2.add(61, 137);

        out1.add(221, 219);
        out2.add(9, 193);
        out1.add(170, 37);
        out2.add(200, 171);
        out1.add(316, 50);
        out2.add(164, 305);

        out1.add(183, 187);
        out2.add(54, 158);
    }

    private void getBrownAndLoweFeatureCentersBinned(List<PairInt> points1,
        List<PairInt> points2) {
        points1.add(new PairInt(128, 87));
        points2.add(new PairInt(32, 82));
        points1.add(new PairInt(150, 68));
        points2.add(new PairInt(56, 66));
        points1.add(new PairInt(148, 79));
        points2.add(new PairInt(53, 76));
        points1.add(new PairInt(165, 85));
        points2.add(new PairInt(66, 83));
        points1.add(new PairInt(144, 67));
        points2.add(new PairInt(52, 64));
    }

    private void getVenturiFeatureCentersBinned(PairIntArray points1,
        PairIntArray points2) {
        points1.add(152, 88);
        points2.add(142, 87);

        points1.add(56, 94);
        points2.add(48, 95);

        points1.add(147, 46);
        points2.add(137, 46);

        points1.add(37, 91);
        points2.add(28, 93);

    }

    private void getVenturiFeatureCentersBinned(List<PairInt> points1,
        List<PairInt> points2) {
        points1.add(new PairInt(152, 88));
        points2.add(new PairInt(142, 87));

        points1.add(new PairInt(56, 94));
        points2.add(new PairInt(48, 95));

        points1.add(new PairInt(147, 46));
        points2.add(new PairInt(137, 46));

        points1.add(new PairInt(37, 91));
        points2.add(new PairInt(28, 93));

    }

    private void getBooksFeatureCentersBinned(PairIntArray points1,
        PairIntArray points2) {
        points1.add(158, 20);
        points2.add(143, 20);

        points1.add(98, 134);
        points2.add(64, 134);

        points1.add(139, 110);
        points2.add(109, 110);

        points1.add(154, 99);
        points2.add(122, 99);

        points1.add(133, 21);
        points2.add(118, 21);
    }

    private void getBooksFeatureCentersBinned(List<PairInt> points1,
        List<PairInt> points2) {

        points1.add(new PairInt(158, 20));
        points2.add(new PairInt(143, 20));

        points1.add(new PairInt(98, 134));
        points2.add(new PairInt(64, 134));

        points1.add(new PairInt(139, 110));
        points2.add(new PairInt(109, 110));

        points1.add(new PairInt(154, 99));
        points2.add(new PairInt(122, 99));

        points1.add(new PairInt(133, 21));
        points2.add(new PairInt(118, 21));

    }

    public static void main(String[] args) {

        FeatureMatcher2Test test = new FeatureMatcher2Test();

        try {
            //test.testCalculateStat2();
            //test.testFeatureMatching();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }

    private void getVenturiFeatureCenters90(PairIntArray out1, PairIntArray out2) {

        out1.add(415, 538);
        out2.add(61, 422);

        // point2 has badly determined orientation
        out1.add(272, 575);
        out2.add(23, 275);

        out1.add(259, 366);
        out2.add(238, 261);

        out1.add(308, 67);
        out2.add(534, 305);

    }

    private void getBooksFeatureCenters90(PairIntArray out1, PairIntArray out2) {
        
        out1.add(7, 57);
        out2.add(581, 8);
        
        out1.add(370, 68);
        out2.add(496, 371);
        
        /* nearby region outside of corner's object changes greatly, so this is 
        a good test point for theta within a contour, but is not matching 
        the entire block because of the differences outside of object due 
        to projection
        
        out1.add(18, 222);
        out2.add(410, 16);
        
        out1.add(449, 163);
        out2.add(413, 449);
        
        out1.add(131, 339);
        out2.add(292, 131);
        */
        
        out1.add(463, 268);
        out2.add(318, 464);
        
        out1.add(547, 334);
        out2.add(215, 547);
        
        out1.add(488, 495);
        out2.add(67, 488);
        
        out1.add(500, 487);
        out2.add(74, 500);
        
        out1.add(215, 413);
        out2.add(216, 215);
        
        out1.add(126, 388);
        out2.add(241, 126);
      
        out1.add(78, 290);
        out2.add(339, 78);
        
        out1.add(92, 308);
        out2.add(323, 92);
        
        // difficult corner:
        out1.add(254, 96);
        out2.add(519, 253);
        
    }
}
