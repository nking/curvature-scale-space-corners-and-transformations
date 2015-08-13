package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.*;
import java.util.Map.Entry;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ShapeMatcherTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public ShapeMatcherTest() {
    }

    public void estCreateNeighborOffsets() throws Exception {

        ShapeMatcher matcher = new ShapeMatcher();

        float[][] offsets = Misc.createNeighborOffsets(2);

        assertTrue(offsets.length == 25);

        Set<PairInt> offsetsSet = new HashSet<PairInt>();

        for (int i = 0; i < offsets.length; ++i) {

            assertTrue(offsets[i].length == 2);

            float x = offsets[i][0];
            float y = offsets[i][1];

            PairInt p = new PairInt(Math.round(x), Math.round(y));

            assertTrue(Math.abs(x) <= 2);
            assertTrue(Math.abs(y) <= 2);

            offsetsSet.add(p);
        }

        assertTrue(offsetsSet.size() == 25);
    }

    public void estRotateOffsets() throws Exception {

        ShapeMatcher matcher = new ShapeMatcher();

        float[][] offsets = Misc.createNeighborOffsets(2);

        Transformer transformer = new Transformer();

        Map<Float, Set<PairInt>> rotOffsetsMap = new HashMap<Float, Set<PairInt>>();

        float rotDeltaInDegrees = 22.5f;

        for (int i = 0; i < 16; ++i) {

            float rotInDegrees = i * rotDeltaInDegrees;

            Set<PairInt> trOffsets = new HashSet<PairInt>();

            float[][] transformed = transformer.transformXY(rotInDegrees,
                offsets);

            //StringBuilder sb = new StringBuilder("rot=")
            //    .append(Float.toString(rotInDegrees)).append(": ");

            for (int j = 0; j < transformed.length; ++j) {

                int x = Math.round(transformed[j][0]);
                int y = Math.round(transformed[j][1]);

                //sb.append(" (").append(Integer.toString(x)).append(",")
                //    .append(Integer.toString(y)).append(") ");

                PairInt p = new PairInt(x, y);
                trOffsets.add(p);
            }

            rotOffsetsMap.put(Float.valueOf(rotInDegrees), trOffsets);

            //log.info(sb.toString());

        }
    }

    protected ImageExt randomImage(int lowValue, int highValue, int width,
        int height) throws NoSuchAlgorithmException {

        ImageExt img = new ImageExt(width, height);

        int range = highValue - lowValue;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1438231301005L;
        log.info("SEED1=" + seed);
        sr.setSeed(seed);

        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                int r = lowValue + sr.nextInt(range);
                int g = lowValue + sr.nextInt(range);
                int b = lowValue + sr.nextInt(range);
                img.setRGB(i, j, r, g, b);
            }
        }

        return img;
    }

    protected void applyRandomOffsets(int maxOffsetSize, ImageExt img)
        throws NoSuchAlgorithmException {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1438231301105L;
        log.info("SEED2=" + seed);
        sr.setSeed(seed);

        int width = img.getWidth();
        int height = img.getHeight();

        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {

                if (sr.nextBoolean()) {
                    int rNoise = sr.nextInt(maxOffsetSize);
                    int gNoise = sr.nextInt(maxOffsetSize);
                    int bNoise = sr.nextInt(maxOffsetSize);
                    int r = img.getR(i, j);
                    int g = img.getG(i, j);
                    int b = img.getB(i, j);
                    if (sr.nextBoolean()) {
                        r -= rNoise;
                    } else {
                        r += rNoise;
                    }
                    if (sr.nextBoolean()) {
                        g -= gNoise;
                    } else {
                        g += gNoise;
                    }
                    if (sr.nextBoolean()) {
                        b -= bNoise;
                    } else {
                        b += bNoise;
                    }
                    img.setRGB(i, j, r, g, b);
                }
            }
        }
    }

        /*
            switch (ds) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getBrownAndLoweFeatureCentersBinned(points1, points2);
                    binFactor = 3;

                    transfomation for images having been binned by factor 3:

                    params=rotationInRadians=6.1807413 rotationInDegrees=354.13039133201386 scale=0.96686685
                        translationX=-81.54607 translationY=-14.233707 originX=0.0 originY=0.0

                    break;
                }
                case 1: {

                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getVenturiFeatureCentersBinned(points1, points2);
                    binFactor = 4;

                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=0.022482984 rotationInDegrees=1.2881800807667518 scale=0.98079264
                        translationX=-8.983139 translationY=4.063554 originX=0.0 originY=0.0

                    break;
                }
                case 2: {

                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getBooksFeatureCentersBinned(points1, points2);

                    binFactor = 4;

                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=6.1616626 rotationInDegrees=353.0372605364883 scale=1.0766312
                        translationX=-23.24295 translationY=-21.994938 originX=0.0 originY=0.0
            */

    public void estFeatureMatching() throws Exception {
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
            runControlledListMatching(fileName1, fileName2, points1, points2);
        }
    }
    
    protected void runControlledListMatching(String fileName1, String fileName2,
        PairIntArray points1, PairIntArray points2) throws Exception {

        ImageProcessor imageProcessor = new ImageProcessor();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        /*
        testing that a known list of matches are found as "matching" by the
        algorithms.
        */

        BinSegmentationHelper helper = new BinSegmentationHelper(fileName1, fileName2);

        if (fileName1.contains("books")) {
            helper.applySteps1(false);
        } else {
            helper.applySteps1(true);
        }

        helper.createSortedCornerRegions();

        ShapeMatcher matcher = new ShapeMatcher();

        int blockHalfWidth = 5;
        boolean useNormalizedIntensities = true;

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

        log.info("params=" + params);

        // NOTE: these are in the trimmed coordinates reference frame
        // which are the same as the images given to Features
        Set<CornerRegion> cornerRegions1 = helper.getCornerRegions1();
        Set<CornerRegion> cornerRegions2 = helper.getCornerRegions2();

        int dither = 1;// 1 is fine for B&L2003

        Features features1 = new Features(gsImg1,
            gXY1, theta1, blockHalfWidth, useNormalizedIntensities);

        Features features2 = new Features(gsImg2,
            gXY2, theta2, blockHalfWidth, useNormalizedIntensities);

        // iterate over the manual list of corners and find the corner regions
        for (int ii = 0; ii < points1.getN(); ++ii) {

            // the point lists are in the original image reference frame so transform them
            final int x1 = (points1.getX(ii) - gXY1.getXRelativeOffset());
            final int y1 = (points1.getY(ii) - gXY1.getYRelativeOffset());
            final int x2 = (points2.getX(ii) - gXY2.getXRelativeOffset());
            final int y2 = (points2.getY(ii) - gXY2.getYRelativeOffset());

            log.info("looking for corner region for (" + points1.getX(ii)
                + "," + points1.getY(ii) + ")");

            Set<CornerRegion> set1 = new HashSet<CornerRegion>();
            for (CornerRegion cr : cornerRegions1) {
                int kMaxIdx = cr.getKMaxIdx();
                int x = cr.getX()[kMaxIdx];
                int y = cr.getY()[kMaxIdx];
                int xDiff = Math.abs(x1 - x);
                int yDiff = Math.abs(y1 - y);
                if (xDiff < 3 && yDiff < 3) {
                    set1.add(cr);
                    try {
                        cr.getRelativeOrientation();
                    } catch(CornerRegion.CornerRegionDegneracyException e) {
                        //log.severe(e.getMessage());
                    }
                }
            }

            Set<CornerRegion> set2 = new HashSet<CornerRegion>();
            for (CornerRegion cr : cornerRegions2) {
                int kMaxIdx = cr.getKMaxIdx();
                int x = cr.getX()[kMaxIdx];
                int y = cr.getY()[kMaxIdx];
                int xDiff = Math.abs(x2 - x);
                int yDiff = Math.abs(y2 - y);
                if (xDiff < 3 && yDiff < 3) {
                    set2.add(cr);
                    try {
                        cr.getRelativeOrientation();
                    } catch(CornerRegion.CornerRegionDegneracyException e) {
                        //log.severe(e.getMessage());
                    }
                }
            }

            log.info(ii + ") nSet1=" + set1.size() + " nSet2=" + set2.size());

            log.info("expected diffRot=90 to 115");

            FeatureComparisonStat best = null;
                        
            FeatureComparisonStat best3 = null;

            int c1 = 0;
            for (CornerRegion cr1 : set1) {

                FeatureComparisonStat[] stats =
                    matcher.findBestMatch(features1, features2, cr1,
                    set2.toArray(new CornerRegion[set2.size()]), dither);

                if (stats.length > 0) {
                    
                    for (FeatureComparisonStat stat : stats) {

                        log.info(ii + ") stat=" + stat.toString());

                        if (matcher.fitIsBetter(best, stat)) {
                            best = stat;
                        }
                        if (matcher.fitIsBetter3(best3, stat)) {
                            best3 = stat;
                        }
                    }
                }
                c1++;
            }
            assertNotNull(best);
            log.info(ii + ") FINAL best=" + best.toString());
            log.info(ii + ") FINAL best3=" + best3.toString());
        }

        // Now search from points1 through all of corners2 to see if best
        // match is what is expected.

        // iterate over the manual list of corners and find the corner regions
        for (int ii = 0; ii < points1.getN(); ++ii) {

            // the point lists are in the original image reference frame so transform them
            final int x1 = (points1.getX(ii) - gXY1.getXRelativeOffset());
            final int y1 = (points1.getY(ii) - gXY1.getYRelativeOffset());
            final int x2 = (points2.getX(ii) - gXY2.getXRelativeOffset());
            final int y2 = (points2.getY(ii) - gXY2.getYRelativeOffset());

            log.info("looking for corner region for (" + points1.getX(ii)
                + "," + points1.getY(ii) + ")");

            Set<CornerRegion> set1 = new HashSet<CornerRegion>();
            for (CornerRegion cr : cornerRegions1) {
                int kMaxIdx = cr.getKMaxIdx();
                int x = cr.getX()[kMaxIdx];
                int y = cr.getY()[kMaxIdx];
                int xDiff = Math.abs(x1 - x);
                int yDiff = Math.abs(y1 - y);
                if (xDiff < 3 && yDiff < 3) {
                    set1.add(cr);
                }
            }

            FeatureComparisonStat best = null;
            FeatureComparisonStat best3 = null;

            for (CornerRegion cr1 : set1) {

                FeatureComparisonStat[] stats
                    = matcher.findBestMatch(features1, features2, cr1,
                      cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]), 
                      dither);

                if (stats.length > 0) {

                    for (FeatureComparisonStat stat : stats) {

                        log.info(ii + ") stat=" + stat.toString());

                        if (matcher.fitIsBetter(best, stat)) {
                            best = stat;
                        }
                        if (matcher.fitIsBetter3(best3, stat)) {
                            best3 = stat;
                        }
                    }
                }
            }

            assertNotNull(best);

            int diffX = best.getImg2Point().getX() - x2;
            int diffY = best.getImg2Point().getY() - y2;
            
            int diffX3 = best3.getImg2Point().getX() - x2;
            int diffY3 = best3.getImg2Point().getY() - y2;

            log.info(ii + ") FINAL diffX,diffY=(" + diffX + "," + diffY
            + ") best for compare against all corners2=" + best.toString());
            
            log.info(ii + ") FINAL diffX,diffY=(" + diffX3 + "," + diffY3
            + ") best3 for compare against all corners2=" + best3.toString());
            
            boolean t1 = (Math.abs(diffX) < 3) && (Math.abs(diffY) < 3);
            boolean t3 = (Math.abs(diffX3) < 3) && (Math.abs(diffY3) < 3);
            assertTrue(t1 || t3);
        }
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

        ImageProcessor imageProcessor = new ImageProcessor();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        /*
        testing that a known list of matches are found as "matching" by the
        algorithms.
        */

        BinSegmentationHelper helper = new BinSegmentationHelper(fileName1, fileName2);

        if (fileName1.contains("books")) {
            helper.applySteps1(false);
        } else {
            helper.applySteps1(true);
        }

        helper.createSortedCornerRegions();

        int blockHalfWidth = 5;
        boolean useNormalizedIntensities = true;

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

        log.info("params=" + params);
        
        Set<CornerRegion> cornerRegions1 = helper.getCornerRegions1();
        Set<CornerRegion> cornerRegions2 = helper.getCornerRegions2();
        
        ShapeMatcher matcher = new ShapeMatcher();
        
        matcher.findSimilarFeatures(gsImg1, gXY1, theta1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2, gXY2, theta2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]));
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

        ShapeMatcherTest test = new ShapeMatcherTest();

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
