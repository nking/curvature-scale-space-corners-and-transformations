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

        float[][] offsets = matcher.createNeighborOffsets();

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

        float[][] offsets = matcher.createNeighborOffsets();

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

    public void testRegionOrientation() throws Exception {

        ImageProcessor imageProcessor = new ImageProcessor();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        /*
        testing that a known list of matches are found as "matching" by the
        algorithms.
        */
        
        PairIntArray points1 = new PairIntArray();
        PairIntArray points2 = new PairIntArray();
        getBrownAndLoweFeatureCenters90(points1, points2);
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        //String fileName1 = "venturi_mountain_j6_0001.png";
        //String fileName2 = "venturi_mountain_j6_0010.png";
        //getVenturiFeatureCenters90(points1, points2);
        
        //binFactor = 3;
        /*
        transfomation for images having been binned by factor 3:

        params=rotationInRadians=6.1807413 rotationInDegrees=354.13039133201386 scale=0.96686685
            translationX=-81.54607 translationY=-14.233707 originX=0.0 originY=0.0
        */
        
        TransformationParameters params = tc.calulateEuclidean(
            points1.getX(0), points1.getY(0), points1.getX(1), points1.getY(1),
            points2.getX(0), points2.getY(0), points2.getX(1), points2.getY(1),
            0, 0);

        log.info("params=" + params);
            
        BinSegmentationHelper helper = new BinSegmentationHelper(fileName1, fileName2);
        
        if (fileName1.contains("brown")) {
            //helper.applySteps0();
            helper.applySteps1();
        } else {
            helper.applySteps1();
        }
        
        helper.createSortedCornerRegions();
        
        ShapeMatcher matcher = new ShapeMatcher();
        
        int blockHalfWidth = 2;
        boolean useNormalizedIntensities = false;
        
        ImageExt img1 = helper.getImage1(); 
        GreyscaleImage gXY1 = helper.getGXY1();
        GreyscaleImage theta1 = helper.getTheta1();
        ImageExt img2 = helper.getImage2(); 
        GreyscaleImage gXY2 = helper.getGXY2();
        GreyscaleImage theta2 = helper.getTheta2();
        float[][] offsets0 = Misc.createNeighborOffsets(2);
        
        Set<CornerRegion> cornerRegions1 = helper.getCornerRegions1();
        Set<CornerRegion> cornerRegions2 = helper.getCornerRegions2();
     
        int dither = 1;
        
        Features features1 = new Features(helper.getGreyscaleImage1(), 
            gXY1, theta1, blockHalfWidth, useNormalizedIntensities);
        
        Features features2 = new Features(helper.getGreyscaleImage2(), 
            gXY2, theta2, blockHalfWidth, useNormalizedIntensities);
        
        
        // iterate over the manual list of corners and find the corner regions
        for (int ii = 0; ii < points1.getN(); ++ii) {

            final int x1 = points1.getX(ii);
            final int y1 = points1.getY(ii);
            final int x2 = points2.getX(ii);
            final int y2 = points2.getY(ii);

            log.info("looking for corner region for (" + x1 + "," + y1 + ")");

            Set<CornerRegion> set1 = new HashSet<CornerRegion>();
            for (CornerRegion cr : cornerRegions1) {
                int kMaxIdx = cr.getKMaxIdx();
                int x = cr.getX()[kMaxIdx];
                int y = cr.getY()[kMaxIdx];
                int xDiff = Math.abs(x1 - x);
                int yDiff = Math.abs(y1 - y);
                if (xDiff < 5 && yDiff < 5) {
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
                if (xDiff < 5 && yDiff < 5) {
                    set2.add(cr);
                    try {
                        cr.getRelativeOrientation();
                    } catch(CornerRegion.CornerRegionDegneracyException e) {
                        log.severe(e.getMessage());
                    }
                }
            }
            
            log.info("nSet1=" + set1.size() + " nSet2=" + set2.size());
            
            log.info("expected diffRot=90 to 115");
            
            FeatureComparisonStat best = null;
            
            int c1 = 0;
            for (CornerRegion cr1 : set1) {
                int c2 = 0;
                for (CornerRegion cr2 : set2) {
                    
                    FeatureComparisonStat stat = null;
                    
                    try {
                        // discard the wrong matches while printing stats of expected
                        // matches
                        float diffRot = AngleUtil.getAngleDifference(
                            cr1.getRelativeOrientationInDegrees(), 
                            cr2.getRelativeOrientationInDegrees());
                        if (Math.abs(Math.abs(diffRot) - 100) > 20) {
                            continue;
                        }
                        
                        stat = matcher.findBestAmongDitheredRotated(
                            features1, features2, cr1, cr2, dither);
                    } catch(CornerRegion.CornerRegionDegneracyException e) {
                        log.severe(e.getMessage());
                    }

                    //TODO: need to use gradient in further fitness function too
                    
                    if (stat != null) {
                        log.info(ii + ") stat=" + stat.toString());
                        
                        if (best == null) {
                            best = stat;
                        } else {
                            if ((best.getSumIntensitySqDiff() >= stat.getSumIntensitySqDiff())
                                && (best.getSumGradientSqDiff() > stat.getSumGradientSqDiff())
                                && (best.getSumThetaDiff() > stat.getSumThetaDiff())
                                ) {
                                best = stat;
                            }
                        }
                    }
                    c2++;
                }
                c1++;
            }
            assertNotNull(best);
            log.info(ii + ") FINAL best=" + best.toString());
        }
       
        // Now search from points1 through all of corners2 to see if best 
        // match is what is expected.
      
        // iterate over the manual list of corners and find the corner regions
        for (int ii = 0; ii < points1.getN(); ++ii) {

            final int x1 = points1.getX(ii);
            final int y1 = points1.getY(ii);
            final int x2 = points2.getX(ii);
            final int y2 = points2.getY(ii);

            log.info("looking for corner region for (" + x1 + "," + y1 + ")");

            Set<CornerRegion> set1 = new HashSet<CornerRegion>();
            for (CornerRegion cr : cornerRegions1) {
                int kMaxIdx = cr.getKMaxIdx();
                int x = cr.getX()[kMaxIdx];
                int y = cr.getY()[kMaxIdx];
                int xDiff = Math.abs(x1 - x);
                int yDiff = Math.abs(y1 - y);
                if (xDiff < 5 && yDiff < 5) {
                    set1.add(cr);
                }
            }
            
            FeatureComparisonStat best = null;
            float bestGradientSSD = Float.MAX_VALUE;
            float secondBestGradientSSD = Float.MAX_VALUE - 1;
                        
            for (CornerRegion cr1 : set1) {                                
                for (CornerRegion cr2 : cornerRegions2) {
                    FeatureComparisonStat stat = null;
                    try {
                        stat = matcher.findBestAmongDitheredRotated(
                            features1, features2, cr1, cr2, dither);
                    } catch(CornerRegion.CornerRegionDegneracyException e) {
                        //log.severe(e.getMessage());
                    }
                    if (stat != null) {
                        log.info(ii + ") stat=" + stat.toString());

                        if (best == null) {
                            best = stat;
                        } else {
                            if ((best.getSumIntensitySqDiff() >= stat.getSumIntensitySqDiff())
                                && (best.getSumGradientSqDiff() > stat.getSumGradientSqDiff())
                                && (best.getSumThetaDiff() > stat.getSumThetaDiff())
                                ) {
                                best = stat;
                            }
                        }
                    }
                }                
            }
            
            int diffX = best.getImg2Point().getX() - x2;
            int diffY = best.getImg2Point().getY() - y2;
            
            assertNotNull(best);
            
            log.info(ii + ") FINAL diffX,diffY=(" + diffX + "," + diffY 
            + ") best for compare against all corners2=" + best.toString());
        }
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
        out1.add(206, 66);
        out2.add(168, 200);
        out1.add(331, 167);
        out2.add(43, 313);
        out1.add(165, 187);
        out2.add(55, 139);
        out1.add(220, 220);
        out2.add(9, 194);
        out1.add(170, 37);
        out2.add(200, 171);
        out1.add(316, 51);
        out2.add(164, 305);
        
        out1.add(165, 186);
        out2.add(55, 139);
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
            test.testRegionOrientation();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }

    private void getVenturiFeatureCenters90(PairIntArray out1, PairIntArray out2) {
        
        out1.add(291, 493);
        out2.add(106, 296);
        
        out1.add(447, 221);
        out2.add(384, 446);
        
        
    }

}
