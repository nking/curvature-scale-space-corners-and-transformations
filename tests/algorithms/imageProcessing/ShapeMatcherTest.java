package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
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

    public void testCreateNeighborOffsets() throws Exception {

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

    public void testRotateOffsets() throws Exception {

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

    public void testCalculateStat() throws Exception {

        ShapeMatcher matcher = new ShapeMatcher();

        // fill a small image randomly from 0 to 255
        // copy that image to img2 then make small random offsets for some
        // pixels.
        ImageExt img1 = randomImage(10, 255, 10, 10);
        ImageExt img2 = (ImageExt)img1.copyImage();
        applyRandomOffsets(2, img2);

        float[][] offsets0 = matcher.createNeighborOffsets();

        FeatureComparisonStat stat = matcher.calculateStat(img1, img2,
            5, 5, 5, 5, offsets0, offsets0);
        assertTrue((stat.getSumSqDiff()/stat.getImg2PointErr()) < 1);

        // copy image1 again, but mult all pixels in img2 by 0.75
        img2 = (ImageExt)img1.copyImage();
        for (int i = 0; i < img2.nPixels; ++i) {
            int r = Math.round(0.75f * img2.getR(i));
            int g = Math.round(0.75f * img2.getG(i));
            int b = Math.round(0.75f * img2.getB(i));
            int col = img2.getCol(i);
            int row = img2.getRow(i);
            img2.setRGB(col, row, r, g, b);
        }
        stat = matcher.calculateStat(img1, img2,
            5, 5, 5, 5, offsets0, offsets0);
        double maxDiff = (255 - 0.75*255);
        assertTrue((stat.getSumSqDiff()/stat.getImg2PointErr()) < 1);

        // copy image 1 and subtract a known amount
        img2 = (ImageExt)img1.copyImage();
        for (int i = 0; i < img2.nPixels; ++i) {
            int r = Math.round(img2.getR(i) - 2);
            int g = Math.round(img2.getG(i) - 2);
            int b = Math.round(img2.getB(i) - 2);
            int col = img2.getCol(i);
            int row = img2.getRow(i);
            img2.setRGB(col, row, r, g, b);
        }
        stat = matcher.calculateStat(img1, img2,
            5, 5, 5, 5, offsets0, offsets0);
        assertTrue((stat.getSumSqDiff()/stat.getImg2PointErr()) < 1);


        // rotate image2 by 67.5 degrees
        // fetch the offsets for transformation
        // compare to the results of previous test
        img2 = (ImageExt)img1.copyImage();
        FeatureComparisonStat stat0 = matcher.calculateStat(img1, img2,
            5, 5, 5, 5, offsets0, offsets0);

        TransformationParameters params = new TransformationParameters();
        params.setOriginX(5);
        params.setOriginY(5);
        params.setRotationInDegrees(67.5f);
        params.setTranslationX(0);
        params.setTranslationY(0);
        params.setScale(1);

        Transformer transformer = new Transformer();
        Image img = transformer.applyTransformation(img1, params, img1.getWidth(),
            img1.getHeight());

        float[][] offset67point5 = transformer.transformXY(67.5f, offsets0);

        stat = matcher.calculateStat(img, img2, 5, 5, 5, 5, offset67point5,
            offsets0);

        assertTrue((stat0.getSumSqDiff()/stat0.getImg2PointErr()) < 1);
    }

    public void testCalculateStat2() throws Exception {

        ImageProcessor imageProcessor = new ImageProcessor();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        // choose regions from Brown & Lowe images to compare
        String fileName1 = null;
        String fileName2 = null;
        ImageExt img1 = null;
        ImageExt img2 = null;
        List<PairInt> points1 = new ArrayList<PairInt>();
        List<PairInt> points2 = new ArrayList<PairInt>();
        int binFactor = 1;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1438231301005L;
        log.info("SEED3=" + seed);
        sr.setSeed(seed);

        for (int ds = 0; ds < 3; ++ds) {

            switch (ds) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    points1 = new ArrayList<PairInt>();
                    points2 = new ArrayList<PairInt>();
                    getBrownAndLoweFeatureCenters(points1, points2);
                    binFactor = 3;
                    /*
                    transfomation for images having been binned by factor 3:

                    params=rotationInRadians=6.1807413 rotationInDegrees=354.13039133201386 scale=0.96686685
                        translationX=-81.54607 translationY=-14.233707 originX=0.0 originY=0.0
                    */
                    break;
                }
                case 1: {

                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    points1 = new ArrayList<PairInt>();
                    points2 = new ArrayList<PairInt>();
                    getVenturiFeatureCenters(points1, points2);
                    
                    binFactor = 4;
                    /*
                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=0.022482984 rotationInDegrees=1.2881800807667518 scale=0.98079264
                        translationX=-8.983139 translationY=4.063554 originX=0.0 originY=0.0
                    */
                    break;
                }
                case 2: {

                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    points1 = new ArrayList<PairInt>();
                    points2 = new ArrayList<PairInt>();
                    getBooksFeatureCenters(points1, points2);
                    binFactor = 4;
                    /*
                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=6.1616626 rotationInDegrees=353.0372605364883 scale=1.0766312
                        translationX=-23.24295 translationY=-21.994938 originX=0.0 originY=0.0
                    */
                    break;
                }
            }

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            img1 = ImageIOHelper.readImageExt(filePath1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            img2 = ImageIOHelper.readImageExt(filePath2);
            img1 = imageProcessor.binImage(img1, binFactor);
            img2 = imageProcessor.binImage(img2, binFactor);

            ShapeMatcher matcher = new ShapeMatcher();

            TransformationParameters params = tc.calulateEuclidean(
                points1.get(0).getX(), points1.get(0).getY(), points1.get(1).getX(), points1.get(1).getY(),
                points2.get(0).getX(), points2.get(0).getY(), points2.get(1).getX(), points2.get(1).getY(),
                0, 0);

            log.info("params=" + params);

            Transformer transformer = new Transformer();

            Set<FeatureComparisonStat> trueMatches = new HashSet<FeatureComparisonStat>();
            Set<FeatureComparisonStat> differentPatches = new HashSet<FeatureComparisonStat>();

            float[][] offsets0 = matcher.createNeighborOffsets();
            float[][] offsetsT = transformer.transformXY(
                Math.round(params.getRotationInDegrees()), offsets0);

            FeatureComparisonStat stat1 = matcher.calculateStat(img1, img2,
                points1.get(0).getX(), points1.get(0).getY(),
                points2.get(0).getX(), points2.get(0).getY(),
                offsetsT, offsets0);
            trueMatches.add(stat1);

            FeatureComparisonStat stat2 = matcher.calculateStat(img1, img2,
                points1.get(1).getX(), points1.get(1).getY(),
                points2.get(1).getX(), points2.get(1).getY(),
                offsetsT, offsets0);
            trueMatches.add(stat2);

            log.info("true match =" + stat1.toString());
            log.info("true match =" + stat2.toString());

            for (int i = 2; i < points1.size(); ++i) {
                FeatureComparisonStat stat = matcher.calculateStat(img1, img2,
                    points1.get(i).getX(), points1.get(i).getY(),
                    points2.get(i).getX(), points2.get(i).getY(),
                    offsetsT, offsets0);
                trueMatches.add(stat);
                log.info("true match =" + stat.toString());
            }

            /*
             INFO: true match =p1=x=117 y=111 p2=x=14 y=105
             avgDiffPix=21.480001 stDevDiffPix=19.142578 avgDivPix=1.0866014 stDevDivPix=0.2703045

             INFO: true match =p1=x=156 y=52 p2=x=64 y=52
             avgDiffPix=15.88 stDevDiffPix=10.870372 avgDivPix=0.9379496 stDevDivPix=0.10855621
            */

            /*
            visit 100 other places randomly in the image avoiding correct match
            and store the stats.

            // offsets from that transformation should not match:
            */

            for (int i = 0; i < 100; i++) {

                int x2 = 7 + sr.nextInt(img2.getWidth() - 12);
                int y2 = 7 + sr.nextInt(img2.getHeight() - 12);
                while ((Math.abs(x2 - points2.get(0).getX()) < 6) &&
                    (Math.abs(y2 - points2.get(0).getY()) < 6)) {
                    x2 = 7 + sr.nextInt(img2.getWidth() - 12);
                    y2 = 7 + sr.nextInt(img2.getHeight() - 12);
                }

                int rotD = sr.nextInt(360);

                float[][] offsetsR = transformer.transformXY(rotD, offsets0);

                FeatureComparisonStat s = matcher.calculateStat(img1, img2,
                    points1.get(0).getX(), points1.get(0).getY(), x2, y2, offsetsR, offsets0);

                differentPatches.add(s);

                float div = s.getSumSqDiff()/s.getImg2PointErr();
                if (div <= 1) {
                    log.info(i + " (expected sumSqDiff/err > 1) for not same blocks ="
                        + s.toString());
                }
            }
        }
    }

    public void testDitherToFindSmallestSqSumDiff() throws Exception {

        ImageProcessor imageProcessor = new ImageProcessor();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        // choose regions from Brown & Lowe images to compare
        String fileName1 = null;
        String fileName2 = null;
        ImageExt img1 = null;
        ImageExt img2 = null;
        PairIntArray points1 = null;
        PairIntArray points2 = null;
        int binFactor = 1;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1438231301005L;
        log.info("SEED3=" + seed);
        sr.setSeed(seed);

        for (int ds = 0; ds < 3; ++ds) {

            switch (ds) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getBrownAndLoweFeatureCenters(points1, points2);
                    binFactor = 3;
                    /*
                    transfomation for images having been binned by factor 3:

                    params=rotationInRadians=6.1807413 rotationInDegrees=354.13039133201386 scale=0.96686685
                        translationX=-81.54607 translationY=-14.233707 originX=0.0 originY=0.0
                    */
                    break;
                }
                case 1: {

                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getVenturiFeatureCenters(points1, points2);
                    binFactor = 4;
                    /*
                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=0.022482984 rotationInDegrees=1.2881800807667518 scale=0.98079264
                        translationX=-8.983139 translationY=4.063554 originX=0.0 originY=0.0
                    */
                    break;
                }
                case 2: {

                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getBooksFeatureCenters(points1, points2);

                    binFactor = 4;
                    /*
                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=6.1616626 rotationInDegrees=353.0372605364883 scale=1.0766312
                        translationX=-23.24295 translationY=-21.994938 originX=0.0 originY=0.0
                    */
                    break;
                }
            }

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            img1 = ImageIOHelper.readImageExt(filePath1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            img2 = ImageIOHelper.readImageExt(filePath2);
            img1 = imageProcessor.binImage(img1, binFactor);
            img2 = imageProcessor.binImage(img2, binFactor);

            ShapeMatcher matcher = new ShapeMatcher();

            TransformationParameters params = tc.calulateEuclidean(
                points1.getX(0), points1.getY(0), points1.getX(1), points1.getY(1),
                points2.getX(0), points2.getY(0), points2.getX(1), points2.getY(1),
                0, 0);

            log.info("params=" + params);


            Transformer transformer = new Transformer();

            int dither = 1;

            float[][] offsets0 = matcher.createNeighborOffsets();
            float[][] offsetsT = transformer.transformXY(
                Math.round(params.getRotationInDegrees()), offsets0);

            if (fileName1.contains("brown_lowe")) {

                /* points for dither test:
                (115, 82)   (23, 76)
                (139, 59)   (48, 56)
                */
                int x1 = 115; int y1 = 82;
                int x2 = 23; int y2 = 76;
                int dither2 = 5;
            
                FeatureComparisonStat stat0 = matcher.ditherToFindSmallestSqSumDiff(
                    img1, img2, x1, y1, x2, y2, offsetsT, offsets0, dither2);

                int deltaX = x1 - stat0.getImg1Point().getX();
                int deltaY = y1 - stat0.getImg1Point().getY();

                log.info("*deltaX=" + deltaX + " deltaY=" + deltaY + " stat=" + stat0);
                assertTrue((deltaX != 0) || (deltaY != 0));

                x1 = 139; y1 = 59;
                x2 = 49; y2 = 56;
                FeatureComparisonStat stat1 = matcher.ditherToFindSmallestSqSumDiff(
                    img1, img2, x1, y1, x2, y2, offsetsT, offsets0, dither2);

                deltaX = x1 - stat1.getImg1Point().getX();
                deltaY = y1 - stat1.getImg1Point().getY();

                log.info("*deltaX=" + deltaX + " deltaY=" + deltaY + " stat=" + stat1);
                assertTrue((deltaX != 0) || (deltaY != 0));

            }

            for (int i = 0; i < points1.getN(); ++i) {

                int x1 = points1.getX(i);
                int y1 = points1.getY(i);

                int x2 = points2.getX(i);
                int y2 = points2.getY(i);

                FeatureComparisonStat stat = matcher.ditherToFindSmallestSqSumDiff(
                    img1, img2, x1, y1, x2, y2, offsetsT, offsets0, dither);

                float div = stat.getSumSqDiff()/stat.getImg2PointErr();
                log.info(stat.toString());
            
                assertTrue(div <= 1);

                if (div > 1) {
                    log.info("ERROR: expected to see div < 1");
                }

                int deltaX = x1 - stat.getImg1Point().getX();
                int deltaY = y1 - stat.getImg1Point().getY();

                log.info("  deltaX=" + deltaX + " deltaY=" + deltaY);

                assertTrue(deltaX == 0);
                assertTrue(deltaY == 0);

                // move center by dx, dy and assert that algorithm finds the true
                // center.  make dither large to test above rounding errors (+-1)
                int dither2 = 2;
                int dx = 2*sr.nextInt(1);
                int dy = 2*sr.nextInt(1);
                if (sr.nextBoolean()) {
                    dx *= -1;
                }
                if (sr.nextBoolean()) {
                    dy *= -1;
                }
                if (dx != 0 || dy != 0) {
                    stat = matcher.ditherToFindSmallestSqSumDiff(
                        img1, img2, (x1 + dx), (y1 + dy), x2, y2, offsetsT, offsets0, 
                        dither2);
                    deltaX = x1 - stat.getImg1Point().getX();
                    deltaY = y1 - stat.getImg1Point().getY();
                    //log.info("stat=" + stat);
                    //log.info("  deltaX=" + deltaX + " deltaY=" + deltaY 
                    //    + " dx=" + dx + " dy=" + dy);
                    assertTrue(Math.abs(deltaX - dx) < 2);
                    assertTrue(Math.abs(deltaY - dy) < 2);
                }
            }
        }
        
    }
    
    public void testFindSimilarFeaturesForRotatedFrames() throws Exception {

        ImageProcessor imageProcessor = new ImageProcessor();

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        // choose regions from Brown & Lowe images to compare
        String fileName1 = null;
        String fileName2 = null;
        ImageExt img1 = null;
        ImageExt img2 = null;
        PairIntArray points1 = null;
        PairIntArray points2 = null;
        int binFactor = 1;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1438231301005L;
        log.info("SEED3=" + seed);
        sr.setSeed(seed);

        for (int ds = 0; ds < 3; ++ds) {

            switch (ds) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getBrownAndLoweFeatureCenters(points1, points2);
                    binFactor = 3;
                    /*
                    transfomation for images having been binned by factor 3:

                    params=rotationInRadians=6.1807413 rotationInDegrees=354.13039133201386 scale=0.96686685
                        translationX=-81.54607 translationY=-14.233707 originX=0.0 originY=0.0
                    */
                    break;
                }
                case 1: {

                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getVenturiFeatureCenters(points1, points2);
                    binFactor = 4;
                    /*
                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=0.022482984 rotationInDegrees=1.2881800807667518 scale=0.98079264
                        translationX=-8.983139 translationY=4.063554 originX=0.0 originY=0.0
                    */
                    break;
                }
                case 2: {

                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    points1 = new PairIntArray();
                    points2 = new PairIntArray();
                    getBooksFeatureCenters(points1, points2);

                    binFactor = 4;
                    /*
                    transfomation for images having been binned by factor binFactor:

                    params=rotationInRadians=6.1616626 rotationInDegrees=353.0372605364883 scale=1.0766312
                        translationX=-23.24295 translationY=-21.994938 originX=0.0 originY=0.0
                    */
                    break;
                }
            }

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            img1 = ImageIOHelper.readImageExt(filePath1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            img2 = ImageIOHelper.readImageExt(filePath2);
            img1 = imageProcessor.binImage(img1, binFactor);
            img2 = imageProcessor.binImage(img2, binFactor);

            ShapeMatcher matcher = new ShapeMatcher();

            TransformationParameters params = tc.calulateEuclidean(
                points1.getX(0), points1.getY(0), points1.getX(1), points1.getY(1),
                points2.getX(0), points2.getY(0), points2.getX(1), points2.getY(1),
                0, 0);

            log.info("params=" + params);

            Map<PairInt, Map<PairInt, Map<Float, FeatureComparisonStat>>> 
                comparisonMap = matcher.findSimilarFeaturesForRotatedFrames(
                img1, img2, points1, points2);
        
            /*
            looking at patterns to determine best solution for all.
            The most distinct best answer(s) may be what can use for filtering
            the map by nearby rotation, then further filtering remaining
            entries using pairwise calculation.
            then further comparing with stat diff sums?
            */

            int p1Count = 0;
            
            for (Entry<PairInt, Map<PairInt, Map<Float, FeatureComparisonStat>>> entry :
                comparisonMap.entrySet()) {

                PairInt p1 = entry.getKey();
                int expectedIndex = -1;
                for (int i = 0; i < points1.getN(); ++i) {
                    int x = points1.getX(i);
                    int y = points1.getY(i);
                    if ((x == p1.getX()) && (y == p1.getY())) {
                        expectedIndex = i;
                        break;
                    }
                }            
                int expectedX2 = points2.getX(expectedIndex);
                int expectedY2 = points2.getY(expectedIndex);
                float expectedRotDeg = params.getRotationInDegrees();

                Map<PairInt, Map<Float, FeatureComparisonStat>> p1PairsMap = entry.getValue();

                FixedSizeSortedVector<FeatureComparisonStat> vec = 
                    new FixedSizeSortedVector<FeatureComparisonStat>(6, 
                        FeatureComparisonStat.class);

                for (PairInt p2 : p1PairsMap.keySet()) {
                    Map<Float, FeatureComparisonStat> p1P2Map = p1PairsMap.get(p2);
                    for (Entry<Float, FeatureComparisonStat> entry3 : p1P2Map.entrySet()) {
                        FeatureComparisonStat stat = entry3.getValue();
                        float div = stat.getSumSqDiff()/stat.getImg2PointErr();
                        if (div <= 1) {
                            vec.add(stat);
                        }
                    }
                }

                //looks like the top ? best always contain the expected point
                boolean foundBest = false;
                int bestIdx = -1;
                FeatureComparisonStat bestStat = null;
                
                for (int i = 0; i < vec.getNumberOfItems(); ++i) {
                    FeatureComparisonStat fs = vec.getArray()[i];
                    int deltaX = Math.abs(fs.getImg2Point().getX() - expectedX2);
                    int deltaY = Math.abs(fs.getImg2Point().getY() - expectedY2);
                    float diffDeg = AngleUtil.getAngleDifference(
                        fs.getImg1RotInDegrees(), expectedRotDeg);
                    if ((deltaX <= 2) && (deltaY <= 2) && (Math.abs(diffDeg) < 22.5f)) {
                        foundBest = true;
                        bestIdx = i;
                        bestStat = fs;
                    }
                }
                
                String str = (bestStat == null) ? 
                    String.format(
                        "%d) best at sorted idx=%d ssd=%.2f er=%.2f rotD=%.2f sep=%d",
                        ds, bestIdx)
                    : String.format(
                        "%d) best at sorted idx=%d ssd=%.2f err=%.2f rotD=%.2f sep=%d",
                        ds, bestIdx, bestStat.getSumSqDiff(), 
                        bestStat.getImg2PointErr(), bestStat.getImg1RotInDegrees(),
                        Math.round(bestStat.getSqDist()));
                //log.info(str);
                if (!foundBest) {
                    // look at p1PairsMap
                    int z = 1;
                }
                
                /*
                Can see from prints below that the SSD is not enough to 
                distinguish the best solution each time.
                
                TODO: add gradient comparisons or implement one of the
                feature descriptors in the literature
                */
                    
                for (int i = 0; i < vec.getNumberOfItems(); ++i) {
                    FeatureComparisonStat fs = vec.getArray()[i];                    
                    String str2 = String.format("%d)   pt(%d) %d %s", ds, p1Count, i, fs.toString());
                    if (i == bestIdx) {
                        str2 =    String.format("%d)***pt(%d) %d %s", ds, p1Count, i, fs.toString());
                    }
                    log.info(str2);
                }
                
                assertTrue(foundBest);
                
                p1Count++;
            }
        
            /*
            analyze total map:
                find for each primary key, the smallest diffSqSum and the smallest diffSqSum/err

                -- For the top <?> of those, 
                   -- are they uniquely matched?
                   -- calculate euclidean transformations from pairs of them
                      -- are the results consistent w/ each other and the found frame rotations?
                         -- if yes, then filter map for rotations near that rotation
                            and remove the points already paired.
                            -- the filter the map to remove any stats matches that
                               are not feasible with the transformation parameters.

                         -- else if not consistent?
            */
        }
    }
    
    protected static void getBrownAndLoweFeatureCenters(PairIntArray out1,
        PairIntArray out2) {
        out1.add(128, 87);
        out2.add(32, 82);
        out1.add(150, 68);
        out2.add(56, 66);
        out1.add(148, 79);
        out2.add(53, 76);
        out1.add(165, 85);
        out2.add(66, 83);
        out1.add(144, 67);
        out2.add(52, 64);
    }

    private void getBrownAndLoweFeatureCenters(List<PairInt> points1,
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
    
    private void getVenturiFeatureCenters(PairIntArray points1, 
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

    private void getVenturiFeatureCenters(List<PairInt> points1, 
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

    private void getBooksFeatureCenters(PairIntArray points1, 
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

    private void getBooksFeatureCenters(List<PairInt> points1, 
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
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }
    
}
