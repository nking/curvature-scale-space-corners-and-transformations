package algorithms.imageProcessing;

import algorithms.imageProcessing.ShapeMatcher.FeatureComparisonStat;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.*;
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
        assertTrue((stat.sumSqDiff/stat.img2PointErr) < 1);
        
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
        assertTrue((stat.sumSqDiff/stat.img2PointErr) < 1);
        
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
        assertTrue((stat.sumSqDiff/stat.img2PointErr) < 1);
        
        
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
        
        assertTrue((stat0.sumSqDiff/stat0.img2PointErr) < 1);       
    }
    
    public void testCalculateStat2() throws Exception {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();
        
        // choose regions from Brown & Lowe images to compare
        String fileName1, fileName2;
        ImageExt img1, img2;
        List<PairInt> points1, points2;
        int binFactor;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1438231301005L;
        log.info("SEED3=" + seed);
        sr.setSeed(seed);


        fileName1 = "brown_lowe_2003_image1.jpg";
        fileName2 = "brown_lowe_2003_image2.jpg";
        points1 = new ArrayList<PairInt>();
        points2 = new ArrayList<PairInt>();
        points1.add(new PairInt(127, 87)); points1.add(new PairInt(150, 68));
        points2.add(new PairInt(32, 82)); points2.add(new PairInt(56, 66));
        binFactor = 3;
        /*
        transfomation for images having been binned by factor 3:
        
        params=rotationInRadians=6.1807413 rotationInDegrees=354.13039133201386 scale=0.96686685
            translationX=-81.54607 translationY=-14.233707 originX=0.0 originY=0.0
        */

        
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
            
            float div = s.sumSqDiff/s.img2PointErr;
            if (div <= 1) {
                log.info(i + " (expected sumSqDiff/err > 1) for not same blocks =" 
                    + s.toString());
            }
        }
                
    }
    
    public static void main(String[] args) {
        
        ShapeMatcherTest test = new ShapeMatcherTest();
        
        try {
            test.testCalculateStat2();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }
}
