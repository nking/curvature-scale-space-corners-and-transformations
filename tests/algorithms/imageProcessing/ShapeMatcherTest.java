package algorithms.imageProcessing;

import algorithms.imageProcessing.ShapeMatcher.FeatureComparisonStat;
import algorithms.util.PairInt;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
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
        //seed = 1438229484875L;
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
        //seed = 1438229484897L;
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
        ImageExt img1 = randomImage(0, 255, 10, 10);
        ImageExt img2 = (ImageExt)img1.copyImage();
        applyRandomOffsets(2, img2); 
        
        float[][] offsets0 = matcher.createNeighborOffsets();
        
        FeatureComparisonStat stat = matcher.calculateStat(img1, img2,
            5, 5, 5, 5, offsets0, offsets0);
        assertTrue(stat.avgDiffPix < 1);
        assertTrue(stat.stDevDiffPix < 1);
        double rch = (2*Math.sqrt(10*10)/25.);
        assertTrue(Math.abs(stat.avgDivPix - 1) <= 3*rch);
        assertTrue(stat.stDevDivPix < Math.sqrt(25*(rch*rch)/24.));
        
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
        assertTrue(stat.avgDiffPix < maxDiff);
        assertTrue(stat.stDevDiffPix < Math.sqrt(25*(maxDiff*maxDiff)/24.));
        assertTrue(Math.abs(stat.avgDivPix - (1./0.75)) < 0.1);
        assertTrue(stat.stDevDivPix < 0.15);
        
        
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
        assertTrue(Math.abs(stat.avgDiffPix - 2) < 0.01);
        assertTrue(stat.stDevDiffPix < 0.03);
        assertTrue(Math.abs(stat.avgDivPix - 1) < 0.33);
        
        
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
        
        assertTrue(Math.abs(stat0.avgDiffPix) < 0.01);
        assertTrue(Math.abs(stat0.avgDiffPix - stat.avgDiffPix) < 0.01);
        assertTrue(Math.abs(stat0.stDevDiffPix) < 0.01);
        assertTrue(Math.abs(stat0.stDevDiffPix - stat.stDevDiffPix) < 0.01);
        assertTrue(Math.abs(stat0.avgDivPix - 1) < 0.01);
        assertTrue(Math.abs(stat0.avgDivPix - stat.avgDivPix) < 0.01);
        assertTrue(Math.abs(stat0.stDevDivPix) < 0.01);
        assertTrue(Math.abs(stat0.stDevDivPix - stat.stDevDivPix) < 0.01);        
    }
}
