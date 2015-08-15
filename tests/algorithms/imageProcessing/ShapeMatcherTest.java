package algorithms.imageProcessing;

import static algorithms.imageProcessing.FeatureMatcherTest.getBrownAndLoweFeatureCenters90;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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

    public void testBlobs() throws Exception {
        String fileName1, fileName2;
        
        for (int i = 0; i < 3; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    break;
                }
                default: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
            }
            
            BinSegmentationHelper helper = new BinSegmentationHelper(fileName1, 
                fileName2);
            
            if (fileName1.contains("books")) {
                helper.applySteps2(false);
            } else {
                helper.applySteps2(true);
            }
            
        }
    }
    
    public void testCreateNeighborOffsets() throws Exception {

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

    public void testRotateOffsets() throws Exception {

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

    public static void main(String[] args) {

        ShapeMatcherTest test = new ShapeMatcherTest();

        try {
            //test.testCalculateStat2();
            //test.testFeatureMatching();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
    }

}
