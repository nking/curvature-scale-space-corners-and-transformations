package algorithms.imageProcessing;

import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.PolynomialFitter;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class RainbowFinderTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public RainbowFinderTest() {
    }
    
    public void testFindRainbowInImage() throws Exception {
        
        String[] fileNames = new String[] {
            "sky_with_rainbow.jpg",
            "sky_with_rainbow_rot45.jpg",
            "sky_with_rainbow2.jpg",
            "brown_lowe_2003_image1.jpg",
        };
        
        for (String fileName : fileNames) {
            
            log.info("fileName=" + fileName);
            
            // revisit infl points.  is there a threshold removing points?
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
     
            SkylineExtractor.setDebugName(fileName);

            ImageHelperForTests helper = new ImageHelperForTests(img1, true);
            
            SkylineExtractor.setDebugName(fileName);
            
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            
            RemovedSets removedSets;
            PairIntArray outputSkyCentroid;
            Set<PairInt> points;
            Set<PairInt> rainbowPoints;
            
            boolean skyIsDarkGrey = false;
            if (fileName.equals("sky_with_rainbow2.jpg")) {
                skyIsDarkGrey = true;
            }
            
            RainbowFinder rFinder = new RainbowFinder();
            
            
            removedSets = skylineExtractor.new RemovedSets();
            outputSkyCentroid = new PairIntArray();
            points = skylineExtractor.extractSkyStarterPoints(
                helper.getTheta(), helper.getGradientXY(), img1, 
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid, 
                removedSets);
        
            rFinder.findRainbowInImage(points, 
                removedSets.getReflectedSunRemoved(), img1, 
                helper.getXOffset(), helper.getYOffset(), 
                helper.getTheta().getWidth(), helper.getTheta().getHeight(), 
                skyIsDarkGrey, removedSets);
            
            log.info(fileName + " rainbowCoeff=" + 
                Arrays.toString(rFinder.getRainbowCoeff())
                + " nPoints=" + rFinder.getRainbowPoints().size()
            );
            
            rainbowPoints = rFinder.getRainbowPoints();
            if (fileName.equals("brown_lowe_2003_image1.jpg")) {
                assertTrue(rainbowPoints.isEmpty());
            } else {
                assertTrue(rainbowPoints.size() > 100);
                assertNotNull(rFinder.getRainbowCoeff());
            }
            
            // spot checks
            if (fileName.equals("sky_with_rainbow2.jpg")) {
                assertTrue(rainbowPoints.contains(new PairInt(70, 237)));
                assertTrue(rainbowPoints.contains(new PairInt(96, 192)));
                assertTrue(rainbowPoints.contains(new PairInt(119, 162)));
                assertTrue(rainbowPoints.contains(new PairInt(520, 168)));
                assertTrue(rainbowPoints.contains(new PairInt(556, 218)));
                assertFalse(rainbowPoints.contains(new PairInt(584, 252)));
            } else if (fileName.equals("sky_with_rainbow.jpg")) {
                assertTrue(rainbowPoints.contains(new PairInt(14, 181)));
                assertTrue(rainbowPoints.contains(new PairInt(85, 164)));
                assertTrue(rainbowPoints.contains(new PairInt(295, 167)));
                assertTrue(rainbowPoints.contains(new PairInt(458, 239)));
                assertTrue(rainbowPoints.contains(new PairInt(556, 326)));
            }
                        
            // ------- rotate images by 90 ---------
            int nRot = 3;
            for (int ii = 0; ii < nRot; ii++) {
                ImageExt img2 = new ImageExt(img1.getHeight(), img1.getWidth());
                for (int col = 0; col < img1.getWidth(); col++) {
                    for (int row = 0; row < img1.getHeight(); row++) {

                        /*
                        | 4 1          1 2 3
                        | 5 2          4 5 6
                        | 6 3
                        */

                        int col2 = img1.getHeight() - 1 - row;
                        int row2 = col;

                        int r = img1.getR(col, row);
                        int g = img1.getG(col, row);
                        int b = img1.getB(col, row);

                        img2.setRGB(col2, row2, r, g, b);
                    }
                }

                helper = new ImageHelperForTests(img2, true);

                SkylineExtractor.setDebugName(fileName);

                skylineExtractor = new SkylineExtractor();

                removedSets = skylineExtractor.new RemovedSets();
                outputSkyCentroid = new PairIntArray();
                points = skylineExtractor.extractSkyStarterPoints(
                    helper.getTheta(), helper.getGradientXY(), img2, 
                    helper.getCannyEdgeFilterSettings(), outputSkyCentroid, 
                    removedSets);

                skyIsDarkGrey = false;
                if (fileName.equals("sky_with_rainbow2.jpg")) {
                    skyIsDarkGrey = true;
                }

                rFinder = new RainbowFinder();

                rFinder.findRainbowInImage(points, 
                    removedSets.getReflectedSunRemoved(), img2, 
                    helper.getXOffset(), helper.getYOffset(), 
                    helper.getTheta().getWidth(), helper.getTheta().getHeight(), 
                    skyIsDarkGrey, removedSets);

                rainbowPoints = rFinder.getRainbowPoints();
                if (fileName.equals("brown_lowe_2003_image1.jpg")) {
                    assertTrue(rainbowPoints.isEmpty());
                } else {
                    assertTrue(rainbowPoints.size() > 100);
                    assertNotNull(rFinder.getRainbowCoeff());
                }

                log.info(fileName + " 90 rotated rainbowCoeff=" + 
                    Arrays.toString(rFinder.getRainbowCoeff()));
                
                img1 = img2;
            }
        }
    }

    public void testGeneratePolynomialPoints() throws Exception {

        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
       
        float[] polynomialCoeff = new float[]{0, 6, -0.5f};

        int noise = 1;
        //GreyscaleImage debug = new GreyscaleImage(20, 20);
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < 100; i++) {
            int x = (int)(sr.nextFloat()*20);
            if (sr.nextBoolean()) {
                x = (sr.nextBoolean()) ? x + noise : x - noise;
            }
            int y = (int)(polynomialCoeff[0] + (x * polynomialCoeff[1]) +
                (x * x * polynomialCoeff[2]));
            if (sr.nextBoolean()) {
                y = (sr.nextBoolean()) ? y + noise : y - noise;
            }
            if ((y >= 0 && y < 20) && (x >= 0 && x < 20)) {
                points.add(new PairInt(x, y)); 
                //debug.setValue(x, y, 255);
            }
        }
        
        //ImageDisplayer.displayImage("parabola", debug);
        
        int n = 10;
        float[] outputX = new float[n];
        float[] outputY = new float[n];
        
        RainbowFinder rFinder = new RainbowFinder();
        rFinder.generatePolynomialPoints(points, polynomialCoeff, outputX, 
            outputY);
        
        System.out.println("outX=" + Arrays.toString(outputX));
        System.out.println("outY=" + Arrays.toString(outputY));
        
        float last = -1;
        for (int i = 0; i < outputX.length; i++) {
            float x = outputX[i];
            assertTrue(x > last);
            last = x;
        }
        
        int minXIdx = MiscMath.findYMinIndex(outputX);
        int maxXIdx = MiscMath.findYMaxIndex(outputX);
        int minYIdx = MiscMath.findYMinIndex(outputY);
        int maxYIdx = MiscMath.findYMaxIndex(outputY);
        
        assertTrue(minXIdx == 0);
        assertTrue(maxXIdx == (n - 1));
        assertTrue((minYIdx == 0) || (minYIdx == (n - 1)));
        assertTrue((maxYIdx == (n/2)) || (maxYIdx == ((n/2) - 1)));
        
        assertTrue(Math.abs(outputY[minYIdx]) < 0.5);
        
        assertTrue(Math.abs(outputY[maxYIdx] - 18) < 0.5);
        
        last = -1;
        for (int i = 0; i <= maxYIdx; i++) {
            float y = outputY[i];
            assertTrue(y > last);
            last = y;
        }
        for (int i = (maxYIdx + 1); i < n; i++) {
            float y = outputY[i];
            assertTrue(y < last);
            last = y;
        }
        
        // ======
        
        float dist = rFinder.maxOfPointMinDistances(points, outputX, outputY);
        
        assertTrue(Math.abs(dist - noise) < 4*noise);
        
        // add a point at a known distance to recover answer too
        points.add(new PairInt(5, 18 + 50));
        
        dist = rFinder.maxOfPointMinDistances(points, outputX, outputY);
        
        assertTrue(Math.abs(dist - 50) < 4*noise);
        
    }
   
    @Test
    public void testPopulatePolygon2() throws Exception {
        
        /*
        if coeff[2] is < 0, the parabola peak is upward (higher y)
        
        the magnitude of the points is in the powers 1/10 of coeff[2]
        */
        
        /*
        PolynomialFitter fitter = new PolynomialFitter();
        float[] c = fitter.solve(new float[]{40, 50, 60, 70, 80}, 
            new float[]{160,175,180,175,160});
        */
        
        RainbowFinder rFinder = new RainbowFinder();
        
        float[] polynomialCoeff = new float[]{
            (float)(2.0644931E-5), 6, -0.05f
        };
        
        int n = 5;
        float[] xPoly = new float[n];
        float[] yPoly = new float[n];

        float xMin = 0;
        int offset = 4;
        for (int i = offset; i < (n + offset); i++) {
            
            float x = 10*(xMin + i);
            float y = polynomialCoeff[0] + (x * polynomialCoeff[1]) +
                (x * x * polynomialCoeff[2]);
            
            xPoly[i - offset] = x;
            yPoly[i - offset] = y;
        }
        
        float dist = 4;
        
        float[] outputXPoly = new float[(2 * n) + 1];
        float[] outputYPoly = new float[outputXPoly.length];

        rFinder.populatePolygon(xPoly, yPoly, dist, outputXPoly, outputYPoly, 
            polynomialCoeff, 300, 300);

        System.out.println("outX=" + Arrays.toString(outputXPoly));
        System.out.println("outY=" + Arrays.toString(outputYPoly));
        
        /*
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 
                100, 100, 200);
            plotter.addPlot(outputXPoly, outputYPoly, xPoly, yPoly, "test hull");
            
            String fileName = plotter.writeFile(Integer.valueOf(100100));
                  
        } catch (IOException e) {
            Logger.getLogger(this.getClass().getName()).severe(e.getMessage());
        }
        */
        int minXIdx = MiscMath.findYMinIndex(outputXPoly);
        int maxXIdx = MiscMath.findYMaxIndex(outputXPoly);
        int minYIdx = MiscMath.findYMinIndex(outputYPoly);
        int maxYIdx = MiscMath.findYMaxIndex(outputYPoly);
        
        assertTrue(minXIdx < maxXIdx);
        
        n = outputXPoly.length;
        
        float last = -10;
        for (int i = 0; i <= maxYIdx; i++) {
            float y = outputYPoly[i];
            assertTrue(y > last);
            last = y;
        }
        for (int i = (maxYIdx + 1); i <= maxXIdx; i++) {
            float y = outputYPoly[i];
            assertTrue(y < last);
            last = y;
        }
        last = -1000;
        for (int i = (maxXIdx + 1); i <= (n - maxYIdx - 2); i++) {
            float y = outputYPoly[i];
            assertTrue(y > last);
            last = y;
        }
        
        for (int i = (n - maxYIdx - 1); i < (n - 1); i++) {
            float y = outputYPoly[i];
            assertTrue(y < last);
            last = y;
        }
        
        assertTrue(outputXPoly[0] == outputXPoly[n - 1]);
        assertTrue(outputYPoly[0] == outputYPoly[n - 1]);
        
    }
    
    public static void main(String[] args) {
        
        try {
            RainbowFinderTest test = new RainbowFinderTest();

            test.testFindRainbowInImage();
            test.testPopulatePolygon2();
            test.testGeneratePolynomialPoints();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
