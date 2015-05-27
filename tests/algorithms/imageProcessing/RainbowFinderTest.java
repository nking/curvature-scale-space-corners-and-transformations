package algorithms.imageProcessing;

import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
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
    
    @Test
    public void testFindRainbowInImage() throws Exception {
        
        String[] fileNames = new String[] {
            "sky_with_rainbow.jpg",
            "sky_with_rainbow2.jpg",
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
            
            RemovedSets removedSets = skylineExtractor.new RemovedSets();
            PairIntArray outputSkyCentroid = new PairIntArray();
            Set<PairInt> points = skylineExtractor.extractSkyStarterPoints(
                helper.getTheta(), helper.getGradientXY(), img1, 
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid, 
                removedSets);
        
            boolean skyIsDarkGrey = false;
            if (fileName.equals("sky_with_rainbow2.jpg")) {
                skyIsDarkGrey = true;
            }
            
            RainbowFinder rFinder = new RainbowFinder();
            
            rFinder.findRainbowInImage(points, 
                removedSets.getReflectedSunRemoved(), img1, 
                helper.getXOffset(), helper.getYOffset(), 
                helper.getTheta().getWidth(), helper.getTheta().getHeight(), 
                skyIsDarkGrey, removedSets);
            
            Set<PairInt> rainbowPoints = rFinder.getRainbowPoints();
            assertTrue(rainbowPoints.size() > 100);
            assertNotNull(rFinder.getRainbowCoeff());
            
            // spot checks
            if (fileName.equals("sky_with_rainbow2.jpg")) {
                assertTrue(rainbowPoints.contains(new PairInt(70, 237)));
                assertTrue(rainbowPoints.contains(new PairInt(96, 192)));
                assertTrue(rainbowPoints.contains(new PairInt(119, 162)));
                assertTrue(rainbowPoints.contains(new PairInt(520, 168)));
                assertTrue(rainbowPoints.contains(new PairInt(556, 218)));
                assertFalse(rainbowPoints.contains(new PairInt(584, 252)));
            } else {
                assertTrue(rainbowPoints.contains(new PairInt(14, 181)));
                assertTrue(rainbowPoints.contains(new PairInt(85, 164)));
                assertTrue(rainbowPoints.contains(new PairInt(295, 167)));
                assertTrue(rainbowPoints.contains(new PairInt(458, 239)));
                assertTrue(rainbowPoints.contains(new PairInt(556, 326)));
            }
        }
    }
   
    public void testPopulatePolygon2() throws Exception {
        
        RainbowFinder rFinder = new RainbowFinder();
        
        /*
        y = x^2 + 3*x
        
        */
        
        float[] polynomialCoeff = new float[] {0.0f, 3.0f, 1.0f};
        
        float[] x = new float[]{-3.5f, -3.f,  -1.5f, 0.f, 0.5f};
        float[] y = new float[]{ 2.f,   0.f, -2.25f, 0.f, 2.0f};
        
        float dist = 1;
        
        float[] outputXPoly = new float[(2 * x.length) + 1];
        float[] outputYPoly = new float[outputXPoly.length];

        rFinder.populatePolygon(x, y, dist, outputXPoly, outputYPoly, 
            polynomialCoeff, 200, 200);

        /*for (int i = 0; i < outputXPoly.length; i++) {
            System.out.println("i=" + i 
                + " (" + outputXPoly[i] + " , " + outputYPoly[i] + ")");
        }*/
        
        int z = 1;
    }
    
    public static void main(String[] args) {
        
        try {
            RainbowFinderTest test = new RainbowFinderTest();

            test.testFindRainbowInImage();
            test.testPopulatePolygon2();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
