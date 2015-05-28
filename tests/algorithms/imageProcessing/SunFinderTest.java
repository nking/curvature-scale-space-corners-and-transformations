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
public class SunFinderTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SunFinderTest() {
    }
    
    public void testFindPhotosphereInImage() throws Exception {
        
        String[] fileNames = new String[] {
            "costa_rica.jpg",
            "arizona-sunrise-1342919937GHz.jpg",
            "stinson_beach.jpg",
            //"stonehenge.jpg"
        };
        
        for (String fileName : fileNames) {
            
            log.info("fileName=" + fileName);
            
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
     
            SkylineExtractor.setDebugName(fileName);

            ImageHelperForTests helper = new ImageHelperForTests(img1, true);
            
            SkylineExtractor.setDebugName(fileName);
            
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            
            boolean skyIsDarkGrey = false;
            
            SunFinder sFinder = new SunFinder();
            
            int xOffset = helper.getXOffset();
            int yOffset = helper.getYOffset();
            
            sFinder.findSunPhotosphere(img1, xOffset, yOffset, 
                skyIsDarkGrey);
            
            if (fileName.equals("costa_rica.jpg") ||
                fileName.equals("arizona-sunrise-1342919937GHz.jpg")
                ) {
                
                log.info(fileName + " sunCoeff=" + 
                    Arrays.toString(sFinder.getSunEllipseCoeff())
                    + " sunPoints=" + sFinder.getSunPoints().size()
                    );
                
                assertTrue(sFinder.getSunPoints().size() > 100);
                
            } else if (fileName.equals("stinson_beach.jpg") ||
                fileName.equals("stonehenge.jpg")
                ) {
                
                assertTrue(sFinder.getSunPoints().isEmpty());
            }
            
            //TODO: add test for sFinder.correctSkylineForSun here or in other methods
                        
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
            
                sFinder = new SunFinder();

                xOffset = helper.getXOffset();
                yOffset = helper.getYOffset();

                sFinder.findSunPhotosphere(img1, xOffset, yOffset, 
                    skyIsDarkGrey);

                if (fileName.equals("costa_rica.jpg") ||
                    fileName.equals("arizona-sunrise-1342919937GHz.jpg")
                    ) {

                    log.info(fileName + " sunCoeff=" + 
                        Arrays.toString(sFinder.getSunEllipseCoeff())
                        + " sunPoints=" + sFinder.getSunPoints().size()
                        );

                    assertTrue(sFinder.getSunPoints().size() > 100);

                } else if (fileName.equals("stinson_beach.jpg") ||
                    fileName.equals("stonehenge.jpg")
                    ) {

                    assertTrue(sFinder.getSunPoints().isEmpty());
                }

                //TODO: add test for sFinder.correctSkylineForSun here or in other methods

                img1 = img2;
            }
        }
    }

    public static void main(String[] args) {
        
        try {
            SunFinderTest test = new SunFinderTest();

            test.testFindPhotosphereInImage();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
