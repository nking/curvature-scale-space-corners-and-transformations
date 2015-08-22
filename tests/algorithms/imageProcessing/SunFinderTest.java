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
            
            int origW = img1.getWidth();
            int origH = img1.getHeight();
     
            for (int nRot = 0; nRot < 4; nRot ++) {
                
                if (nRot > 0) {
                    
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
                    img1 = img2;
                }
                
                ImageHelperForTests helper = new ImageHelperForTests(img1, true);
            
                SkylineExtractor skylineExtractor = new SkylineExtractor();
                SkylineExtractor.setDebugName(fileName);
                
                RemovedSets removedSets;
                PairIntArray outputSkyCentroid;
                Set<PairInt> skyPoints;

                int xOffset = helper.getXOffset();
                int yOffset = helper.getYOffset();
                
                // --- find sky points needed for 2nd method tested in SunFinder ---
                removedSets = skylineExtractor.new RemovedSets();
                outputSkyCentroid = new PairIntArray();
                skyPoints = skylineExtractor.extractSkyStarterPoints(
                    helper.getTheta(), helper.getGradientXY(), img1,
                    helper.getCannyEdgeFilterSettings(), outputSkyCentroid,
                    removedSets);

                GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
                    img1, xOffset, yOffset);

                boolean skyIsDarkGrey = skylineExtractor.skyIsDarkGrey(allSkyColor);

                SunFinder sFinder = new SunFinder();

                sFinder.findSunPhotosphere(img1, xOffset, yOffset,
                    skyIsDarkGrey);
            
                if (fileName.equals("costa_rica.jpg")
                    || fileName.equals("arizona-sunrise-1342919937GHz.jpg")) {

                    log.info(fileName + " sunCoeff="
                        + Arrays.toString(sFinder.getSunEllipseCoeff())
                        + " sunPoints=" + sFinder.getSunPoints().size()
                    );

                    assertTrue(sFinder.getSunPoints().size() > 100);

                } else if (fileName.equals("stinson_beach.jpg")
                    || fileName.equals("stonehenge.jpg")) {

                    assertTrue(sFinder.getSunPoints().isEmpty());
                }
                
                if (!sFinder.getSunPoints().isEmpty()) {
                    //---- test for sFinder.correctSkylineForSun here -----
                    int nSkyPointsBeforeFindClouds = skyPoints.size();
                    skylineExtractor.findClouds(skyPoints, new HashSet<PairInt>(),
                        img1, helper.getTheta());
                    skylineExtractor.addEmbeddedIfSimilarToSky(skyPoints, img1,
                        xOffset, yOffset, removedSets);
                    Set<PairInt> exclude = new HashSet<PairInt>();
                    exclude.addAll(removedSets.getHighContrastRemoved());
                    int nAdded = skylineExtractor.addImageBoundaryEmbeddedSkyIfSimilar(
                        skyPoints, exclude, img1, xOffset, yOffset, removedSets,
                        !sFinder.getSunPoints().isEmpty());

                    sFinder.correctSkylineForSun(skyPoints, img1, xOffset,
                        yOffset, helper.getGradientXY());

                    // spot checks
                    if (fileName.equals("arizona-sunrise-1342919937GHz.jpg")
                        || fileName.equals("costa_rica.jpg")) {

                        int x0 = 357;
                        int y0 = 373;
                        int x1 = 340;
                        int y1 = 372;

                        if (fileName.equals("costa_rica.jpg")) {
                            x0 = 60;
                            y0 = 263;
                            x1 = 66;
                            y1 = 266;
                        }
                        
                        int h = origH;
                        int w = origW;

                        // rotate the coords by 90 degrees for each ii of nRot
                        for (int iii = 0; iii < nRot; iii++) {
                                                        
                            int col2 = h - 1 - y0;
                            int row2 = x0;
                            x0 = col2;
                            y0 = row2;

                            col2 = h - 1 - y1;
                            row2 = x1;
                            x1 = col2;
                            y1 = row2;
                            
                            int tmp = h;
                            h = w;
                            w = tmp;
                        }

                        log.info(String.format(
                            "checking for point (%d,%d) in nRot %d of file %s",
                            x0, y0, nRot, fileName));

                        assertFalse(skyPoints.contains(new PairInt(x0, y0)));

                        log.info(String.format(
                            "checking for point (%d,%d) in nRot %d of file %s",
                            x1, y1, nRot, fileName));

                        assertFalse(skyPoints.contains(new PairInt(x1, y1)));
                    }
                }
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
