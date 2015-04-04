package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import static org.junit.Assert.fail;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class SkylineExtractorTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SkylineExtractorTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
   
    @Test
    public void testSkyline() throws Exception {
        
        String[] fileNames = new String[] {
            /*"brown_lowe_2003_image1.jpg",
            //"brown_lowe_2003_image1_rot.jpg",
            //"brown_lowe_2003_image2.jpg",
            "venturi_mountain_j6_0001.png",
            //"venturi_mountain_j6_0010.png",
            "seattle.jpg",
            "arches.jpg",
            "stinson_beach.jpg",
            "cloudy_san_jose.jpg",
            "30.jpg",
            "sky_with_rainbow.jpg",
            "sky_with_rainbow2.jpg",
            "stonehenge.jpg",
            "norwegian_mtn_range.jpg",
            "halfdome.jpg",
            "costa_rica.jpg",
            "new-mexico-sunrise_w725_h490.jpg",
            "arizona-sunrise-1342919937GHz.jpg",
            "arches_sun_01.jpg", */
            "stlouis_arch.jpg", 
            //"contrail.jpg"
        };
        
        for (String fileName : fileNames) {
            
            log.info("fileName=" + fileName);
            
            // revisit infl points.  is there a threshold removing points?
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            Image img1 = ImageIOHelper.readImage(filePath1);
            int image1Width = img1.getWidth();
            int image1Height = img1.getHeight();
      
            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);
            detector.useOutdoorModeAndExtractSkyline();
            detector.findCorners();
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            
            Image originalColorImage = ImageIOHelper.readImage(filePath1);
            
            //------------------------
            Set<PairInt> skyPoints = new HashSet<PairInt>();
            int binFactor = skylineExtractor.determineBinFactorForSkyMask(
                detector.getTheta().getNPixels());
            
            RemovedSets removedSets = 
                skylineExtractor.filterAndExtractSkyFromGradient(
                originalColorImage, detector.getTheta(), 
                detector.getGradientXY(), binFactor, skyPoints);

            PerimeterFinder perimeterFinder = new PerimeterFinder();
            int[] skyRowMinMax = new int[2];
            Map<Integer, PairInt> skyRowColRange = perimeterFinder.find(
                skyPoints, skyRowMinMax);
        
            skylineExtractor.rightAndLowerDownSizingSkyPointCorrections(
                skyPoints, binFactor, 
                skyRowColRange, skyRowMinMax, originalColorImage,
                detector.getTheta().getWidth(), detector.getTheta().getHeight(),
                detector.getTheta().getXRelativeOffset(), 
                detector.getTheta().getYRelativeOffset());

debugPlot(skyPoints, originalColorImage, 
detector.getTheta().getXRelativeOffset(), detector.getTheta().getYRelativeOffset(), 
"skylinetest_" + fileName + "_after_downsize_corrections");            
            
            GreyscaleImage mask = detector.getGradientXY().createWithDimensions();
            mask.fill(1);
            for (PairInt p : skyPoints) {
                int x = p.getX();
                int y = p.getY();
                mask.setValue(x, y, 0);
            }
            
            int nSkyPointsBeforeFindClouds = skyPoints.size();
        
            Map<Integer, PixelColors> pixelColorsMap = new HashMap<Integer, PixelColors>();

            Map<PairInt, Set<PixelColors>> skyColorsMap = new HashMap<PairInt, Set<PixelColors>>();

            //TODO: this will be revised when have narrowed down which color spaces
            // to use.
            skylineExtractor.populatePixelColorMaps(
                skyPoints, originalColorImage, 
                mask, pixelColorsMap,
                skyColorsMap);

            skylineExtractor.findClouds(
                skyPoints, new HashSet<PairInt>(), originalColorImage, 
                mask,
                pixelColorsMap, skyColorsMap);

debugPlot(skyPoints, originalColorImage, 
detector.getTheta().getXRelativeOffset(), detector.getTheta().getYRelativeOffset(), 
"skylinetest_" + fileName + "_after_findClouds");

            GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
                originalColorImage, detector.getTheta().getXRelativeOffset(), 
                detector.getTheta().getYRelativeOffset());

            boolean skyIsDarkGrey = skylineExtractor.skyIsDarkGrey(allSkyColor);

            Set<PairInt> sunPoints = new HashSet<PairInt>();
            double[] ellipFitParams = 
                skylineExtractor.findSunConnectedToSkyPoints(skyPoints, 
                removedSets.getReflectedSunRemoved(), originalColorImage, 
                detector.getTheta().getXRelativeOffset(), 
                detector.getTheta().getYRelativeOffset(), skyIsDarkGrey, sunPoints);

            //TODO: adjust this:
            if ((ellipFitParams != null) && 
                (
                ((ellipFitParams[2]/ellipFitParams[3])> 7)
                //TODO: reconsider this rule for sun on edge of image
                || ((ellipFitParams[0] < 0) || (ellipFitParams[1] < 0))
                )
                ) {
                sunPoints.clear();
            }

            // should not see sun and rainbow in same image
            Set<PairInt> rainbowPoints = sunPoints.isEmpty() ?
                skylineExtractor.findRainbowPoints(skyPoints, 
                    removedSets.getReflectedSunRemoved(), 
                    originalColorImage, 
                    detector.getTheta().getXRelativeOffset(), 
                    detector.getTheta().getYRelativeOffset(),
                    skyIsDarkGrey) :
                new HashSet<PairInt>();

            // find remaining embedded points that look like high contrast clouds
            // (uses a color range to try to avoid including objects that aren't sky)
            Set<PairInt> embeddedPoints = skylineExtractor.findEmbeddedNonPoints(
                skyPoints, sunPoints, rainbowPoints);

debugPlot(embeddedPoints, originalColorImage, 
detector.getTheta().getXRelativeOffset(), detector.getTheta().getYRelativeOffset(), 
"skylinetest_" + fileName + "_embedded");
            
            //--------------------------
        }
    }
    
    public static void main(String[] args) {
        
        try {
            SkylineExtractorTest test = new SkylineExtractorTest();

            test.testSkyline();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    static int count=0;
    private void debugPlot(Set<PairInt> skyPoints, Image originalColorImage, 
        int xRelativeOffset, int yRelativeOffset, String label) {
        
        Image clr = originalColorImage.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(skyPoints, xRelativeOffset, 
                yRelativeOffset, clr, 0, 0, 255);

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + label + "_" + count + ".png", clr);

            //ImageDisplayer.displayImage(label, clr);
            
        } catch (IOException e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }

}
