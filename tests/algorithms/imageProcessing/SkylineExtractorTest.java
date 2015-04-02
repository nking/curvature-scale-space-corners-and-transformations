package algorithms.imageProcessing;

import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import static org.junit.Assert.fail;

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
   
    private void testSkyline() throws Exception {
        
        String[] fileNames = new String[] {
            "brown_lowe_2003_image1.jpg",
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
            "arches_sun_01.jpg", 
            "stlouis_arch.jpg", 
            "contrail.jpg"
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
            Set<PairInt> skyPoints = new HashSet<PairInt>();
            int binFactor = skylineExtractor.determineBinFactorForSkyMask(
                detector.getTheta().getNPixels());
            
            RemovedSets removedSets = 
                skylineExtractor.filterAndExtractSkyFromGradient(
                originalColorImage, detector.getTheta(), 
                detector.getGradientXY(), binFactor, skyPoints);
            
            
            Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);

            ImageIOHelper.addToImage(skyPoints, 0, 0, image1, 255, 0, 0);

            String dirPath = ResourceFinder.findDirectory("bin");
            String outFilePath = dirPath + "/skylinetest_01" + 
                fileName + ".png";

            ImageIOHelper.writeOutputImage(outFilePath, image1);
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

}
