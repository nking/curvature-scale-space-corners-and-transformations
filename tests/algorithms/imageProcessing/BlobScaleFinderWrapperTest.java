package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;
import junit.framework.Test;
import junit.framework.TestSuite;

/**
 *
 * @author nichole
 */
public class BlobScaleFinderWrapperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public BlobScaleFinderWrapperTest() {
    }

    public void test0() throws Exception {

        boolean rotate = true;
        
        String fileName1, fileName2;

        for (int i = 0; i < 4; ++i) {
            //fileName1 = "valve_gaussian.png";
            //fileName2 = "valve_gaussian.png";
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
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
                default: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    rotate = false;
                    break;
                }
            }

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img1Orig = ImageIOHelper.readImageExt(filePath1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            ImageExt img2Orig = ImageIOHelper.readImageExt(filePath2);
            
            if (rotate) {
                TransformationParameters params90 = new TransformationParameters();
                params90.setRotationInDegrees(90);
                params90.setOriginX(0);
                params90.setOriginY(0);
                params90.setTranslationX(0);
                params90.setTranslationY(img1Orig.getWidth() - 1);
                Transformer transformer = new Transformer();
                img1Orig = (ImageExt) transformer.applyTransformation(img1Orig, 
                    params90, img1Orig.getHeight(), img1Orig.getWidth());
                int z = 1;
            }

            log.info("fileName=" + fileName1 + ", " + fileName2);
            
            BlobScaleFinderWrapper scaleFinder = null;
            
            if (i == 3) {
                scaleFinder = new BlobScaleFinderWrapper(img1Orig, img2Orig, true);
            } else {
                scaleFinder = new BlobScaleFinderWrapper(img1Orig, img2Orig);
            }
            
            scaleFinder.setToDebug();

            TransformationParameters params = scaleFinder.calculateScale();

            assertNotNull(params);
            
            log.info("FINAL PARAMS for " + fileName1 + " " + params.toString());

            if (fileName1.contains("brown_lowe_2003_image1")) {
                // one portion of image scale is ~ 0.9
                assertTrue(Math.abs(params.getScale() - 1) < 0.15);
            } else {
                assertTrue(Math.abs(params.getScale() - 1) < 0.12);
            }
            
            scaleFinder = null;
            System.gc();
        }
    }
    
    /**
     * Test suite
     * @return static Test
    */
    public static Test suite(){
        return new TestSuite(BlobScaleFinderWrapperTest.class);
    }

    /**
     * Set up a Junit test runner
     * @param args Not used.
    */
    public static void main(String[] args) {

        junit.textui.TestRunner.run(suite());
    }
}
