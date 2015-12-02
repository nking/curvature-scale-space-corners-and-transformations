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

        boolean rotate = false;
        boolean useBinning = true;

        String fileName1, fileName2;

        for (int i = 0; i < 5; ++i) {
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
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    useBinning = true;
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    useBinning = true;
                    rotate = false;
                    break;
                }
            }

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img1Orig = ImageIOHelper.readImageExt(filePath1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            ImageExt img2Orig = ImageIOHelper.readImageExt(filePath2);
            
            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);
            
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
            
            if (useBinning) {
                scaleFinder = new BlobScaleFinderWrapper(img1Orig, img2Orig, 
                    true, fileName1Root);
            } else {
                scaleFinder = new BlobScaleFinderWrapper(img1Orig, img2Orig, 
                    fileName1Root);
            }
            
            TransformationParameters params = scaleFinder.calculateScale();

            assertNotNull(params);
            
            log.info("FINAL PARAMS for " + fileName1 + " " + params.toString());
            
            assertTrue(Math.abs(params.getScale() - 1) < 0.2);
            
            scaleFinder = null;
            System.gc();
        }
    }
    
    /**
     * Test suite
     * @return static Test
    */
    public static Test suite() {
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
