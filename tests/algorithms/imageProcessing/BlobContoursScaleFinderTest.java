package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BlobContoursScaleFinderTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public void testFindScaleAdapMeans() throws Exception {

        boolean rotate = true;
        
        String fileName1, fileName2;

        for (int i = 1; i < 3;/*3;*/ ++i) {
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

            log.info("test for " + fileName1 + ", " + fileName2);
            
            BlobScaleFinderWrapper scaleFinder = new BlobScaleFinderWrapper(
                img1Orig, img2Orig);
            
            scaleFinder.setToDebug();

            TransformationParameters params = scaleFinder.calculateScale();

            assertNotNull(params);
            
            log.info("FINAL PARAMS for " + fileName1 + " " + params.toString());

            assertTrue(Math.abs(params.getScale() - 1) < 0.1);
        }
    }
    
    public void estFindScale() throws Exception {

        String fileName1, fileName2;

        for (int i = 0; i < 1;/*3;*/ ++i) {
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

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img1Orig = ImageIOHelper.readImageExt(filePath1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            ImageExt img2Orig = ImageIOHelper.readImageExt(filePath2);

            // this uses CannyEdge filter to make gradient image:  true is for outdoorMode
            ImageHelperForTests helper1 = new ImageHelperForTests(img1Orig, true);

            // this uses CannyEdge filter to make gradient image:  true is for outdoorMode
            ImageHelperForTests helper2 = new ImageHelperForTests(img2Orig, true);

            if (!fileName1.contains("books")) {
                // perform sky subtraction
                SkylineExtractor skylineExtractor = new SkylineExtractor();
                // sky are the zeros in this:
                img1Orig = skylineExtractor.createSkyMaskedImage(
                    helper1.getTheta(), helper1.getGradientXY(), img1Orig,
                    helper1.getCannyEdgeFilterSettings());

                skylineExtractor = new SkylineExtractor();
                img2Orig = skylineExtractor.createSkyMaskedImage(
                    helper2.getTheta(), helper2.getGradientXY(), img2Orig,
                    helper2.getCannyEdgeFilterSettings());
            }

            BlobScaleFinderWrapper scaleFinder = new BlobScaleFinderWrapper(
                img1Orig, img2Orig);
            
            scaleFinder.setToDebug();

            TransformationParameters params = scaleFinder.calculateScale0();

            assertNotNull(params);

            assertTrue(Math.abs(params.getScale() - 1) < 0.1);
        }

    }
}
