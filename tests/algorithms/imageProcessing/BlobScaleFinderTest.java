package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BlobScaleFinderTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public void testFindScale() throws Exception {

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

            BlobScaleFinder scaleFinder = new BlobScaleFinder();
            
            scaleFinder.setToDebug();

            TransformationParameters params = scaleFinder.calculateScale(
                img1Orig, img2Orig);

            assertNotNull(params);

            assertTrue(Math.abs(params.getScale() - 1) < 0.1);
        }

    }

    public void estExtractBoundsOfBlobs() throws Exception {

        int width = 10;
        int height = 10;
        int xCenter = 50;
        int yCenter = 50;
        Set<PairInt> rectangle10 = DataForTests.getRectangle(width, height,
            xCenter, yCenter);

        Set<PairInt> rectangle5 = DataForTests.getRectangle(2*width, height,
            2*xCenter, 2*yCenter);

        List<Set<PairInt>> blobs = new ArrayList<Set<PairInt>>();
        blobs.add(rectangle10);
        blobs.add(rectangle5);

        BlobScaleFinder blobFinder = new BlobScaleFinder();

        List<PairIntArray> blobPerimeters = new ArrayList<PairIntArray>();

        int dimension = (2 * xCenter) + xCenter;
        
        // make a fake image for the edge guide
        GreyscaleImage img = new GreyscaleImage(dimension, dimension);
        Arrays.fill(img.getValues(), 200);
        for (Set<PairInt> blobSets : blobs) {
            for (PairInt p : blobSets) {
                img.setValue(p.getX(), p.getY(), 50);
            }
        }

 //NOTE: update test.  expansion of blob?       
        
        blobFinder.extractBoundsOfBlobs(img, blobs, blobPerimeters, dimension,
            dimension, true);

        assertTrue(blobs.size() == 2);
        assertTrue(blobPerimeters.size() == 2);

        // the extracted bounds are tested in EdgeExtractorForBlobBorderTest

    }

    public void estSumIntensity() throws Exception {

        int width = 10;
        int height = 10;
        int xCenter = 10;
        int yCenter = 10;
        Set<PairInt> rectangle10 = DataForTests.getRectangle(width, height,
            xCenter, yCenter);

        boolean toggle = true;
        GreyscaleImage img = new GreyscaleImage(xCenter*2, yCenter*2);
        for (PairInt p : rectangle10) {
            if (toggle) {
                img.setValue(p.getX(), p.getY(), 100 + 3);
            } else {
                img.setValue(p.getX(), p.getY(), 100 - 3);
            }
            toggle = !toggle;
        }
        // average = 100,  each of 36 will be +3 or -3 so sum should be 0

        BlobScaleFinder blobFinder = new BlobScaleFinder();

        double result = blobFinder.sumIntensity(img, rectangle10);

        int expected = rectangle10.size()*100;

        log.info("expected=" + expected + " result=" + result);

        assertTrue(Math.abs(result - expected) < 0.1);

    }

}
