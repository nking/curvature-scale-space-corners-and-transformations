package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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

            BlobScaleFinder scaleFinder = new BlobScaleFinder();
            
            TransformationParameters params = scaleFinder.calculateScale(
                img1Orig, img2Orig);
            
            assertNotNull(params);
            
            assertTrue(Math.abs(params.getScale() - 1) < 0.1);
        }
                
    }

    public void testExtractBoundsOfBlobs() throws Exception {
                
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
        
        blobFinder.extractBoundsOfBlobs(blobs, blobPerimeters, dimension, 
            dimension, true);
        
        assertTrue(blobs.size() == 2);
        assertTrue(blobPerimeters.size() == 2);
        
        // the extracted bounds are tested in EdgeExtractorForBlobBorderTest
        
        // --------- filterBlobsByFirstList -----
        List<Set<PairInt>> blobs0 = new ArrayList<Set<PairInt>>();
        blobs0.add(rectangle10);
        int tolerance = 10;
        
        List<Set<PairInt>> filtered = blobFinder.filterBlobsByFirstList(
            blobs0, blobs, tolerance);
        
        assertTrue(filtered.size() == 1);
        
        assertTrue(filtered.get(0).equals(rectangle10));
    }
    
    public void testExtractBlobsFromSegmentedImage() throws Exception {
        
        int width = 10;
        int height = 10;
        int xCenter = 10; 
        int yCenter = 10;
        Set<PairInt> rectangle10 = DataForTests.getRectangle(width, height, 
            xCenter, yCenter);
        
        GreyscaleImage img = new GreyscaleImage(xCenter*2, yCenter*2);
        for (PairInt p : rectangle10) {
            img.setValue(p.getX(), p.getY(), 100);
        }
        
        BlobScaleFinder blobFinder = new BlobScaleFinder();
        
        List<Set<PairInt>> outputBlobs = new ArrayList<Set<PairInt>>();
        
        int smallestGroupLimit = 10;
        int largestGroupLimit = 5000;
        
        blobFinder.extractBlobsFromSegmentedImage(2, img, outputBlobs,
            smallestGroupLimit, largestGroupLimit);
        
        assertTrue(outputBlobs.size() == 2);
        
        boolean foundRectangle10 = false;
        for (Set<PairInt> p : outputBlobs) {
            if (p.equals(rectangle10)) {
                foundRectangle10 = true;
                break;
            }
        }
        assertTrue(foundRectangle10);
      
        //----------------
        outputBlobs.clear();
        List<PairIntArray> outputBounds = new ArrayList<PairIntArray>();
        
        blobFinder.extractBlobsFromSegmentedImage(2, img, outputBlobs, 
            outputBounds, smallestGroupLimit, largestGroupLimit);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        for (int x = (xCenter - (width/2)); x < (xCenter + (width/2)); ++x) {
            expected.add(new PairInt(x, 5));
            expected.add(new PairInt(x, 14));
        }
        for (int y = (yCenter - (height/2)); y < (yCenter + (height/2)); ++y) {
            expected.add(new PairInt(5, y));
            expected.add(new PairInt(14, y));
        }
        assertTrue(expected.size() == 36);   
        
        int idx = -1;
        for (int i = 0; i < outputBounds.size(); ++i) {
            for (int j = 0; j < outputBounds.get(i).getN(); ++j) {
                int x = outputBounds.get(i).getX(j);
                int y = outputBounds.get(i).getY(j);
                PairInt p = new PairInt(x, y);
                if (expected.contains(p)) {
                    idx = i;
                    break;
                }
            }
        }
        assertTrue(idx > -1);
        
        PairIntArray pai = outputBounds.get(idx);
        for (int j = 0; j < pai.getN(); ++j) {
            int x = pai.getX(j);
            int y = pai.getY(j);
            PairInt p = new PairInt(x, y);
            assertTrue(expected.remove(p));
        }
        
        assertTrue(expected.isEmpty());
    }
    
    public void testSumIntensity() throws Exception {
        
        boolean doNormalize = false;
        
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
        
        int result = blobFinder.sumIntensity(img, rectangle10, doNormalize);
        
        int expected = doNormalize ? 0 : rectangle10.size()*100;
        
        log.info("expected=" + expected + " result=" + result);
        
        assertTrue(result == expected);
        
    }
    
    public void testFilterBlobsByFeatures() throws Exception {
        
        int width = 10;
        int height = 10;
        int xCenter = 10; 
        int yCenter = 10;
        Set<PairInt> rectangle10 = DataForTests.getRectangle(width, height, 
            xCenter, yCenter);
        
        width = 15;
        height = 15;
        xCenter = 50; 
        yCenter = 50;
        Set<PairInt> rectangle50 = DataForTests.getRectangle(width, height, 
            xCenter, yCenter);
        
        /*
        filling rectangle10 and rectangle50 with values scaling from 100 to 136
        in img1.
        
        Doing the same in img2 but adding 50 to each value.
        */
        
        boolean doNormalize = false;
        
        int add = doNormalize ? 50 : 0;
        
        GreyscaleImage img1 = new GreyscaleImage(xCenter*2, yCenter*2);
        int count = 0;
        boolean toggle = true;
        for (PairInt p : rectangle10) {
            if (toggle) {
                img1.setValue(p.getX(), p.getY(), 100 + count);
            } else {
                img1.setValue(p.getX(), p.getY(), 100 - count);
            }
            toggle = !toggle;
            count++;
        }
        count = 0;
        for (PairInt p : rectangle50) {
            img1.setValue(p.getX(), p.getY(), 100 + count);
            count++;
        }
        
        GreyscaleImage img2 = new GreyscaleImage(xCenter*2, yCenter*2);
        count = 0;
        toggle = true;
        for (PairInt p : rectangle10) {
            if (toggle) {
                img2.setValue(p.getX(), p.getY(), 100 + add + count);
            } else {
                img2.setValue(p.getX(), p.getY(), 100 + add - count);
            }
            toggle = !toggle;
            count++;
        }
        count = 0;
        for (PairInt p : rectangle50) {
            img2.setValue(p.getX(), p.getY(), 100 + add + count);
            count++;
        }
        
        List<Set<PairInt>> blobs1 = new ArrayList<Set<PairInt>>();
        blobs1.add(rectangle10);
        blobs1.add(rectangle50);
        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        blobs2.add(rectangle10);
        blobs2.add(rectangle50);
                
        BlobScaleFinder blobFinder = new BlobScaleFinder();
        Map<Integer, List<Integer>> blobPairs = blobFinder.filterBlobsByFeatures(
            img1, img2, blobs1, blobs2, doNormalize);
        
        assertTrue(blobPairs.size() == 2);
        assertTrue(blobPairs.get(0).size() == 1);
        assertTrue(blobPairs.get(0).get(0).intValue() == 0);
        assertTrue(blobPairs.get(1).size() == 1);
        assertTrue(blobPairs.get(1).get(0).intValue() == 1);
        
    }
    
}
