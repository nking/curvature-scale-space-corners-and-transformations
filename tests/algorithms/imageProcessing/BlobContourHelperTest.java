package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BlobContourHelperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public BlobContourHelperTest() {
    }

    public void test1() throws Exception {
        
        boolean useBinned = false;
        
        String filePath = ResourceFinder.findFileInTestResources("blox.gif");
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        BlobPerimeterHelper bph = new BlobPerimeterHelper(img, "blox");
        bph.increaseLargestGroupLimit(100000);
        
        assertEquals(1, bph.getBinFactor(true));
        
        bph.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        
        bph.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        
        BlobContourHelper bch = new BlobContourHelper(bph, "blox");
       
        List<List<CurvatureScaleSpaceContour>> contoursList1 = 
            bch.generatePerimeterContours(SegmentationType.GREYSCALE_KMPP, 
                useBinned);
        
        List<List<CurvatureScaleSpaceContour>> contoursList2 = 
            bch.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY, 
                useBinned);
        
        assertNotNull(contoursList1);
        assertNotNull(contoursList2);
        
        assertTrue(contoursList1.size() >= 5);
        assertTrue(contoursList2.size() >= 5);
        
        int nNonZero1 = bch.sumPointsOfInterest(
            SegmentationType.GREYSCALE_KMPP, useBinned);
        int nNonZero2 = bch.sumPointsOfInterest(
            SegmentationType.COLOR_POLARCIEXY, useBinned);
        
        assertTrue(nNonZero1 >= 2*5);
        assertTrue(nNonZero2 >= 2*5);
    }

    public void test2() throws Exception {
        
        boolean useBinned = true;
        
        String filePath = ResourceFinder.findFileInTestResources("valve_gaussian.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        BlobPerimeterHelper bph = new BlobPerimeterHelper(img, "valve");
        bph.increaseLargestGroupLimit(100000);
                
        bph.createBinnedGreyscaleImage(300);
        
        assertEquals(3, bph.getBinFactor(true));
        
        bph.applySegmentation(SegmentationType.GREYSCALE_KMPP, useBinned);
        
        bph.applySegmentation(SegmentationType.COLOR_POLARCIEXY, useBinned);
        
        BlobContourHelper bch = new BlobContourHelper(bph, "valve");
       
        List<List<CurvatureScaleSpaceContour>> contoursList1 = 
            bch.generatePerimeterContours(SegmentationType.GREYSCALE_KMPP, 
                useBinned);
        
        List<List<CurvatureScaleSpaceContour>> contoursList2 = 
            bch.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY, 
                useBinned);
        
        assertNotNull(contoursList1);
        assertNotNull(contoursList2);
        
        assertTrue(contoursList1.size() >= 5);
        assertTrue(contoursList2.size() >= 10);
        
        int nNonZero1 = bch.sumPointsOfInterest(
            SegmentationType.GREYSCALE_KMPP, useBinned);
        int nNonZero2 = bch.sumPointsOfInterest(
            SegmentationType.COLOR_POLARCIEXY, useBinned);
        
        assertTrue(nNonZero1 >= 2*5);
        assertTrue(nNonZero2 >= 2*10);
    }

    public void test3() throws Exception {
        
        boolean useBinned = false;
        
        String filePath = ResourceFinder.findFileInTestResources("blox.gif");
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        BlobPerimeterHelper bph = new BlobPerimeterHelper(img, "blox3_1");
        bph.increaseLargestGroupLimit(100000);
        
        assertEquals(1, bph.getBinFactor(true));
        
        bph.applySegmentation(SegmentationType.ADAPTIVE_MEAN, useBinned);
        
        bph.applySegmentation(SegmentationType.COLOR_POLARCIEXY_LARGE, useBinned);
        
        BlobContourHelper bch = new BlobContourHelper(bph, "blox3_2");
       
        List<List<CurvatureScaleSpaceContour>> contoursList1 = 
            bch.generatePerimeterContours(SegmentationType.ADAPTIVE_MEAN, 
                useBinned);
        
        List<List<CurvatureScaleSpaceContour>> contoursList2 = 
            bch.generatePerimeterContours(SegmentationType.COLOR_POLARCIEXY_LARGE, 
                useBinned);
        
        assertNotNull(contoursList1);
        assertNotNull(contoursList2);
        
        //assertTrue(contoursList1.size() >= 5);
        //assertTrue(contoursList2.size() >= 5);
        
        int nNonZero1 = bch.sumPointsOfInterest(
            SegmentationType.ADAPTIVE_MEAN, useBinned);
        int nNonZero2 = bch.sumPointsOfInterest(
            SegmentationType.COLOR_POLARCIEXY_LARGE, useBinned);
        
        //assertTrue(nNonZero1 >= 2*5);
        //assertTrue(nNonZero2 >= 2*5);
    }

}
