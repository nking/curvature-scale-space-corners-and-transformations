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
public class BlobPerimeterHelperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public BlobPerimeterHelperTest() {
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
        
        GreyscaleImage segImg = bph.getSegmentationImage(
            SegmentationType.GREYSCALE_KMPP);
        
        assertNotNull(segImg);
        
        segImg = bph.getSegmentationImage(
            SegmentationType.COLOR_POLARCIEXY);
        
        assertNotNull(segImg);
        
        assertTrue(bph.getSmallestGroupLimit() > 1);
        assertTrue(bph.getLargestGroupLimit() > 1);
        assertTrue(bph.getSmallestGroupLimitBinned() > 1);
        assertTrue(bph.getLargestGroupLimitBinned() > 1);
        
        assertNotNull(bph.getGreyscaleImage());
        assertNotNull(bph.getImage());
        
        List<Set<PairInt>> blobs = bph.getBlobs(SegmentationType.GREYSCALE_KMPP, 
            useBinned);
        
        assertNotNull(blobs);
        
        List<Set<PairInt>> blobs2 = bph.getBlobs(SegmentationType.COLOR_POLARCIEXY, 
            useBinned);
        
        assertNotNull(blobs2);
        
        assertTrue(blobs2.size() >= 10);
   
        List<PairIntArray> perimeterList = bph.getBlobPerimeters(
            SegmentationType.GREYSCALE_KMPP, useBinned);
        
        List<PairIntArray> perimeterList2 = bph.getBlobPerimeters(
            SegmentationType.COLOR_POLARCIEXY, useBinned);
        
        assertEquals(blobs.size(), perimeterList.size());
        
        assertEquals(blobs2.size(), perimeterList2.size());
       
    }

}
