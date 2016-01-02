package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class BlobsAndPerimetersTest extends TestCase {
    
    public BlobsAndPerimetersTest() {
    }

    public void test0() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources("two_circles_color2.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        SegmentedImageHelper imgHelper = new SegmentedImageHelper(img);
        SegmentationType type = SegmentationType.COLOR_POLARCIEXY;
        imgHelper.applySegmentation(type);
        boolean useBinned = false;
        
        boolean filterOutImageBoundaryBlobs = true;
        boolean filterOutZeroValuePixels = false;
        List<Set<PairInt>> blobs = 
            BlobsAndPerimeters.extractBlobsFromSegmentedImage(
            imgHelper, type, useBinned, filterOutImageBoundaryBlobs,
            filterOutZeroValuePixels);
                
        // 1 of the blobs is touching image boundaries so is filtered out
        assertEquals(2, blobs.size());

        List<PairIntArray> perimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
            imgHelper, type, blobs, useBinned, false);
        
        assertEquals(2, perimeters.size());
    
        BlobsAndPerimeters.removeRedundantBlobs(blobs);
        
        assertEquals(2, blobs.size());
   
        int w = 100;
        int h = 100;
        
        for (Set<PairInt> blob : blobs) {
            
            int n0 = blob.size();
            
            BlobsAndPerimeters.growRadius(blob, w, h);
            
            int n1 = blob.size();
            
            assertTrue(n1 > n0);
            
            BlobsAndPerimeters.shrinkRadius(blob);
            
            int n2 = blob.size();
            
            assertTrue(n2 < n1);
        }
        
    }

}
