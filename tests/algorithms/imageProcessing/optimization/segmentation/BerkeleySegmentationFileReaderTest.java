package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BerkeleySegmentationFileReaderTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public BerkeleySegmentationFileReaderTest() {
    }

    public void test0() throws IOException, Exception {
        
        boolean enable = false;
        
        if (!enable) {
            return;
        }
       
        String fileRoot = "241004";
        
        String dir = ResourceFinder.findTestResourcesDirectory() +
            "/berkeleySegSubset";
            
        String imgFilePath = dir + "/" + fileRoot + "/" + fileRoot + ".jpg";
        
        String segFilePath = dir + "/" + fileRoot + "/" + "241004_1103.seg";
        
        log.info("fileName=" + fileRoot);
                        
        ImageExt img = ImageIOHelper.readImageExt(imgFilePath);

        MiscDebug.writeImage(img, fileRoot);
        
        BerkeleySegmentationFileReader reader = new BerkeleySegmentationFileReader();
        
        List<Set<PairInt>> sets = reader.readFile(segFilePath);
        
        assertNotNull(sets);
        
        assertEquals(17, sets.size());
        
        ImageIOHelper.addAlternatingColorPointSetsToImage(sets, 0, 0, 0, img);
        
        MiscDebug.writeImage(img, fileRoot + "_seg_");
        
        int[][] xyCenN = reader.readCentroids(segFilePath);
        assertEquals(3, xyCenN.length);
        assertEquals(17, xyCenN[0].length);
        assertEquals(17, xyCenN[1].length);
        assertEquals(17, xyCenN[2].length);
    }
}
