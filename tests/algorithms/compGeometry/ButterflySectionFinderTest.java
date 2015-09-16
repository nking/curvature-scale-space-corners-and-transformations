package algorithms.compGeometry;

import algorithms.imageProcessing.Image;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ButterflySectionFinderTest extends TestCase {
    
    public void testFindButterflySections() throws Exception {
        
        String fileName = "blob_butterfly_01.dat";
        //String fileName = "blob_not_butterfly_01.dat";

        String filePath = ResourceFinder.findDirectory("testresources") + "/" + fileName;

        PairIntArray closedCurve = Misc.deserializePairIntArray(filePath);
        
        Image img = new Image(512, 512);
        MiscDebug.writeImage(closedCurve, img, 0, "butterfly");

        assertNotNull(closedCurve);
        
        assertTrue(closedCurve.getN() > 0);
        
        ButterflySectionFinder finder = new ButterflySectionFinder();
        
        List<Set<PairInt>> sections = finder.findButterflySections(closedCurve);
        
        assertTrue(sections.size() == 1);
        
        Set<PairInt> bPoints = sections.get(0);
        
        assertTrue(bPoints.size() == 10);
        
        assertTrue(bPoints.contains(new PairInt(374, 263)));
        assertTrue(bPoints.contains(new PairInt(374, 261)));
        assertTrue(bPoints.contains(new PairInt(375, 263)));
        assertTrue(bPoints.contains(new PairInt(375, 262)));
        assertTrue(bPoints.contains(new PairInt(376, 263)));
        assertTrue(bPoints.contains(new PairInt(376, 262)));
        assertTrue(bPoints.contains(new PairInt(377, 263)));
        assertTrue(bPoints.contains(new PairInt(377, 262)));
        assertTrue(bPoints.contains(new PairInt(378, 264)));
        assertTrue(bPoints.contains(new PairInt(378, 262)));
        
    }
    
}
