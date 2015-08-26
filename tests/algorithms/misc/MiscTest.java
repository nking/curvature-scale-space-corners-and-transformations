package algorithms.misc;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class MiscTest extends TestCase {
    
    public MiscTest() {
    }

    public void testPersistToFile() throws Exception {
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < 10; ++i) {
            points.add(new PairInt(i, i));
        }
        
        String fileName = "blob_100_150.dat";
        
        String outPath = Misc.persistToFile(fileName, points);
        
        assertTrue((new File(outPath)).exists());
        
        String filePath = ResourceFinder.findDirectory("bin") + "/" + fileName;
        
        Set<PairInt> pointsR = Misc.deserializeSetPairInt(filePath);

        assertTrue(pointsR.size() == points.size());
        
        for (PairInt p : pointsR) {
            assertTrue(points.remove(p));
        }
        
        assertTrue(points.isEmpty());
    }    
    
}
