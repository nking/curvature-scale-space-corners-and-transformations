package algorithms.misc;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
        
        Set<PairInt> points2 = new HashSet<PairInt>(points);
        
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
        
        PairIntArray pointsR2 = Misc.deserializePairIntArray(filePath);

        assertTrue(pointsR2.getN() == points2.size());
        
        for (int i = 0; i < pointsR2.getN(); ++i) {
            PairInt p2 = new PairInt(pointsR2.getX(i), pointsR2.getY(i));
            assertTrue(points2.remove(p2));
        }
        
        assertTrue(points2.isEmpty());
    }    
    
    public void testReverse() {
        
        List<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < 4; ++i) {
            Integer index = Integer.valueOf(i);
            list.add(index);
        }
        Misc.<Integer>reverse(list);
        
        for (int i = 0; i < 4; ++i) {
            Integer index = list.get(i);
            assertEquals(3 - i, index.intValue());
        }
    }
}
