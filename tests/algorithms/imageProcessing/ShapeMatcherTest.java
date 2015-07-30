package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ShapeMatcherTest extends TestCase {
    
    public ShapeMatcherTest() {
    }

    public void testCreateNeighborOffsets() throws Exception {
        
        ShapeMatcher matcher = new ShapeMatcher();
        
        float[][] offsets = matcher.createNeighborOffsets();
        
        assertTrue(offsets.length == 25);
        
        Set<PairInt> offsetsSet = new HashSet<PairInt>();
        
        for (int i = 0; i < offsets.length; ++i) {
            
            assertTrue(offsets[i].length == 2);
            
            float x = offsets[i][0];
            float y = offsets[i][1];
            
            PairInt p = new PairInt(Math.round(x), Math.round(y));
            
            assertTrue(Math.abs(x) <= 2);
            assertTrue(Math.abs(y) <= 2);
            
            offsetsSet.add(p);
        }
        
        assertTrue(offsetsSet.size() == 25);
    }
}
