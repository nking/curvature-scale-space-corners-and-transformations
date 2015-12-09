package algorithms.compGeometry;

import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FurthestPairTest extends TestCase {
    
    public FurthestPairTest() {
    }
    
    /*
            @ 3
            |      
    @       |      @
    @--------------@
    3       |      3
            |
            |
            @ 4
    */
    public void test0() throws Exception {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        points.add(new PairInt(3, 0));
        points.add(new PairInt(3, 1));
        points.add(new PairInt(0, 3));
        points.add(new PairInt(-3, 1));
        points.add(new PairInt(-3, 0));
        points.add(new PairInt(0, -4));
               
        FurthestPair furthestP = new FurthestPair();
        
        int maxDistSq = Integer.MIN_VALUE;
        PairInt p0 = null;
        PairInt p1 = null;
        for (PairInt p : points) {
            for (PairInt p2 : points) {
                if (p.equals(p2)) {
                    continue;
                }
                int dSq = furthestP.distanceSq(p, p2);
                if (dSq > maxDistSq) {
                    maxDistSq = dSq;
                    p0 = p;
                    p1 = p2;
                }
            }
        }
        
        PairInt[] result = furthestP.find(points);
        
        assertNotNull(result);
        assertTrue(result.length == 2);
        
        PairInt expectedP0 = new PairInt(0, 3);
        PairInt expectedP1 = new PairInt(0, -4);
        
        assertTrue(expectedP0.equals(result[0]) || expectedP0.equals(result[1]));
        assertTrue(expectedP1.equals(result[0]) || expectedP1.equals(result[1]));
        assertFalse(result[0].equals(result[1]));
        
    }
    
    public void test2() throws Exception {
        
    }
}
