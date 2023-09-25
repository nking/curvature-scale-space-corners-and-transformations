package algorithms.compGeometry;

import algorithms.util.PairIntWithIndex;
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
public class ClosestPairBetweenSetsTest extends TestCase {
    
    public ClosestPairBetweenSetsTest() {
    }

    public void testFindClosestPair() {
        
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        
        for (int i = 0; i < 10; ++i) {
            set1.add(new PairInt(i, i));
        }
        //9,9  12,12
        for (int i = 12; i < 24; ++i) {
            set2.add(new PairInt(i, i));
        }
        
        ClosestPairBetweenSets cp = new ClosestPairBetweenSets();
                
        ClosestPairBetweenSets.ClosestPairInt result = cp.findClosestPair(set1, set2);
        
        PairIntWithIndex p1 = result.point0;
        PairIntWithIndex p2 = result.point1;
        
        if (p1.getPixIndex() == 1) {
            assertEquals(9, p1.getX());
            assertEquals(9, p1.getY());
            
            assertEquals(2, p2.getPixIndex());
            assertEquals(12, p2.getX());
            assertEquals(12, p2.getY());
        } else {
            assertEquals(1, p2.getPixIndex());
            assertEquals(9, p2.getX());
            assertEquals(9, p2.getY());
            
            assertEquals(2, p1.getPixIndex());
            assertEquals(12, p1.getX());
            assertEquals(12, p1.getY());
        }
    }


    public void testBruteForceMinDistance() {
        
        List<PairIntWithIndex> p = new ArrayList<PairIntWithIndex>();
        
        for (int i = 8; i < 10; ++i) {
            p.add(new PairIntWithIndex(i, i, 1));
        }
        //9,9  12,12
        for (int i = 12; i < 14; ++i) {
            p.add(new PairIntWithIndex(i, i, 2));
        }
        
        ClosestPairBetweenSets instance = new ClosestPairBetweenSets();
        ClosestPairBetweenSets.ClosestPairInt result = instance.bruteForceMinDistance(p);
        
        PairIntWithIndex p1 = result.point0;
        PairIntWithIndex p2 = result.point1;
        
        if (p1.getPixIndex() == 1) {
            assertEquals(9, p1.getX());
            assertEquals(9, p1.getY());
            
            assertEquals(2, p2.getPixIndex());
            assertEquals(12, p2.getX());
            assertEquals(12, p2.getY());
        } else {
            assertEquals(1, p2.getPixIndex());
            assertEquals(9, p2.getX());
            assertEquals(9, p2.getY());
            
            assertEquals(2, p1.getPixIndex());
            assertEquals(12, p1.getX());
            assertEquals(12, p1.getY());
        }
    }
    
}
