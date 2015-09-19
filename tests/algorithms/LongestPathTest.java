package algorithms;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class LongestPathTest extends TestCase {
    
    public LongestPathTest() {
    }
    
    public void testFindMaxCostPath() throws Exception {
        
        Set<PairInt> points = getData1();
        
        PairInt src = new PairInt(163, 243);
        Set<PairInt> potentialDest = new HashSet<PairInt>();
        potentialDest.add(new PairInt(164, 244));
        
        LongestPath longestPath = new LongestPath();
        
        Set<PairInt> path = longestPath.findMaxCostPath(points, src, potentialDest);
        
        assertTrue(path.size() == points.size());
        
    }

    private Set<PairInt> getData1() {
        
        Set<PairInt> points = new HashSet<PairInt>();
        points.add(new PairInt(163, 244));
        points.add(new PairInt(163, 243));
        points.add(new PairInt(162, 243));
        points.add(new PairInt(162, 242));
        points.add(new PairInt(162, 241));
        points.add(new PairInt(161, 242));
        points.add(new PairInt(160, 242));
        points.add(new PairInt(161, 240));
        points.add(new PairInt(160, 240));
        points.add(new PairInt(159, 240));
        points.add(new PairInt(159, 241));
        
        return points;
    }
}
