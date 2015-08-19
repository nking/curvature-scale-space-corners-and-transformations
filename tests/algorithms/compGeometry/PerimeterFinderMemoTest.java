package algorithms.compGeometry;

import algorithms.util.BitVectorRepresentation;
import algorithms.util.PairInt;
import algorithms.util.PathStep;
import algorithms.util.VeryLongBitString;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PerimeterFinderMemoTest extends TestCase {
    
    public PerimeterFinderMemoTest() {
    }

    public void testStoreNoSolutionState() {
        
        Set<PairInt> points = getSet0();
        
        BitVectorRepresentation bvr = new BitVectorRepresentation(points);
        
        PairInt currentXY = points.iterator().next();
        Set<PairInt> availableMoves = new HashSet<PairInt>();
        
        for (int i = (currentXY.getX() - 1); i <= (currentXY.getX() + 1); ++i) {
            for (int j = (currentXY.getY() - 1); j <= (currentXY.getY() + 1); ++j) {
                PairInt p = new PairInt(i, j);
                if (p.equals(currentXY)) {
                    continue;
                }
                availableMoves.add(p);
            }
        }
        Set<PairInt> visited = new HashSet<PairInt>();
        int nVisited = 10;
        for (PairInt p : points) {
            if (p.equals(currentXY) || availableMoves.contains(p)) {
                continue;
            }
            visited.add(p);
            if (visited.size() == nVisited) {
                break;
            }
        }
        
        VeryLongBitString totalVisitedBitVector = bvr.createBitstring(visited);
        
        PathStep step = new PathStep(currentXY, availableMoves, 
            totalVisitedBitVector);
        
        LinkedList<PathStep> path = new LinkedList<PathStep>();
        path.add(step);
        
        PerimeterFinderMemo memo = new PerimeterFinderMemo();
        
        memo.storeNoSolutionState(path);
        
        PathStep stepCopy = step.copy();
        
        assertTrue(memo.isANoSolutionState(stepCopy));
        
        availableMoves.remove(availableMoves.iterator().next());
        PathStep step2 = new PathStep(currentXY, availableMoves, 
            bvr.createBitstring(visited));
        assertFalse(memo.isANoSolutionState(step2));
        
        visited.remove(visited.iterator().next());
        step2 = new PathStep(currentXY, availableMoves, 
            bvr.createBitstring(visited));
        assertFalse(memo.isANoSolutionState(step2));
    }

    private Set<PairInt> getSet0() {
        
        Set<PairInt> points = new HashSet<PairInt>();

        for (int x = 10; x < 50; ++x) {
            for (int y = 10; y < 50; ++y) {
                points.add(new PairInt(x, y));
            }
        }
        
        return points;
    }
    
}
