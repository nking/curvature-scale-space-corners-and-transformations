package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class HopcroftKarpTest extends TestCase {
    
    public HopcroftKarpTest() {
    }

    public void test0() {
        
        GraphWithoutWeights g = getTestGraph();
        
        int[] expectedM = getTestGraphExpectedMatching();
        
        HopcroftKarp hk = new HopcroftKarp();
        
        int[] m = hk.hopcroftKarpV0(g);
                        
        assertEquals(expectedM.length, m.length);
        
        for (int i = 0; i < expectedM.length; ++i) {        
            assertEquals(expectedM[i], m[i]);
        }
        
        System.out.println(Arrays.toString(m));
        
        MinCostUnbalancedAssignment bipartite2 = 
            new MinCostUnbalancedAssignment();
        
        int s = Math.min(g.getNLeft(), g.getNRight());
        
       
    }
    
    private GraphWithoutWeights getTestGraph() {
        
        GraphWithoutWeights g = new GraphWithoutWeights(5, 5);
        
        int n1 = g.getNLeft();
        int n2 = g.getNRight();
        
        Map<Integer, Set<Integer>> adjMap = g.getAdjacencyMap();
        for (int i = 0; i < n1; ++i) {
            Integer index = Integer.valueOf(i);
            adjMap.put(index, new HashSet<Integer>());
        }
        adjMap.get(Integer.valueOf(0))
            .add(Integer.valueOf(0));
        adjMap.get(Integer.valueOf(0))
            .add(Integer.valueOf(1));
        
        adjMap.get(Integer.valueOf(1))
            .add(Integer.valueOf(0));
        adjMap.get(Integer.valueOf(1))
            .add(Integer.valueOf(4));
        
        adjMap.get(Integer.valueOf(2))
            .add(Integer.valueOf(2));
        adjMap.get(Integer.valueOf(2))
            .add(Integer.valueOf(3));
        
        adjMap.get(Integer.valueOf(3))
            .add(Integer.valueOf(0));
        adjMap.get(Integer.valueOf(3))
            .add(Integer.valueOf(4));
        
        adjMap.get(Integer.valueOf(4))
            .add(Integer.valueOf(0));
        adjMap.get(Integer.valueOf(4))
            .add(Integer.valueOf(3));
        
        if (true) {
            return g;
        }
        // connections in opposite direction
        adjMap.get(Integer.valueOf(0))
            .add(Integer.valueOf(0));
        adjMap.get(Integer.valueOf(1))
            .add(Integer.valueOf(0));
        
        adjMap.get(Integer.valueOf(0))
            .add(Integer.valueOf(1));
        adjMap.get(Integer.valueOf(4))
            .add(Integer.valueOf(1));
        
        adjMap.get(Integer.valueOf(2))
            .add(Integer.valueOf(2));
        adjMap.get(Integer.valueOf(3))
            .add(Integer.valueOf(2));
        
        adjMap.get(Integer.valueOf(0))
            .add(Integer.valueOf(3));
        adjMap.get(Integer.valueOf(4))
            .add(Integer.valueOf(3));
        
        adjMap.get(Integer.valueOf(0))
            .add(Integer.valueOf(4));
        adjMap.get(Integer.valueOf(3))
            .add(Integer.valueOf(4));
        
        return g;
    }
    
    private int[] getTestGraphExpectedMatching() {
    
        int[] m = new int[5];
        m[0] = 1;
        m[1] = 4;
        m[2] = 2;
        m[3] = 0;
        m[4] = 3;
        
        return m;
    }
}
