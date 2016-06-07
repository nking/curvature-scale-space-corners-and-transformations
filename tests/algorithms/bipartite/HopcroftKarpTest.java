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
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class HopcroftKarpTest extends TestCase {
    
    Logger log = Logger.getLogger(this.getClass().getName());
    
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
        
        log.info(Arrays.toString(m));
        
        MinCostUnbalancedAssignment bipartite2 = 
            new MinCostUnbalancedAssignment();
        
        int s = expectedM.length;
        Graph g2 = getTestGraph0();
        Map<Integer, Integer> m2 = new HashMap<Integer, Integer>();
        m2 = bipartite2.hopcroftKarp(g2, s);
        for (Entry<Integer, Integer> entry : m2.entrySet()) {
            log.info(entry.getKey() + ":" +
                entry.getValue());
        }
        log.info("bfs/dfs hk-> " + Arrays.toString(m));
    
        /*
        0 1
        1 4  <-- or 1,0
        2 2
        3 0  <-- or 3,4
        4 3
        */
        
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
        
        return g;
    }
    
    private Graph getTestGraph0() {
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        weights.put(new PairInt(0, 0), Integer.valueOf(2));
        weights.put(new PairInt(0, 1), Integer.valueOf(3));
        
        weights.put(new PairInt(1, 0), Integer.valueOf(2));
        weights.put(new PairInt(1, 4), Integer.valueOf(1));
        
        weights.put(new PairInt(2, 2), Integer.valueOf(2));
        weights.put(new PairInt(2, 3), Integer.valueOf(2));
        
        weights.put(new PairInt(3, 0), Integer.valueOf(1));
        weights.put(new PairInt(3, 4), Integer.valueOf(3));
    
        weights.put(new PairInt(4, 0), Integer.valueOf(1));
        weights.put(new PairInt(4, 3), Integer.valueOf(2));
                
        /*
        0 1
        1 4  <-- or 1,0
        2 2
        3 0  <-- or 3,4
        4 3
        */
        
        Graph g = new Graph(6, 5, weights, true);

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
