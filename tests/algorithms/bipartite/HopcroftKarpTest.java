package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;

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
        TIntIntMap m2 = new TIntIntHashMap();
        m2 = bipartite2.hopcroftKarp(g2, s);
        
        TIntIntIterator iter = m2.iterator();
        for (int i = m2.size(); i-- > 0;) {
            iter.advance();
            log.info(iter.key() + ":" + iter.value());
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
        
        TIntObjectMap<TIntSet> adjMap = g.getAdjacencyMap();
        for (int i = 0; i < n1; ++i) {
            adjMap.put(i, new TIntHashSet());
        }
        adjMap.get(0).add(0);
        adjMap.get(0).add(1);
        
        adjMap.get(1).add(0);
        adjMap.get(1).add(4);
        
        adjMap.get(2).add(2);
        adjMap.get(2).add(3);
        
        adjMap.get(3).add(0);
        adjMap.get(3).add(4);
        
        adjMap.get(4).add(0);
        adjMap.get(4).add(3);
        
        return g;
    }
    
    private Graph getTestGraph0() {
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
        
        weights.put(new PairInt(0, 0), 2);
        weights.put(new PairInt(0, 1), 3);
        
        weights.put(new PairInt(1, 0), 2);
        weights.put(new PairInt(1, 4), 1);
        
        weights.put(new PairInt(2, 2), 2);
        weights.put(new PairInt(2, 3), 2);
        
        weights.put(new PairInt(3, 0), 1);
        weights.put(new PairInt(3, 4), 3);
    
        weights.put(new PairInt(4, 0), 1);
        weights.put(new PairInt(4, 3), 2);
                
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
