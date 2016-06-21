package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MinCostUnbalancedAssignmentTest extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public MinCostUnbalancedAssignmentTest() {
    }
    
    public void testCreateResidualDigraph() {
        
        // test graphs on pg 21
        Graph g = getTestGraph(false);
        
        TIntIntMap m = new TIntIntHashMap();
        m.put(3, 0);
        m.put(1, 2);
        m.put(4, 3);
        
        ResidualDigraph rM = new ResidualDigraph(g, m);
        
        assertEquals(6, rM.getForwardLinksRM().size());
        assertEquals(3, rM.getBackwardLinksRM().size());
       
        
        assertTrue(rM.getForwardLinksRM().get(0) 
            .contains(0));
        assertTrue(rM.getForwardLinksRM().get(0) 
            .contains(2));
        assertTrue(rM.getForwardLinksRM().get(1) 
            .contains(1));
        assertTrue(rM.getForwardLinksRM().get(2) 
            .contains(3));
        assertTrue(rM.getForwardLinksRM().get(3) 
            .contains(1));
        assertTrue(rM.getForwardLinksRM().get(3) 
            .contains(4));
        assertTrue(rM.getForwardLinksRM().get(4) 
            .contains(4));
        assertTrue(rM.getForwardLinksRM().get(5) 
            .contains(3));
        
        assertEquals(rM.getBackwardLinksRM().get(0), 3);
        assertEquals(rM.getBackwardLinksRM().get(2), 1);
        assertEquals(rM.getBackwardLinksRM().get(3), 4);
        
        // ---------------------------------------------------
        /*
        DoubleLinkedCircularList[] forest = 
            bipartite.buildForest(rM);
                
        int k = 0;
        for (DoubleLinkedCircularList list : forest) {
            if (list == null) {
                continue;
            }
            // traverse w/ getLeft for FIFO order
            HeapNode node = list.getSentinel();
            int i = 0;
            while (i < list.getNumberOfNodes()) {
                node = node.getRight();
                List<PathNode> path = bipartite.extractNodes(node);
                log.info("*k=[" + k + "] i=[" + i + "]");
                for (PathNode node1 : path) {
                    log.info("  node=" + node1.toString());
                }
                ++i;
            }
            ++k;
        }
        */
        
    }

    private Graph getTestGraph(boolean createSourceAndSinkEdges) {
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
         
        weights.put(new PairInt(0, 0), 1);
        weights.put(new PairInt(0, 2), 1);
        
        weights.put(new PairInt(1, 1), 2);
        weights.put(new PairInt(1, 2), 1);
        
        weights.put(new PairInt(2, 3), 2);
        
        weights.put(new PairInt(3, 0), 1);
        weights.put(new PairInt(3, 1), 2);
        weights.put(new PairInt(3, 4), 3);
    
        weights.put(new PairInt(4, 3), 2);
        weights.put(new PairInt(4, 4), 1);
        
        weights.put(new PairInt(5, 3), 2);
        
        Graph g = new Graph(6, 5, weights, createSourceAndSinkEdges);

        return g;
    }
}
