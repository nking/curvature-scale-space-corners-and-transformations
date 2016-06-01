package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
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
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        m.put(Integer.valueOf(3), Integer.valueOf(0));
        m.put(Integer.valueOf(1), Integer.valueOf(2));
        m.put(Integer.valueOf(4), Integer.valueOf(3));
        
        ResidualDigraph rM = new ResidualDigraph(g, m);
        
        assertEquals(6, rM.getForwardLinksRM().size());
        assertEquals(3, rM.getBackwardLinksRM().size());
       
        for (int i = 0; i < 6; ++i) {
            Integer x = Integer.valueOf(i);
            assertTrue(rM.getLeftRM().contains(x));
            
        }
        for (int i = 0; i < 5; ++i) {
            Integer y = Integer.valueOf(i);
            assertTrue(rM.getRightRM().contains(y));        
        }
        
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(0)) 
            .contains(Integer.valueOf(0)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(0)) 
            .contains(Integer.valueOf(2)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(1)) 
            .contains(Integer.valueOf(1)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(2)) 
            .contains(Integer.valueOf(3)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(3)) 
            .contains(Integer.valueOf(1)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(3)) 
            .contains(Integer.valueOf(4)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(4)) 
            .contains(Integer.valueOf(4)));
        assertTrue(rM.getForwardLinksRM().get(Integer.valueOf(5)) 
            .contains(Integer.valueOf(3)));
        
        assertTrue(rM.getBackwardLinksRM().get(Integer.valueOf(0))
            .equals(Integer.valueOf(3)));
        assertTrue(rM.getBackwardLinksRM().get(Integer.valueOf(2))
            .equals(Integer.valueOf(1)));
        assertTrue(rM.getBackwardLinksRM().get(Integer.valueOf(3))
            .equals(Integer.valueOf(4)));
               
        
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
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
                
        weights.put(new PairInt(0, 0), Integer.valueOf(1));
        weights.put(new PairInt(0, 2), Integer.valueOf(1));
        
        weights.put(new PairInt(1, 1), Integer.valueOf(2));
        weights.put(new PairInt(1, 2), Integer.valueOf(1));
        
        weights.put(new PairInt(2, 3), Integer.valueOf(2));
        
        weights.put(new PairInt(3, 0), Integer.valueOf(1));
        weights.put(new PairInt(3, 1), Integer.valueOf(2));
        weights.put(new PairInt(3, 4), Integer.valueOf(3));
    
        weights.put(new PairInt(4, 3), Integer.valueOf(2));
        weights.put(new PairInt(4, 4), Integer.valueOf(1));
        
        weights.put(new PairInt(5, 3), Integer.valueOf(2));
        
        Graph g = new Graph(6, 5, weights, createSourceAndSinkEdges);

        return g;
    }
}
