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
        
        Graph g = new Graph(6, 5);
        
        g.getEdgeWeights().put(new PairInt(0, 0), Integer.valueOf(1));
        g.getEdgeWeights().put(new PairInt(0, 2), Integer.valueOf(1));
        
        g.getEdgeWeights().put(new PairInt(1, 1), Integer.valueOf(2));
        g.getEdgeWeights().put(new PairInt(1, 2), Integer.valueOf(1));
        
        g.getEdgeWeights().put(new PairInt(2, 3), Integer.valueOf(2));
        
        g.getEdgeWeights().put(new PairInt(3, 0), Integer.valueOf(1));
        g.getEdgeWeights().put(new PairInt(3, 1), Integer.valueOf(2));
        g.getEdgeWeights().put(new PairInt(3, 4), Integer.valueOf(3));
    
        g.getEdgeWeights().put(new PairInt(4, 3), Integer.valueOf(2));
        g.getEdgeWeights().put(new PairInt(4, 4), Integer.valueOf(1));
        
        g.getEdgeWeights().put(new PairInt(5, 3), Integer.valueOf(2));
         
        ResidualDigraph expectedDigrph = new ResidualDigraph();
        for (int i = 0; i < 6; ++i) {
            Integer x = Integer.valueOf(i);
            expectedDigrph.getLeftRM().add(x);
            expectedDigrph.getForwardLinksRM().put(x, 
                new HashSet<Integer>());
        }
        for (int i = 0; i < 5; ++i) {
            Integer y = Integer.valueOf(i);
            expectedDigrph.getRightRM().add(y);
        }
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(0)) 
            .add(Integer.valueOf(0));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(0)) 
            .add(Integer.valueOf(2));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(1)) 
            .add(Integer.valueOf(1));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(2)) 
            .add(Integer.valueOf(3));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(3)) 
            .add(Integer.valueOf(1));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(3)) 
            .add(Integer.valueOf(4));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(4)) 
            .add(Integer.valueOf(4));
        expectedDigrph.getForwardLinksRM().get(Integer.valueOf(5)) 
            .add(Integer.valueOf(3));
        
        expectedDigrph.getBackwardLinksRM().put(Integer.valueOf(0), 
            Integer.valueOf(3));
        expectedDigrph.getBackwardLinksRM().put(Integer.valueOf(2), 
            Integer.valueOf(1));
        expectedDigrph.getBackwardLinksRM().put(Integer.valueOf(3), 
            Integer.valueOf(4));
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        m.put(Integer.valueOf(3), Integer.valueOf(0));
        m.put(Integer.valueOf(1), Integer.valueOf(2));
        m.put(Integer.valueOf(4), Integer.valueOf(3));
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        ResidualDigraph rM = bipartite.createResidualGraph(g, m);
        
        assertEquals(expectedDigrph.getLeftRM().size(),
            rM.getLeftRM().size());
        
        assertEquals(expectedDigrph.getRightRM().size(),
            rM.getRightRM().size());
        
        assertEquals(expectedDigrph.getForwardLinksRM().size(),
            rM.getForwardLinksRM().size());
        
        assertEquals(expectedDigrph.getBackwardLinksRM().size(),
            rM.getBackwardLinksRM().size());
        
        for (Integer x : expectedDigrph.getLeftRM()) {
            assertTrue(rM.getLeftRM().contains(x));
        }
        for (Integer y : expectedDigrph.getRightRM()) {
            assertTrue(rM.getRightRM().contains(y));
        }
        for (Entry<Integer, Set<Integer>> entry : 
            expectedDigrph.getForwardLinksRM().entrySet()) {
            
            Integer key = entry.getKey();
            Set<Integer> links = entry.getValue();
            
            assertTrue(rM.getForwardLinksRM().containsKey(key));
            
            assertEquals(rM.getForwardLinksRM().get(key).size(),
                links.size());
            
            for (Integer y : links) {
                assertTrue(rM.getForwardLinksRM().get(key).contains(y));
            }
        }
        
        for (Entry<Integer, Integer> entry : 
            expectedDigrph.getBackwardLinksRM().entrySet()) {
            
            Integer key = entry.getKey();
            Integer value = entry.getValue();
            
            assertTrue(rM.getBackwardLinksRM().containsKey(key));
            
            assertEquals(rM.getBackwardLinksRM().get(key),
                value);
        }
        
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
        // ------------------------------------
        int s = Math.min(g.getNLeft(), g.getNRight());
        Map<Integer, Integer> aMatching = 
        //    bipartite.hopcroftKarp(g);
            bipartite.hopcroftKarp(g, rM, s);
        for (Entry<Integer, Integer> entry : aMatching.entrySet()) {
            log.info("matched left " + entry.getKey() + " to right " +
                entry.getValue());
        }
    }

}
