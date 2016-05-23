package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.Graph;
import algorithms.bipartite.MinCostUnbalancedAssignment.ResidualDigraph;
import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
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
        
        Graph g = new Graph();
        for (int i = 0; i < 6; ++i) {
            Integer x = Integer.valueOf(i);
            g.leftG.add(x);
        }
        
        for (int i = 0; i < 5; ++i) {
            Integer y = Integer.valueOf(i);
            g.rightG.add(y);
        }
        
        g.edgeWeights.put(new PairInt(0, 0), Integer.valueOf(1));
        g.edgeWeights.put(new PairInt(0, 2), Integer.valueOf(1));
        
        g.edgeWeights.put(new PairInt(1, 1), Integer.valueOf(2));
        g.edgeWeights.put(new PairInt(1, 2), Integer.valueOf(1));
        
        g.edgeWeights.put(new PairInt(2, 3), Integer.valueOf(2));
        
        g.edgeWeights.put(new PairInt(3, 0), Integer.valueOf(1));
        g.edgeWeights.put(new PairInt(3, 1), Integer.valueOf(2));
        g.edgeWeights.put(new PairInt(3, 4), Integer.valueOf(3));
    
        g.edgeWeights.put(new PairInt(4, 3), Integer.valueOf(2));
        g.edgeWeights.put(new PairInt(4, 4), Integer.valueOf(1));
        
        g.edgeWeights.put(new PairInt(5, 3), Integer.valueOf(2));
         
        ResidualDigraph expectedDigrph = new ResidualDigraph();
        for (int i = 0; i < 6; ++i) {
            Integer x = Integer.valueOf(i);
            expectedDigrph.leftRM.add(x);
            expectedDigrph.forwardLinksRM.put(x, 
                new HashSet<Integer>());
        }
        for (int i = 0; i < 5; ++i) {
            Integer y = Integer.valueOf(i);
            expectedDigrph.rightRM.add(y);
        }
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(0)) 
            .add(Integer.valueOf(0));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(0)) 
            .add(Integer.valueOf(2));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(1)) 
            .add(Integer.valueOf(1));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(2)) 
            .add(Integer.valueOf(3));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(3)) 
            .add(Integer.valueOf(1));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(3)) 
            .add(Integer.valueOf(4));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(4)) 
            .add(Integer.valueOf(4));
        expectedDigrph.forwardLinksRM.get(Integer.valueOf(5)) 
            .add(Integer.valueOf(3));
        
        expectedDigrph.backwardLinksRM.put(Integer.valueOf(0), 
            Integer.valueOf(3));
        expectedDigrph.backwardLinksRM.put(Integer.valueOf(2), 
            Integer.valueOf(1));
        expectedDigrph.backwardLinksRM.put(Integer.valueOf(3), 
            Integer.valueOf(4));
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        m.put(Integer.valueOf(3), Integer.valueOf(0));
        m.put(Integer.valueOf(1), Integer.valueOf(2));
        m.put(Integer.valueOf(4), Integer.valueOf(3));
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        ResidualDigraph rM = bipartite.createResidualGraph(g, m);
        
        assertEquals(expectedDigrph.leftRM.size(),
            rM.leftRM.size());
        
        assertEquals(expectedDigrph.rightRM.size(),
            rM.rightRM.size());
        
        assertEquals(expectedDigrph.forwardLinksRM.size(),
            rM.forwardLinksRM.size());
        
        assertEquals(expectedDigrph.backwardLinksRM.size(),
            rM.backwardLinksRM.size());
        
        for (Integer x : expectedDigrph.leftRM) {
            assertTrue(rM.leftRM.contains(x));
        }
        for (Integer y : expectedDigrph.rightRM) {
            assertTrue(rM.rightRM.contains(y));
        }
        for (Entry<Integer, Set<Integer>> entry : 
            expectedDigrph.forwardLinksRM.entrySet()) {
            
            Integer key = entry.getKey();
            Set<Integer> links = entry.getValue();
            
            assertTrue(rM.forwardLinksRM.containsKey(key));
            
            assertEquals(rM.forwardLinksRM.get(key).size(),
                links.size());
            
            for (Integer y : links) {
                assertTrue(rM.forwardLinksRM.get(key).contains(y));
            }
        }
        
        for (Entry<Integer, Integer> entry : 
            expectedDigrph.backwardLinksRM.entrySet()) {
            
            Integer key = entry.getKey();
            Integer value = entry.getValue();
            
            assertTrue(rM.backwardLinksRM.containsKey(key));
            
            assertEquals(rM.backwardLinksRM.get(key),
                value);
        }
        
        DoubleLinkedCircularList[] forest = 
            bipartite.buildForest(rM);
        
        int k = 0;
        for (DoubleLinkedCircularList list : forest) {
            if (list == null) {
                continue;
            }
            
            HeapNode node = list.getSentinel();
            for (int i = 0; i < list.getNumberOfNodes(); ++i) {
                node = node.getRight();
                String nodeType = 
                    node.getClass().getSimpleName().contains("Left") ?
                    " Left" : " Right";
                log.info("k=[" + k + "] key=" + node.getKey()
                    + nodeType + " index=" + node.getData());
            }
            ++k;
        }
    }

}
