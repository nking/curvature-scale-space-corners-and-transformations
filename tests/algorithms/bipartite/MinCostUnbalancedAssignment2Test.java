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
public class MinCostUnbalancedAssignment2Test extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public MinCostUnbalancedAssignment2Test() {
    }
    
    public void testRefine() {
        
        // test graphs on pg 49
        Graph g = getTestGraph();
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        m.put(Integer.valueOf(1), Integer.valueOf(2));
        m.put(Integer.valueOf(3), Integer.valueOf(0));
        m.put(Integer.valueOf(4), Integer.valueOf(3));
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        FlowNetwork gFlow = new FlowNetwork(g, m);

        assertEquals(1.0f, gFlow.getSourceToLeftFlow(1));
        assertEquals(1.0f, gFlow.getSourceToLeftFlow(3));
        assertEquals(1.0f, gFlow.getSourceToLeftFlow(4));
        
        assertEquals(0.0f, gFlow.getSourceToLeftFlow(0));
        assertEquals(0.0f, gFlow.getSourceToLeftFlow(2));
        assertEquals(0.0f, gFlow.getSourceToLeftFlow(5));
        
        assertEquals(1.0f, gFlow.getRightToSinkFlow(0));
        assertEquals(1.0f, gFlow.getRightToSinkFlow(2));
        assertEquals(1.0f, gFlow.getRightToSinkFlow(3));
        
        assertEquals(0.0f, gFlow.getRightToSinkFlow(1));
        assertEquals(0.0f, gFlow.getRightToSinkFlow(4));
        
        int q = 2;
        int s = m.size();
        double eps = Math.pow(q, 3) * gFlow.getMaxC();
        bipartite.refine(gFlow, s, (float)eps, q);
        
        int z = 1;
    }

    private Graph getTestGraph() {
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
                
        weights.put(new PairInt(0, 0), Integer.valueOf(2));
        weights.put(new PairInt(0, 2), Integer.valueOf(3));
        
        weights.put(new PairInt(1, 1), Integer.valueOf(2));
        weights.put(new PairInt(1, 2), Integer.valueOf(1));
        
        weights.put(new PairInt(2, 3), Integer.valueOf(2));
        
        weights.put(new PairInt(3, 0), Integer.valueOf(1));
        weights.put(new PairInt(3, 1), Integer.valueOf(2));
        weights.put(new PairInt(3, 4), Integer.valueOf(3));
    
        weights.put(new PairInt(4, 3), Integer.valueOf(1));
        weights.put(new PairInt(4, 4), Integer.valueOf(2));
        
        weights.put(new PairInt(5, 3), Integer.valueOf(2));
        
        Graph g = new Graph(6, 5, weights, true);

        return g;
    }
}
