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
public class MinCostUnbalancedAssignment3Test extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public MinCostUnbalancedAssignment3Test() {
    }
    
    public void testFlow() {
        
        // test graphs on pg 49
        Graph g = getTestGraph();
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        Map<Integer, Integer> m = bipartite.flowAssign(g);
        
        int z = 1;
    }

    private Graph getTestGraph() {
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        /*
        for best cost, w/o regard to maximizing nMatches:
        1  2
        3  0
        4  3
        
        unmatched:
        0,2,5     1,4
        ------------------
        
        for maximizing number of matches, then min cost:
        results are same as Hopcroft-Karp for this example
        which is a matching size of 5:
        5  3   or  2 3 is equivalent
        4  4
        3  0
        1  1
        0  2
        
        -----
        so far, the full algorithm progresses to this, but can find
        no further augmenting paths:
          0  1
          1  1
          3  4
          4  3  <--- this one is min-cost, optimal for the pair
                     but it prevents 5 3 
        */
                
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
