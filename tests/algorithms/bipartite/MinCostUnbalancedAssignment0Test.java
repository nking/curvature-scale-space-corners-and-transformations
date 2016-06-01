package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MinCostUnbalancedAssignment0Test extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public MinCostUnbalancedAssignment0Test() {
    }
    
    public void testBuildForest() {
        
        // test graphs on pg 21
        Graph g = getTestGraph();
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        ResidualDigraph rM = new ResidualDigraph(g, m);
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        // ------------------------------------
        int s = Math.min(g.getNLeft(), g.getNRight());
        Map<Integer, Integer> aMatching = 
            bipartite.hopcroftKarp(g, rM, s);
        for (Entry<Integer, Integer> entry : aMatching.entrySet()) {
            log.info("matched left " + entry.getKey() + " to right " +
                entry.getValue());
        }
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
