package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
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
        
        TIntIntMap m = new TIntIntHashMap();
                
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        // ------------------------------------
        int s = Math.min(g.getNLeft(), g.getNRight());
        TIntIntMap aMatching = bipartite.hopcroftKarp(g, s);

        TIntIntIterator iter = aMatching.iterator();
        
        for (int i = aMatching.size(); i-- > 0;) {
            iter.advance();
            log.info("matched left " + 
                iter.key() + " to right " +
                iter.value());
        }
    }
    
    private Graph getTestGraph() {
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
                
        weights.put(new PairInt(0, 0), 2);
        weights.put(new PairInt(0, 2), 3);
        
        weights.put(new PairInt(1, 1), 2);
        weights.put(new PairInt(1, 2), 1);
        
        weights.put(new PairInt(2, 3), 2);
        
        weights.put(new PairInt(3, 0), 1);
        weights.put(new PairInt(3, 1), 2);
        weights.put(new PairInt(3, 4), 3);
    
        weights.put(new PairInt(4, 3), 1);
        weights.put(new PairInt(4, 4), 2);
        
        weights.put(new PairInt(5, 3), 2);
        
        Graph g = new Graph(6, 5, weights, true);

        return g;
    }
}
