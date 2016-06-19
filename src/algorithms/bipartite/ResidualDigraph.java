package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class ResidualDigraph {
    
    /**
     * left (==X) vertices in graph
     */
    private final int nLeft;

    /**
     * right (==Y) vertices in graph
     */
    private final int nRight;

    /**
     * links X to Y (that is, left to right). 
     * These are "idle" arcs, f=0, in the
     * flow network N_G. They correspond to "unmarried" in the matched M graph.
     */
    private Map<Integer, Set<Integer>> forwardLinksRM
        = new HashMap<Integer, Set<Integer>>();

    /**
     * links Y to X (that is, right to left). These are "saturated" arcs, f=1,
     * in the flow network N_G. They correspond to "married" in the matched M
     * graph.
     */
    private Map<Integer, Integer> backwardLinksRM
        = new HashMap<Integer, Integer>();

    public ResidualDigraph(Graph g, Map<Integer, Integer> m) {                
        
        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
            
        for (Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {
            
            PairInt p = entry.getKey();
            
            Integer x = Integer.valueOf(p.getX());
            Integer y = Integer.valueOf(p.getY());
                 
            Set<Integer> ys = forwardLinksRM.get(x);
            if (ys == null) {
                ys = new HashSet<Integer>();
                forwardLinksRM.put(x, ys);
            }
            ys.add(y);
        }
        
        for (Entry<Integer, Integer> entry : m.entrySet()) {
            Integer x = entry.getKey();
            Integer y = entry.getValue();
                        
            backwardLinksRM.put(y, x);
            
            forwardLinksRM.get(x).remove(y);
        }
    }
            
    public ResidualDigraph(GraphWithoutWeights g) {                
    
        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
        
        for (Entry<Integer, Set<Integer>> entry : 
            g.getAdjacencyMap().entrySet()) {
            
            Integer uIndex = entry.getKey();
            
            Set<Integer> vIndexes = entry.getValue();
            for (Integer v : vIndexes) {            
                Set<Integer> ys = forwardLinksRM.get(uIndex);
                if (ys == null) {
                    ys = new HashSet<Integer>();
                    forwardLinksRM.put(uIndex, ys);
                }
                ys.add(v);
            }
        }
    }
 
    public int countOfForwardBipartiteLinks() {
        int n = 0;
        for (Entry<Integer, Set<Integer>> entry : forwardLinksRM.entrySet()) {
           n += entry.getValue().size();
        }
        return n;
    }
    
    /**
     * @return nLeft number of left nodes
     */
    public int getNLeft() {
        return nLeft;
    }

   /**
     * @return number of right nodes
     */
    public int getNRight() {
        return nRight;
    }
    
    /**
     * @return the forwardLinksRM
     */
    public Map<Integer, Set<Integer>> getForwardLinksRM() {
        return forwardLinksRM;
    }

    /**
     * @return the backwardLinksRM
     */
    public Map<Integer, Integer> getBackwardLinksRM() {
        return backwardLinksRM;
    }

    public Map<Integer, Integer> extractMatchings() {

        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        for (Entry<Integer, Integer> entry : backwardLinksRM.entrySet()) {
            Integer rightIndex = entry.getKey();
            Integer leftIndex = entry.getValue();
            m.put(leftIndex, rightIndex);
        }
        
        return m;
    }

}
