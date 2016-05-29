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
     * nodes that correspond to left (==X) vertices in G
     */
    private Set<Integer> leftRM = new HashSet<Integer>();

    /**
     * nodes that correspond to right (==Y) vertices in G
     */
    private Set<Integer> rightRM = new HashSet<Integer>();

    /**
     * links X to Y (that is, left to right). These are "idle" arcs, f=0, in the
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
        
        for (int i = 0; i < g.getNLeft(); ++i) {
           getLeftRM().add(Integer.valueOf(i));
        }
        
        for (int i = 0; i < g.getNRight(); ++i) {
            getRightRM().add(Integer.valueOf(i));
        }
        
        for (Entry<PairInt, Integer> entry : g.getEdgeWeights().entrySet()) {
            
            PairInt p = entry.getKey();
            
            Integer x = Integer.valueOf(p.getX());
            Integer y = Integer.valueOf(p.getY());
            
            Set<Integer> ys = getForwardLinksRM().get(x);
            if (ys == null) {
                ys = new HashSet<Integer>();
                getForwardLinksRM().put(x, ys);
            }
            
            ys.add(y);
        }
        
        for (Entry<Integer, Integer> entry : m.entrySet()) {
            Integer x = entry.getKey();
            Integer y = entry.getValue();
            
            if (getForwardLinksRM().containsKey(x)) {
                getForwardLinksRM().get(x).remove(y);
            }
            
            getBackwardLinksRM().put(y, x);
        }
    }
            
    public ResidualDigraph(GraphWithoutWeights g) {                
    
        for (int i = 0; i < g.getNLeft(); ++i) {
           getLeftRM().add(Integer.valueOf(i));
        }
        
        for (int i = 0; i < g.getNRight(); ++i) {
            getRightRM().add(Integer.valueOf(i));
        }
        
        for (Entry<Integer, Set<Integer>> entry : 
            g.getAdjacencyMap().entrySet()) {
            
            Integer uIndex = entry.getKey();
            
            Set<Integer> vIndexes = entry.getValue();
            for (Integer v : vIndexes) {            
                Set<Integer> ys = getForwardLinksRM().get(uIndex);
                if (ys == null) {
                    ys = new HashSet<Integer>();
                    getForwardLinksRM().put(uIndex, ys);
                }
                ys.add(v);
            }
        }
    }
 
    /**
     * @return the leftRM
     */
    public Set<Integer> getLeftRM() {
        return leftRM;
    }

    /**
     * @return the rightRM
     */
    public Set<Integer> getRightRM() {
        return rightRM;
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
