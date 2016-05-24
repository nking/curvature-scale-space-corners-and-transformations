package algorithms.bipartite;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
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
}
