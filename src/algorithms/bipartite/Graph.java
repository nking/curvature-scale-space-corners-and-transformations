package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class Graph {
    
    /**
     * left (==X) vertices in graph
     */
    private Set<Integer> leftG = new HashSet<Integer>();

    /**
     * right (==Y) vertices in graph
     */
    private Set<Integer> rightG = new HashSet<Integer>();

    /**
     * map of edge weights with key = pairint of left index and right index and
     * value being the edge weight.
     */
    private Map<PairInt, Integer> edgeWeights
        = new HashMap<PairInt, Integer>();

    /**
     * @return the leftG
     */
    public Set<Integer> getLeftG() {
        return leftG;
    }

    /**
     * @return the rightG
     */
    public Set<Integer> getRightG() {
        return rightG;
    }

    /**
     * @return the edgeWeights
     */
    public Map<PairInt, Integer> getEdgeWeights() {
        return edgeWeights;
    }
}
