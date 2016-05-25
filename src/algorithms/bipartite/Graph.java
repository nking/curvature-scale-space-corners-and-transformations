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
    private final int nLeft;

    /**
     * right (==Y) vertices in graph
     */
    private final int nRight;

    /**
     * map of edge weights with key = pairint of left index and right index and
     * value being the edge weight.
     */
    private Map<PairInt, Integer> edgeWeights
        = new HashMap<PairInt, Integer>();

    public Graph(int nLeftVertices, int nRightVertices) {
        this.nLeft = nLeftVertices;
        this.nRight = nRightVertices;
    }
    
    /**
     * @return the number of left vertices
     */
    public int getNLeft() {
        return nLeft;
    }

    /**
     * @return the number of right vertices
     */
    public int getNRight() {
        return nRight;
    }

    /**
     * @return the edgeWeights
     */
    public Map<PairInt, Integer> getEdgeWeights() {
        return edgeWeights;
    }
}
