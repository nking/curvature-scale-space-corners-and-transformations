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
public class GraphWithoutWeights {
    
    /**
     * left (==X) vertices in graph
     */
    private final int nLeft;

    /**
     * right (==Y) vertices in graph
     */
    private final int nRight;

    private Map<Integer, Set<Integer>> adjacencyMap
        = new HashMap<Integer, Set<Integer>>();

    public GraphWithoutWeights(int nLeftVertices, int nRightVertices) {
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
     * @return the adjacencyMap
     */
    public Map<Integer, Set<Integer>> getAdjacencyMap() {
        return adjacencyMap;
    }
}
