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
    
    public GraphWithoutWeights(Graph g) {
        
        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
                
        for (Entry<PairInt, Integer> entry : 
            g.getEdgeWeights().entrySet()) {
            
            Integer index1 = entry.getKey().getX();
            Integer index2 = entry.getKey().getY();
            
            Set<Integer> indexes2 = adjacencyMap.get(index1);
            if (indexes2 == null) {
                indexes2 = new HashSet<Integer>();
                adjacencyMap.put(index1, indexes2);
            }
            indexes2.add(index2);
        }
        
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
