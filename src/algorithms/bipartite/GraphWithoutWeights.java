package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

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

    private TIntObjectMap<TIntSet> adjacencyMap = new
        TIntObjectHashMap<TIntSet>();

    public GraphWithoutWeights(int nLeftVertices, 
        int nRightVertices) {
        this.nLeft = nLeftVertices;
        this.nRight = nRightVertices;
    }
    
    public GraphWithoutWeights(Graph g) {
        
        this.nLeft = g.getNLeft();
        this.nRight = g.getNRight();
                
        TObjectIntIterator<PairInt> iter = 
            g.getEdgeWeights().iterator();

        for (int i = g.getEdgeWeights().size(); i-- > 0;) {

            iter.advance();

            PairInt p = iter.key();

            int idx1 = p.getX();
            int idx2 = p.getY();

            TIntSet indexes2 = adjacencyMap.get(idx1);
            if (indexes2 == null) {
                indexes2 = new TIntHashSet();
                adjacencyMap.put(idx1, indexes2);
            }
            indexes2.add(idx2);
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
    public TIntObjectMap<TIntSet> getAdjacencyMap() {
        return adjacencyMap;
    }
}
