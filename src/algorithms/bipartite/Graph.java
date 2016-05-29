package algorithms.bipartite;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * class to hold a bipartite weighted graph.  Note that
 * source and sink nodes are created and edges from 
 * the source to left vertexes are created and assigned
 * a cost of zero
 * and edges from the right vertexes to the sink node
 * are created and assigned a cost of zero.
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
     * <pre>
     * |-
     * </pre>
     */
    private final int sourceNode;
    
    /**
     * <pre>
     * -|
     * </pre>
     */
    private int sinkNode;
    
    /**
     * map of edge weights with key = pairint of left index and right index and
     * value being the edge weight.
     */
    private Map<PairInt, Integer> edgeWeights
        = new HashMap<PairInt, Integer>();

    private Map<PairInt, Integer> sourceEdgeWeights = null;
    private Map<PairInt, Integer> sinkEdgeWeights = null;

    public Graph(int nLeftVertices, int nRightVertices,
        Map<PairInt, Integer> theEdgeWeights, 
        boolean createSourceAndSinkEdges) {
        
        this.nLeft = nLeftVertices;
        this.nRight = nRightVertices;
        
        edgeWeights.putAll(theEdgeWeights);
        
        /*
        For each vertex x in X, there is a left-dummy arc |- -> x, 
           directed from the source node |- to the node x. 
           The per-unit cost of a left-dummy arc is zero: 
           c(|-, x) := 0. 
        For each vertex y in Y , there is a right-dummy arc y -> -|, 
           directed from the node 
           y to the sink node -| and of cost zero: 
           c(y, -|) := 0.            
        */
        if (createSourceAndSinkEdges) {
            this.sourceNode = nLeft;
            this.sinkNode = nRight;
            sourceEdgeWeights = new HashMap<PairInt, Integer>();
            sinkEdgeWeights = new HashMap<PairInt, Integer>();
            for (int i = 0; i < nLeft; ++i) {
                PairInt p = new PairInt(sourceNode, i);
                sourceEdgeWeights.put(p, 0);
            }
            for (int i = 0; i < nLeft; ++i) {
                PairInt p = new PairInt(i, sinkNode);
                 sinkEdgeWeights.put(p, 0);
            }
        } else {
            this.sourceNode = -1;
            this.sinkNode = -1;
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
     * @return the edgeWeights
     */
    public Map<PairInt, Integer> getEdgeWeights() {
        return edgeWeights;
    }

    /**
     * @return the sourceEdgeWeights
     */
    public Map<PairInt, Integer> getSourceEdgeWeights() {
        return sourceEdgeWeights;
    }

    /**
     * @return the sinkEdgeWeights
     */
    public Map<PairInt, Integer> getSinkEdgeWeights() {
        return sinkEdgeWeights;
    }

    /**
     * @return the sourceNode
     */
    public int getSourceNode() {
        return sourceNode;
    }

    /**
     * @return the sinkNode
     */
    public int getSinkNode() {
        return sinkNode;
    }
    
    Graph copyToCreateSourceSink() {
    
        if (sourceNode == -1) {
            Graph g = new Graph(nLeft, nRight, edgeWeights,
                true);
            return g;
        }
    
        throw new IllegalStateException("currently, the "
            + "method is for use on graphs that were"
            + " constructed without source and sink nodes");
    }
}
