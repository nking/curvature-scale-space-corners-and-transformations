package algorithms.shortestPath;

import algorithms.graphs.TopologicalSort;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * given a weighted directed graph with weight function, the shortest path from
 * u to v is the path whose sum of weights is the minimum.
 * 
 * The DAG shortest path accepts a DAG as input, so does not have negative
 * weight cycles.
 * 
 * The runtime complexity is O(V + E).
 * That's in contrast to the Bellman-Ford single source shortest path 
 * algorithm which is O(V*E), but can handle negative weight cycles by returning 
 * false at the end of the algorithm.
 * 
 * implemented from pseudocode from Cormen et al. "Introduction to Algorithms".
 * 
 * from Cormen et al. :
 * useful for determining critical paths in PERT chart analysis.
 * edges hold jobs to be perform and their weights are the estimated time
 * to complete the jobs.
 * If edge (u,v) enters vertex v, and edge (v,x) leaves v, then job (u,v) must
 * be performed prior to (v,x).
 *   *A critical path is the longest time to perform an ordered sequence of jobs
 * and that is the longest path through the DAG.   The weight on a critical
 * path is the minimum time to perform all jobs.
 *   Can calculate the critical path by either:
 *     - negating the edge weights and run this code on that, or
 *     - or modify a version of this code to use -infinity in the initialization
 *       method and change '>' to '<' in the relax method.
 * 
 * @author nichole
 */
public class DAGShortestPath {
    
    /**
     * directed acyclic graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     */
    protected SimpleLinkedListNode[] g = null;
    
    /* edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
    */
    protected TIntIntMap[] w = null;
    
    protected int[] dist = null;

    protected int[] predecessor = null;   
    
    protected int src = -1;
    
    public DAGShortestPath() {
    }
    
    /**
     * find the single shortest path in dAG with edge weights w starting from s.
     * @param dAG directed acyclic graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     * Note that all vertexes, including edge vertexes, must be present as an
     * index of the array dAG, i.e. all vertexes must have numerical value 
     * less than dAG.length.
     * @param weights the edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
     * @param sourceVertex the source vertex index
     */
    public void find(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
        
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }
        
        /*
        topologically sort the vertices of G
        initialize single source (g, s)
        for each vertex u, taken in topologically sorted order
            do for each vertex v in adj[u]
               do relax(u, v, w)
        */ 
        
        init(dAG, weights, sourceVertex);
        
        find();
    }
        
    private void find() {
       
        TopologicalSort ts = new TopologicalSort(g);
        
        int[] sortedVertexes = ts.sort();
                
        for (int u : sortedVertexes) {
            
            TIntIntMap uWeights = w[u];
            
            if (uWeights == null) {
                continue;
            }
            
            SimpleLinkedListNode next = g[u];
            
            while (next != null && next.getKey() != -1) {
            
                int v = next.getKey();
                
                if (!uWeights.containsKey(v)) {
                    throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
                }
                int wUV = uWeights.get(v);                
                
                int dUPlusWUV = wUV;
                if (dist[u] == Integer.MAX_VALUE) {
                    dUPlusWUV = Integer.MAX_VALUE;
                } else {
                    dUPlusWUV += dist[u];
                }
                                
                if (dist[v] > dUPlusWUV) {
                    dist[v] = dUPlusWUV;
                    predecessor[v] = u;
                }
                
                next = next.getNext();
            }
        }
    }
    
    private void init(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
        
        g = dAG.clone();
        for (int i = 0; i < dAG.length; ++i) {
            g[i] = new SimpleLinkedListNode(dAG[i]);
        }
        w = weights.clone();
        for (int i = 0; i < weights.length; ++i) {
            if (weights[i] != null) {
                w[i] = new TIntIntHashMap(weights[i]);
            }
        }
        src = sourceVertex;
    
        dist = new int[g.length];
        predecessor = new int[g.length];
        
        Arrays.fill(dist, Integer.MAX_VALUE);
        Arrays.fill(predecessor, -1);
        
        dist[src] = 0;
    }
    
    /**
     * get shortest path from source to destIndex
     * @param destVertex
     * @return 
     */
    public int[] getShortestPathToVertex(int destVertex) {
        if (destVertex < 0 || destVertex >= g.length) {
            throw new IllegalArgumentException("destIndex cannot be null");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        
        int[] p = new int[g.length];
        p[p.length - 1] = destVertex;
                
        for (int i = p.length - 2; i > -1; --i) {
            if (destVertex == src) {
                int len = p.length - 1 - i;
                int[] t = new int[len];
                System.arraycopy(p, i + 1, t, 0, len);
                return t;
            } else if (destVertex == -1) {
                throw new IllegalStateException("path did not complete correctly");
            }
            p[i] = predecessor[destVertex];
            destVertex = p[i];
        }
        
        if (p[0] != src) {
            throw new IllegalStateException("path did not complete correctly for destIndex");
        }
        
        return p;
    }
    
    public int getSumOfPath(int[] vertexes) {
        if (vertexes == null) {
            throw new IllegalArgumentException("vertexes cannot be null");
        }
        if (dist == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        int sum = 0;
        for (int idx : vertexes) {
            sum += dist[idx];
        }
        return sum;
    }
    
}
