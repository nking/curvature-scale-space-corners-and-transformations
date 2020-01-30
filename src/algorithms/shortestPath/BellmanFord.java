package algorithms.shortestPath;

import algorithms.misc.MiscMath;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import java.util.Arrays;

/**
 * solves the single source paths problem.
 * 
 * It can handle negative weights, but if there is a negative cycle, it returns false.
 * It solves the problem using dynamic programming (= dividing the problem into
 * sub-problems and combining them into the solution).
 *
 * It is non-greedy solution that uses the relax function to find the minimum
 * edge for a node.  relax is invoked |V| - 1 times.
 * 
 * Runtime complexity is O(|V||E|).
 * 
 * implemented following Cormen et al. "Introduction to Algorithms"
 * 
 * @author nichole
 */
public class BellmanFord {
    
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
    
    private int sentinel = Integer.MAX_VALUE;
    
    // this is recalculated in constructor
    private int maxValue = sentinel - 1;
    
    /**
     * find the single shortest paths in dAG with edge weights w starting from s.
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
     * @return returns false if a negative cycle is present, else returns true 
     * and the results are usable.
     */
    public boolean find(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
        
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }
        
        init(dAG, weights, sourceVertex);
        
        return find();
    }
        
    private boolean find() {
        /*
        init
        for i=1 to |V[G]| - 1
           for each edge (u,v) in E[G]
              relax(u,v,w)
        for each edge (u,v) in E[G]
           if d(v) > d(u) + w(u,v)
              return false
        return true
        */
        for (int i = 0; i < (g.length - 1); ++i) {
            for (int u = 0; u < g.length; ++u) {
                
                TIntIntMap uWeights = w[u];

                if (uWeights == null) {
                    continue;
                }

                SimpleLinkedListNode vNode = g[u];

                while (vNode != null && vNode.getKey() != -1) {

                    int v = vNode.getKey();

                    if (!uWeights.containsKey(v)) {
                        throw new IllegalStateException("no weight found for edge "
                                + u + " to " + v);
                    }
                    int wUV = uWeights.get(v);

                    relax(u, v, wUV);

                    vNode = vNode.getNext();
                }
            }
        }
        
        return checkForNegativeCycle();
    }
    
    private void init(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
        
        g = dAG.clone();
        for (int i = 0; i < g.length; ++i) {
            g[i] = new SimpleLinkedListNode(g[i]);
        }
        w = Arrays.copyOf(weights, weights.length);
        src = sourceVertex;
        
        maxValue = calcUpperLimitKeyValue();
        
        sentinel = maxValue + 1;
    
        dist = new int[g.length];
        predecessor = new int[g.length];
        
        Arrays.fill(dist, sentinel);
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
        int u, v;
        for (int i = 1; i < vertexes.length; ++i) {
            u = vertexes[i - 1];
            v = vertexes[i];
            
            TIntIntMap uWeights = w[u];
            
            if (!uWeights.containsKey(v)) {
                throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
            }
            
            sum += uWeights.get(v);
        }
        return sum;
    }

    private int calcUpperLimitKeyValue() {
        
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < w.length; ++i) {
            TIntIntMap map = w[i];
            TIntIntIterator iter = map.iterator();
            int value;
            for (int j = 0; j < map.size(); ++j) {
                iter.advance();
                value = iter.value();
                if (value > max) {
                    max = value;
                }
            }
        }
        
        // if # of bits in max * g.length is larger than max_integer,
        //   use max_integer, else use it
        int b1 = MiscMath.numberOfBits(max);
        int b2 = MiscMath.numberOfBits(g.length);
        if ((b1 + b2) > 31) {
            max = Integer.MAX_VALUE - 1;
        } else {
            max *= g.length;
        }
        
        return max;
    }
    
    private void relax(int u, int v, int wUV) {
        int dUPlusWUV = wUV;
        if (dist[u] == sentinel) {
            dUPlusWUV = sentinel;
        } else {
            dUPlusWUV += dist[u];
        }

        if (dist[v] > dUPlusWUV) {
            dist[v] = dUPlusWUV;
            predecessor[v] = u;
        }
    }

    private boolean checkForNegativeCycle() {
        
        for (int u = 0; u < g.length; ++u) {

            TIntIntMap uWeights = w[u];

            if (uWeights == null) {
                continue;
            }

            SimpleLinkedListNode vNode = g[u];

            while (vNode != null && vNode.getKey() != -1) {

                int v = vNode.getKey();

                if (!uWeights.containsKey(v)) {
                    throw new IllegalStateException("no weight found for edge "
                            + u + " to " + v);
                }
                int wUV = uWeights.get(v);

                int dUPlusWUV = wUV;
                if (dist[u] == sentinel) {
                    dUPlusWUV = sentinel;
                } else {
                    dUPlusWUV += dist[u];
                }

                if (dist[v] > dUPlusWUV) {
                    return false;
                }

                vNode = vNode.getNext();
            }
        }
        
        return true;
    }
}
