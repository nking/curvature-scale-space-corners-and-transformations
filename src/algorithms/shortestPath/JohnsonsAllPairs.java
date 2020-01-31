package algorithms.shortestPath;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;

/**
 *
 * all-pairs shortest path algorithm for sparse graphs.
 * finds shortest paths between every pair of vertices v, v' in the graph.
 * there can be negative edge weights in the input graph, but no negative weight
 * cycles.
 * 
 * runtime complexity is O(|V|^2 + |V||E|)
 * 
 * implemented following Cormen et al. "Introduction to Algorithms".
 * 
 * @author nichole
 */
public class JohnsonsAllPairs {
    
    /**
     * graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     */
    protected SimpleLinkedListNode[] g = null;
    
    /* edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
    */
    protected TIntIntMap[] w = null;
    
    protected int[][] dist = null;

    protected int[][] predecessor = null;   
        
    private int sentinel = Integer.MAX_VALUE;
    
    public JohnsonsAllPairs() {
    }
    
    /*
    compute G' where V[G'] = V[G] union S, and E[G'] = E[G] + (s,v),
       and w(s,v)=0, for all v in V[G].
       (add new node q to graph and set zero weight edges for it to all other nodes)
    If Bellman-Ford(G', w, s) = false,
       exit with message that a negative cycle was found in G.
       (use Bellman-Ford starting from q to get h(v) which will be used as a minimum weight)
    for each vertex v in V[G']
       set h(v) = delta(s,v) computed by Bellman-Ford.
    for each edge (u,v) in E[G']
       do w_hat(u,v) = w(u,v) + h(u) - h(v)
       (reweight the graph so it can handle negative weights)
    for each vertex u in V[G]
       run Dijkstras(G, w_hat, u) to compute delta_hat(u,v) for all v in V[G]
       for each vertex v in V[G]
          do d(u,v) = delta_hat(u, v) + h(v) - h(u)
    return D
    
    where D = d_i_j and d_i_j = delta(i,j)
    */
    
    /**
     * find the shortest paths between all pairs of vertices in an edge-weighted, directed graph. 
     * @param graph directed weighted graph where the main index of the object array
     * is the vertex number, and each vertex has a linked list of outgoing 
     * edges connecting to the vertexes held in each linked list node key.
     * Note that all vertexes, including edge vertexes, must be present as an
     * index of the array graph, i.e. all vertexes must have numerical value 
     * less than dAG.length.
     * @param weights the edge weights in format of outer array index being the
     * first vertex and the map within having a key being the 2nd vertex of the
     * edge with value being the edge weight.
     * @return returns false if a negative cycle is present, else returns true 
     * and the results are usable.
     */
    public boolean find(SimpleLinkedListNode[] graph, TIntIntMap[] weights) {
        
        if (graph == null || graph.length == 0) {
            throw new IllegalArgumentException("graph cannot be null");
        }
        
        init(graph, weights);
        
        return find();
    }
    
    private void init(SimpleLinkedListNode[] graph, TIntIntMap[] weights) {
        
        g = graph.clone();
        for (int i = 0; i < g.length; ++i) {
            if (graph[i] != null) {
                g[i] = new SimpleLinkedListNode(graph[i]);
            }
        }
        w = weights.clone();
        for (int i = 0; i < weights.length; ++i) {
            if (weights[i] != null) {
                w[i] = new TIntIntHashMap(weights[i]);
            }
        }
                
        dist = new int[g.length][g.length];
        predecessor = new int[g.length][g.length];
        
        for (int i = 0; i < g.length; ++i) {
            dist[i] = new int[g.length];
            Arrays.fill(dist[i], sentinel);
            predecessor[i] = new int[g.length];
            Arrays.fill(predecessor[i], -1);
        }
                
        //dist[src] = 0;
    }
    
    private boolean find() {
        
        G g2 = addZeroWeightNode();
        
        BellmanFord bf = new BellmanFord(g2.g, g2.w, g2.newNode);
        boolean noNegativeWgtCycle = bf.find();
        if (!noNegativeWgtCycle) {
            return noNegativeWgtCycle;
        }
        
        // re-weight the graph:
        
        int[] hv = Arrays.copyOf(bf.dist, bf.dist.length);
        
        for (int u = 0; u < g2.g.length; ++u) {
            
            TIntIntMap uWeights = g2.w[u];
           
            if (uWeights == null) {
                continue;
            }
            
            SimpleLinkedListNode vNode = g2.g[u];
            
            while (vNode != null && vNode.getKey() != -1) {
            
                int v = vNode.getKey();
                
                if (!uWeights.containsKey(v)) {
                    throw new IllegalStateException("no weight found for edge " 
                    + u + " to " + v);
                }
                int wUV = uWeights.get(v) + hv[u] - hv[v];
                
                uWeights.put(v, wUV);
                
                vNode = vNode.getNext();
            }
        }
        
        /*
        for each vertex u in V[G]
           run Dijkstras(G, w_hat, u) to compute delta_hat(u,v) for all v in V[G]
           for each vertex v in V[G]
              do d(u,v) = delta_hat(u, v) + h(v) - h(u)
        */
        for (int u = 0; u < g.length; ++u) {
            
            Dijkstras dijkstras = new Dijkstras(g, g2.w, u);
            
            dijkstras.find();
                            
            SimpleLinkedListNode v2Node = g[u];

            while (v2Node != null && v2Node.getKey() != -1) {

                int v2 = v2Node.getKey();

                //do d(u,v) = delta_hat(u, v) + h(v) - h(u)
                dist[u][v2] = dijkstras.dist[v2] + hv[v2] - hv[u];

                System.arraycopy(dijkstras.predecessor, 0, predecessor[u], 
                    0, dijkstras.predecessor.length);

                v2Node = v2Node.getNext();
            }
        }
        
        return true;
    }
    
    private G addZeroWeightNode() {
        
        int nV = g.length;
        
        // add an extra node to the graph and connect it to all vertices with
        //    an outgoing edge of weight 0
        
        SimpleLinkedListNode[] g2 = new SimpleLinkedListNode[nV + 1];
        for (int i = 0; i < g.length; ++i) {
            g2[i] = new SimpleLinkedListNode(g[i]);
        }
        g2[nV] = new SimpleLinkedListNode();
        for (int i = 0; i < g.length; ++i) {
            g2[nV].insert(i);
        }
        
        TIntIntMap[] w2 = new TIntIntMap[nV + 1];
        for (int i = 0; i < w.length; ++i) {
            if (w[i] != null) {
                w2[i] = new TIntIntHashMap();
                TIntIntIterator iter = w[i].iterator();
                for (int j = 0; j < w[i].size(); ++j) {
                    iter.advance();
                    w2[i].put(iter.key(), iter.value());
                }
            }
        }
        w2[nV] = new TIntIntHashMap();
        for (int i = 0; i < w.length; ++i) {
            w2[nV].put(i, 0);
        }
        
        G gPrime = new G();
        gPrime.g = g2;
        gPrime.w = w2;
        gPrime.newNode = nV;
        
        return gPrime;
    }
    
    /**
     * get shortest path from source to destIndex
     * @param srcVertex
     * @param destVertex
     * @return 
     */
    public int[] getShortestPathToVertex(int srcVertex, int destVertex) {
        if (destVertex < 0 || destVertex >= g.length) {
            throw new IllegalArgumentException("destIndex cannot be null");
        }
        if (predecessor == null) {
            throw new IllegalStateException("find must be run before this method can be used");
        }
        
        int[] p = new int[g.length];
        p[p.length - 1] = destVertex;
                
        for (int i = p.length - 2; i > -1; --i) {
            if (destVertex == srcVertex) {
                int len = p.length - 1 - i;
                int[] t = new int[len];
                System.arraycopy(p, i + 1, t, 0, len);
                return t;
            } else if (destVertex == -1) {
                throw new IllegalStateException("path did not complete correctly");
            }
            p[i] = predecessor[srcVertex][destVertex];
            destVertex = p[i];
        }
        
        if (p[0] != srcVertex) {
            throw new IllegalStateException("path did not complete correctly for "
                + " srcVertex=" + srcVertex + " destIndex=" + destVertex);
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

    public static class G {
        public SimpleLinkedListNode[] g;
        public TIntIntMap[] w;
        public int newNode;
    }
  
}
