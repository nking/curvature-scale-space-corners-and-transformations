package algorithms.search;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.MiscMath;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
   Uniform Cost Search.
 * given a weighted directed graph with weight function, finds the greedy 
 * but optimal best-first shortest path with smaller number of items in the queue than
 * Dijkstras. 
   Uniform Cost Search is a variant of Best-First-Search.
  
   NOTE: to save more space, could refactor to use bit-vectors for state variables.
  
 * All edge weights must be non-negative.
 * 
 * The runtime complexity is O(V + E) due to use of a YFastTrie as the min-priority heap.
 * Note that if the heap wrapper has to choose a Fibonacci instead due to
 * memory constraints, the runtime complexity is O(V*log_2V + E) instead.

   implemented from
   "Position Paper: Dijkstra’s Algorithm versus Uniform Cost Search or a Case Against Dijkstra’s Algorithm"
    by Felner
    Proceedings, The Fourth International Symposium on Combinatorial Search (SoCS-2011)
    https://www.aaai.org/ocs/index.php/SOCS/SOCS11/paper/viewFile/4017/4357

 * 
 * @author nichole
 */
public class UniformCostSearch {
    
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
    
    protected int[] visited = null;

    protected int[] predecessor = null;   
    
    protected int src = -1;

    protected int dest = -1;
    
    private int sentinel = Integer.MAX_VALUE;
        
    // this is recalculated in constructor
    private int maxValue = sentinel - 1;

    // key is cost of path so far plus the edge weight
    protected MinHeapForRT2012 heap = null;

    // refs to nodes internal to heap for decrease key operations
    protected HeapNode[] nodes = null;
    
    private Logger log = Logger.getLogger(getClass().getSimpleName());
    
    private Level logLevel = Level.FINEST;
    
    /**
     *
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
     * @param destVertex the destination vertex index
     * 
     */
    public UniformCostSearch(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex,
        int destVertex) {
        
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }
        
        /*
        initialize single source (g, s)
          initial node is sourceVertex and has cost 0
          initialize the priority queue and add the source node
          initialize visited array to -1
          note: the cost function is g(n), the sum of the weights of the 
            edges from the source node to node n along the shortest currently 
            known path. g(v) = g(u) + w(u, v).
        loop do
          if priority queue empty, return failure
          node = queue.extractMin
          if (visited[goal] is completed) done
          set visited[node] to explored
          for each v adjacent to node
            g(v) = g(u) + w(u,v)
            if (v is not explored and is not in queue)
              insert v into queue
            else if child is in queue and dist[child] > path_cost
              decrease key of v to path_cost
        */ 

        init(dAG, weights, sourceVertex, destVertex);        
    }
        
    /**
     *  find the single shortest path in dAG with edge weights w starting from s.
     */
    public void find() {
                               
        while (heap.getNumberOfNodes() > 0) {

            HeapNode uNode = heap.extractMin();
            
            int u = ((Integer)uNode.getData()).intValue();

            log.log(logLevel, "u: " + toString(u));
            
            if (visited[dest] == 2) {
                log.log(logLevel, "exit heap.n=" + heap.getNumberOfNodes());
                return;
            }

            visited[u] = 2;
            
            // null the entry in nodes so it isn't used in decrease key
            nodes[u] = null;
            
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
            
                int dUPlusWUV = uWeights.get(v) + dist[u];
                
                log.log(logLevel, "  v: " + toString(v) + " dist=" + dUPlusWUV);
                                
                if (visited[v] == 0) {
                    visited[v] = 1;
                    log.log(logLevel, "  add to min-heap v=" + v);
                    dist[v] = dUPlusWUV;
                    int key = dist[v];
                    HeapNode node = new HeapNode(key);
                    node.setData(Integer.valueOf(v));
                    heap.insert(node);
                    nodes[v] = node;
                    predecessor[v] = u;
                } else {
                    if (dist[v] > dUPlusWUV) {
                        log.log(logLevel, "    decrease key to " + dUPlusWUV);
                        dist[v] = dUPlusWUV;
                        predecessor[v] = u;
                        heap.decreaseKey(nodes[v], dUPlusWUV);
                    }
                }
                
                vNode = vNode.getNext();
            }
        }
    }
    
    private void init(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex,
        int destVertex) {
        
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
        dest = destVertex;
        
        maxValue = calcUpperLimitKeyValue();
        
        sentinel = maxValue + 1;
    
        dist = new int[g.length];
        predecessor = new int[g.length];
        visited = new int[g.length];
        
        Arrays.fill(dist, sentinel);
        Arrays.fill(predecessor, -1);
        Arrays.fill(visited, 0);
        
        dist[src] = 0;
                
        initHeap();
    }
    
    private void initHeap() {
        
        int n = g.length;
                
        int nBits = (int)Math.ceil(Math.log(maxValue/Math.log(2)));
        
        //int maxValue, int approxN, int maxNumberOfBits
        heap = new MinHeapForRT2012(sentinel, n, nBits);

        nodes = new HeapNode[n];

        int key = dist[src];

        HeapNode node = new HeapNode(key);
        node.setData(Integer.valueOf(src));

        heap.insert(node);

        nodes[src] = node;
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
        
        log.log(logLevel, "    dist[]=" + Arrays.toString(dist));
        
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
            
            if (uWeights == null || !uWeights.containsKey(v)) {
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
            if (map == null) {
                continue;
            }
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
    
    private String toString(int u) {
        StringBuffer sb = new StringBuffer();
        sb.append("node=").append(u).append(": visited=").append(visited[u]).
               append(", dist=").append(dist[u])
                .append(", prev=").append(predecessor[u])
                .append(" is in queue=").append(!(nodes[u] == null));
        return sb.toString();
    }
}
