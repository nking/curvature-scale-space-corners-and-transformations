package algorithms.shortestPath;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.MiscMath;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;

/**
 * given a weighted directed graph with weight function, the shortest path from
 * u to v is the path whose sum of weights is the minimum.
 * 
 * All edge weights must be non-negative.
 * 
 * The runtime complexity is O(V + E) due to use of a YFastTrie as the min-priority heap.
 * Note that if the heap wrapper has to choose a Fibonacci instead due to
 * memory constraints, the runtime complexity is O(V*log_2V + E) instead.
 * 
 * implemented from pseudocode from Cormen et al. "Introduction to Algorithms".
 * 
 * from Cormen et al. :
 * 
 * @author nichole
 */
public class Dijkstras {
    
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

    // key is total estimate from srcIdx to destIdx for the given refIdx
    //    (that is the distance from srcIdx to refIdx + refIdx + heuristic)
    protected MinHeapForRT2012 heap = null;

    // refs to nodes internal to heap for decrease key operations
    protected HeapNode[] nodes = null;
    
    public Dijkstras() {
    }
    
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
     */
    public Dijkstras(SimpleLinkedListNode[] dAG, TIntIntMap[] weights, int sourceVertex) {
        
        if (dAG == null || dAG.length == 0) {
            throw new IllegalArgumentException("dAG cannot be null");
        }
        if (sourceVertex < 0 || sourceVertex >= dAG.length) {
            throw new IllegalArgumentException("sourceIndex cannot be null");
        }
        
        /*
        initialize single source (g, s)
        s=0
        add vertexes to heap
        while heap is not empty
          do u = extract-min
            s = s union u
            for each vertex v in adj[u]
               do relax(u, v, w)
        */ 
        
        init(dAG, weights, sourceVertex);        
    }
        
    /**
     *  find the single shortest path in dAG with edge weights w starting from s.
     */
    public void find() {
                               
        while (heap.getNumberOfNodes() > 0) {

            HeapNode uNode = heap.extractMin();
            
            int u = ((Integer)uNode.getData()).intValue();
            
            TIntIntMap uWeights = w[u];
            
            // null the entry in nodes so it isn't used in decrease key
            nodes[u] = null;
            
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
                if (nodes[v] == null) {
                    vNode = vNode.getNext();                   
                    continue;
                }
                int wUV = uWeights.get(v);
            
                int dUPlusWUV = wUV;
                if (dist[u] == sentinel) {
                    dUPlusWUV = sentinel;
                } else {
                    dUPlusWUV += dist[u];
                }
                                
                if (dist[v] > dUPlusWUV) {
                    dist[v] = dUPlusWUV;
                    predecessor[v] = u;
                    heap.decreaseKey(nodes[v], dUPlusWUV);
                }
                
                vNode = vNode.getNext();
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
        
        maxValue = calcUpperLimitKeyValue();
        
        sentinel = maxValue + 1;
    
        dist = new int[g.length];
        predecessor = new int[g.length];
        
        Arrays.fill(dist, sentinel);
        Arrays.fill(predecessor, -1);
        
        dist[src] = 0;
                
        initHeap();
    }
    
    private void initHeap() {
        
        int n = g.length;
                
        int nBits = (int)Math.ceil(Math.log(maxValue/Math.log(2)));
        
        //int maxValue, int approxN, int maxNumberOfBits
        heap = new MinHeapForRT2012(sentinel, n, nBits);

        nodes = new HeapNode[n];

        // initialize all except the source node as having infinite distance
        for (int i = 0; i < g.length; ++i) {

            if (g[i] == null) {
                continue;
            }
            
            int key = dist[i];

            HeapNode node = new HeapNode(key);
            node.setData(Integer.valueOf(i));

            heap.insert(node);

            nodes[i] = node;
        }            
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
    
}
