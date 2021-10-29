package algorithms.mst;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.heapsAndPQs.HeapNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Stack;

/**
 * minimum spanning tree is the subset of edges in a weighted undirected graph
 * that connect all vertexes for a total minimum cost (sum of edge weights).
 * 
 *
 * Implemented from pseudo code in Cormen et al. Introduction to Algorithms and
 * from http://en.wikipedia.org/wiki/Prim's_algorithm.
   Useful also was
  http://www.i-programmer.info/projects/61-algorithms/534-minimum-spanning-tree.html?start=1
  
 * Time complexity for different implementations:
 *
 *     Minimum edge weight data structure	Time complexity (total)
 *     ----------------------------------   -----------------------
 *     adjacency matrix, searching          O(N^2)
 *     binary heap and adjacency list       O((N + E) lg2 N) = O(E lg2 N)
 *     Fibonacci heap and adjacency list	O(E + N lg2 N)
 *     YFastTrie and adjacency list         O(E + V)
 *     
 * this implementation uses a YFastTrie min priority queue and adjacency list.
 * 
 * @author nichole
 */
public class PrimsMST {

    private int[] prev = null;
    
    private TIntObjectMap<TIntIntMap> adjCostMap;
    
    /**
     * 
     * @param adjCostMap key=vertex1 index, 
     *   value=map with key = vertex2 index and
     *   value = cost of edge between vertex1 and vertex2.  Note that the 
     *   algorithm assumes the map key values are 0 through the number of vertexes
     *   without gaps.
     * @param maximumWeightInGraph the maximum value of any weight in the graph.
     * This sets the word size of the YFastTrie used as the min priority q which
     * is used by default if the VM has enough memory (a large number of items
     * requires more memory).  If the YFastTrie is expected to consume more 
     * memory than available, this class will use a Fibonacci Heap instead.
     * 
     */
    public void calculateMinimumSpanningTree(
        final TIntObjectMap<TIntIntMap> adjCostMap, 
        int maximumWeightInGraph) {

        this.adjCostMap = adjCostMap;
        
        int nVertexes = adjCostMap.size();
        
        boolean[] inQ = new boolean[nVertexes];
        Arrays.fill(inQ, true);
        prev = new int[nVertexes];
        Arrays.fill(prev, -1);
        
        int sentinel = (maximumWeightInGraph < Integer.MAX_VALUE) ? 
            (maximumWeightInGraph + 1) : Integer.MAX_VALUE;
        
        int maxNumberOfBits = (int)Math.ceil(Math.log(sentinel)/Math.log(2));
        
        MinHeapForRT2012 heap = new MinHeapForRT2012(sentinel,
            nVertexes, maxNumberOfBits);
   
        // extra data structure needed to hold references to the nodes to be
        //     able to read the nodes still in the queue.
        List<HeapNode> nodes = new ArrayList<HeapNode>();

        // initialize heap by adding all nodes
        for (int i = 0; i < nVertexes; i++) {
            HeapNode v = new HeapNode();
            if (i == 0) {
                v.setKey(0);
            } else {
                v.setKey(sentinel);
            }
            // i is the index in nodes list in inQ array
            v.setData(Integer.valueOf(i)); 
            heap.insert(v);
            nodes.add(v);
        }
        
        //O(|V|)
        while (heap.getNumberOfNodes() > 0) {

            // O(log_2 log_2(w_bits)) or O(log_2(|V|))
            HeapNode u = heap.extractMin(); 
           
            int uIdx = ((Integer)u.getData()).intValue();
            inQ[uIdx] = false;
            
            TIntIntMap adjMap0 = adjCostMap.get(uIdx);
            if (adjMap0 == null) {
                continue;
            }
            
            TIntIntIterator iter = adjMap0.iterator();
            for (int i = 0; i < adjMap0.size(); ++i) {
                iter.advance();                 
                int vIdx = iter.key();
                int cost = iter.value();
                long distV = nodes.get(vIdx).getKey();
               
                if (inQ[vIdx] && (cost < distV)) {
                    prev[vIdx] = uIdx;
                    heap.decreaseKey(nodes.get(vIdx), cost); 
                }
            }
        }
        
        //System.out.println(Arrays.toString(prev));
    }
    
    public int[] getPredeccessorArray() {
        if (prev == null) {
            return null;
        }
        return Arrays.copyOf(prev, prev.length);
    }
    
    public int findRoot() {
        // in prev array, since all vertices are connected, all except one
        // should have a predecessor.
        int root = -1;
        for (int i = 0; i < prev.length; ++i) {
            if (prev[i] == -1) {
                root = i;
                break;
            }
        }
        return root;
    }
    
    public Map<Integer, LinkedList<Integer>> makeTreeFromPrev() {
        Map<Integer, LinkedList<Integer>> tree = new HashMap<Integer, LinkedList<Integer>>();
        int parent;
        LinkedList<Integer> children;
        for (int child = 0; child < prev.length; ++child) {
            parent = prev[child];
            if (parent == -1) {
                continue;
            }
            children = tree.get(parent);
            if (children == null) {
                children = new LinkedList<Integer>();
                tree.put(parent, children);
            }
            children.add(child);
        }
        return tree;
    }
    
    /**
     * walk the tree in prev as a pre-order traversal and return the indexes
     * of the nodes in that order.
     * The pre-order traversal visits subtrees of root, then left, then right.
     * @return 
     */
    public TIntList getPreorderIndexes() {
        if (prev == null) {
            return null;
        }
        int root = findRoot();
        return getPreorderIndexes(root);
    }
    private TIntList getPreorderIndexes(int root) {
        TIntList out = new TIntArrayList();
        Map<Integer, LinkedList<Integer>> tree = makeTreeFromPrev();
        int node = root;
        int pNode;
        Stack<Integer> stack = new Stack<Integer>();
        stack.add(node);
        // removing the tree nodes as they are used
        LinkedList<Integer> children;
        while (!stack.isEmpty() || node != -1) {
            if (node != -1) {
                out.add(node);
                stack.push(node);
                pNode = node;
                children = tree.get(pNode);
                if (children != null) {
                    node = children.removeFirst();
                    if (children.isEmpty()) {
                        tree.remove(pNode);
                    }
                } else {
                    node = - 1;
                }
            } else {
                node = stack.pop();// discard as it's already in out
                // get next node from tree.
                // since it might be n-ary tree, will keep taking from front of list
                //    though there's no reason to choose first over last
                pNode = node;
                children = tree.get(pNode);
                if (children != null) {
                    node = children.removeFirst();
                    if (children.isEmpty()) {
                        tree.remove(pNode);
                    }
                } else {
                    node = - 1;
                }
            }
        }
        return out;
    }
    
    /**
     * 
     * @param adjCostMap adjacency map with cost.  key=index1, value = map
     * with key=index2 and value=cost.
     * @return maximum cost found in adjCostMap
     */
    public static int maxEdgeCost(TIntObjectMap<TIntIntMap> adjCostMap) {
        int max = Integer.MIN_VALUE;
        TIntObjectIterator<TIntIntMap> iter = adjCostMap.iterator();
        int i, j, idx1, idx2, c;
        TIntIntIterator iter2;
        TIntIntMap map;
        for (i = 0; i < adjCostMap.size(); ++i) {
            iter.advance();
            idx1 = iter.key();
            map = iter.value();
            iter2 = map.iterator();            
            for (j = 0; j < map.size(); ++j) {
                iter2.advance();
                idx2 = iter2.key();
                c = iter2.value();
                if (c > max) {
                    max = c;
                }
            }
        }
        return max;
    }    
}
