package algorithms.mst;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.heapsAndPQs.HeapNode;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

/**
 * Minimum spanning tree is the minimal network that spans all nodes in a tree
 * and has the smallest cost (sum of edges).
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
 * this implementation uses a Fibonacci heap and adjacency list
 *      
 * @author nichole
 */
public class PrimsMSTFaster {

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

        // initialize heap
        for (int i = 0; i < nVertexes; i++) {
            HeapNode v = new HeapNode();
            if (i == 0) {
                v.setKey(0);
            } else {
                v.setKey(sentinel);
            }
            v.setData(Integer.valueOf(i));
            heap.insert(v);
            nodes.add(v);
        }
        
        while (heap.getNumberOfNodes() > 0) {

            HeapNode u = heap.extractMin(); 
           
            Integer uIndex = (Integer)u.getData();
            int uIdx = uIndex.intValue();
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
    
    public int[] getPreOrderWalkOfTree() {
        
        //pre-order is
        //root, left subtree, right subtree
        //given the top node as the starter
        //level 0            [0]
        //level 1      [1]            [6]
        //level 2   [2]   [3]       [7] [8]
        //level 3       [4][5]
        // process node sees 0,1,2,3,4,5,6,7,8 

        TIntObjectMap<TIntList> nodeMap = 
            createReverseMap();
                
        int count = 0;
        int[] walk = new int[prev.length];
        
        Integer node = Integer.valueOf(0);
        
        TIntSet inW = new TIntHashSet();
        
        Stack<Integer> stack = new Stack<Integer>();
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                //process node
                if (!inW.contains(node.intValue())) {
                    walk[count] = node.intValue();
                    count++;
                    inW.add(node.intValue());
                }
                stack.push(node);
                int origIdx = node.intValue();
                TIntList children = nodeMap.get(node.intValue());
                if (children == null) {
                    node = null;
                } else {
                    node = children.get(0);
                    //NOTE: an expensive delete. could change structure
                    boolean rm = children.remove(node.intValue());
                    assert(rm);
                }
                if ((children != null) && children.isEmpty()) {
                    nodeMap.remove(origIdx);
                }
            } else {
                node = stack.pop();
                int origIdx = node.intValue();
                TIntList children = nodeMap.get(node.intValue());
                if (children == null) {
                    node = null;
                } else {
                    node = children.get(0);
                    boolean rm = children.remove(node.intValue());
                    assert(rm);
                }
                if ((children != null) && children.isEmpty()) {
                    nodeMap.remove(origIdx);
                }
            }
        }
        
        return walk;
    }
    
    protected LinkedList<Integer> getPreOrderWalkOfTreeWithMarkers(
        TIntObjectMap<TIntList> nodeMap) {
        
        //pre-order is
        //root, left subtree, right subtree
        //given the top node as the starter
        //level 0            [0]
        //level 1      [1]            [6]
        //level 2   [2]   [3]       [7] [8]
        //level 3       [4][5]
        // process node sees 0,1,2,-1, 3,4,-1, 5, -1, 6,7,-1, 8, -1 
        // (-1's added where no children)
      
        TIntSet added = new TIntHashSet();
        
        // key = node, map = children
        TIntObjectMap<LinkedList<Integer>> cMap 
            = new TIntObjectHashMap<LinkedList<Integer>>();
        
        LinkedList<Integer> walk = new LinkedList<Integer>();
                
        Integer node = Integer.valueOf(0);
                
        Stack<Integer> stack = new Stack<Integer>();
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                int idx = node.intValue();
                if (!added.contains(idx)) {
                    walk.add(idx);
                    added.add(idx);
                    
                    TIntList c = nodeMap.get(idx);
                    if (c != null) {
                        LinkedList<Integer> cL = new LinkedList<Integer>();
                        cMap.put(idx, cL);
                        TIntIterator iter = c.iterator();
                        while (iter.hasNext()) {
                            cL.add(Integer.valueOf(iter.next()));
                        }
                    } else {
                        walk.add(-1);
                    }
                }
                stack.push(node);
                if (!cMap.containsKey(idx)) {
                    node = null;
                } else {
                    LinkedList<Integer> cL = cMap.get(idx);
                    node = cL.removeFirst();
                    if (cL.isEmpty()) {
                        cMap.remove(idx);
                    }
                }
            } else {
                // add a marker
                walk.add(-1);
                node = stack.pop();
                int idx = node.intValue();
                if (!cMap.containsKey(idx)) {
                    node = null;
                } else {
                    LinkedList<Integer> cL = cMap.get(idx);
                    node = cL.removeFirst();
                    if (cL.isEmpty()) {
                        cMap.remove(idx);
                    }
                }
            }
        }
        
        return walk;
    }
    
    protected LinkedList<Integer> getPostOrderWalkOfTreeWithMarkers(
        TIntObjectMap<TIntList> nodeMap) {
        
        //post-order traversal:  left subtree, right subtree, root
        // given the top node as the starter
        //level 0            [0]
        //level 1      [1]            [6]
        //level 2   [2]   [3]       [7] [8]
        //level 3       [4][5]
        // process node sees  -1, 2,-1,4,-1,5,3,1,-1,7,-1,8,6,0
        //  (-1's added where there were no children)
                
        ArrayDeque<Integer> children = new ArrayDeque<Integer>();
        
        LinkedList<Integer> walk = new LinkedList<Integer>();
               
        Stack<Integer> stack = new Stack<Integer>();
        Stack<Integer> stack2 = new Stack<Integer>();
        
        Integer node = Integer.valueOf(0);
        
        stack.push(node);
            
        while (!stack.isEmpty()) {
            node = stack.pop();
            stack2.push(node);
            int idx = node.intValue();
            TIntList c = nodeMap.get(idx);
            if (c != null) {
                TIntIterator iter = c.iterator();
                while (iter.hasNext()) {
                    stack.push(iter.next());
                }
            } else {
                stack2.push(Integer.valueOf(-1));
            }
        }
        
        // remove first -1
        if (!stack2.isEmpty()) {
            stack2.pop();
        }

        while (!stack2.isEmpty()) {
            node = stack2.pop();
            walk.add(node.intValue());
        }
        
        return walk;
    }
    
    public TIntObjectMap<TIntList> createReverseMap() {
    
        TIntObjectMap<TIntList> revPrevMap 
            = new TIntObjectHashMap<TIntList>();
        
        for (int i = 0; i < prev.length; ++i) {
            int parentIdx = prev[i];
            if (parentIdx == -1) {
                continue;
            }
            TIntList indexes = revPrevMap.get(parentIdx);
            if (indexes == null) {
                indexes = new TIntArrayList();
                revPrevMap.put(parentIdx, indexes);
            }
            indexes.add(i);
        }
                
        return revPrevMap;
    }
    
}
