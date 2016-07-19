package algorithms.mst;

import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
 *     
 * this implementation uses a Fibonacci heap and adjacency list
 *      
 * @author nichole
 */
public class PrimsMST {

    private int[] prev = null;
    
    private TIntObjectMap<Set<PairInt>> adjCostMap;
    
    /**
     * 
     * @param nVertexes
     * @param adjCostMap key=vertex index, 
     *   value=set of pairint where each pairint has 
     *   x = adjacent vertex and y = cost of edge.
     * @return 
     */
    public void calculateMinimumSpanningTree(
        final int nVertexes, final TIntObjectMap<Set<PairInt>>
            adjCostMap) {

        this.adjCostMap = adjCostMap;
        
        boolean[] inQ = new boolean[nVertexes];
        Arrays.fill(inQ, true);
        prev = new int[nVertexes];
        Arrays.fill(prev, -1);
        
        Heap heap = new Heap();
        
        List<HeapNode> nodes = new ArrayList<HeapNode>();

        // initialize heap
        for (int i = 0; i < nVertexes; i++) {
        	HeapNode v = new HeapNode();
        	if (i == 0) {
                v.setKey(0);
            } else {
                v.setKey(Integer.MAX_VALUE);
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
            
            Set<PairInt> adjSet = adjCostMap.get(uIdx);
            if (adjSet == null) {
                continue;
            }
            
            for (PairInt indexCost : adjSet) {
                int vIdx = indexCost.getX();
                int cost = indexCost.getY();
                long distV = nodes.get(vIdx).getKey();
               
                if (inQ[vIdx] && (cost < distV)) {
                    prev[vIdx] = uIndex.intValue();
                    heap.decreaseKey(nodes.get(vIdx), cost); 
                }
            }
        }
        
        //System.out.println(Arrays.toString(prev));
    }
    
    public int[] getPrecessorArray() {
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
        //level 1     [1]           [4]
        //level 2   [2] [3]       [5] [6]
        // process node sees 0,1,2,3,4,5,6
        
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
    
    private TIntObjectMap<TIntList> createReverseMap() {
    
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
