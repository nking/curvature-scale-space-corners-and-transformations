package algorithms.mst;

import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
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
            inQ[uIndex.intValue()] = false;
            
            Set<PairInt> adjSet = adjCostMap.get(uIndex.intValue());
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
    
    private TIntObjectMap<TIntSet> createReverseMap(
        int[] pi) {
    
        TIntObjectMap<TIntSet> revPrevMap 
            = new TIntObjectHashMap<TIntSet>();
        for (int i = 0; i < pi.length; ++i) {
            int idx = pi[i];
            if (idx == -1) {
                continue;
            }
            TIntSet indexes = revPrevMap.get(idx);
            if (indexes == null) {
                indexes = new TIntHashSet();
                revPrevMap.put(idx, indexes);
            }
            indexes.add(i);
        }
        return revPrevMap;
    }
    
    public HeapNode extractResultsAsTree() {
        
        // can either construct a tree from children and parent
        // pairs in prev
        // or can return a list of edges,
        
        // make a reverse map of prev
        TIntObjectMap<TIntSet> revPrevMap 
            = createReverseMap(prev);
        
        // here, the keys are the indexes
        HeapNode mst = new HeapNode();
        mst.setKey(0);
        
        Set<HeapNode> set = new HashSet<HeapNode>();
        HeapNode[] tmp = new HeapNode[prev.length];
        ArrayDeque<HeapNode> queue = new ArrayDeque<HeapNode>();                
        
        for (int i = 0; i < prev.length; i++) {
            HeapNode node;
            if (i != 0) {
                node = new HeapNode();
                node.setKey(i);
            } else {
                node = mst;
            }
            set.add(node);
            tmp[i] = node;
        }

        // we start w/ parent node = root  with the first read of pendingQueue
        HeapNode currentParentNode = null;
        queue.add(mst);
        
        while (!set.isEmpty()) {
            while (!queue.isEmpty()) {
                currentParentNode = queue.poll();
                boolean removed = set.remove(currentParentNode);
                
                int currentParentIdx = (int)currentParentNode.getKey();

                //revPrevMap has key = prev[i], values=set{i}
                TIntSet childrenIndexes = 
                    revPrevMap.get(currentParentIdx);
                if (childrenIndexes == null) {
                    continue;
                }
                TIntIterator iter = childrenIndexes.iterator();
                while (iter.hasNext()) {
                    int childIdx = iter.next();                
                    HeapNode currentChildNode = tmp[childIdx];
                    currentParentNode.addChild(currentChildNode);
                    queue.add(currentChildNode);
                    removed = set.remove(currentChildNode);
                }
            }
            // because the mst has leaves, there may be nodes left in the map still:
            if (!set.isEmpty()) {
                
                currentParentNode = set.iterator().next();
                int currentParentIdx = (int)currentParentNode.getKey();
                boolean removed = set.remove(currentParentNode);
                
                TIntSet childrenIndexes = 
                    revPrevMap.get(currentParentIdx);
                if (childrenIndexes == null) {
                    continue;
                }
                TIntIterator iter2 = childrenIndexes.iterator();
                while (iter2.hasNext()) {
                    int childIdx = iter2.next();
                    HeapNode currentChildNode = tmp[childIdx];
                    currentParentNode.addChild(currentChildNode);
                    queue.add(currentChildNode);
                    removed = set.remove(currentChildNode);
                }
            }
        }

        return mst;
    }

}
