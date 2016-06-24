package algorithms.mst;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.Stack;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PrimsMSTTest extends TestCase {
    
    public PrimsMSTTest() {
    }
    
    /**
     * test created from Cormen et al. Chap 24, MST, Fig 23.4.
     * @throws Exception 
     */
    public void test0() throws Exception {
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4    |          2     \         |    9
         *  [a]       11      [i]         4     14      [e]
         *       8    |     7                    |   10
         *         \ [h] / -- 1  --[g]-- 2 -- \ [f]/
         *   
         * letter  index        
         *    a      0
         *    b      1
         *    c      2
         *    d      3
         *    e      4
         *    f      5
         *    g      6
         *    h      7
         *    i      8
         */

        int nVertexes = 9;
        TIntObjectMap<Set<PairInt>> adjCostMap 
            = new TIntObjectHashMap<Set<PairInt>>();
        
        adjCostMap.put(0, new HashSet<PairInt>());
        adjCostMap.get(0).add(new PairInt(1, 4));
        adjCostMap.get(0).add(new PairInt(7, 8));
        
        adjCostMap.put(1, new HashSet<PairInt>());
        adjCostMap.get(1).add(new PairInt(0, 4));
        adjCostMap.get(1).add(new PairInt(7, 11));
        
        adjCostMap.put(2, new HashSet<PairInt>());
        adjCostMap.get(2).add(new PairInt(3, 7));
        adjCostMap.get(2).add(new PairInt(5, 4));
        adjCostMap.get(2).add(new PairInt(8, 2));
        
        adjCostMap.put(3, new HashSet<PairInt>());
        adjCostMap.get(3).add(new PairInt(2, 7));
        adjCostMap.get(3).add(new PairInt(4, 9));
        adjCostMap.get(3).add(new PairInt(5, 14));
        
        adjCostMap.put(4, new HashSet<PairInt>());
        adjCostMap.get(4).add(new PairInt(5, 10));
        adjCostMap.get(4).add(new PairInt(3, 9));
        
        adjCostMap.put(5, new HashSet<PairInt>());
        adjCostMap.get(5).add(new PairInt(4, 10));
        adjCostMap.get(5).add(new PairInt(2, 4));
        adjCostMap.get(5).add(new PairInt(6, 2));
        
        adjCostMap.put(6, new HashSet<PairInt>());
        adjCostMap.get(6).add(new PairInt(5, 2));
        adjCostMap.get(6).add(new PairInt(7, 1));
        
        adjCostMap.put(7, new HashSet<PairInt>());
        adjCostMap.get(7).add(new PairInt(0, 8));        
        adjCostMap.get(7).add(new PairInt(1, 11));
        adjCostMap.get(7).add(new PairInt(8, 7));
        adjCostMap.get(7).add(new PairInt(6, 1));
      
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(
            nVertexes, adjCostMap);                
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4               2     \              9
         *  [a]              [i]          4            [e]
         *       8                                  
         *         \ [h]   -- 1  --[g]-- 2 -- \ [f]
         */
        /*
        0 1  
        2 3  
        2 5  
        2 8  
        3 4  
        5 6  
        6 7  
        7 0  
        */
        
        int[] predecessorArray = prims.getPrecessorArray();
        assertTrue(Arrays.equals(
            new int[]{-1, 0, 5, 2, 3, 6, 7, 0, 2}, 
            predecessorArray));
        
        HeapNode mst = prims.extractResultsAsTree();
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4               2     \              9
         *  [a]              [i]          4            [e]
         *       8                                  
         *         \ [h]   -- 1  --[g]-- 2 -- \ [f]
        */
        
        // traverse the tree to find these from root=0 
        /*
        0 1, 0 7
        7 6
        6 5
        5 2
        2 8, 2 3
        3 4 
        */
        assertEquals(0, mst.getKey());
        assertEquals(2, mst.getNumberOfChildren());
        
        // map of expected parent and child keys
        TIntObjectMap<TIntSet> expected 
            = new TIntObjectHashMap<TIntSet>();
        expected.put(0, new TIntHashSet());
        expected.get(0).add(1);
        expected.get(0).add(7);
        
        expected.put(7, new TIntHashSet());
        expected.get(7).add(6);
        
        expected.put(6, new TIntHashSet());
        expected.get(6).add(5);
        
        expected.put(5, new TIntHashSet());
        expected.get(5).add(2);
        
        expected.put(2, new TIntHashSet());
        expected.get(2).add(8);
        expected.get(2).add(3);
        
        expected.put(3, new TIntHashSet());
        expected.get(3).add(4);
         
        // traverse the tree using pre-order traversal
        //  root then left to right children
        //level 0            [0]
        //level 1     [1]           [4]
        //level 2   [2] [3]       [5] [6]
        //process node sees 0,1,2,3,4,5,6
        /*
        node=0, s=() p=()
           get children and add to s
        node=1, s=(1,4) p=(0)
           get children and add to s
        node=2, s=(2,3,1,4) p=(0,1)
           set node = null 
           visits else block
        node=null, s=(2,3,1,4) p=(0,1,2)
           pop from s node=2, but see it in p, so keep popping, node=3
        node=3, s=(1,4) p=(0,1,2)
           set node = null, visits else block
        node=null, s=(1,4) p=(0,1,2,3)
           pop from s, node=1, in p, so node=4
        node=4, s=(), p=(0,1,2,3)
           get children and add to s
        node=5, s=(5,6), p=(0,1,2,3,4)
           set node=null, visit else block
        node=null, s=(5,6), p=(0,1,2,3,4,5)
           pop from s, node=5, but see it in p, so, pop node=6
        node=6, s=(), p=(0,1,2,3,4,5)
                p=(0,1,2,3,4,5,6)
        */       
        
        HeapNode node = mst;
        HeapNode parentNode = null;
        Stack<HeapNode> stack = new Stack<HeapNode>();
        TIntSet visited = new TIntHashSet();
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {

                // process node: add to visited and 
                // remove all node->children values in expectedMap
                int index = (int)node.getKey();
                
                visited.add(index);
                
                // this is node = node.getLeft()
                DoubleLinkedCircularList children =
                    node.getChildren();
                if (children != null && (children.getNumberOfNodes() > 0)) {
                    // put all the children into stack
                    // and remove from expected map
                    HeapNode nodeC = children.getSentinel().getRight();
                    node = nodeC;
                    int indexC = (int)nodeC.getKey();
                    
                    assertTrue(expected.containsKey(index));
                    assertTrue(expected.get(index).remove(indexC));
                    if (expected.get(index).isEmpty()) {
                        expected.remove(index);
                    }
                    
                    for (int i = 0; i < (children.getNumberOfNodes() - 1); ++i) {
                        nodeC = nodeC.getRight();
                        stack.add(nodeC);
                        
                        indexC = (int)nodeC.getKey();
                       
                        assertTrue(expected.containsKey(index));
                        assertTrue(expected.get(index).remove(indexC));
                        if (expected.get(index).isEmpty()) {
                            expected.remove(index);
                        }
                    }
                } else {
                    node = null;
                }
            } else {
                node = stack.pop();
                int index = (int)node.getKey();
                while (visited.contains(index)) {
                    node = stack.pop();
                    index = (int)node.getKey();
                }
            }
        }
        
        assertEquals(9, visited.size());
        
        assertEquals(0, expected.size());
    }

    public void test1() throws Exception {
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4    |          2     \         |    9
         *  [a]       11      [i]         4     14      [e]
         *       8    |     7                    |   10
         *         \ [h] / -- 1  --[g]-- 2 -- \ [f]/
         *   
         * letter  index        
         *    a      0
         *    b      1
         *    c      2
         *    d      3
         *    e      4
         *    f      5
         *    g      6
         *    h      7
         *    i      8
        
         *         / [1]           /[2]-- 7  -- [3]\
         *       4    |          2     \         |    9
         *  [0]       11      [8]         4     14      [4]
         *       8    |     7                    |   10
         *         \ [7] / -- 1  --[6]-- 2 -- \ [5]/
         */

        int nVertexes = 9;
        TIntObjectMap<Set<PairInt>> adjCostMap 
            = new TIntObjectHashMap<Set<PairInt>>();
        
        adjCostMap.put(0, new HashSet<PairInt>());
        adjCostMap.get(0).add(new PairInt(1, 4));
        adjCostMap.get(0).add(new PairInt(7, 8));
        
        adjCostMap.put(2, new HashSet<PairInt>());
        adjCostMap.get(2).add(new PairInt(3, 7));
        adjCostMap.get(2).add(new PairInt(8, 2));
        
        adjCostMap.put(3, new HashSet<PairInt>());
        adjCostMap.get(3).add(new PairInt(4, 9));
        
        adjCostMap.put(5, new HashSet<PairInt>());
        adjCostMap.get(5).add(new PairInt(2, 4));
        
        adjCostMap.put(6, new HashSet<PairInt>());
        adjCostMap.get(6).add(new PairInt(5, 2));
        
        adjCostMap.put(7, new HashSet<PairInt>());
        adjCostMap.get(7).add(new PairInt(6, 1));
      
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(
            nVertexes, adjCostMap);                
        
        /*         / [1]           /[2]-- 7  -- [3]\
         *       4               2     \              9
         *  [0]              [8]          4            [4]
         *       8
         *         \ [7]   -- 1  --[6]-- 2 -- \ [5]
        
        [-1, 0, 5, 2, 3, 6, 7, 0, 2]
             *  *  *  *  *  *  *  
         */
     
        int[] predecessorArray = prims.getPrecessorArray();
       
        assertTrue(Arrays.equals(
            new int[]{-1, 0, 5, 2, 3, 6, 7, 0, 2}, 
            predecessorArray));
                
    }
}
