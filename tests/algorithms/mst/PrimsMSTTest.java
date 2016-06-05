package algorithms.mst;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
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
        Map<Integer, Set<PairInt>> adjCostMap 
            = new HashMap<Integer, Set<PairInt>>();
        
        adjCostMap.put(Integer.valueOf(0), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(0)).add(new PairInt(1, 4));
        adjCostMap.get(Integer.valueOf(0)).add(new PairInt(7, 8));
        
        adjCostMap.put(Integer.valueOf(1), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(1)).add(new PairInt(0, 4));
        adjCostMap.get(Integer.valueOf(1)).add(new PairInt(7, 11));
        
        adjCostMap.put(Integer.valueOf(2), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(2)).add(new PairInt(3, 7));
        adjCostMap.get(Integer.valueOf(2)).add(new PairInt(5, 4));
        adjCostMap.get(Integer.valueOf(2)).add(new PairInt(8, 2));
        
        adjCostMap.put(Integer.valueOf(3), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(3)).add(new PairInt(2, 7));
        adjCostMap.get(Integer.valueOf(3)).add(new PairInt(4, 9));
        adjCostMap.get(Integer.valueOf(3)).add(new PairInt(5, 14));
        
        adjCostMap.put(Integer.valueOf(4), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(4)).add(new PairInt(5, 10));
        adjCostMap.get(Integer.valueOf(4)).add(new PairInt(3, 9));
        
        adjCostMap.put(Integer.valueOf(5), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(5)).add(new PairInt(4, 10));
        adjCostMap.get(Integer.valueOf(5)).add(new PairInt(2, 4));
        adjCostMap.get(Integer.valueOf(5)).add(new PairInt(6, 2));
        
        adjCostMap.put(Integer.valueOf(6), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(6)).add(new PairInt(5, 2));
        adjCostMap.get(Integer.valueOf(6)).add(new PairInt(7, 1));
        
        adjCostMap.put(Integer.valueOf(7), new HashSet<PairInt>());
        adjCostMap.get(Integer.valueOf(7)).add(new PairInt(0, 8));        
        adjCostMap.get(Integer.valueOf(7)).add(new PairInt(1, 11));
        adjCostMap.get(Integer.valueOf(7)).add(new PairInt(8, 7));
        adjCostMap.get(Integer.valueOf(7)).add(new PairInt(6, 1));
      
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
        Map<Integer, Set<Integer>> expected = 
            new HashMap<Integer, Set<Integer>>();
        expected.put(Integer.valueOf(0), new HashSet<Integer>());
        expected.get(Integer.valueOf(0)).add(Integer.valueOf(1));
        expected.get(Integer.valueOf(0)).add(Integer.valueOf(7));
        
        expected.put(Integer.valueOf(7), new HashSet<Integer>());
        expected.get(Integer.valueOf(7)).add(Integer.valueOf(6));
        
        expected.put(Integer.valueOf(6), new HashSet<Integer>());
        expected.get(Integer.valueOf(6)).add(Integer.valueOf(5));
        
        expected.put(Integer.valueOf(5), new HashSet<Integer>());
        expected.get(Integer.valueOf(5)).add(Integer.valueOf(2));
        
        expected.put(Integer.valueOf(2), new HashSet<Integer>());
        expected.get(Integer.valueOf(2)).add(Integer.valueOf(8));
        expected.get(Integer.valueOf(2)).add(Integer.valueOf(3));
        
        expected.put(Integer.valueOf(3), new HashSet<Integer>());
        expected.get(Integer.valueOf(3)).add(Integer.valueOf(4));
         
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
        Set<Integer> visited = new HashSet<Integer>();
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {

                // process node: add to visited and 
                // remove all node->children values in expectedMap
                Integer index = Integer.valueOf((int)node.getKey());
                
                visited.add(index);
                
                // this is node = node.getLeft()
                DoubleLinkedCircularList children =
                    node.getChildren();
                if (children != null && (children.getNumberOfNodes() > 0)) {
                    // put all the children into stack
                    // and remove from expected map
                    HeapNode nodeC = children.getSentinel().getRight();
                    node = nodeC;
                    Integer indexC = Integer.valueOf((int)nodeC.getKey());
                    
                    assertTrue(expected.containsKey(index));
                    assertTrue(expected.get(index).remove(indexC));
                    if (expected.get(index).isEmpty()) {
                        expected.remove(index);
                    }
                    
                    for (int i = 0; i < (children.getNumberOfNodes() - 1); ++i) {
                        nodeC = nodeC.getRight();
                        stack.add(nodeC);
                        
                        indexC = Integer.valueOf((int)nodeC.getKey());
                       
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
                Integer index = Integer.valueOf((int)node.getKey());
                while (visited.contains(index)) {
                    node = stack.pop();
                    index = Integer.valueOf((int)node.getKey());
                }
            }
        }
        
        assertEquals(9, visited.size());
        
        assertEquals(0, expected.size());
    }

}
