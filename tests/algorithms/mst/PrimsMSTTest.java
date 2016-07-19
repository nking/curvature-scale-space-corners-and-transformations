package algorithms.mst;

import algorithms.util.PairInt;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
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
        
        int[] treeWalk = prims.getPreOrderWalkOfTree();
        
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
