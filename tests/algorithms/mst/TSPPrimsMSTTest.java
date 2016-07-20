package algorithms.mst;

import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class TSPPrimsMSTTest extends TestCase {
    
    public void test0() {
        
        /* Test of traveling salesman approximate tour. 
         * 
         * From Cormen et al. chap 35 Approximation Algorithms, 35.2 Traveling-Salesman problem, Fig 35.2
         * 
         * 
         *   6  |---|---|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   5  |---|--[A]--|--[D]--|---|---|----
         *      |   |   |   |   |   |   |   |
         *   4  |---|---|---|---|--[E]--|---|----
         *      |   |   |   |   |   |   |   |
         *   3  |---|--[B]--|--[F]--|--[G]--|----
         *      |   |   |   |   |   |   |   |
         *   2  |--[C]--|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   1  |---|---|--[H]--|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   0  |---|---|---|---|---|-------|----
         *      0   1   2   3   4   5   6   7
        
         TSP-approx:
        
          L is pre-order walk of mst prims tree
        H is the tour of those visits in that order

        an MST is
                  A0
            B1         D3
           C2 H7       E4
                      F5 G6
        
        Prims here is giving:
                    0
               1        3
            2         4
                    5  6
                  7

        pre-order is
        root, left subtree, right subtree
        // given the top node as the starter
        //level 0            [0]
        //level 1     [1]           [4]
        //level 2   [2] [3]       [5] [6]
        // process node sees 0,1,2,3,4,5,6

        then pre-order trversal as a tour:
        a, b, c, b, h, b, a, d, e, f, e, g, e, d, a
        removing vertexes already encountered:
        a, b, c,  h, d, e, f, g,

        a pre-order for prims mst here is:
                   0
               1        3
            2         4
                    5  6
                  7
        
         0, 1, 2, 3, 4, 5, 7, 6
        
        In contrast, an optimal cost, non-crossing:
        a, b, c,  h, f, g, e, d
        */
                
        List<PairInt> points = new ArrayList<PairInt>();
        points.add(new PairInt(2, 5));
        points.add(new PairInt(2, 3));
        points.add(new PairInt(1, 2));
        points.add(new PairInt(4, 5));
        points.add(new PairInt(5, 4));
        points.add(new PairInt(4, 3));
        points.add(new PairInt(6, 3));
        points.add(new PairInt(3, 1));
        
        int[] expected = new int[9];
        expected[0] = 0;
        expected[1] = 1;
        expected[2] = 2;
        expected[3] = 7;
        expected[4] = 3;
        expected[5] = 4;
        expected[6] = 5;
        expected[7] = 6;
        expected[8] = 0;
        /*A B C D E F G H
          0 1 2 3 4 5 6 7
        */
         
        TIntObjectMap<Set<PairInt>>
            adjCostMap = new TIntObjectHashMap<Set<PairInt>>();
        
        for (int i = 0; i < points.size(); ++i) {
            int x1 = points.get(i).getX();
            int y1 = points.get(i).getY();            
            Set<PairInt> set1 = adjCostMap.get(i);
            if (set1 == null) {
                set1 = new HashSet<PairInt>();
                adjCostMap.put(i, set1);
            }
            
            for (int j = 0; j < points.size(); ++j) {
                if (i == j) {
                    continue;
                }
                int x2 = points.get(j).getX();
                int y2 = points.get(j).getY();
                
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                int dist = diffX * diffX + diffY * diffY;
//System.out.println("i=" + i + " j=" + j + "  dist=" + dist);                
                set1.add(new PairInt(j, dist));
            }
        }
        
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(points.size(), 
            adjCostMap);
        int[] walk = prims.getPreOrderWalkOfTree();
        System.out.println("prims walk=" +
            Arrays.toString(walk));
        assertTrue(Arrays.equals(
            new int[]{0, 1, 2, 3, 4, 5, 7, 6}, walk));
        
        TSPPrimsMST tsp = new TSPPrimsMST();
        int[] tour = tsp.approxTSPTour(points.size(), adjCostMap);
        
        assertEquals(expected.length, tour.length);
        for (int i = 0; i < tour.length; ++i) {
            int tourIdx = tour[i];
            assertEquals(expected[i], tourIdx);
        }
    }
}
