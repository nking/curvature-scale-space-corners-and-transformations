package algorithms.shortestPath;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class AStarTest extends TestCase {
    
    public AStarTest() {
    }
    
    public void testSearch0() throws Exception {
        /*
         *     5
         *
         *     4 [8]-----------[6]
         *        |             |
         *     3  |       [3]   |
         *        |       /     |
         *     2 [7]   [2]-[4]--.
         *       /    /
         *     1[5] [1]
         *      | /                     [9]
         *-1  [0]/  1   2   3   4   5
         *
         *    -1
         */

        PairInt[] points = new PairInt[10];
        points[0] = new PairInt(0, 0);
        points[1] = new PairInt(1, 1);
        points[2] = new PairInt(2, 2);
        points[3] = new PairInt(3, 3);
        points[4] = new PairInt(3, 2);
        points[5] = new PairInt(0, 1);
        points[6] = new PairInt(4, 4);
        points[7] = new PairInt(0, 2);
        points[8] = new PairInt(0, 4);
        points[9] = new PairInt(6, 0);
        
        List<LinkedList<Integer>> adjList = new ArrayList<LinkedList<Integer>>();
        
        // index 0
        LinkedList<Integer> list = new LinkedList<Integer>();
        list.add(Integer.valueOf(1));
        list.add(Integer.valueOf(5));
        adjList.add(list);
        
        // index 1
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(0));
        list.add(Integer.valueOf(2));
        adjList.add(list);
        
        // index 2
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(1));
        list.add(Integer.valueOf(3));
        list.add(Integer.valueOf(4));
        adjList.add(list);
        
        // index 3
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(2));
        adjList.add(list);
        
        // index 4
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(2));
        list.add(Integer.valueOf(6));
        adjList.add(list);
        
        // index 5
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(0));
        list.add(Integer.valueOf(7));
        adjList.add(list);
        
        // index 6
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(4));
        list.add(Integer.valueOf(8));
        adjList.add(list);
        
        // index 7
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(5));
        list.add(Integer.valueOf(8));
        adjList.add(list);
        
        // index 8
        list = new LinkedList<Integer>();
        list.add(Integer.valueOf(6));
        list.add(Integer.valueOf(7));
        adjList.add(list);
        
        // index 9
        list = new LinkedList<Integer>();
        adjList.add(list);
        
        int srcIndx = 0; 
        int destIndx = 3;
        AStar aStar = new AStar(points, adjList, srcIndx, destIndx);
        
        int[] indexes = aStar.search();
        assertTrue(indexes.length == 4);
        assertEquals(0, indexes[0]);
        assertEquals(1, indexes[1]);
        assertEquals(2, indexes[2]);
        assertEquals(3, indexes[3]);
        
        // ----
        srcIndx = 3; 
        destIndx = 0;
        aStar = new AStar(points, adjList, srcIndx, destIndx);
        indexes = aStar.search();
        assertTrue(indexes.length == 4);
        assertEquals(3, indexes[0]);
        assertEquals(2, indexes[1]);
        assertEquals(1, indexes[2]);
        assertEquals(0, indexes[3]);
        
        srcIndx = 0; 
        destIndx = 9;
        aStar = new AStar(points, adjList, srcIndx, destIndx);
        indexes = aStar.search();
        assertNull(indexes);
    }
         /*
         *     5
         *
         *     4 [8]-----------[6]
         *        |             |
         *     3  |       [3]   |
         *        |       /     |
         *     2 [7]   [2]-[4]--.
         *       /    /
         *     1[5] [1]
         *      | /                     [9]
         *-1  [0]/  1   2   3   4   5
         *
         *    -1
         */
}
