package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TopologicalSortTest extends TestCase {

    public void testSort() {
        
        /*
        from Cormen et al. Fig 22.7
        
        11/16 undershorts     socks 17/18    
               |         \     |               watch 9/10
               V          V    V
          12/15 pants ----> shoes 13/14
               |
               |     shirt  1/8
               V     /   V
          6/7 belt <    tie  2/5
                  \      |
                   \     V
                    > jacket 3/4
        
        // socks, undershorts, pants, shoes, watch, shirt, belt, tie, jacket
        */
        Map<String, Integer> objs = new HashMap<String, Integer>();
        objs.put("shirt", 0);      // belt, tie
        objs.put("belt", 1);       // jacket
        objs.put("tie", 2);        // jacket
        objs.put("jacket", 3);     //
        objs.put("watch", 4);      //
        objs.put("undershorts", 5);//  pants, shoes
        objs.put("pants", 6);      //  shoes
        objs.put("shoes", 7);      //
        objs.put("socks", 8);      //  shoes
        
        SimpleLinkedListNode[] connected = new SimpleLinkedListNode[9];
        for (int i = 0; i < connected.length; ++i) {
            connected[i] = new SimpleLinkedListNode();
        }
        connected[objs.get("shirt")].insert(objs.get("belt")); connected[objs.get("shirt")].insert(objs.get("tie"));
        connected[objs.get("belt")].insert(objs.get("jacket"));
        connected[objs.get("tie")].insert(objs.get("jacket"));
        connected[objs.get("undershorts")].insert(objs.get("pants")); connected[objs.get("undershorts")].insert(objs.get("shoes"));
        connected[objs.get("pants")].insert(objs.get("shoes"));
        connected[objs.get("socks")].insert(objs.get("shoes"));
        // socks, undershorts, pants, shoes, watch, shirt, belt, tie, jacket
        int[] expResult = new int[]{objs.get("socks"), objs.get("undershorts"),
            objs.get("pants"), objs.get("shoes"), objs.get("watch"),
            objs.get("shirt"), objs.get("belt"), objs.get("tie"),
            objs.get("jacket")
        };

        TopologicalSort ts = new TopologicalSort(connected);
        
        int[] result = ts.sort();
        
        //System.out.println("expected = " + Arrays.toString(expResult));
        //System.out.println("result = " + Arrays.toString(result));
        
        assertTrue(Arrays.equals(expResult, result));
    }
    
    /**
     * Test of sort method, of class TopologicalSort.
     */
    public void estSort_SimpleDAG() {

        System.out.println("testSort_SimpleDAG");
        
        // constructing tests from MIT open courseware
        // network_optimization/MIT15_082JF10_av03.pdf
        
        SimpleLinkedListNode[] connected = new SimpleLinkedListNode[8];

        /*              <5> --------> <0>
         *            >  |       >     |
         *         .     V   .         V
         *      <4> --> <1>           <7> ---> <2>
         *      >  .                 >        .>
         *     .     .            .      .
         *    .         >      .   .
         *  <6> ------> <3> .
         */
        /*
        connected[0] = new SimpleLinkedListNode();
        connected[0].insert(7);

        connected[1] = new SimpleLinkedListNode();
        connected[1].insert(0);

        connected[2] = new SimpleLinkedListNode();

        connected[3] = new SimpleLinkedListNode();
        connected[3].insert(2);
        connected[3].insert(7);

        connected[4] = new SimpleLinkedListNode();
        connected[4].insert(1);
        connected[4].insert(3);
        connected[4].insert(5);

        connected[5] = new SimpleLinkedListNode();
        connected[5].insert(0);
        connected[5].insert(1);

        connected[6] = new SimpleLinkedListNode();
        connected[6].insert(3);
        connected[6].insert(4);

        connected[7] = new SimpleLinkedListNode();
        connected[7].insert(2);

        int[] expResult = new int[]{6, 4, 3, 5, 1, 0, 7, 2};

        TopologicalSort ts = new TopologicalSort(connected);
        
        int[] result = ts.sort();

        assertTrue(Arrays.equals(expResult, result));
        */
    }
   
    
    public void estSort2() {

        System.out.println("testSort2");
        
        // constructing test from Cormen et al.'s "Introduction to Algorithms"
        boolean[][] connected = null;

        /*    *0 \         *3
         *    ||    \      ||        *8
         *    \/       \-> \/
         *    *1 --------> *4
         *    ||
         *    ||    /*5
         *    ||   / ||
         *    \/<-/  ||
         *    *2     \/
         *      \    *6
         *       \   ||
         *       \   ||
         *        \->\/
         *           *7
         */
        connected = new boolean[9][];
        for (int i = 0; i < connected.length; i++) {
            connected[i] = new boolean[connected.length];
        }
        connected[0][1] = true;
        connected[0][4] = true;
        connected[1][2] = true;
        connected[1][4] = true;
        connected[2][7] = true;
        connected[3][4] = true;
        connected[5][2] = true;
        connected[5][6] = true;
        connected[6][7] = true;

        
        /*    *0 \         *3
         *    ||    \      ||        *8
         *    \/       \-> \/
         *    *1 --------> *4
         *    ||
         *    ||    /*5
         *    ||   / ||
         *    \/<-/  ||
         *    *2     \/
         *      \    *6
         *       \   ||
         *       \   ||
         *        \->\/
         *           *7

           sorted=[8, 5, 6, 3,  0,  4,  1, 2,  7]
                d=[1, 2, 3, 11, 7, 13, 14, 4, 17]
                f=[9, 8, 5, 11, 7, 15, 14, 4, 17]
                   0  1  2   3  4   5   6  7   8
                p=[i, 0, 1,  i, 1,  i,  5, 2,  i]

        This solution:
        0   1   3   4    5    2    6    7    8
        --->-------->    ---->---------->
        ------------>    --------->----->
                 --->              
           
                      
        Book solution:
        3   0   1   4    8    5    2    6    7
        ------------>         --------->----->
            --->---->         ---->---------->
            -------->
                                        
        So, they are equivalent in finding the arrangement with fewest overlapping edges,
        but the solution here is optimized to finish faster (if one thinks of the edges
        as start and finish times of a scheduled process, for example).
        
         */

        /*SimpleDAG dag = new SimpleDAG(connected);
        
        TopologicalSort ts = new TopologicalSort();
        
        int[] result = ts.sort(dag);
        
        int[] currentSolution = new int[]{0,   1,   3,   4,    5,    2,    6,    7,    8};
        int[] bookSolution = new int[]{3, 0, 1, 4, 8, 5, 2, 6, 7};

        
        System.out.println("sorted=" + Arrays.toString(result));
        System.out.println("currentSolution=" + Arrays.toString(currentSolution));
        System.out.println("bookSolution=" + Arrays.toString(bookSolution));

        assertTrue(Arrays.equals(currentSolution, result));
        */
    }
}
