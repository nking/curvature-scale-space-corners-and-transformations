package algorithms.search;

import algorithms.util.LinkedListCostNode;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.list.TIntList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BeamSearchTest extends TestCase {
    
    public BeamSearchTest() {
    }
    
    public void testSearch() {
        int k = 2;
        
        /*
                        0
           1                      2
         3    4            5         6         7
       8  9  10  11     12 13 14   15 16 17  18 19 20
        */
        SimpleLinkedListNode[] adjList = new SimpleLinkedListNode[21];
        for (int i = 0; i < 21; i++) {
            adjList[i] = new SimpleLinkedListNode(i);
        }
        adjList[0].insert(1);
        adjList[0].insert(2);
        
        adjList[1].insert(3);
        adjList[1].insert(4);
        
        adjList[2].insert(5);
        adjList[2].insert(6);
        adjList[2].insert(7);
        
        adjList[3].insert(8);
        adjList[3].insert(9);
        
        adjList[4].insert(10);
        adjList[4].insert(11);
        
        adjList[5].insert(12);
        adjList[5].insert(13);
        adjList[5].insert(14);
        
        adjList[6].insert(15);
        adjList[6].insert(16);
        adjList[6].insert(17);
        
        adjList[7].insert(18);
        adjList[7].insert(19);
        adjList[7].insert(20);
        
        BeamSearch b = new BeamSearch(adjList, 0, k);
        TIntList searched = b.search();
        
        assertEquals(21-6, searched.size());
    }
    
    public void testSearch2() {
        int k = 2;
        
        // making cost of middle edges from 5,6,7 higher.
        // same for edge from 2 to 6
        
        int c = 2;
        int cHigh = c*10;
        /*
                        0
            1                      2
          3    4            5         6         7
        8  9  10  11     12 13 14   15 16 17  18 19 20   
        */
        LinkedListCostNode[] adjList = new LinkedListCostNode[21];
        for (int i = 0; i < 21; i++) {
            adjList[i] = new LinkedListCostNode(i);
        }
        adjList[0].insert(1, c);
        adjList[0].insert(2, c);
        
        adjList[1].insert(3, c);
        adjList[1].insert(4, c);
        
        adjList[2].insert(5, c);
        adjList[2].insert(6, cHigh);
        adjList[2].insert(7, c);
        
        adjList[3].insert(8, c);
        adjList[3].insert(9, c);
        
        adjList[4].insert(10, c);
        adjList[4].insert(11, c);
        
        adjList[5].insert(12, c);
        adjList[5].insert(13, cHigh);
        adjList[5].insert(14, c);
        
        adjList[6].insert(15, c);
        adjList[6].insert(16, cHigh);
        adjList[6].insert(17, c);
        
        adjList[7].insert(18, c);
        adjList[7].insert(19, cHigh);
        adjList[7].insert(20, c);
        
        BeamSearch b = new BeamSearch(adjList, 0, k);
        TIntList searched = b.search();
        
        
        assertEquals(21-6, searched.size());
        
        //6,15,16,17,13,19
        int[] notPresent = new int[]{6,15,16,17,13,19};
        
        TIntSet present = new TIntHashSet(searched);
        for (int t : notPresent) {
            assertFalse(present.contains(t));
        }
    }
}
