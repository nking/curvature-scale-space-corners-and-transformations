package algorithms.search;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

/**
 *
 * @author nichole
 */
public class UniformCostSearchTest extends TestCase {
    
    public UniformCostSearchTest() {
    }
 
    public void testGreedy() {
        /*
           edges   weight
             s0  a1    1
             s0  d4    5
             a1  b2    3
             b2  c3    1
             b2  g5    1
             c3  g5    1
             d4  g5    1
           s  a  b  c  d  g
           0  1  2  3  4  5
        */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[6];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
         }
         g[0].insert(1); g[0].insert(4);
         g[1].insert(2);
         g[2].insert(3); g[2].insert(5);
         g[3].insert(5);
         g[4].insert(5); 
         TIntIntMap[] w = new TIntIntMap[6];
         w[0] = new TIntIntHashMap(); w[0].put(1, 1); w[0].put(4, 5);
         w[1] = new TIntIntHashMap(); w[1].put(2, 3); 
         w[2] = new TIntIntHashMap(); w[2].put(3, 1); w[2].put(5, 1);
         w[3] = new TIntIntHashMap(); w[3].put(5, 1);
         w[4] = new TIntIntHashMap(); w[4].put(5, 1);
         
         int dest = 5;
         UniformCostSearch srch = new UniformCostSearch(g, w, 0, dest);
         srch.find();
         
         int[] p;
         int dist;
         // s0  a1  b2  g5
         
         p = srch.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 1, 2, 5}, p));
         dist = srch.getSumOfPath(p);
         assertEquals(5, dist);
    }
    
    public void testOptimal() {
        /*
           edges   weight
             s0  a1    1
             s0  d4    2
             a1  b2    3
             b2  c3    1
             b2  g5    1
             c3  g5    1
             d4  g5    1
           s  a  b  c  d  g
           0  1  2  3  4  5
        */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[6];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
         }
         g[0].insert(1); g[0].insert(4);
         g[1].insert(2);
         g[2].insert(3); g[2].insert(5);
         g[3].insert(5);
         g[4].insert(5); 
         TIntIntMap[] w = new TIntIntMap[6];
         w[0] = new TIntIntHashMap(); w[0].put(1, 1); w[0].put(4, 2);
         w[1] = new TIntIntHashMap(); w[1].put(2, 3); 
         w[2] = new TIntIntHashMap(); w[2].put(3, 1); w[2].put(5, 1);
         w[3] = new TIntIntHashMap(); w[3].put(5, 1);
         w[4] = new TIntIntHashMap(); w[4].put(5, 1);
         
         int dest = 5;
         UniformCostSearch srch = new UniformCostSearch(g, w, 0, dest);
         srch.find();
         
         int[] p;
         int dist;
         
         p = srch.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 4, 5}, p));
         dist = srch.getSumOfPath(p);
         assertEquals(3, dist);
    }
}
