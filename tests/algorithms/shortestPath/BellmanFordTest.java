package algorithms.shortestPath;

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
public class BellmanFordTest extends TestCase {
    
    public BellmanFordTest() {
    }
    
    public void test0() {
        
        /*
        from Cormen at al. Fig 24.4
        
             edge     weight
             s0 t1    6   
             s0 y3    7
             t1 x2    5
             t1 y3    8
             t1 z4    -4
             x2 t1    -2
             y3 x2    -3    
             y3 z4    9     
             z4 s0    2
             z4 x2    7
        
           s  t  x  y  z
           0  1  2  3  4
        */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[5];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
         }
         g[0].insert(1); g[0].insert(3);
         g[1].insert(2); g[1].insert(3); g[1].insert(4);
         g[2].insert(1);
         g[3].insert(2); g[3].insert(4);
         g[4].insert(0); g[4].insert(2);
         TIntIntMap[] w = new TIntIntMap[5];
         w[0] = new TIntIntHashMap(); w[0].put(1,6); w[0].put(3,7);
         w[1] = new TIntIntHashMap(); w[1].put(2,5); w[1].put(3,8); w[1].put(4,-4);
         w[2] = new TIntIntHashMap(); w[2].put(1,-2);
         w[3] = new TIntIntHashMap(); w[3].put(2,-3); w[3].put(4,9);
         w[4] = new TIntIntHashMap(); w[4].put(0,2); w[4].put(2,7);
         
         int src = 0;
         
         BellmanFord bf = new BellmanFord(g, w, src);
         boolean noNegativeCycles = bf.find();
         assertTrue(noNegativeCycles);
         
         int[] p;
         int dist, dest;
         // s  t  x  y  z
         // 0  1  2  3  4
         
         dest = 1;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 2, 1}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(2, dist);
         
         dest = 2;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 2}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(4, dist);
         
         dest = 3;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(7, dist);
         
         dest = 4;
         p = bf.getShortestPathToVertex(dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 2, 1, 4}, p));
         dist = bf.getSumOfPath(p);
         assertEquals(-2, dist);
    }
}
