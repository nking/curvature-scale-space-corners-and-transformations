package algorithms.shortestPath;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DAGShortestPathTest extends TestCase {
    
    public DAGShortestPathTest() {
    }
    
     public void test0() {
        
        // Fig 24.5 from Cormen et al. "Intro to Alg"
        /*
         * r=0, s=1, t=2, x=3, y=4, z=5
         *
         *  <0> --5--> <1> --6--> <3>
         *    \         |        >| \
         *      3       2      /      \
         *       \      |    7         \
         *        \    \/  /             \     
         *         --> <2> -- 4--><4>     1
         *                 \                \
         *                       \           \
         *                           2        \/
         *                               \-> <5>
         */
        
         SimpleLinkedListNode[] g = new SimpleLinkedListNode[6];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
         }
         g[0].insert(1); g[0].insert(2);
         g[1].insert(2); g[1].insert(3);
         g[2].insert(3); g[2].insert(4); g[2].insert(5);
         g[3].insert(5);
         TIntIntMap[] w = new TIntIntMap[6];
         w[0] = new TIntIntHashMap(); w[0].put(1, 5); w[0].put(2, 3);
         w[1] = new TIntIntHashMap(); w[1].put(2, 2); w[1].put(3, 6);
         w[2] = new TIntIntHashMap(); w[2].put(3, 7); w[2].put(4, 4); w[2].put(5, 2);
         w[3] = new TIntIntHashMap(); w[3].put(5, 1);
         
         int src = 0;
         int dest = 4;
         DAGShortestPath sp = new DAGShortestPath();
         sp.find(g, w, src);
         
         int[] path = sp.getShortestPathToVertex(dest);
         assertTrue(path.length == 3);
         assertTrue(path[0] == src);
         assertTrue(path[1] == 2);
         assertTrue(path[2] == dest);
         
         dest = 5;
         path = sp.getShortestPathToVertex(dest);
         assertTrue(path.length == 3);
         assertTrue(path[0] == src);
         assertTrue(path[1] == 2);
         assertTrue(path[2] == dest);
     }
     
}
