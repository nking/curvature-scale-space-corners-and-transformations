package algorithms.shortestPath;

import algorithms.shortestPath.JohnsonsAllPairs.G;
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
public class JohnsonsAllPairsTest extends TestCase {
    
    public JohnsonsAllPairsTest() {
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp(); 
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    private G createDataset() {
                /*
        from Cormen at al. Fig 25.6
        
             edge     weight
             s0 t1    3   
             s0 x2    8
             s0 y3    -4
             t1 y3    7    
             t1 z4    1
             x2 t1    4   
             y3 z4    6    
             z4 s0    2 
             z4 x2    -5
        
           s  t  x  y  z
           0  1  2  3  4
        */
        SimpleLinkedListNode[] g = new SimpleLinkedListNode[5];
         for (int i = 0; i < g.length; ++i) {
             g[i] = new SimpleLinkedListNode();
         }
         g[0].insert(1); g[0].insert(2); g[0].insert(3);
         g[1].insert(3); g[1].insert(4);
         g[2].insert(1);
         g[3].insert(4);
         g[4].insert(0);g[4].insert(2);
         TIntIntMap[] w = new TIntIntMap[5];
         w[0] = new TIntIntHashMap(); w[0].put(1,3); w[0].put(2,8); w[0].put(3,-4);
         w[1] = new TIntIntHashMap(); w[1].put(3,7); w[1].put(4,1); 
         w[2] = new TIntIntHashMap(); w[2].put(1,4);
         w[3] = new TIntIntHashMap(); w[3].put(4,6);
         w[4] = new TIntIntHashMap(); w[4].put(0,2); w[4].put(2,-5);

         G graph = new G();
         graph.g = g;
         graph.w = w;
         graph.newNode = -1;
         
         return graph;
    }
    
    public void testInit() {
        
        G graph = createDataset();
        
        SimpleLinkedListNode[] g = graph.g;
        TIntIntMap[] w = graph.w;
        
        assertEquals(5, g.length);
        assertEquals(g.length, w.length);
        
        for (int u = 0; u < g.length; ++u) {            
            SimpleLinkedListNode vNode = g[u];
            TIntIntMap wUVMap = w[u];
            while (vNode != null && vNode.getKey() != -1) {
                int v = vNode.getKey();
                assertTrue(wUVMap.containsKey(v));
                
                vNode = vNode.getNext();
            }
        }
    }
    
    public void test0() {
        
        G graph = createDataset();
                  
         JohnsonsAllPairs j = new JohnsonsAllPairs(graph.g, graph.w);
         boolean noNegativeCycles = j.find();
         assertTrue(noNegativeCycles);
         
         int[] p;
         int dist, dest, src, expDist;
         // s  t  x  y  z
         // 0  1  2  3  4
         
         src = 0;
         dest = 3;
         expDist = -4;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{0, 3}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         src = 0;
         dest = 1;
         expDist = 1;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{0, 3, 4, 2, 1}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         src = 1;
         dest = 2;
         expDist = -4;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{1, 4, 2}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         src = 1;
         dest = 3;
         expDist = -1;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{1, 4, 0, 3}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         src = 2;
         dest = 3;
         expDist = 3;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{2, 1, 4, 0, 3}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         src = 4;
         dest = 1;
         expDist = -1;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{4, 2, 1}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         src = 3;
         dest = 1;
         expDist = 5;
         p = j.getShortestPathToVertex(src, dest);
         assertTrue(Arrays.equals(new int[]{3, 4, 2, 1}, p));
         dist = j.getSumOfPath(p);
         assertEquals(expDist, dist);
         
         // s  t  x  y  z
         // 0  1  2  3  4
         
    }
}
