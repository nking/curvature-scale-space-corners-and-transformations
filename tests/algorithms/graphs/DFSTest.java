package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DFSTest extends TestCase {

    public DFSTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testDFS() {

        /*
         *   u --+ v     w
         *   |   + |    /|
         *   |  /  |   / |
         *   + /   | +   | +..
         *   |/    +     +....
         *   x +-- y     z
         *
         *   where u=0, v=1, w=2, x=3, y=4, z=5
         */
        SimpleLinkedListNode[] directedEdges = new SimpleLinkedListNode[6];
        directedEdges[0] = new SimpleLinkedListNode(1); // edge from 0:1
        directedEdges[0].insert(3);                   // edge 0:3
        directedEdges[1] = new SimpleLinkedListNode(4); // edge 1:4
        directedEdges[2] = new SimpleLinkedListNode(4); // edge 2:4
        directedEdges[2].insert(5); // edge 2:5
        directedEdges[3] = new SimpleLinkedListNode(1); // edge 3:1
        directedEdges[4] = new SimpleLinkedListNode(3); // edge 4:3
        directedEdges[5] = new SimpleLinkedListNode(5); // edge 5:5

        DFS dfs = new DFS(directedEdges);
        dfs.walk();
        
        int[] dIdxs = dfs.getOrderedBeginIndexes();
        int[] fIdxs = dfs.getOrderedEndIndexes();
        
        System.out.println("dTimes=" + Arrays.toString(dfs.getTd()));
        System.out.println("fTimes=" + Arrays.toString(dfs.getTf()));
        System.out.println("dIndexes=" + Arrays.toString(dIdxs));
        System.out.println("fIndexes=" + Arrays.toString(fIdxs));
                
        assertTrue(Arrays.equals(new int[]{0, 3, 1, 4, 2, 5}, dIdxs));
        assertTrue(Arrays.equals(new int[]{4, 1, 3, 0, 5, 2}, fIdxs));
        
        System.out.println("ITERATIVE:");
        
        DFSIterative dfs2 = new DFSIterative();
        dfs2.walk(directedEdges);
        
        dIdxs = dfs2.getOrderedBeginIndexes();
        fIdxs = dfs2.getOrderedEndIndexes();
        
        System.out.println("dTimes=" + Arrays.toString(dfs2.getTd()));
        System.out.println("fTimes=" + Arrays.toString(dfs2.getTf()));
        System.out.println("dIndexes=" + Arrays.toString(dIdxs));
        System.out.println("fIndexes=" + Arrays.toString(fIdxs));
                
        assertTrue(Arrays.equals(new int[]{0, 3, 1, 4, 2, 5}, dIdxs));
        assertTrue(Arrays.equals(new int[]{4, 1, 3, 0, 5, 2}, fIdxs));
        
    }

    public void testWalk2() {
        System.out.println("testWalk2");
        
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
        
        // jacket, tie, belt, shirt, watch, shoes, pants, undershorts, socks
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
        
        // jacket, tie, belt, shirt, watch, shoes, pants, undershorts, socks
        int[] expected = new int[] {
            objs.get("jacket"), objs.get("tie"), objs.get("belt"),
            objs.get("shirt"), objs.get("watch"), objs.get("shoes"),
            objs.get("pants"), objs.get("undershorts"), objs.get("socks"),
        };
        
        int[] dIdxs;
        int[] fIdxs;

        System.out.println("recursive:");
        DFS dfs = new DFS(connected);
        dfs.walk();
        
        dIdxs = dfs.getOrderedBeginIndexes();
        fIdxs = dfs.getOrderedEndIndexes();
        
        System.out.println("dTimes=" + Arrays.toString(dfs.getTd()));
        System.out.println("fTimes=" + Arrays.toString(dfs.getTf()));
        System.out.println("dIndexes=" + Arrays.toString(dIdxs));
        System.out.println("fIndexes=" + Arrays.toString(fIdxs));
        
        assertTrue(Arrays.equals(expected, fIdxs));
        
        System.out.println("iterative:");
        DFSIterative dfs2 = new DFSIterative();
        dfs2.walk(connected);
        
        dIdxs = dfs2.getOrderedBeginIndexes();
        fIdxs = dfs2.getOrderedEndIndexes();
        
        System.out.println("dTimes=" + Arrays.toString(dfs2.getTd()));
        System.out.println("fTimes=" + Arrays.toString(dfs2.getTf()));
        System.out.println("dIndexes=" + Arrays.toString(dIdxs));
        System.out.println("fIndexes=" + Arrays.toString(fIdxs));
                
        assertTrue(Arrays.equals(expected, fIdxs));
    }
}
