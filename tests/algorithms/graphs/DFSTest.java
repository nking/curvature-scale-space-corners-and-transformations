package algorithms.graphs;

import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
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
    
    public void test0() {
        Logger.getLogger(getClass().getSimpleName()).info("test0");
        
        // constructing test from Cormen et al.'s "Introduction to Algorithms"
        SimpleLinkedListNode[] connected = new SimpleLinkedListNode[9];
        
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
        for (int i = 0; i < connected.length; i++) {
            connected[i] = new SimpleLinkedListNode();
        }
        connected[0].insert(1);
        connected[0].insert(4);
        connected[1].insert(2);
        connected[1].insert(4);
        connected[2].insert(7);
        connected[3].insert(4);
        connected[5].insert(2);
        connected[5].insert(6);
        connected[6].insert(7);
        
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

        Book solution for indexes of tf (==fIndexes) reversed:
        3   0   1   4    8    5    2    6    7
        ------------>         --------->----->
            --->---->         ---->---------->
            -------->
           
        result:
            dIndexes=[0, 4, 1, 2, 7, 3, 5, 6, 8]
            fIndexes=[4, 7, 2, 1, 0, 3, 6, 5, 8]
                prev=[-1, 0, 1, -1, 0, -1, 5, 2, -1]
        NOTE: if were tracking prev of prev etc to top-most
              parent, the value at prev[7] would become 0 instead of 2
              and prev[2] value would become 0 instead of 1
        */
        
        //DFS dfs = new DFS(connected);
        //dfs.walk();
        //DFSIterative dfs = new DFSIterative();
        //dfs.walk(connected);
        DFSWithIndependentSets dfs = new DFSWithIndependentSets();
        dfs.walk(connected);
        
        int[] dIdxs = dfs.getOrderedBeginIndexes();
        int[] fIdxs = dfs.getOrderedEndIndexes();
        int[] prev = dfs.getPredecessorIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("prev=" + Arrays.toString(prev));
        
        Logger.getLogger(getClass().getSimpleName()).info(dfs.printIndependentSetsInTF());
    }

    public void testDFS() {
        
        Logger.getLogger(getClass().getSimpleName()).info("testDFS");
        
        /*   0     1     2
         *   u --+ v     w
         *   |   + |    /|
         *   |  /  |   / |
         *   + /   | +   | +..
         *   |/    +     +....
         *   x +-- y     z
         *   3     4     5
         *   where u=0, v=1, w=2, x=3, y=4, z=5
        
            4, 1, 3, 0, 
            5, 2
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

        int[] dIdxs;
        int[] fIdxs;
        
        DFS dfs = new DFS(directedEdges);
        dfs.walk();
        
        dIdxs = dfs.getOrderedBeginIndexes();
        fIdxs = dfs.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
                
        assertTrue(Arrays.equals(new int[]{0, 3, 1, 4, 2, 5}, dIdxs));
        assertTrue(Arrays.equals(new int[]{4, 1, 3, 0, 5, 2}, fIdxs));
        
        Logger.getLogger(getClass().getSimpleName()).info("ITERATIVE:");
        
        DFSIterative dfs2 = new DFSIterative();
        dfs2.walk(directedEdges);
        
        dIdxs = dfs2.getOrderedBeginIndexes();
        fIdxs = dfs2.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs2.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs2.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
                
        assertTrue(Arrays.equals(new int[]{0, 3, 1, 4, 2, 5}, dIdxs));
        assertTrue(Arrays.equals(new int[]{4, 1, 3, 0, 5, 2}, fIdxs));
        
        /*
        4, 1, 3, 0, 
            5, 2
        */
        TIntSet expectedNodes0 = new TIntHashSet();
        TIntSet expectedNodes1 = new TIntHashSet();
        expectedNodes0.add(4); expectedNodes0.add(1); expectedNodes0.add(3); expectedNodes0.add(0);
        expectedNodes1.add(5); expectedNodes1.add(2);
        
        DFSWithIndependentSets dfs3 = new DFSWithIndependentSets();
        dfs3.walk(directedEdges);
        dIdxs = dfs3.getOrderedBeginIndexes();
        fIdxs = dfs3.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs3.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs3.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
        Logger.getLogger(getClass().getSimpleName()).info(dfs3.printIndependentSetsInTF());
        
    }

    public void testWalk2() {
        Logger.getLogger(getClass().getSimpleName()).info("testWalk2");
        
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
             3       2   1     0      4       7      6        5         8
        
        dTimes=[1, 2, 3, 6, 9, 11, 12, 14, 17]
        fTimes=[4, 5, 7, 8, 10, 13, 15, 16, 18]
        dIndexes=[0, 2, 3, 1, 4, 5, 7, 6, 8]
        fIndexes=[3, 2, 1, 0, 4, 7, 6, 5, 8]
        parent=8; set=8,
        parent=5; set=7, 6, 5,
        parent=4; set=4,
        parent=0; set=3, 2, 1, 0,
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
        connected[objs.get("pants")].insert(objs.get("shoes")); connected[objs.get("pants")].insert(objs.get("belt"));
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

        Logger.getLogger(getClass().getSimpleName()).info("recursive:");
        DFS dfs = new DFS(connected);
        dfs.walk();
        
        dIdxs = dfs.getOrderedBeginIndexes();
        fIdxs = dfs.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
        
        assertTrue(Arrays.equals(expected, fIdxs));
        
        Logger.getLogger(getClass().getSimpleName()).info("iterative:");
        DFSIterative dfs2 = new DFSIterative();
        dfs2.walk(connected);
        
        dIdxs = dfs2.getOrderedBeginIndexes();
        fIdxs = dfs2.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs2.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs2.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
                
        assertTrue(Arrays.equals(expected, fIdxs));
        
        DFSWithIndependentSets dfs3 = new DFSWithIndependentSets();
        dfs3.walk(connected);
        dIdxs = dfs3.getOrderedBeginIndexes();
        fIdxs = dfs3.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs3.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs3.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
        Logger.getLogger(getClass().getSimpleName()).info(dfs3.printIndependentSetsInTF());
    }
}
