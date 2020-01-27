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

    boolean testSort_SimpleDAG = true;
    boolean test0 = true;
    boolean testDFS = true;
    boolean testWalk2 = true;
            
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
        if (!test0) {
            return;
        }

        Logger.getLogger(getClass().getSimpleName()).info("test0");

        for (int j = 0; j < 4; ++j) {

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
            
            4
            7210
            3
            65
            8
             */
            int[] dIdxs;
            int[] fIdxs;

            if (j == 0) {
                DFS dfs = new DFS(connected);
                dfs.walk();
                dIdxs = dfs.getOrderedBeginIndexes();
                fIdxs = dfs.getOrderedEndIndexes();
            } else if (j == 1) {

                DFSIterative dfs2 = new DFSIterative();
                dfs2.walk(connected);

                dIdxs = dfs2.getOrderedBeginIndexes();
                fIdxs = dfs2.getOrderedEndIndexes();
            } else if (j == 2) {
                DFSWithIndependentSets dfs3 = new DFSWithIndependentSets();
                dfs3.walk(connected);

                dIdxs = dfs3.getOrderedBeginIndexes();
                fIdxs = dfs3.getOrderedEndIndexes();

                Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs3.getTd()));
                Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs3.getTf()));
                Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
                Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));

                Logger.getLogger(getClass().getSimpleName()).info(dfs3.printIndependentSets());

                assertEquals(2, dfs3.getIndependentSets().size());
            } else if (j == 3) {
                DFSIterativeWithIndependentSets dfs4 = new DFSIterativeWithIndependentSets();
                dfs4.walk(connected);
                assertEquals(2, dfs4.getIndependentSets().size());
            }
        }
    }

    public void testDFS() {
        if (!testDFS) {
            return;
        }
        
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
        for (int j = 0; j < 4; ++j) {
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

            TIntSet expectedNodes0 = new TIntHashSet();
            TIntSet expectedNodes1 = new TIntHashSet();
            expectedNodes0.add(4);
            expectedNodes0.add(1);
            expectedNodes0.add(3);
            expectedNodes0.add(0);
            expectedNodes1.add(5);
            expectedNodes1.add(2);

            if (j == 0) {
                DFSWithIndependentSets dfs3 = new DFSWithIndependentSets();
                dfs3.walk(directedEdges);
                dIdxs = dfs3.getOrderedBeginIndexes();
                fIdxs = dfs3.getOrderedEndIndexes();

                Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs3.getTd()));
                Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs3.getTf()));
                Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
                Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
                Logger.getLogger(getClass().getSimpleName()).info(dfs3.printIndependentSets());

                assertEquals(1, dfs3.getIndependentSets().size());

            } else if (j == 1) {
                DFSIterativeWithIndependentSets dfs4 = new DFSIterativeWithIndependentSets();
                dfs4.walk(directedEdges);
                assertEquals(1, dfs4.getIndependentSets().size());

            } else if (j == 2) {
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

            } else if (j == 3) {
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
            }
        }
    }

    public void testWalk2() {
        if (!testWalk2) {
            return;
        }
        
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
        
        */
        
        for (int j = 0; j < 4; ++j) {
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

            if (j == 0) {
                DFSWithIndependentSets dfs3 = new DFSWithIndependentSets();
                dfs3.walk(connected);
                dIdxs = dfs3.getOrderedBeginIndexes();
                fIdxs = dfs3.getOrderedEndIndexes();

                Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs3.getTd()));
                Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs3.getTf()));
                Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
                Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
                Logger.getLogger(getClass().getSimpleName()).info(dfs3.printIndependentSets());
                
                assertEquals(2, dfs3.getIndependentSets().size());

            } else if (j == 1) {

                DFSIterativeWithIndependentSets dfs4 = new DFSIterativeWithIndependentSets();
                dfs4.walk(connected);
                assertEquals(2, dfs4.getIndependentSets().size());
            } else if (j == 2) {
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

            } else if (j == 3) {

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
            }
        }
    }
    
    public void testSort_SimpleDAG() {
        if (!testSort_SimpleDAG) {
            return;
        }

        Logger.getLogger(getClass().getSimpleName()).info("testSort_SimpleDAG");
        
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
        expected=[6, 4, 3, 5, 1, 0, 7, 2]
        result=  [6, 4, 5, 3, 1, 0, 7, 2]
        6 4 5
        3
        1
        0 7 2
        
        dTimes=[1, 2, 3, 7, 9, 11, 12, 15]
        fTimes=[4, 5, 6, 8, 10, 13, 14, 16]
        dIndexes=[0, 7, 2, 1, 3, 4, 5, 6]
        fIndexes=[2, 7, 0, 1, 3, 5, 4, 6]
        parent=5; set=7, 6, 5, 4, 3, 2, 1, 0, 
        */
        
        for (int i = 0; i < 8; ++i) {
            connected[i] = new SimpleLinkedListNode();
        }
        
        connected[0].insert(7);

        connected[1].insert(0);

        connected[3].insert(2);
        connected[3].insert(7);

        connected[4].insert(1);
        connected[4].insert(5);
        connected[4].insert(3);

        connected[5].insert(0);
        connected[5].insert(1);

        connected[6].insert(3);
        connected[6].insert(4);

        connected[7].insert(2);

           
        DFSWithIndependentSets dfs3 = new DFSWithIndependentSets();
        dfs3.walk(connected);
        int[] dIdxs = dfs3.getOrderedBeginIndexes();
        int[] fIdxs = dfs3.getOrderedEndIndexes();
        
        Logger.getLogger(getClass().getSimpleName()).info("dTimes=" + Arrays.toString(dfs3.getTd()));
        Logger.getLogger(getClass().getSimpleName()).info("fTimes=" + Arrays.toString(dfs3.getTf()));
        Logger.getLogger(getClass().getSimpleName()).info("dIndexes=" + Arrays.toString(dIdxs));
        Logger.getLogger(getClass().getSimpleName()).info("fIndexes=" + Arrays.toString(fIdxs));
        Logger.getLogger(getClass().getSimpleName()).info(dfs3.printIndependentSets());
                
        assertEquals(1, dfs3.getIndependentSets().size());
    
        DFSIterativeWithIndependentSets dfs4 = new DFSIterativeWithIndependentSets();
        dfs4.walk(connected);
        assertEquals(1, dfs4.getIndependentSets().size());
    }
}
