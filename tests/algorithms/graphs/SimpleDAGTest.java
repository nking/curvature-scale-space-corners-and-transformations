package algorithms.graphs;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SimpleDAGTest extends TestCase {

    public SimpleDAGTest(String testName) {
        super(testName);
    }

    // constructing tests from MIT open courseware
    // network_optimization/MIT15_082JF10_av03.pdf

    @Override
    protected void setUp() throws Exception {

        super.setUp();

    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    private boolean[][] getAdjacencyList1() {
        
        /*              <5> --------> <0>
         *            >  |       >     |
         *         .     V   .         V
         *      <4> --> <1>           <7> ---> <2>
         *      >  .                 >        .>
         *     .     .            .      .
         *    .         >      .   .
         *  <6> ------> <3> .
         */
        
        boolean[][] connected = new boolean[8][];

        connected[0] = new boolean[8];
        connected[0][7] = true;

        connected[1] = new boolean[1];
        connected[1][0] = true;

        connected[2] = new boolean[0];

        connected[3] = new boolean[8];
        connected[3][2] = true;
        connected[3][7] = true;

        connected[4] = new boolean[6];
        connected[4][1] = true;
        connected[4][3] = true;
        connected[4][5] = true;

        connected[5] = new boolean[2];
        connected[5][0] = true;
        connected[5][1] = true;

        connected[6] = new boolean[5];
        connected[6][3] = true;
        connected[6][4] = true;

        connected[7] = new boolean[3];
        connected[7][2] = true;
        
        return connected;
    }
    
    private LinkedListNode[] getExpectedResult() {
        
        LinkedListNode[] expResult = new LinkedListNode[8];
        expResult[0] = new LinkedListNode(0);
        expResult[0].insertOutgoing(7);
        expResult[0].insertIncoming(1);
        expResult[0].insertIncoming(5);

        expResult[1] = new LinkedListNode(1);
        expResult[1].insertIncoming(4);
        expResult[1].insertIncoming(5);
        expResult[1].insertOutgoing(0);

        expResult[2] = new LinkedListNode(2);
        expResult[2].insertIncoming(3);
        expResult[2].insertIncoming(7);

        expResult[3] = new LinkedListNode(3);
        expResult[3].insertOutgoing(2);
        expResult[3].insertOutgoing(7);
        expResult[3].insertIncoming(4);
        expResult[3].insertIncoming(6);

        expResult[4] = new LinkedListNode(4);
        expResult[4].insertOutgoing(1);
        expResult[4].insertOutgoing(3);
        expResult[4].insertOutgoing(5);
        expResult[4].insertIncoming(6);

        expResult[5] = new LinkedListNode(5);
        expResult[5].insertOutgoing(0);
        expResult[5].insertOutgoing(1);
        expResult[5].insertIncoming(4);

        expResult[6] = new LinkedListNode(6);
        expResult[6].insertOutgoing(3);
        expResult[6].insertOutgoing(4);

        expResult[7] = new LinkedListNode(7);
        expResult[7].insertOutgoing(2);
        expResult[7].insertIncoming(0);
        expResult[7].insertIncoming(3);
        
        assertExpectedResult(expResult);
        
        return expResult;
    }

    private void assertExpectedResult(LinkedListNode[] expResult) {
        
        assertTrue(expResult[0].nOutgoing == 1);
        assertTrue(expResult[1].nOutgoing == 1);
        assertTrue(expResult[2].nOutgoing == 0);
        assertTrue(expResult[3].nOutgoing == 2);
        assertTrue(expResult[4].nOutgoing == 3);
        assertTrue(expResult[5].nOutgoing == 2);
        assertTrue(expResult[6].nOutgoing == 2);
        assertTrue(expResult[7].nOutgoing == 1);

        assertTrue(expResult[0].outgoing[0] == 7);

        assertTrue(expResult[1].outgoing[0] == 0);

        assertTrue(expResult[3].outgoing[0] == 2);
        assertTrue(expResult[3].outgoing[1] == 7);

        assertTrue(expResult[4].outgoing[0] == 1);
        assertTrue(expResult[4].outgoing[1] == 3);
        assertTrue(expResult[4].outgoing[2] == 5);

        assertTrue(expResult[5].outgoing[0] == 0);
        assertTrue(expResult[5].outgoing[1] == 1);

        assertTrue(expResult[6].outgoing[0] == 3);
        assertTrue(expResult[6].outgoing[1] == 4);

        assertTrue(expResult[7].outgoing[0] == 2);
        
        /*              <5> --------> <0>
         *            >  |       >     |
         *         .     V   .         V
         *      <4> --> <1>           <7> ---> <2>
         *      >  .                 >        .>
         *     .     .            .      .
         *    .         >      .   .
         *  <6> ------> <3> .
         */

        assertTrue(expResult[0].nIncoming == 2);
        assertTrue(expResult[1].nIncoming == 2);
        assertTrue(expResult[2].nIncoming == 2);
        assertTrue(expResult[3].nIncoming == 2);
        assertTrue(expResult[4].nIncoming == 1);
        assertTrue(expResult[5].nIncoming == 1);
        assertTrue(expResult[6].nIncoming == 0);
        assertTrue(expResult[7].nIncoming == 2);

        assertTrue(expResult[0].incoming[0] == 1);
        assertTrue(expResult[0].incoming[1] == 5);

        assertTrue(expResult[1].incoming[0] == 4);
        assertTrue(expResult[1].incoming[1] == 5);

        assertTrue(expResult[2].incoming[0] == 3);
        assertTrue(expResult[2].incoming[1] == 7);

        assertTrue(expResult[3].incoming[0] == 4);
        assertTrue(expResult[3].incoming[1] == 6);

        assertTrue(expResult[4].incoming[0] == 6);

        assertTrue(expResult[5].incoming[0] == 4);

        assertTrue(expResult[7].incoming[0] == 0);
        assertTrue(expResult[7].incoming[1] == 3);
    }

    public void testGetNodes1() {

        System.out.println("testGetNodes1");
        
        boolean[][] connected = getAdjacencyList1();

        SimpleDAG dag = new SimpleDAG(connected);

        LinkedListNode[] expResult = getExpectedResult();
        
        /*              <5> --------> <0>
         *            >  |       >     |
         *         .     V   .         V
         *      <4> --> <1>           <7> ---> <2>
         *      >  .                 >        .>
         *     .     .            .      .
         *    .         >      .   .
         *  <6> ------> <3> .
         *
         */

        LinkedListNode[] nodes = dag.getNodes();

        assertTrue(nodes.length == expResult.length);

        for (int i = 0; i < nodes.length; i++) {

            LinkedListNode e = expResult[i];
            LinkedListNode r = nodes[i];

            assertTrue(e.key == i);
            assertTrue(r.key == i);

            assertTrue(e.nIncoming == r.nIncoming);
            assertTrue(e.nOutgoing == r.nOutgoing);

            for (int ii = 0; ii < e.nIncoming; ii++) {
                assertTrue(e.incoming[ii] == r.incoming[ii]);
            }

            for (int ii = 0; ii < e.nOutgoing; ii++) {
                assertTrue(e.outgoing[ii] == r.outgoing[ii]);
            }
        }
    }

    /**
     * Test of removeNode method, of class SimpleDAG.
     */
    public void testRemoveNode() {

        System.out.println("removeNode");
        
        boolean[][] connected = getAdjacencyList1();

        SimpleDAG dag = new SimpleDAG(connected);

        for (int i = 0; i < dag.getNodes().length; i++) {

            dag.removeNode(i);

            assertNull(dag.getNodes()[i]);
        }

        assertTrue(dag.isEmpty());
    }
    
}
