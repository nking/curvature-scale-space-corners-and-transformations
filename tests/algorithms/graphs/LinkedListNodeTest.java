package algorithms.graphs;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LinkedListNodeTest extends TestCase {

    // constructing tests from MIT open courseware
    // network_optimization/MIT15_082JF10_av03.pdf

    public LinkedListNodeTest(String testName) {
        super(testName);
    }

    // constructing tests from MIT open courseware
    // network_optimization/MIT15_082JF10_av03.pdf
    protected boolean[][] connected = null;

    @Override
    protected void setUp() throws Exception {

        super.setUp();

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
        connected = new boolean[8][];

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
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of insertIncoming method, of class LinkedListNode.
     */
    public void testInsertIncoming() {

        System.out.println("testInsertIncoming");

        LinkedListNode node0 = new LinkedListNode(0);
        node0.insertOutgoing(2);

        assertTrue(node0.key == 0);
        assertTrue(node0.nOutgoing == 1);
        assertTrue(node0.outgoing[0] == 2);

        node0.insertIncoming(5);
        node0.insertIncoming(7);
        assertTrue(node0.nIncoming == 2);
        assertTrue(node0.incoming[0] == 5);
        assertTrue(node0.incoming[1] == 7);


        node0.removeOutgoing(3);
        assertTrue(node0.nOutgoing == 1);
        assertTrue(node0.outgoing[0] == 2);

        node0.removeOutgoing(2);
        assertTrue(node0.nOutgoing == 0);


        node0.removeIncoming(3);
        assertTrue(node0.nIncoming == 2);
        assertTrue(node0.incoming[0] == 5);
        assertTrue(node0.incoming[1] == 7);

        node0.removeIncoming(5);
        assertTrue(node0.nIncoming == 1);
        assertTrue(node0.incoming[0] == 7);

        node0.removeIncoming(7);
        assertTrue(node0.nIncoming == 0);
    }
}
