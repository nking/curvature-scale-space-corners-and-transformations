package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PointLinkedListNodeTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public PointLinkedListNodeTest(String testName) {
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

    /**
     * Test of insert method, of class PointLinkedListNode.
     */
    public void testInsert() {
        log.info("insert");

        int key = 1234;
        int key2 = 5678;

        PointLinkedListNode node = new PointLinkedListNode();

        node.insert(key, key, key);
        assertTrue(node.secondPointIndex == key);
        assertTrue(node.next == null);

        node.insert(key2, key2, key2);
        assertTrue(node.secondPointIndex == key);
        assertTrue(node.next.secondPointIndex == key2);
        assertTrue(node.next.next == null);

    }

    /**
     * Test of insertIfDoesNotAlreadyExist method, of class PointLinkedListNode.
     */
    public void testInsertIfDoesNotAlreadyExist() {
        log.info("insertIfDoesNotAlreadyExist");

        int key = 1234;
        int key2 = 5678;

        PointLinkedListNode node = new PointLinkedListNode();

        node.insertIfDoesNotAlreadyExist(key, key, key);
        assertTrue(node.secondPointIndex == key);
        assertTrue(node.next == null);

        node.insertIfDoesNotAlreadyExist(key2, key2, key2);
        assertTrue(node.secondPointIndex == key);
        assertTrue(node.next.secondPointIndex == key2);
        assertTrue(node.next.next == null);

        node.insertIfDoesNotAlreadyExist(key, key, key);
        assertTrue(node.secondPointIndex == key);
        assertTrue(node.next.secondPointIndex == key2);
        assertTrue(node.next.next == null);
    }

    /**
     * Test of delete method, of class PointLinkedListNode.
     */
    public void testDelete_PointLinkedListNode() {

        log.info("delete");

        int key = 1234;
        int key2 = 5678;

        PointLinkedListNode node = new PointLinkedListNode();

        node.insert(key, key, key);
        node.insert(key2, key2, key2);

        node.delete(key);

        assertTrue(node.secondPointIndex == key2);
        assertTrue(node.next == null);

        node.delete(key2);

        assertTrue(node.secondPointIndex == -1);
    }

    /**
     * Test of search method, of class PointLinkedListNode.
     */
    public void testSearch() {

        log.info("search");

        int key = 1234;
        int key2 = 5678;

        PointLinkedListNode node = new PointLinkedListNode();

        node.insert(key, key, key);
        node.insert(key2, key2, key2);

        PointLinkedListNode rNode = node.search(key);
        PointLinkedListNode rNode2 = node.search(key2);

        assertNotNull(rNode);
        assertTrue(rNode.secondPointIndex == key);

        assertNotNull(rNode2);
        assertTrue(rNode2.secondPointIndex == key2);
    }
}
