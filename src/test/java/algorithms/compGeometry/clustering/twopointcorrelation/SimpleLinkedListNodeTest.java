package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SimpleLinkedListNodeTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());
    
    public SimpleLinkedListNodeTest(String testName) {
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
     * Test of insert method, of class SimpleLinkedListNode.
     */
    public void testInsert() {
        log.info("insert");

        int key = 1234;
        int key2 = 5678;

        SimpleLinkedListNode node = new SimpleLinkedListNode();

        node.insert(key);
        assertTrue(node.key == key);
        assertTrue(node.next == null);

        node.insert(key2);
        assertTrue(node.key == key);
        assertTrue(node.next.key == key2);
        assertTrue(node.next.next == null);

    }

    /**
     * Test of insertIfDoesNotAlreadyExist method, of class SimpleLinkedListNode.
     */
    public void testInsertIfDoesNotAlreadyExist() {
        log.info("insertIfDoesNotAlreadyExist");

        int key = 1234;
        int key2 = 5678;

        SimpleLinkedListNode node = new SimpleLinkedListNode();

        node.insertIfDoesNotAlreadyExist(key);
        assertTrue(node.key == key);
        assertTrue(node.next == null);

        node.insertIfDoesNotAlreadyExist(key2);
        assertTrue(node.key == key);
        assertTrue(node.next.key == key2);
        assertTrue(node.next.next == null);

        node.insertIfDoesNotAlreadyExist(key);
        assertTrue(node.key == key);
        assertTrue(node.next.key == key2);
        assertTrue(node.next.next == null);
    }

    /**
     * Test of delete method, of class SimpleLinkedListNode.
     */
    public void testDelete_SimpleLinkedListNode() {

        log.info("delete");

        int key = 1234;
        int key2 = 5678;

        SimpleLinkedListNode node = new SimpleLinkedListNode();

        node.insert(key);
        node.insert(key2);

        node.delete(key);

        assertTrue(node.key == key2);
        assertTrue(node.next == null);

        node.delete(key2);

        assertTrue(node.key == -1);
    }

    /**
     * Test of search method, of class SimpleLinkedListNode.
     */
    public void testSearch() {

        log.info("search");

        int key = 1234;
        int key2 = 5678;

        SimpleLinkedListNode node = new SimpleLinkedListNode();

        node.insert(key);
        node.insert(key2);

        SimpleLinkedListNode rNode = node.search(key);
        SimpleLinkedListNode rNode2 = node.search(key2);

        assertNotNull(rNode);
        assertTrue(rNode.key == key);

        assertNotNull(rNode2);
        assertTrue(rNode2.key == key2);
    }
}
