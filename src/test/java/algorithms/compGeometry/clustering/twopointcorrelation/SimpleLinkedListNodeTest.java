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
        assertTrue(node.key == key2);
        assertTrue(node.next.key == key);
        assertTrue(node.next.next == null);
        
        assertNull(node.search(9876));
        
        boolean threwException = false;
        try {
            node.insert(-1);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            node.insertIfDoesNotAlreadyExist(-1);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);

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
        assertTrue(node.key == key2);
        assertTrue(node.next.key == key);
        assertTrue(node.next.next == null);

        node.insertIfDoesNotAlreadyExist(key);
        assertTrue(node.key == key2);
        assertTrue(node.next.key == key);
        assertTrue(node.next.next == null);
        
        int key3 = 12;
        
        SimpleLinkedListNode r = node.insertIfDoesNotAlreadyExist(key3);
        assertNotNull(r);
        
        r = node.insertIfDoesNotAlreadyExist(key3);
        assertNull(r);
        
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
        
        node = new SimpleLinkedListNode();
        node.delete(new SimpleLinkedListNode(2));
        
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.delete(key);
        assertTrue(node.key == -1);
        assertTrue(node.next == null);
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.delete(new SimpleLinkedListNode(key));
        assertTrue(node.key == -1);
        assertTrue(node.next == null);
        
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.insert(key2);
        node.delete(key2);
        assertTrue(node.key == key);
        assertTrue(node.next == null);
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.insert(key2);
        node.delete(new SimpleLinkedListNode(key2));
        assertTrue(node.key == key);
        assertTrue(node.next == null);
        
        
        int key3 = 12;
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.insert(key2);
        node.insert(key3);
        node.delete(key2);
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key);
        node.delete(-1);
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key);
        node.delete(key2);
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key);
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.insert(key2);
        node.insert(key3);
        node.delete(new SimpleLinkedListNode(key2));
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key);
        node.delete(-1);
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key);
        node.delete(new SimpleLinkedListNode(key2));
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key);
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.insert(key2);
        node.insert(key3);
        node.delete(new SimpleLinkedListNode(key));
        assertTrue(node.key == key3);
        assertTrue(node.next.key == key2);
        
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
        
        
        node = new SimpleLinkedListNode();
        node.insert(key);
        node.insert(key2);
        assertTrue(node.contains(key));
    }
   
    public void testGetKeys() {

        log.info("testGetKeys");

        int key = 1234;
        int key2 = 5678;

        SimpleLinkedListNode node = new SimpleLinkedListNode();

        node.insert(key);
        node.insert(key2);
        
        assertTrue(node.getNumberOfKeys() == 2);
        
        int[] keys = node.getKeys();
        assertTrue(keys.length == 2);
        assertTrue(keys[0] == key2);
        assertTrue(keys[1] == key);
        
        
        node = new SimpleLinkedListNode();
        keys = node.getKeys();
        assertTrue(keys.length == 0);
        assertTrue(node.getNumberOfKeys() == 0);
        
        assertFalse(node.equals(new Object()));
        
        node = new SimpleLinkedListNode();
        for (int i = 0; i < 12; i++) {
            node.insert(i);
        }
        keys = node.getKeys();
        assertTrue(keys.length == 12);
    }
    
}
