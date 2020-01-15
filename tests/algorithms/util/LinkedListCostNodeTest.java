package algorithms.util;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * 
 * @author nichole
 */
public class LinkedListCostNodeTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());
    
    public LinkedListCostNodeTest(String testName) {
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
     * Test of insert method, of class LinkedListCostNode.
     */
    public void testInsert() {
        log.info("insert");

        int key = 1234;
        int key2 = 5678;

        LinkedListCostNode node = new LinkedListCostNode();

        node.insert(key);
        assertTrue(node.key == key);
        assertTrue(node.getCost() == LinkedListCostNode.DEFAULT_COST);
        assertTrue(node.next == null);

        node.insert(key2);
        assertTrue(node.key == key2);
        assertTrue(node.getCost() == LinkedListCostNode.DEFAULT_COST);
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
     * Test of insertIfDoesNotAlreadyExist method, of class LinkedListCostNode.
     */
    public void testInsertIfDoesNotAlreadyExist() {
        log.info("insertIfDoesNotAlreadyExist");

        int key = 1234;
        int key2 = 5678;
        
        int cost = 1;
        int cost2 = 2;
       

        LinkedListCostNode node = new LinkedListCostNode();

        node.insertIfDoesNotAlreadyExist(key, cost);
        assertTrue(node.key == key);
        assertTrue(node.cost == cost);
        assertTrue(node.next == null);

        node.insertIfDoesNotAlreadyExist(key2, cost2);
        assertTrue(node.key == key2);
        assertTrue(node.cost == cost2);
        assertTrue(node.next.key == key);
        assertTrue(node.next.next == null);

        node.insertIfDoesNotAlreadyExist(key);
        assertTrue(node.key == key2);
        assertTrue(node.cost == cost2);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.getNext()).getCost() == cost);
        assertTrue(node.next.next == null);
        
        int key3 = 12;
        
        LinkedListCostNode r = (LinkedListCostNode)node.insertIfDoesNotAlreadyExist(key3);
        assertNotNull(r);
        
        assertNull(node.insertIfDoesNotAlreadyExist(key3));
        
    }

    /**
     * Test of delete method, of class SimpleLinkedListNode.
     */
    public void testDelete_LinkedListCostNode() {

        log.info("delete");

        int key = 1234;
        int key2 = 5678;
        
        int cost = 1;
        int cost2 = 2;

        LinkedListCostNode node = new LinkedListCostNode();

        node.insert(key, cost);
        node.insert(key2, cost2);

        node.delete(key);

        assertTrue(node.key == key2);
        assertTrue(node.cost == cost2);
        assertTrue(node.next == null);

        node.delete(key2);

        assertTrue(node.key == -1);
        
        node = new LinkedListCostNode();
        node.delete(new LinkedListCostNode(2, Integer.MAX_VALUE));
        
        
        node = new LinkedListCostNode();
        node.insert(key);
        node.delete(key);
        assertTrue(node.key == -1);
        assertTrue(node.cost == LinkedListCostNode.DEFAULT_COST);
        assertTrue(node.next == null);
        
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.delete(new LinkedListCostNode(key));
        assertTrue(node.key == -1);
        assertTrue(node.cost == LinkedListCostNode.DEFAULT_COST);
        assertTrue(node.next == null);
        
        
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.insert(key2, cost2);
        node.delete(key2);
        assertTrue(node.key == key);
        assertTrue(node.cost == cost);
        assertTrue(node.next == null);
        
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.insert(key2, cost2);
        node.delete(new LinkedListCostNode(key2));
        assertTrue(node.key == key);
        assertTrue(node.cost == cost);
        assertTrue(node.next == null);
        
        int key3 = 12;
        int cost3 = 3;
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.insert(key2, cost2);
        node.insert(key3, cost3);
        node.delete(key2);
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.next).cost == cost);
        node.delete(-1);
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.next).cost == cost);
        node.delete(key2);
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.next).cost == cost);
        
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.insert(key2, cost2);
        node.insert(key3, cost3);
        
        node.delete(new LinkedListCostNode(key2));
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.next).cost == cost);
        
        node.delete(-1);
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.next).cost == cost);
        
        node.delete(new LinkedListCostNode(key2));
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key);
        assertTrue(((LinkedListCostNode)node.next).cost == cost);
        
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.insert(key2, cost2);
        node.insert(key3, cost3);
        node.delete(new LinkedListCostNode(key));
        assertTrue(node.key == key3);
        assertTrue(node.cost == cost3);
        assertTrue(node.next.key == key2);
        assertTrue(((LinkedListCostNode)node.next).cost == cost2);
        
    }

    /**
     * Test of search method, of class LinkedListCostNode.
     */
    public void testSearch() {

        log.info("search");

        int key = 1234;
        int key2 = 5678;
        
        int cost = 1;
        int cost2 = 2;

        LinkedListCostNode node = new LinkedListCostNode();

        node.insert(key, cost);
        node.insert(key2, cost2);

        LinkedListCostNode rNode = (LinkedListCostNode)node.search(key);
        LinkedListCostNode rNode2 = (LinkedListCostNode)node.search(key2);

        assertTrue(rNode.key == key);
        assertTrue(rNode.cost == cost);

        assertTrue(rNode2.key == key2);
        assertTrue(rNode2.cost == cost2);
        
        
        node = new LinkedListCostNode();
        node.insert(key, cost);
        node.insert(key2, cost2);
        assertTrue(node.contains(key));
        assertTrue(node.contains(key2));
    }
   
    public void testGetKeys() {

        log.info("testGetKeys");

        int key = 1234;
        int key2 = 5678;
        
        int cost = 1;
        int cost2 = 2;

        LinkedListCostNode node = new LinkedListCostNode();

        node.insert(key, cost);
        node.insert(key2, cost2);
        
        assertTrue(node.getNumberOfKeys() == 2);
        
        int[] keys = node.getKeys();
        assertTrue(keys.length == 2);
        assertTrue(keys[0] == key2);
        assertTrue(keys[1] == key);
        
        
        node = new LinkedListCostNode();
        keys = node.getKeys();
        assertTrue(keys.length == 0);
        assertTrue(node.getNumberOfKeys() == 0);
        
        assertFalse(node.equals(new Object()));
        
        node = new LinkedListCostNode();
        for (int i = 0; i < 12; i++) {
            node.insert(i);
        }
        keys = node.getKeys();
        assertTrue(keys.length == 12);
    }
        
}
