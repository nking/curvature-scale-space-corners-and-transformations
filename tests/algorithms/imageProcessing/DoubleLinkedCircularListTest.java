package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DoubleLinkedCircularListTest extends TestCase {

    public DoubleLinkedCircularListTest(String testName) {
        super(testName);
    }

    public void testInsert() {

        System.out.println("insert");

        int key = 1234;
        int key2 = 5678;
        
        String data1 = "data1";
        String data2 = "data2";

        HeapNode node = new HeapNode(key);
        node.setData(data1);
        
        HeapNode node2 = new HeapNode(key2);
        node2.setData(data2);

        DoubleLinkedCircularList list = new DoubleLinkedCircularList();

        list.insert(node);
        assertTrue(list.getSentinel().getRight().getKey() == key);
        assertTrue(list.getSentinel().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(list.getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
        
        assertEquals(list.getSentinel().getRight().getData(), data1);
        
        list.insert(node2);
        assertTrue(list.getSentinel().getRight().getKey() == key2);
        assertTrue(list.getSentinel().getRight().getRight().getKey() == key);
        assertTrue(list.getSentinel().getRight().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(list.getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
        
        assertEquals(list.getSentinel().getRight().getData(), data2);
    }

    /**
     * Test of delete method, of class DoubleLinkedCircularList.
     */
    public void testDelete_HeapNode() {

        System.out.println("delete");

        int key = 1234;
        int key2 = 5678;

        String data1 = "data1";
        String data2 = "data2";
        
        HeapNode node = new HeapNode(key);
        HeapNode node2 = new HeapNode(key2);
        node.setData(data1);
        node2.setData(data2);

        DoubleLinkedCircularList list = new DoubleLinkedCircularList();

        list.insert(node);
        list.insert(node2);

        list.remove(node);

        HeapNode rNode = list.search(key);
        assertNull(rNode);

        HeapNode rNode2 = list.search(key2);
        assertNotNull(rNode2);
        assertTrue(node2.equals(rNode2));
        assertEquals(rNode2.getData(), data2);

        assertTrue(list.getSentinel().getRight().getKey() == key2);
        assertTrue(list.getSentinel().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(list.getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
    }

    /**
     * Test of delete method, of class DoubleLinkedCircularList.
     */
    public void testDelete_int() {
        int key = 1234;
        int key2 = 5678;
        int key3 = 90;

        String data1 = "data1";
        String data2 = "data2";
        String data3 = "data3";
        
        HeapNode node = new HeapNode(key);
        HeapNode node2 = new HeapNode(key2);
        HeapNode node3 = new HeapNode(key3);
        node.setData(data1);
        node2.setData(data2);
        node3.setData(data3);

        DoubleLinkedCircularList list = new DoubleLinkedCircularList();

        list.insert(node);
        list.insert(node2);
        list.insert(node3);

        assertTrue(list.remove(key));
        assertTrue(list.remove(key2));

        HeapNode rNode = list.search(key);
        assertNull(rNode);

        HeapNode rNode3 = list.search(key3);
        assertNotNull(rNode3);
        assertEquals(rNode3.getData(), data3);
        assertTrue(node3.equals(rNode3));

        assertTrue(list.getSentinel().getRight().getKey() == key3);
        assertTrue(list.getSentinel().getRight().getRight().getKey() == DoubleLinkedCircularList.sentinelKey);
        assertTrue(list.getSentinel().getRight().getLeft().getKey() == DoubleLinkedCircularList.sentinelKey);
    }

    public void testSearch() {
        System.out.println("search");

        int key = 1234;
        int key2 = 5678;

        String data1 = "data1";
        String data2 = "data2";
        
        HeapNode node = new HeapNode(key);
        HeapNode node2 = new HeapNode(key2);
        node.setData(data1);
        node2.setData(data2);
        
        DoubleLinkedCircularList list = new DoubleLinkedCircularList();

        list.insert(node);

        list.insert(node2);

        HeapNode rNode = list.search(key);
        HeapNode rNode2 = list.search(key2);

        assertNotNull(rNode);
        assertTrue(node.equals(rNode));
        assertEquals(rNode.getData(), data1);

        assertNotNull(rNode2);
        assertTrue(node2.equals(rNode2));
        assertEquals(rNode2.getData(), data2);
    }
}
