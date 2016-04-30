package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
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

        /*
        
        sentinel -> 1234 -> [ sentinel.right ]
        sentinel <- 1234 <- [ sentinel.right ]
        
        sentinel -> 5678 -> 1234 -> [ sentinel.right ]
        sentinel <- 5678 <- 1234 <- [ sentinel.right ]
        
        HeapNode rightOfSentinel = sentinel.getRight();
        node.setRight(rightOfSentinel);
        rightOfSentinel.setLeft(node);
        sentinel.setRight(node);
        node.setLeft(sentinel);
        number++;
        */
        
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
    
    public void testInsertAfter() {
        
        /*
        sentinel -> 2nd inserted -> 1st inserted -> [ sentinel.right ]
        sentinel <- 2nd inserted <- 1st inserted <- [ sentinel.right ]
        *
        * subsequent traversal by FIFO should use :
        *    sentinel.getRight().getLeft() and proceed left for n=number items
        * 
        * subsequent traversal by LIFO should use :
        *    sentinel.getRight() and proceed right for n=number items
        */
        
        DoubleLinkedCircularList dlcl = new DoubleLinkedCircularList();
        
        assertEquals(0, dlcl.getNumberOfNodes());
        
        List<HeapNode> nodes = new ArrayList<HeapNode>();
        List<Long> expectedFIFOKeys = new ArrayList<Long>();
        for (int i = 0; i < 10; i += 2) {
            HeapNode node = new HeapNode(i);
            dlcl.insert(node);
            nodes.add(node);
            expectedFIFOKeys.add(Long.valueOf(i));
        }
        
        assertEquals(nodes.size(), dlcl.getNumberOfNodes());
        assertEquals(expectedFIFOKeys.size(), dlcl.getNumberOfNodes());
        
        // check LIFO traversal has expected keys
        
        HeapNode node = dlcl.getSentinel();
        for (int i = 0; i < dlcl.getNumberOfNodes(); ++i) {
            node = node.getLeft();
            assertEquals(expectedFIFOKeys.get(i).longValue(), node.getKey());
        }

        HeapNode insertNode = new HeapNode(1);
        dlcl.insertAfter(nodes.get(0), insertNode);
        expectedFIFOKeys.add(1, Long.valueOf(1));
        
        assertEquals(expectedFIFOKeys.size(), dlcl.getNumberOfNodes());
        node = dlcl.getSentinel();
        for (int i = 0; i < dlcl.getNumberOfNodes(); ++i) {
            node = node.getLeft();
            assertEquals(expectedFIFOKeys.get(i).longValue(), node.getKey());
        }
        
    }
    
    public void testInsertRemoveGetNumberOfNodes() {
        
        Set<HeapNode> nodes = new HashSet<HeapNode>();
        
        DoubleLinkedCircularList linkedList = new DoubleLinkedCircularList();
        
        for (int i = 0; i < 100000; ++i) {
            HeapNode node = new HeapNode(i);
            linkedList.insert(node);
            nodes.add(node);
        }
        
        int nRemaining = nodes.size();
        for (HeapNode node : nodes) {
            
            assertEquals(nRemaining, linkedList.getNumberOfNodes());
            
            linkedList.remove(node);
            
            nRemaining--;
        }
        
        assertEquals(0, nRemaining);
        
        assertEquals(nRemaining, linkedList.getNumberOfNodes());
    }
    
    public void testInsertRemoveGetNumberOfNodes2() {
        
        Set<HeapNode> nodes = new HashSet<HeapNode>();
        
        DoubleLinkedCircularList linkedList = new DoubleLinkedCircularList();
        
        for (int i = 0; i < 1000; ++i) {
            HeapNode node = new HeapNode(i);
            linkedList.insert(node);
            nodes.add(node);            
        }
        
        Set<HeapNode> removed = new HashSet<HeapNode>();
        int nRemaining = nodes.size();
        for (HeapNode node : nodes) {
            
            assertEquals(nRemaining, linkedList.getNumberOfNodes());
            
            linkedList.remove(node);
            
            removed.add(node);
                        
            nRemaining--;
            
            if (nRemaining == (nodes.size()/2)) {
                break;
            }
        }
        nodes.removeAll(removed);        
        removed.clear();
        
        Random r = new Random();
        for (int i = 0; i < 10; ++i) {
            
            HeapNode node = new HeapNode(i + 1000000);
            
            if ((i & 1) == 1) {
                
                // choose an existing node from nodes
                int nIter = r.nextInt(100) + 1;
                int n = 0;
                Iterator<HeapNode> iter = nodes.iterator();
                HeapNode existingNode = null;
                while (n < nIter) {
                    existingNode = iter.next();
                    n++;
                }

                linkedList.insertAfter(existingNode, node);
                
            } else {
                
                linkedList.insert(node);
            }
            
            nodes.add(node);            
        }

        nRemaining = nodes.size();
        for (HeapNode node : nodes) {
            
            assertEquals(nRemaining, linkedList.getNumberOfNodes());
            
            linkedList.remove(node);
            
            removed.add(node);
                        
            nRemaining--;
        }
              
        assertEquals(0, nRemaining);
        
        assertEquals(nRemaining, linkedList.getNumberOfNodes());
    }

}
