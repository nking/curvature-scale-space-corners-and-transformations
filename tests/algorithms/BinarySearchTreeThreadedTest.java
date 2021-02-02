package algorithms;

import algorithms.heapsAndPQs.HeapNode;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

/**
 *
 * @author nichole
 */
public class BinarySearchTreeThreadedTest extends TestCase {
    
    public BinarySearchTreeThreadedTest() {
    }
    
    //http://adtinfo.org/libavl.html/Inserting-into-a-TBST.html 
    public void testInsert() {
        
        BinarySearchTreeThreaded<HeapNode> bst 
            = new BinarySearchTreeThreaded<HeapNode>();
        
        HeapNode node6 = new HeapNode();
        node6.setKey(6);
        HeapNode node4 = new HeapNode();
        node4.setKey(4);
        HeapNode node5 = new HeapNode();
        node5.setKey(5);
        HeapNode node2 = new HeapNode();
        node2.setKey(2);
        HeapNode node1 = new HeapNode();
        node1.setKey(1);
        HeapNode node3 = new HeapNode();
        node3.setKey(3);
        
        bst.insert(node6);
        bst.insert(node4);
        bst.insert(node5);
        bst.insert(node2);
        bst.insert(node1);
        
        assertEquals(6, bst.root.getKey());
        assertEquals(4, bst.root.getLeft().getKey());
        assertEquals(2, bst.root.getLeft()
            .getLeft().getKey());
        
        assertEquals(4, bst.root.getLeft()
            .getLeft().getRight().getKey());
        
        bst.insert(node3);
        
        assertEquals(3, bst.root.getLeft()
            .getLeft().getRight().getKey());
    }
    
    public void test0() throws NoSuchAlgorithmException {

        /*
        if use n = 6, here is the bst structure:
        
              r 0
       0                 1
                     0            2
                              1            5
                                      4        null
                                   3     5
                                 2   4
                               1  5

        */

        BinarySearchTreeThreaded<HeapNode> bst 
            = new BinarySearchTreeThreaded<HeapNode>();
        
        int n = 6;
        
        HeapNode[] nodes = new HeapNode[2*n];
        
        for (int i = 0; i < n/2; ++i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
            assertNotNull(bst.search(node));
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
            assertNotNull(bst.search(node));
        }
        
        assertEquals(n, bst.getNumberOfNodes());
        
        assertEquals(n - 1, bst.maximum().getKey());
        
        assertEquals(0L, bst.minimum().getKey());
        
        for (int i = 0; i < n; ++i) {
            HeapNode node = bst.search(nodes[i]);
            assertEquals((long)i, node.getKey());
            assertTrue(bst.threadMap.containsKey(node));
            assertNotNull(bst.search(node));
            
            /*
            smallest element in the tree with key greater
            than node.key.            
            */
            if (i < (n - 1)) {
                HeapNode next = bst.successor(node);
                assertEquals((long)(i + 1), next.getKey());
            }
            
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            if (i > 0) {
                HeapNode prev = bst.predecessor(node);
                assertEquals((long)(i - 1), prev.getKey());
            }
            
        }
        
        // remove some nodes randomly
        Set<HeapNode> rm = new HashSet<HeapNode>();
        Set<Integer> rmInt = new HashSet<Integer>();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1465940667831L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        rm.add(nodes[0]);
        bst.delete(nodes[0]);
        HeapNode nod = bst.search(nodes[0]);
        assertNull(nod);
        assertFalse(bst.threadMap.containsKey(nodes[0]));
        rmInt.add(Integer.valueOf(0));
        
        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            HeapNode r = nodes[idx];
            if (!rm.contains(r)) {
                assertTrue(bst.threadMap.containsKey(r));
                bst.delete(r);
                rm.add(r);
                rmInt.add(Integer.valueOf((int)r.getKey()));
                assertFalse(bst.threadMap.containsKey(r));
                assertNull(bst.search(r));
            }
        }
        
        assertEquals((n - rm.size()), bst.getNumberOfNodes());
        
        boolean minChecked = false;
        
        for (int i = 0; i < n; ++i) {
            HeapNode node0 = nodes[i];
            HeapNode node = bst.search(node0);
            if (rm.contains(node0)) {
                assertNull(node);
            } else {
                if (!minChecked) {
                    HeapNode minNode = bst.minimum();
                    assertEquals((long)i, minNode.getKey());
                    minChecked = true;
                }
                assertEquals((long)i, node.getKey());            
                if (i < (n - 1)) {
                    HeapNode next = bst.successor(node);
                    int expected = i + 1;
                    while (expected < n) {
                        if (rmInt.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals((long)expected, next.getKey());
                            break;
                        }
                    }
                }
            }
        }
            
        // ==== then add n more nodes and repeat assertions
        for (int i = n; i < 2*n; ++i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
            assertNotNull(bst.search(node));
        }

        assertEquals((n - rm.size()) + n, bst.getNumberOfNodes());

        Integer maxExpected = Integer.valueOf(2*n - 1);
        while (rmInt.contains(maxExpected)) {
            maxExpected = Integer.valueOf(maxExpected.intValue() - 1);
        }
        
        assertEquals(maxExpected.intValue(), 
            bst.maximum().getKey());

        for (int i = n; i < 2*n; ++i) {
            HeapNode node = bst.search(nodes[i]);
            assertEquals((long)i, node.getKey());
            if (i < (n - 1)) {
                HeapNode next = bst.successor(node);
                assertEquals((long)(i + 1), next.getKey());
            }
        }

        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            HeapNode r = nodes[idx];
            if (!rm.contains(r)) {
                rm.add(r);
                bst.delete(r);
                assertNull(bst.search(r));
            }
        }

        assertEquals((2*n - rm.size()), bst.getNumberOfNodes());

        minChecked = false;

        for (int i = 0; i < 2*n; ++i) {
            HeapNode node0 = nodes[i];
            HeapNode node = bst.search(node0);
            if (rm.contains(node0)) {
                assertNull(node);
            } else {
                if (!minChecked) {
                    HeapNode minNode = bst.minimum();
                    assertEquals((long)i, minNode.getKey());
                    minChecked = true;
                }                
                assertEquals((long)i, node.getKey());

                if (i < (n - 1)) {
                    HeapNode next = bst.successor(node);
                    int expected = i + 1;
                    while (expected < n) {
                        if (rmInt.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals((long)expected, next.getKey());
                            break;
                        }
                    }
                }
            }
        }
    }

}
