package algorithms;

import algorithms.imageProcessing.HeapNode;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BinarySearchTreeTest extends TestCase {
    
    public BinarySearchTreeTest() {
    }
   
    public void test0() throws Exception {
        
        int n = 100;
        
        BinarySearchTree<HeapNode> bst = new 
            BinarySearchTree<HeapNode>();
        
        HeapNode[] nodes = new HeapNode[2*n];
        
        for (int i = 0; i < n/2; ++i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
        }
        
        assertEquals(n, bst.getNumberOfNodes());
        
        assertEquals(n - 1, bst.maximum().getKey());
        
        assertEquals(0L, bst.minimum().getKey());
        
        for (int i = 0; i < n; ++i) {
            HeapNode node = bst.search(nodes[i]);
            assertEquals((long)i, node.getKey());
            if (i < (n - 1)) {
                HeapNode next = bst.successor(node);
                assertEquals((long)(i + 1), next.getKey());
            }
        }
        
        // remove some nodes randomly
        Set<HeapNode> rm = new HashSet<HeapNode>();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1465900260715L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        rm.add(nodes[0]);
        bst.delete(nodes[0]);
                
        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            HeapNode r = nodes[idx];
            if (!rm.contains(r)) {
                rm.add(r);
                bst.delete(r);
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
                        if (rm.contains(nodes[expected])) {
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
        }
        
        assertEquals((n - rm.size()) + n, bst.getNumberOfNodes());
        
        assertEquals(199L, bst.maximum().getKey());
        
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
                        if (rm.contains(nodes[expected])) {
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
