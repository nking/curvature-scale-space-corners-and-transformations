package algorithms;

import algorithms.imageProcessing.HeapNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class BinarySearchTreeTest extends TestCase {
    
    public BinarySearchTreeTest() {
    }
    
    public void testDelete() {
        
        /*
        {//DEBUG
                    if (r.getKey() == 62) {
                        System.out.println("debug:");
                        HeapNode[] dns = bst.printPreOrderTraversalIterative(bst.root);
                        
                        for (HeapNode dn : dns) {
                            System.out.println("node=" + dn);
                        }
                        int z = 0;
                    }
                }
        */
        
        BinarySearchTree<HeapNode> bst = new 
            BinarySearchTree<HeapNode>();
        
        TIntObjectMap<HeapNode> nodeMap = new TIntObjectHashMap<HeapNode>();
        
        
        int[] ins = new int[]{15, 5, 3, 12, 10, 13, 6, 7, 16, 20, 18, 23};
        for (int key : ins) {
            HeapNode node = new HeapNode();
            node.setKey(key);
            nodeMap.put(key, node);
            bst.insert(node);
            assertNotNull(bst.search(node));
        }
        
        {//DEBUG
            System.out.println("debug:");
            HeapNode[] dns = bst.printPreOrderTraversalIterative(bst.root);
            for (HeapNode dn : dns) {
                System.out.println("node=" + dn);
            }
        }
        
        assertEquals(12, bst.getNumberOfNodes());

        bst.delete(nodeMap.get(5));
        
        assertEquals(11, bst.getNumberOfNodes());
        
        HeapNode r = bst.search(nodeMap.get(5));
        
        /*{//DEBUG
            System.out.println("debug:");
            HeapNode[] dns = bst.printPreOrderTraversalIterative(bst.root);
            for (HeapNode dn : dns) {
                System.out.println("node=" + dn);
            }
        }*/
        
        assertNull(r);
    }

    public void test0() throws Exception {

        int n = 100;

        BinarySearchTree<HeapNode> bst = new BinarySearchTree<HeapNode>();

        HeapNode[] nodes = new HeapNode[2 * n];

        int count = 0;
        for (int i = 0; i < n / 2; ++i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
            assertNotNull(bst.search(node));
            count++;
            assertEquals(count, bst.getNumberOfNodes());
        }
        for (int i = (n - 1); i >= (n / 2); --i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
            assertNotNull(bst.search(node));
            count++;
            assertEquals(count, bst.getNumberOfNodes());
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
        seed = 1490050571981L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        assertEquals(n - rm.size(), bst.getNumberOfNodes());
        
        rm.add(nodes[0]);
        bst.delete(nodes[0]);
        
        assertEquals(n - rm.size(), bst.getNumberOfNodes());
                
        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            HeapNode r = nodes[idx];
            if (!rm.contains(r)) {
                rm.add(r);
                
                bst.delete(r);
        
                assertEquals(n - rm.size(), bst.getNumberOfNodes());
                
                assertNull(bst.search(r));
            }
        }
        
        assertEquals((n - rm.size()), bst.getNumberOfNodes());
        
        boolean minChecked = false;
        
        for (int i = 0; i < n; ++i) {
            HeapNode node0 = nodes[i];
            HeapNode node = bst.search(node0);
            if (!rm.contains(node0)) {
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
        
        count = bst.getNumberOfNodes();
       
        // ==== then add n more nodes and repeat assertions
        for (int i = n; i < 2*n; ++i) {
            HeapNode node = new HeapNode();
            node.setKey(i);
            nodes[i] = node;
            bst.insert(node);
            assertNotNull(bst.search(node));
            count++;
            assertEquals(count, bst.getNumberOfNodes());
        }
        
        int n2 = n * 2;
        
        assertEquals(n2 - rm.size(), bst.getNumberOfNodes());
        
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
                
                assertEquals(n2 - rm.size(), bst.getNumberOfNodes());
                
                HeapNode node2 = bst.search(r);
                if (node2 != null) {
                    System.out.println("debug:");
                    HeapNode[] dns = bst.printPreOrderTraversalIterative(bst.root);
                    for (HeapNode dn : dns) {
                        System.out.println("node=" + dn);
                    }
                    System.out.println("r=" + r.toString() 
                        + " node2=" + node2);
                    bst.delete(r);
                }
                assertNull(node2);
            }
        }
        
        assertEquals((n2 - rm.size()), bst.getNumberOfNodes());
        
        minChecked = false;
        
        for (int i = 0; i < 2*n; ++i) {
            HeapNode node0 = nodes[i];
            HeapNode node = bst.search(node0);
            if (!rm.contains(node0)) {
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
