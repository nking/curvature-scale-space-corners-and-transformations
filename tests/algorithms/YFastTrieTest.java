package algorithms;

import algorithms.misc.Misc;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import thirdparty.ods.*;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class YFastTrieTest extends TestCase {
    
    public YFastTrieTest() {
    }

    public void est00() {
        
        System.out.println("test00");
        
        int w = 4;
		
        YFastTrie bt = new YFastTrie(w);
        
        assertEquals(0, bt.size());
  
        System.out.println("add 2 0010");
        
        int node2 = 2;
        int node3 = 3;
        int node4 = 4;
        int node5 = 5;
        int node6 = 6;
        
        boolean added = bt.add(node2);
        assertTrue(added);
        
        assertEquals(1, bt.size());
                
        System.out.println("add 3 0011");
        
        added = bt.add(node3);
        assertTrue(added);
                
        assertEquals(2, bt.size());
        
        System.out.println("add 5 0101");
         
        added = bt.add(node5);
        assertTrue(added);
                
        assertEquals(3, bt.size());
        
        System.out.println("add 4 0100");
       
        added = bt.add(node4);
        assertTrue(added);
       
        assertEquals(4, bt.size());
        
        assertEquals(node3, bt.find(node3));
        assertEquals(node2, bt.find(node2));
        assertEquals(node4, bt.find(node4));
        assertEquals(node5, bt.find(node5));
        
        int next = bt.successor(node4);
        assertEquals(node5, next);
        
        next = bt.successor(0);
        assertEquals(node2, next);
        
        next = bt.successor(node2);
        assertEquals(node3, next);
        
        next = bt.successor(node3);
        assertEquals(node4, next);
        
        next = bt.successor(node4);
        assertEquals(node5, next);
        
        int prev;
       
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(node5, prev);
        
        prev = bt.predecessor(node6);
        assertEquals(node5, prev);
       
        prev = bt.predecessor(node5);
        assertEquals(node4, prev);
        
        prev = bt.predecessor(node4);
        assertEquals(node3, prev);
        
        prev = bt.predecessor(node3);
        assertEquals(node2, prev);
        
        prev = bt.predecessor(node2);
        assertEquals(-1, prev);
        
        assertEquals(node2, bt.minimum());
        assertEquals(node5, bt.maximum());
        
        assertTrue(bt.remove(node3));
        assertEquals(3, bt.size());
        assertEquals(-1, bt.find(node3));
        assertEquals(node2, bt.find(node2));
        assertEquals(node4, bt.find(node4));
        assertEquals(node5, bt.find(node5));
        
        added = bt.add(node3);
        assertTrue(added);
        assertEquals(node3, bt.find(node3));
        
        assertTrue(bt.remove(node5));
        assertEquals(-1, bt.find(node5));
        
        int node4r = bt.extractMaximum();
        assertEquals(node4, node4r);
        assertEquals(-1, bt.find(node4));
        
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(node3, prev);
        
        assertEquals(node2, bt.minimum());
        assertEquals(node3, bt.maximum());
        
        int node2r = bt.extractMinimum();
        assertEquals(node2, node2r);
        assertEquals(-1, bt.find(node2));
    }
       
    public void est0() {
        
        System.out.println("test0");
        
        int w = 4;
		
        YFastTrie bt = new YFastTrie(w);
        
        assertEquals(0, bt.size());
  
        System.out.println("add 2 0010");
        
        boolean added = bt.add(2);
        assertTrue(added);
        
        assertEquals(1, bt.size());
                
        System.out.println("add 3 0011");
        
        added = bt.add(3);
        assertTrue(added);
                
        assertEquals(2, bt.size());
        
        System.out.println("add 5 0101");
         
        added = bt.add(5);
        assertTrue(added);
                
        assertEquals(3, bt.size());
        
        System.out.println("add 4 0100");
       
        added = bt.add(4);
        assertTrue(added);
       
        assertEquals(4, bt.size());
        
        //bt.debugNodes();
        
        /*
        add "2", 0010
                                               r
                               0(pad)                   null
                       0(pad)       
                           1                  
                          0            
        -----------
        add "2", 0010
        add "3", 0011
                                               r
                               0(pad)                   null
                       0(pad)       
                           1                  
                         0   1                  
        ------------
        add "2", 0010
        add "3", 0011
        add "5", 0101
                                               r
                               0(pad)                   null
                       0(pad)              1     
                           1           0    
                         0   1           1
        ------------
        add "2", 0010
        add "3", 0011
        add "5", 0101
        add "4", 0100
                                               r
                               0(pad)                   null
                       0(pad)              1     
                           1           0    
                         0   1       0   1
        */
        
        assertEquals(3, bt.find(3));
        assertEquals(2, bt.find(2));
        assertEquals(4, bt.find(4));
        assertEquals(5, bt.find(5));
        
        int next = bt.successor(4);
        assertEquals(5, next);
        
        next = bt.successor(0);
        assertEquals(2, next);
        
        next = bt.successor(2);
        assertEquals(3, next);
        
        next = bt.successor(3);
        assertEquals(4, next);
        
        next = bt.successor(4);
        assertEquals(5, next);
        
        int prev;
       
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(5, prev);
        
        prev = bt.predecessor(6);
        assertEquals(5, prev);
       
        prev = bt.predecessor(5);
        assertEquals(4, prev);
        
        prev = bt.predecessor(4);
        assertEquals(3, prev);
        
        prev = bt.predecessor(3);
        assertEquals(2, prev);
        
        prev = bt.predecessor(2);
        assertEquals(-1, prev);
        
        assertEquals(2, bt.minimum());
        assertEquals(5, bt.maximum());
        
        assertTrue(bt.remove(3));
        assertEquals(3, bt.size());
        assertEquals(-1, bt.find(3));
        assertEquals(2, bt.find(2));
        assertEquals(4, bt.find(4));
        assertEquals(5, bt.find(5));

        assertEquals(-1, bt.find(0));
        
        added = bt.add(3);
        assertTrue(added);
        assertEquals(3, bt.find(3));
        
        assertTrue(bt.remove(5));
        assertEquals(-1, bt.find(5));
        assertTrue(bt.remove(4));
        assertEquals(-1, bt.find(4));
        
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(3, prev);
        
        assertEquals(2, bt.minimum());
        assertEquals(3, bt.maximum());
    }
        
    public void est1() throws Exception {
    
        System.out.println("test1");
        
        YFastTrie bt = new YFastTrie();
        
        int n = 100;
        
        List<Integer> nodes = new ArrayList<Integer>(2*n);
        
        for (int i = 0; i < n/2; ++i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(i);
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(i);
        }
        
        assertEquals(n, bt.size());
                
        assertEquals(n - 1, bt.maximum());
        
        assertEquals(0L, bt.minimum());
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            int foundIndex = bt.find(index);
            assertTrue(foundIndex > -1);
            assertEquals(index.intValue(), foundIndex);
            
            /*
            smallest element in the tree with key greater
            than node.key.            
            */
            
            if (index.intValue() < (n - 1)) {
                int next = bt.successor(index);
                assertEquals((index.intValue() + 1), next);
            }
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            if (index.intValue() > 0) {
                int prev = bt.predecessor(index);
                assertEquals((index.intValue() - 1), prev);
            }
        }
        
        // remove some nodes randomly
        Set<Integer> rm = new HashSet<Integer>();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1465940667831L;
        sr.setSeed(seed);
        System.out.println("SEED=" + seed);
        
        bt.remove(nodes.get(0));
        int nod = bt.find(nodes.get(0));
        assertEquals(-1, nod);
        rm.add(nodes.get(0));
        
        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                assertTrue(bt.remove(r));
                rm.add(r);
            }
        }
        
        assertEquals((n - rm.size()), bt.size());
                
        boolean minChecked = false;
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            int foundIndex = bt.find(index);
            if (rm.contains(index)) {
                assertEquals(-1, foundIndex);
            } else {
                if (!minChecked) {
                    int min = bt.minimum();
                    assertEquals(index.intValue(), min);
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex);
                if (index.intValue() < (n - 1)) {
                    int next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, next);
                            break;
                        }
                    }
                }
            }
        }

        // ==== then add n more nodes and repeat assertions
        for (int i = n; i < 2*n; ++i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(i);
        }

        assertEquals((n - rm.size()) + n, bt.size());

        Integer maxExpected = Integer.valueOf(2*n - 1);
        while (rm.contains(maxExpected)) {
            maxExpected = Integer.valueOf(maxExpected.intValue() - 1);
        }
        
        assertEquals(maxExpected.intValue(), bt.maximum());
        
        for (int i = n; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            int foundIndex = bt.find(index);
            assertEquals(index.intValue(), foundIndex);
            if (index.intValue() < (n - 1)) {
                int next = bt.successor(index);
                assertEquals(index.intValue() + 1, next);
            }
        }

        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                rm.add(r);
                bt.remove(r);
                assertEquals(-1, bt.find(r));
            }
        }
        
        assertEquals((2*n - rm.size()), bt.size());

        minChecked = false;
        
        for (int i = 0; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            int foundIndex  = bt.find(index);
            if (rm.contains(index)) {
                assertEquals(-1, foundIndex);
            } else {  
                if (!minChecked) {
                    int min = bt.minimum();
                    assertEquals(index.intValue(), min);
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex);

                if (index.intValue() < (n - 1)) {
                    int next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, next);
                            break;
                        }
                    }
                }
            }
        }
    }
    
     public void test_random_large() {
        
        System.out.println("test_random_large");
        
        Random rng = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        seed = 1498089496315L;
        rng.setSeed(seed);
        System.out.println("test_random_large SEED=" + seed);
        
        int n = 100000;
        
        int w = 30;
        
        int maxC = (1 << w) - 1;
		
        YFastTrie bt = new YFastTrie(w);
     
        TIntList list = new TIntArrayList();
     
        for (int i = 0; i < n; ++i) {
            int v = rng.nextInt(maxC);
            //System.out.println("i=" + i + " v=" + v);
            bt.add(v);
            list.add(v);
        }
        
        list.sort();
        
        assertEquals(n, bt.size());
        
        int i = 0;
        while (bt.size() > 0) {
            
            int bv = bt.extractMinimum();
            int v = list.get(i);
            
            assertEquals(v, bv);
            
            ++i;
        }
     }
    
}
