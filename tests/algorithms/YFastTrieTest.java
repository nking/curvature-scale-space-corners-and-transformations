package algorithms;

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

    public void test0() {
        
        System.out.println("test0");
        
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
        
    public void est1() throws Exception {
    
        //TODO: update these for YFastTrie
        
        System.out.println("test1");
        
		XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
		
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
		
        XFastTrie<XFastTrieNode<Integer>, Integer> bt
            = new XFastTrie<XFastTrieNode<Integer>, Integer>(node, it);
        
        int n = 100;
        
        List<Integer> nodes = new ArrayList<Integer>(2*n);
        
        for (int i = 0; i < n/2; ++i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            Integer ii = Integer.valueOf(i);
            nodes.add(ii);
            bt.add(ii);
        }
        
        assertEquals(n, bt.size());
                
        assertEquals(n - 1, bt.maximum().intValue());
        
        assertEquals(0L, bt.minimum().intValue());
        
        for (int i = 0; i < n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            assertNotNull(foundIndex);
            assertEquals(index.intValue(), foundIndex.intValue());
            
            /*
            smallest element in the tree with key greater
            than node.key.            
            */
            
            if (index.intValue() < (n - 1)) {
                Integer next = bt.successor(index);
                assertEquals((index.intValue() + 1), 
                    it.intValue(next));
            }
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            if (index.intValue() > 0) {
                Integer prev = bt.predecessor(index);
                assertEquals((index.intValue() - 1), 
                    it.intValue(prev));
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
        Integer nod = bt.find(nodes.get(0));
        assertNull(nod);
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
            Integer foundIndex = bt.find(index);
            if (rm.contains(index)) {
                assertNull(foundIndex);
            } else {
                if (!minChecked) {
                    Integer min = bt.minimum();
                    assertEquals(index.intValue(), min.intValue());
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex.intValue());
                if (index.intValue() < (n - 1)) {
                    Integer next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, it.intValue(next));
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
            bt.add(ii);
        }

        assertEquals((n - rm.size()) + n, bt.size());

        Integer maxExpected = Integer.valueOf(2*n - 1);
        while (rm.contains(maxExpected)) {
            maxExpected = Integer.valueOf(maxExpected.intValue() - 1);
        }
        
        assertEquals(maxExpected.intValue(), bt.maximum().intValue());
        
        for (int i = n; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex = bt.find(index);
            assertEquals(index.intValue(), foundIndex.intValue());
            if (index.intValue() < (n - 1)) {
                Integer next = bt.successor(index);
                assertEquals(index.intValue() + 1, 
                    it.intValue(next));
            }
        }

        for (int i = 0; i < n/2; ++i) {
            int idx = sr.nextInt(n);
            Integer r = nodes.get(idx);
            if (!rm.contains(r)) {
                rm.add(r);
                bt.remove(r);
                assertNull(bt.find(r));
            }
        }
        
        assertEquals((2*n - rm.size()), bt.size());

        minChecked = false;
        
        for (int i = 0; i < 2*n; ++i) {
            Integer index = nodes.get(i);
            Integer foundIndex  = bt.find(index);
            if (rm.contains(index)) {
                assertNull(foundIndex);
            } else {  
                if (!minChecked) {
                    Integer min = bt.minimum();
                    assertEquals(index.intValue(), min.intValue());
                    minChecked = true;
                }
                assertEquals(index.intValue(), foundIndex.intValue());

                if (index.intValue() < (n - 1)) {
                    Integer next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, it.intValue(next));
                            break;
                        }
                    }
                }
            }
        }
    }

    /*
    class I implements Integerizer<Integer> { 
			public int intValue(Integer i) { return i; }
		}
		int n = 50;
		YFastTrie<Integer> t = new YFastTrie<Integer>(new I());
		Random rand = new Random(0);
		System.out.println("Adding: ");
		for (int i = 0; i < n; i++) {
			int x = rand.nextInt(100*n);
			System.out.print(x + ((i < n - 1) ? "," : ""));
			t.add(x);
			// t.checkIt();
		}
		System.out.println();
		for (Pair<Integer> p : t.xft) {
			System.out.print(p.x + ",");
		}
		System.out.println();
		for (Pair<Integer> p : t.xft) {
			System.out.print(p.t + ",");
		}

		System.out.println(t);
		System.out.print("Searching: ");
		for (int i = 0; i < n; i++) {
			int x = rand.nextInt(100*n);
			System.out.print(x + "=>" + t.find(x) + ",");
		}
		System.out.println();
		System.out.println(t);
		System.out.print("Removing: ");
		for (int i = 0; i < n/2; i++) {
			Integer x = t.find(rand.nextInt(100*n));
			if (x != null) {
				System.out.print(x + ((i < n/2-1) ? "," : ""));
				System.out.flush();
				t.remove(x);
			}
			// t.checkIt();
		}
		System.out.println();
		System.out.println("Size = " + t.size());
		System.out.println(t);
		System.out.print("Searching: ");
		for (int i = 0; i < n; i++) {
			int x = rand.nextInt(100*n);
			System.out.print(x + "=>" + t.find(x) + ",");
		}
		System.out.println();
		System.out.println("done");
    */
    
    
}
