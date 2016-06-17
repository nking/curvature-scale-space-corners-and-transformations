package thirdparty.ods;

import algorithms.imageProcessing.HeapNode;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class XFastTrieTest extends TestCase {
    
    public XFastTrieTest() {
    }

    public void test0() {
        
        System.out.println("test0");
        
		XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
		
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        int w = 4;
		
        XFastTrie<XFastTrieNode<Integer>, Integer> bt
            = new XFastTrie<XFastTrieNode<Integer>, Integer>(node, it, w);
        
        assertEquals(0, bt.size());
  
        System.out.println("add 2 0010");
        
        boolean added = bt.add(Integer.valueOf(2));
        assertTrue(added);
        
        assertEquals(1, bt.size());
                
        System.out.println("add 3 0011");
        
        added = bt.add(Integer.valueOf(3));
        assertTrue(added);
                
        assertEquals(2, bt.size());
        
        System.out.println("add 5 0101");
         
        added = bt.add(Integer.valueOf(5));
        assertTrue(added);
                
        assertEquals(3, bt.size());
        
        System.out.println("add 4 0100");
       
        added = bt.add(Integer.valueOf(4));
        assertTrue(added);
       
        assertEquals(4, bt.size());
        
        bt.debugNodes();
        
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
        
        assertEquals(3, bt.find(Integer.valueOf(3)).intValue());
        assertEquals(2, bt.find(Integer.valueOf(2)).intValue());
        assertEquals(4, bt.find(Integer.valueOf(4)).intValue());
        assertEquals(5, bt.find(Integer.valueOf(5)).intValue());
        
        XFastTrieNode<Integer> next = bt.successor(4);
        assertEquals(5, it.intValue(next.x));
        
        next = bt.successor(0);
        assertEquals(2, it.intValue(next.x));
        
        next = bt.successor(2);
        assertEquals(3, it.intValue(next.x));
        
        next = bt.successor(3);
        assertEquals(4, it.intValue(next.x));
        
        next = bt.successor(4);
        assertEquals(5, it.intValue(next.x));
        
        XFastTrieNode<Integer> prev;
       
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(5, it.intValue(prev.x));
        
        prev = bt.predecessor(6);
        assertEquals(5, it.intValue(prev.x));
       
        prev = bt.predecessor(5);
        assertEquals(4, it.intValue(prev.x));
        
        prev = bt.predecessor(4);
        assertEquals(3, it.intValue(prev.x));
        
        prev = bt.predecessor(3);
        assertEquals(2, it.intValue(prev.x));
        
        prev = bt.predecessor(2);
        assertNull(prev);
        
        assertEquals(2, bt.minimum().intValue());
        assertEquals(5, bt.maximum().intValue());
        
        assertTrue(bt.remove(3));
        assertEquals(3, bt.size());
        assertNull(bt.find(Integer.valueOf(3)));
        assertEquals(2, bt.find(Integer.valueOf(2)).intValue());
        assertEquals(4, bt.find(Integer.valueOf(4)).intValue());
        assertEquals(5, bt.find(Integer.valueOf(5)).intValue());

        assertNull(bt.find(Integer.valueOf(0)));
        
        added = bt.add(Integer.valueOf(3));
        assertTrue(added);
        assertEquals(3, bt.find(Integer.valueOf(3)).intValue());
        
        assertTrue(bt.remove(5));
        assertNull(bt.find(Integer.valueOf(5)));
        assertTrue(bt.remove(4));
        assertNull(bt.find(Integer.valueOf(4)));
        
        prev = bt.predecessor((1<<(w-1)) - 1);
        assertEquals(3, it.intValue(prev.x));
        
        assertEquals(2, bt.minimum().intValue());
        assertEquals(3, bt.maximum().intValue());
    }
        
    public void test1() throws Exception {
    
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
                XFastTrieNode<Integer> next = bt.successor(index);
                assertEquals((index.intValue() + 1), 
                    it.intValue(next.x));
            }
            /*
            the largest element in the tree with key smaller 
            than node.key
            */
            if (index.intValue() > 0) {
                XFastTrieNode<Integer> prev 
                    = bt.predecessor(index);
                assertEquals((index.intValue() - 1), 
                    it.intValue(prev.x));
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
                    XFastTrieNode<Integer> next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, it.intValue(next.x));
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
                XFastTrieNode<Integer> next = bt.successor(index);
                assertEquals(index.intValue() + 1, 
                    it.intValue(next.x));
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
                    XFastTrieNode<Integer> next = bt.successor(index);
                    int expected = index.intValue() + 1;
                    while (expected < n) {
                        if (rm.contains(Integer.valueOf(expected))) {
                            ++expected;
                        } else {
                            assertEquals(expected, it.intValue(next.x));
                            break;
                        }
                    }
                }
            }
        }
    }
    
    public void estMain() {
        System.out.println("main");
        
        int n = 20;
		
		XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
		
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
		
        BinaryTrie<XFastTrieNode<Integer>, Integer> t
            = new BinaryTrie<XFastTrieNode<Integer>, Integer>(node, it);
        
		System.out.println(t.getClass());
		Random rand = new Random(0);
		System.out.println("Adding: ");
		for (int i = 0; i < n; i++) {
			int x = rand.nextInt(100*n);
			System.out.print(x + ((i < n - 1) ? "," : ""));
			t.add(x);
			//t.checkIt();
		}
		System.out.println();
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
			//t.checkIt();
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
    }
    
}
