package thirdparty.ods;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/*
The class is adapted from the open datastructures source code
http://opendatastructures.org/ods-java.pdf

"The source code available there is released under a Creative Commons
Attribution license, meaning that anyone is free to share: to copy, distribute
and transmit the work; and to remix: to adapt the work, including
the right to make commercial use of the work. The only condition on
these rights is attribution: you must acknowledge that the derived work
contains code and/or text from opendatastructures.org.
https://github.com/patmorin/ods
*/
@SuppressWarnings("unchecked")
public class XFastTrie<S extends XFastTrieNode<T>, T> 
	extends BinaryTrie<S, T> {

	/**
	 * The hash tables used to store prefixes
	 */
    protected final List<Map<Integer, S>> t;
	
	public XFastTrie(S sampleNode, Integerizer<T> it)  {
		super(sampleNode, it);
        t = new ArrayList<Map<Integer, S>>();
		S nil = (S)new XFastTrieNode<T>();
		nil.prefix = Integer.MIN_VALUE;
		for (int i = 0; i <= w; i++) {
            t.add(new HashMap<Integer, S>());
		}
		t.get(0).put(Integer.valueOf(0), r);
	}
    
    public XFastTrie(S sampleNode, Integerizer<T> it,
        int smallerWordSize)  {
		super(sampleNode, it, smallerWordSize);
        t = new ArrayList<Map<Integer, S>>();
		S nil = (S)new XFastTrieNode<T>();
		nil.prefix = Integer.MIN_VALUE;
		for (int i = 0; i <= w; i++) {
            t.add(new HashMap<Integer, S>());
		}
		t.get(0).put(Integer.valueOf(0), r);
	}

	@SuppressWarnings("unchecked")
	public XFastTrie(Integerizer<T> it)  {
		this((S)new XFastTrieNode<T>(), it);
	}
	
	public boolean add(T x) {
		if (super.add(x)) {
			int i, c = 0, ix = it.intValue(x);
			S u = (S)r.child[(ix>>>w-1) & 1];
			for (i = 1; i <= w; i++) {
				u.prefix = ix >>> w-i;
				t.get(i).put(Integer.valueOf(u.prefix), u);
				c = (ix >>> w-i-1) & 1;
				u = (u.child[c] != null) ? (S)u.child[c] : null;
			}
			return true;
		}
		return false;
	}
	
	public boolean remove(T x) {
		// 1 - find leaf, u, containing x
		int i = 0, c, ix = it.intValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can remove is " + maxC);
        }
		S u = r;
        S v;		
        int l = 0, h = w+1;
        int prefix = -1;
        // binary search over range w
		while (h-l > 1) {
			i = (l+h)/2;
			prefix = ix >>> w-i;
            v = t.get(i).get(Integer.valueOf(prefix));
			if (v == null) {
				h = i;
			} else {
				u = v;
				l = i;
			}
		}
       
		// 2 - remove u from linked list
		S pred = (u.child[prev] != null) ?
            (S)u.child[prev] : null;   // predecessor
		S succ = (u.child[next] != null) ?
            (S)u.child[next] : null;   // successor
		pred.child[next] = succ;
		succ.child[prev] = pred;
		u.child[next] = u.child[prev] = null;
		S w = u;
		// 3 - delete nodes on path to u
		while (w != r && w.child[left] == null && w.child[right] == null) {
			if (w == w.parent.child[left]) {
				w.parent.child[left] = null;
            } else { // u == u.parent.child[right] 
				w.parent.child[right] = null;
            }
            prefix = w.prefix;
			t.get(i--).remove(Integer.valueOf(w.prefix));
			w = (w.parent != null) ? (S)w.parent : null;
		}
		// 4 - update jump pointers
		w.jump = (w.child[left] == null) ? succ : pred;
		w = (w.parent != null) ? (S)w.parent : null;
		while (w != null) {
			if (w.jump == u)
				w.jump = (w.child[left] == null) ? succ : pred;
			w = (w.parent != null) ? (S)w.parent : null;
		}
		n--;
		return true;
	}

    /**
     * find node with key ix.
     * runtime complexity is O(1)
     * @param ix
     * @return 
     */
	protected S findNode(int ix) {
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can search for is " + maxC);
        }
		S q = t.get(w).get(Integer.valueOf(ix));
        return q;
	}
	
    /**
     * find node key, with key x.
     * runtime complexity is O(1).
     * @param x
     * @return 
     */
	public T find(T x) {
        
        int ix = it.intValue(x);
		S q = findNode(ix);
        if (q == null) {
            return null;
        }
        
        return q.x;
    }
    
    /**
	 * Find the node that contains the successor of x.
	 * runtime complexity is O(log_2(w)).
     * @param x
	 * @return The node after the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
	public S successor(T x) {
        int ix = it.intValue(x);
        return successor(ix);
    }
    
    /**
	 * Find the node that contains the successor of x.
	 * runtime complexity is O(log_2(w)).
     * @param ix
	 * @return The node after the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
	protected S successor(int ix) {
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can search for is " + maxC);
        }
		// find lowest node that is an ancestor of ix
		int l = 0, h = w+1;
		S v, u = r;
        int prefix;
        // binary search over range w
		while (h-l > 1) {
			int i = (l+h)/2;
			prefix = ix >>> w-i;
            v = t.get(i).get(Integer.valueOf(prefix));
			if (v == null) {
				h = i;
			} else {
				u = v;
				l = i;
			}
		}
		BinaryTrieNode<T> successor;
		if (l == w) {
            successor = u.child[next];
        } else {
		    int c = (ix >>> w-l-1) & 1;
            successor = (c == 0) ? u.jump : u.jump.child[1];
        }
		return (successor != null) ? (S)successor : null;
	}
    
    /**
	 * Find the node that contains the predecessor of x.
	 * runtime complexity is O(log2(w)).
     * @param x
	 * @return The node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
	public S predecessor(T x) {
		int ix = it.intValue(x);
        return predecessor(ix);
    }
    
    /**
	 * Find the node that contains the predecessor of x.
	 * runtime complexity is O(log2(w)).
     * @param ix
	 * @return The node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
	protected S predecessor(int ix) {
		
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can search for is " + maxC);
        }
		// find lowest node that is an ancestor of ix
		int l = 0, h = w+1;
		S v, u = r;
        int prefix = -1;
        // binary search over range w
		while (h-l > 1) {
			int i = (l+h)/2;
			prefix = ix >>> w-i;
            v = t.get(i).get(Integer.valueOf(prefix));
			if (v == null) {
				h = i;
			} else {
				u = v;
				l = i;
			}
		}
        
        if (l == w) {
            if (u.child[prev] == null) {
                return null;
            }
            if (u.child[prev].x == null) {
                return null;
            }
            return (S)u.child[prev]; 
        }
                         
        int c = (ix >>> w-l-1) & 1;
        
        if (c == 1 && u.jump != null) {
            return (S)u.jump;
        }
        
        //TODO: revisit this.  below never occurs
        
        XFastTrieNode<T> pred;	
        if (u.jump.child[1] == null) {
            pred = null;
        } else {
            pred = (XFastTrieNode<T>) u.jump.child[1];
        }
		return (pred != null) ? (S)pred : null;
	}

	public void clear() {
		super.clear();
		for (Map<Integer, S> m : t) 
			m.clear();
	}
    
    public T minimum() {
        if (t.get(w).containsKey(Integer.valueOf(0))) {
            return t.get(w).get(Integer.valueOf(0)).x;
        }
        S min = successor(0);
        return (min != null) ? min.x : null;
    }
    
    public T maximum() {
        if (t.get(w).containsKey(Integer.valueOf(maxC))) {
            return t.get(w).get(Integer.valueOf(maxC)).x;
        }
        S max = predecessor(maxC);
        return (max != null) ? max.x : null;
    }
	
    /**
     * print out the dummy link keys then print the
     * root tree in level traversal
     */
    void debugNodes() {
        Set<Integer> dummyHashCodes = new HashSet<Integer>();
        S node = dummy;
        //System.out.println("dummy.hashCode=" + dummy.hashCode());
        System.out.print("\ndummy=");
        do {
            int dhc = node.hashCode();
            System.out.print(node.x + ", ");
            dummyHashCodes.add(Integer.valueOf(dhc));
            node = (S)node.child[1];
        } while (!node.equals(dummy));
        
        System.out.println();
        
        if (r == null) {
            return;
        }
        
        for (int i = 1; i <= w; ++i) {
            System.out.println("level=" + i);
            Map<Integer, S> nodesMap = t.get(i);
            for (Entry<Integer, S> entry : nodesMap.entrySet()) {
                S nodeL = entry.getValue();
                System.out.println(nodeL.toString2());
            }
        }
    }
    
}
