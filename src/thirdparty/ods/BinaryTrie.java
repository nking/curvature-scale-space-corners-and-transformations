package thirdparty.ods;

/*
The class is adapted from the open datastructures source code
http://opendatastructures.org/ods-java.pdf

"The source code available there is released under a Creative Commons
Attribution license, meaning that anyone is free to share: to copy, distribute
and transmit the work; and to remix: to adapt the work, including
the right to make commercial use of the work. The only condition on
these rights is attribution: you must acknowledge that the derived work
contains code and/or text from opendatastructures.org.
http://github.com/patmorin/ods
*/
import java.util.ArrayDeque;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

@SuppressWarnings({"unchecked"})
public class BinaryTrie<S extends BinaryTrieNode<T>, T> 
    implements SSet<T> {
	
	protected static final int prev = 0;
	protected static final int next = 1;
	protected static final int left = 0;
	protected static final int right = 1;
	
	protected int w = 32;
    protected final int maxC;
	
	/**
	* The root node
	*/
	protected final S r;
	
	/**
	 * The dummy node in the doubly-linked list
	 */
	protected final S dummy;
	
	/**
	 * For converting elements of type T into integers
	 */
	protected final Integerizer<T> it;

	/**
	 * The number of elements stored in the trie
	 */
	int n;

	/**
	 * To make a node factory
	 */
	protected final S sampleNode;
	
	/**
	 * Allocate a new node.  if S is an extension
     * of BinaryTrieNode,
     * make sure the class is public and if it's an
     * inner class, it has to be static too.
	 * @return
	 */
	protected S newNode() {
		try {
            S u = (S)sampleNode.getClass().newInstance();
			u.parent = u.child[0] = u.child[1] = null;
			return u;
		} catch (Exception e) {
            throw new UnsupportedOperationException(
               "sampleNode constructor is not reachable "
               + "as a public merhod with a no arguments constructor");
		}
	}

    /**
     * 
     * @param sampleNode a node instance that is used for 
     * its class type when creating other nodes, such as the 
     * root and linked-list sentinel nodes.
     * @param it class to provide the inner node key which
     * prefixes are extracted from
     */
	public BinaryTrie(S sampleNode, Integerizer<T> it) {
		this.sampleNode = sampleNode;
		this.dummy = newNode();
		dummy.child[prev] = dummy.child[next] = dummy;
		this.r = newNode();
		r.jump = dummy;
		this.it = it;
		n = 0;
        maxC = (1 << 31) - 1;
	}
	
    /**
     * 
     * @param sampleNode a node instance that is used for 
     * its class type when creating other nodes, such as the 
     * root and linked-list sentinel nodes.
     * @param it class to provide inner node keys which the
     * prefixes are extracted from
     * @param smallerWordSize when a word size smaller than
     * 32 bits is known a priori, can set this.
     */
	public BinaryTrie(S sampleNode, Integerizer<T> it,
        int smallerWordSize) {
        if (smallerWordSize < 32 && smallerWordSize > 1) {
            this.w = smallerWordSize;
        } else {
            throw new IllegalStateException("smallerWordSize "
                + " shoulw be greater than 1 and less than 32");
        }
        maxC = (1 << (w - 1)) - 1;
		this.sampleNode = sampleNode;
		this.dummy = newNode();
		dummy.child[prev] = dummy.child[next] = dummy;
		this.r = newNode();
		r.jump = dummy;
		this.it = it;
		n = 0;
	}
	
    /**
     * add x to the tree.  runtime complexity is O(w).
     * @param x
     * @return 
     */
	public boolean add(T x) {
		int i, c = 0, ix = it.intValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can add is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		// 1 - search for ix until falling out of the trie
		for (i = 0; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			if (u.child[c] == null) break;
			u = u.child[c];
		}		
		if (i == w) return false; // already contains x - abort
		BinaryTrieNode<T> pred = (c == right) ? u.jump : u.jump.child[0];
		u.jump = null;  // u will have two children shortly
		// 2 - add path to ix
		for (; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			u.child[c] = newNode();
			u.child[c].parent = u;
			u = u.child[c];
		}
		u.x = x;
		// 3 - add u to linked list
		u.child[prev] = pred;
		u.child[next] = pred.child[next];
		u.child[prev].child[next] = u;
		u.child[next].child[prev] = u;
		// 4 - walk back up, updating jump pointers
		BinaryTrieNode<T> v = u.parent;
		while (v != null) {
			if ((v.child[left] == null 
	        	&& (v.jump == null ||
                it.intValue((T)v.jump.x) > ix))
			|| (v.child[right] == null 
	    		&& (v.jump == null || 
                (v.jump.x != null && it.intValue((T)v.jump.x) < ix))
                )) {
				v.jump = u;
            }
			v = v.parent;
		}
        
		n++;
        
		return true;
	}

    public T find(T x) {
        int ix = it.intValue(x);
        T v = findValue(ix);
        if (v == null) {
            return v;
        } else if (it.intValue(v) == ix) {
            return v;
        }
        return null;
    }
    
	protected T findValue(int ix) {
		int i, c = 0;
		if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can find is " + maxC);
        }
        BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		if (i == w) return u.x;  // found it
        if (c == 0) {
            u = u.jump;
        } else {
            u = u.jump.child[next]; 
        }
		return u == dummy ? null : u.x;
	}

	public boolean remove(T x) {
		// 1 - find leaf, u, containing x
		int i, c, ix = it.intValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can remove is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			if (u.child[c] == null) return false;
			u = u.child[c];
		}
		// 2 - remove u from linked list
		u.child[prev].child[next] = u.child[next];
		u.child[next].child[prev] = u.child[prev];
		BinaryTrieNode<T> v = u;
		// 3 - delete nodes on path to u
		for (i = w-1; i >= 0; i--) {
			c = (ix >>> w-i-1) & 1;
			v = v.parent;
			v.child[c] = null;
			if (v.child[1-c] != null) break;
		}
		// 4 - update jump pointers
		c = (ix >>> w-i-1) & 1;
		v.jump = u.child[1-c];
		v = v.parent;
		i--;
		for (; i >= 0; i--) {
			c = (ix >>> w-i-1) & 1;
			if (v.jump == u) {
				v.jump = u.child[1-c];
            }
			v = v.parent;
		}
		n--;
		return true;
	}

	/**
	 * Part of SSet interface, but not really relevant here
	 * TODO: We can still implement this
	 */
	public Comparator<? super T> comparator() {
		return new Comparator<T>() {
			public int compare(T a, T b) {
				return it.intValue(a) - it.intValue(b);
			}
		};
	}

	public T findGE(T x) {
        int ix = it.intValue(x);
		return findValue(ix);
	}

	/**
	 * Find the node that contains the predecessor of x.
	 * runtime complexity is O(w).
     * @param x
	 * @return The node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
	public S predecessor(T x) {
		int i, c = 0, ix = it.intValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value argument is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		BinaryTrieNode<T> pred;
		if (i == w) pred = u.child[prev];
		else pred = (c == 1) ? u.jump : u.jump.child[0]; 
		return (pred != null) ? (S)pred : null;
	}
    
    /**
	 * Find the node that contains the successor of x.
	 * runtime complexity is O(w).
     * @param x
	 * @return The node after the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
	public S successor(T x) {
		int i, c = 0, ix = it.intValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value argument is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		BinaryTrieNode<T> successor;
		if (i == w) successor = u.child[next];
		else successor = (c == 0) ? u.jump : u.jump.child[1]; 
		return (successor != null) ? (S)successor : null;
	}
	
    @Override
	public T findLT(T x) {
		S pred = predecessor(x);
		return (pred == dummy) ? null : ((S)(pred.child[next])).x;
	}

	/**
	 * This is just a simple linked-list iterator
	 * @author morin
	 *
	 */
	protected class TrieIterator implements Iterator<T> {
		
        protected BinaryTrieNode<T> p;
		
		public TrieIterator(S p) {
			this.p = p;
		}
		
		public boolean hasNext() {
			return !p.equals(dummy);
		}

        @Override
		public T next() {
			T x = p.x;
			p = p.child[1];
			return x;
		}
		
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}

	public Iterator<T> iterator(T x) {
		BinaryTrieNode<T> pred = predecessor(x);
		return new TrieIterator((S)pred.child[next]);
	}

	public Iterator<T> iterator() {
		return new TrieIterator((S)dummy.child[next]);
	}
	
	public int size() {
		return n;
	}
	
	public void clear() {
		n = 0;
		r.child[0] = r.child[1] = null;
		r.jump = dummy;
		dummy.child[0] = dummy.child[1] = dummy;
	}
	
    /**
     * print out the dummy link keys then print the
     * root tree in level  traversal
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
        
        int dhc = dummy.hashCode();
        
        Deque<S> q0 = new ArrayDeque<S>();
        Deque<S> q1 = new ArrayDeque<S>();
        q0.offer(r);
        
        int count = 0;
        boolean skip = true;
        while(!q0.isEmpty() && count < w) {
            System.out.println("count=" + count);
            while(!q0.isEmpty()) {
                node = q0.poll();
                if (!node.equals(r)) {
                    System.out.println(node.toString2());
                }
                if (node.child[0] != null) {
                    int hc = node.child[0].hashCode();
                    if (count < w) {
                        q1.offer((S) node.child[0]);
                    }
                }
                if (node.child[1] != null) {
                    int hc = node.child[1].hashCode();
                    if (count < w) {
                        q1.offer((S) node.child[1]);
                    }
                }
            }
            if (!skip) {
                count++;
            } else {
                skip = false;
            }
            q0.addAll(q1);
            q1.clear();
        }        
    }
}
