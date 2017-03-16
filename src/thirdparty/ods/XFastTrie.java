package thirdparty.ods;

import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;

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
    protected final List<TIntObjectHashMap<S>> t;
	
	public XFastTrie(S sampleNode, Integerizer<T> it)  {
		super(sampleNode, it);
        t = new ArrayList<TIntObjectHashMap<S>>();
		S nil = (S)new XFastTrieNode<T>();
		nil.prefix = Integer.MIN_VALUE;
		for (int i = 0; i <= w; i++) {
            t.add(new TIntObjectHashMap<S>());
		}
		t.get(0).put(0, r);
	}
    
    public XFastTrie(S sampleNode, Integerizer<T> it,
        int smallerWordSize)  {
		super(sampleNode, it, smallerWordSize);
        t = new ArrayList<TIntObjectHashMap<S>>();
		S nil = (S)new XFastTrieNode<T>();
		nil.prefix = Integer.MIN_VALUE;
		for (int i = 0; i <= w; i++) {
            t.add(new TIntObjectHashMap<S>());
		}
		t.get(0).put(0, r);
	}

	@SuppressWarnings("unchecked")
	public XFastTrie(Integerizer<T> it)  {
		this((S)new XFastTrieNode<T>(), it);
	}
	
    /**
     * runtime complexity is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor, else is 32
     * and l is the prefix tree already filled leading
     * up to the value x.
     * @param x
     * @return 
     */
	public boolean add(T x) {
        final int ix = it.intValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can add is " + maxC
               + " . ix=" + ix);
        }
        S u = r;
        S v;
        int i;
        int l = 0, h = w+1;
        int prefix = -1;
        // binary search over range w;  rt is < O(lg_2(w))
		while (h-l > 1) {
			i = (l+h)/2;
			prefix = ix >>> w-i;
            
            v = t.get(i).get(prefix);
			if (v == null) {
				h = i;
			} else {
				u = v;
				l = i;
			}
		}
        
        if (l == w) return false; // already contains x - abort
       
        int c = (ix >>> w-l-1) & 1;
        
        S pred = (c == right) ? (S)u.jump : (S)u.jump.child[0];
		u.jump = null;  // u will have two children shortly
        i = l;
        // 2 - add path to ix. rt is O(w-l)
		for (; i < w; i++) {
			c = (ix >>> w-i-1) & 1;
			u.child[c] = newNode();
			u.child[c].parent = u;
			u = (S) u.child[c];
		}
		u.x = x;
		// 3 - add u to linked list
		u.child[prev] = pred;
		u.child[next] = pred.child[next];
		u.child[prev].child[next] = u;
		u.child[next].child[prev] = u;
		// 4 - walk back up, updating jump pointers
	    v = (u.parent != null) ? (S)u.parent : null;
		while (v != null) {
			if ((v.child[left] == null 
	        	&& (v.jump == null ||
                it.intValue(v.jump.x) > ix))
			|| (v.child[right] == null 
	    		&& (v.jump == null || 
                (v.jump.x != null && it.intValue(v.jump.x) < ix))
                )) {
				v.jump = u;
            }
			v = (v.parent != null) ? (S)v.parent : null;
		}
        
		n++;
       
        u = (S) r.child[(ix >>> w - 1) & 1];
        for (i = 1; i <= w; i++) {
            u.prefix = ix >>> w - i;
            t.get(i).put(u.prefix, u);
            c = (ix >>> w - i - 1) & 1;
            u = (u.child[c] != null) ? (S) u.child[c] : null;
        }
        return true;
	}
	
    /**
     * runtime complexity is O(log_2(w)) + O(w-l)
     * where w is the number of
     * bits set in the constructor, else is 32
     * and l is the prefix tree already filled leading
     * up to the value x.
     * 
     * @param x
     * @return 
     */
    @Override
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
            v = t.get(i).get(prefix);
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
			t.get(i--).remove(w.prefix);
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
		S q = t.get(w).get(ix);
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
	 * Find the key of the node that contains the successor of x.
	 * runtime complexity is O(log_2(w)) where w is the number of
     * bits set in the constructor, else is 32.
     * @param x
	 * @return The node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
    @Override
	public T successor(T x) {
        S q = successorNode(x);
        if (q != null) {
            return q.x;
        }
        return null;
    }
   
	protected T successor(int ix) {
        S q = successorNode(ix);
        if (q != null) {
            return q.x;
        }
        return null;
    }
    
	protected S successorNode(T x) {
        int ix = it.intValue(x);
        return successorNode(ix);
    }
    
	protected S successorNode(int ix) {
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
            v = t.get(i).get(prefix);
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
	 * Find the key of the node that contains the predecessor of x.
	 * runtime complexity is O(log_2(w)) where w is the number of
     * bits set in the constructor, else is 32.
     * @param x
	 * @return The node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
    @Override
	public T predecessor(T x) {
        S q = predecessorNode(x);
        if (q != null) {
            // root node will return null here too
            return q.x;
        }
        return null;
    }
    
	protected T predecessor(int ix) {
        S q = predecessorNode(ix);
        if (q != null) {
            return q.x;
        }
        return null;
    }
    
	protected S predecessorNode(T x) {
        int ix = it.intValue(x);
        return predecessorNode(ix);
    }
    
	protected S predecessorNode(int ix) {
		
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
            v = t.get(i).get(prefix);
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
                
        XFastTrieNode<T> pred;	
        if (u.jump.child[0] == null) {
            pred = null;
        } else {
            pred = (XFastTrieNode<T>) u.jump.child[0];
        }
		return (pred != null) ? (S)pred : null;
	}

	public void clear() {
		super.clear();
		for (TIntObjectHashMap<S> m : t) 
			m.clear();
	}
    
    /**
     * find the key of the minimum value node.
     * runtime complexity is O(log_2(w)) where w is the number of
     * bits set in the constructor, else is 32.
     * @return 
     */
    @Override
    public T minimum() {
        if (t.get(w).containsKey(0)) {
            return t.get(w).get(0).x;
        }
        return successor(0);
    }
    
    /**
     * find the key of the minimum value node.
     * runtime complexity is O(log_2(w)) where w is the number of
     * bits set in the constructor, else is 32.
     * @return 
     */
    @Override
    public T maximum() {
        if (t.get(w).containsKey(maxC)) {
            return t.get(w).get(maxC).x;
        }
        return predecessor(maxC);
    }
	
    /**
     * print out the dummy link keys then print the
     * root tree in level traversal
     */
    void debugNodes() {
        TIntSet dummyHashCodes = new TIntHashSet();
        S node = dummy;
        //System.out.println("dummy.hashCode=" + dummy.hashCode());
        System.out.print("\ndummy=");
        do {
            int dhc = node.hashCode();
            System.out.print(node.x + ", ");
            dummyHashCodes.add(dhc);
            node = (S)node.child[1];
        } while (!node.equals(dummy));
        
        System.out.println();
        
        if (r == null) {
            return;
        }
        
        for (int i = 1; i <= w; ++i) {
            System.out.println("level=" + i);
            TIntObjectHashMap<S> nodesMap = t.get(i);
           
            TIntObjectIterator<S> iter = nodesMap.iterator();
            for (int ii = nodesMap.size(); ii-- > 0;) {
                iter.advance();
                S nodeL = iter.value();
                System.out.println(nodeL.toString2());
            }
        }
    }
    
}
