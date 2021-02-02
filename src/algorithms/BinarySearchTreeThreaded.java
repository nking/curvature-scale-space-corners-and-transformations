package algorithms;

import algorithms.heapsAndPQs.HeapNode;
import java.lang.reflect.Array;
import java.util.HashMap;
import java.util.Map;

/**
 * Building upon the BinarySearchTree by adding
 * double threading.
 * 
 * adapted from the following: 
 * http://adtinfo.org/libavl.html/TBST-Data-Types.html
 * The GNU libavl 2.0.2 has the following license:
 * http://adtinfo.org/libavl.html/Code-License.html
 * which is 
 * GNU General Public License as
   published by the Free Software Foundation; either version 2.
 * 
 * TODO:
 * The class is not completely tested yet.  Need to
 * test more thoroughly for integrity of the bst 
 * after root node deletions.
 * 
 * @author nichole
 */
public class BinarySearchTreeThreaded<T extends HeapNode> {

    //TODO: note, have made a snall change which needs
    // a careful look and possibly changes at compariions.
    // the change was to allow nodes with the same key
    // to be inserted (and no new nodes are created internally)
    
    /**
     * tracks whether a node's left and right links 
     * are child pointers or 
     * threads for left and for right.
     */
    protected Map<T, Integer> threadMap = new HashMap<T, Integer>();
    
    protected int n = 0;
    
    protected T root = null;
    
    public BinarySearchTreeThreaded() {
    }
    
    /*
    if the left subtree is empty, 
        LLINK points to the in-order predecessor
    if the right subtree is empty, 
        RLINK points to the in-order successor.
    */
    private void setChildLinksState(T node, boolean leftIsAChild,
        boolean rightIsAChild) {
        /*
        using set bit operations for states:
        0 = both are unset.
        1 = if set, left is left (child) node, else is
                     in-order predecessor
        2 = if set, right is right (child), else is 
                     in-order successor
        where in=order is left subtree, root, right subtree
        */
        Integer v = threadMap.get(node);
        int vsets = (v == null) ? 0 : v.intValue();
        if (leftIsAChild) {            
            vsets |= (1 << 1);
        } else {
            // unset bit 1
            vsets &= ~(1 << 1);
        }
        if (rightIsAChild) {
            vsets |= (1 << 2);
        } else {
            // unset bit 1
            vsets &= ~(1 << 2);
        }
        threadMap.put(node, Integer.valueOf(vsets));
    }
    
    private void setLeftIsAChild(T node, boolean leftIsAChild) {
        /*
        using set bit operations for states:
        0 = both are unset.
        1 = if set, left is left (child) node, else is
                     in-order predecessor
        2 = if set, right is right (child), else is 
                     in-order successor
        where in-order is left subtree, root, right subtree
        */
        Integer v = threadMap.get(node);
        int vsets = (v == null) ? 0 : v.intValue();
        if (leftIsAChild) {            
            vsets |= (1 << 1);
        } else {
            // unset
            vsets &= ~(1 << 1);
        }
        threadMap.put(node, Integer.valueOf(vsets));
    }
    
    private void setRightIsAChild(T node, boolean rightIsAChild) {
        /*
        using set bit operations for states:
        0 = both are unset.
        1 = if set, left is left (child) node, else is
                     in-order predecessor
        2 = if set, right is right (child), else is 
                     in-order successor
        where in-order is left subtree, root, right subtree
        */
        Integer v = threadMap.get(node);
        int vsets = (v == null) ? 0 : v.intValue();
        
        if (rightIsAChild) {
            vsets |= (1 << 2);
        } else {
            // unset
            vsets &= ~(1 << 2);
        }
        threadMap.put(node, Integer.valueOf(vsets));
    }
    
    private boolean hasALeftChild(T node) {
        Integer v = threadMap.get(node);
        assert(v != null);
        return (v.intValue() & (1 << 1)) != 0;
    }
    private boolean hasALeftChild(Integer v) {
        assert(v != null);
        return (v.intValue() & (1 << 1)) != 0;
    }
    private boolean hasARightChild(T node) {
        Integer v = threadMap.get(node);
        assert(v != null);
        return (v.intValue() & (1 << 2)) != 0;
    }
    private boolean hasARightChild(Integer v) {
        assert(v != null);
        return (v.intValue() & (1 << 2)) != 0;
    }
    
    //------- editing below for the threaded tree ------
    public void insert(T insertNode) {
        tblProbe(insertNode);               
    }
    
    @SuppressWarnings({"unchecked"})
    private void tblProbe(T item) {

        //http://adtinfo.org/libavl.html/Inserting-into-a-TBST.html#256
        
        assert(item != null);
        
        threadMap.put(item, Integer.valueOf(0));
               
        // insertion point
        T p = root;
        // new node
        T nn = null;
        int dir = 0;
        
        if (root != null) {
            
            while (true) {
                if (item.getKey() == p.getKey()) {
                    dir = 0;
                    break;
                } else if (item.getKey() > p.getKey()) {
                    dir = 1;
                } else {
                    dir = 0;
                }
                if (dir == 0) {
                    if (!hasALeftChild(p)) {
                        break;
                    }
                    p = (T)p.getLeft();
                } else {
                    if (!hasARightChild(p)) {
                        break;
                    }
                    p = (T)p.getRight();
                }
            }
        } else {
            p = root;
            dir = 0;
        }
       
        //TODO: check this
        boolean isRoot = false;
        if (p == null) {
            root = item;
            p = root;
            isRoot = true;
            // set item.right to null if default ever changes
        }
       
        setChildLinksState(item, false, false);
        
        if (dir == 0) {
            item.setLeft(p.getLeft());
            if (!isRoot) {
                setLeftIsAChild(p, true);
                item.setRight(p);
            }
            p.setLeft(item);
        } else {
            item.setRight(p.getRight());
            if (!isRoot) {
                setRightIsAChild(p, true);
                item.setLeft(p);
            }
            p.setRight(item);
        }
        
        n++;
        
        assert(threadMap.size() == n);
    }
    
    @SuppressWarnings({"unchecked"})
    public T search(T item) {
        
        if (root == null) {
            return null;
        }
        if (item == null) {
            return null;
        }
        
        T current = root;
        
        while (true) {
            
            int dir = 0;
            if (current.getKey() == item.getKey()) {
                return current;
            } 
            if (item.getKey() > current.getKey()) {
                dir = 1;
            }
            if (dir == 0) {
                if (hasALeftChild(current)) {
                    current = (T)current.getLeft();
                } else {
                    return null;
                }
            } else {
                if (hasARightChild(current)) {
                    current = (T)current.getRight();
                } else {
                    return null;
                }
            }
        }
    }
    
    /**
     * finds node with equal key and its parent if there is one
     * as []{foundNode, parentOfFoundNode}, else null if not found
     * @param node
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    private T[] searchForNodeAndParent(T item, int[] outDir) {
        
        if (root == null) {
            return null;
        }
        
        T parent = null;
        T current = root;
        outDir[0] = 0;
        
        // item is the node to delete
        // p is the found node to delete
        // q is the parent
        
        while (true) {
    
            //TODO: might want to consider changes to this
            // if need to return the identical node
            if (current.getKey() == item.getKey()) {
                break;
            } 
            
            if (item.getKey() > current.getKey()) {
                outDir[0] = 1;
            } else {
                outDir[0] = 0;
            }            
            
            parent = current;
            
            if (outDir[0] == 0) {
                if (!hasALeftChild(current)) {
                    return null;
                }
                current = (T)current.getLeft();
            } else {
                if (!hasARightChild(current)) {
                    return null;
                }
                current = (T)current.getRight();
            }
        }
            
        T[] results = (T[]) Array.newInstance(T.getType(), 2);
        results[0] = current;
        results[1] = parent;
        return results;
    }
    
    @SuppressWarnings({"unchecked"})
    public T delete(T node) {
        
        if (root == null) {
            return null;
        }
        
        Integer v = threadMap.get(node);
        
        if (v == null) {
            return null;
        }
        
        int[] dirA = new int[1];
                
        T[] foundAndParent = searchForNodeAndParent(node, dirA);
        
        if (foundAndParent == null) {
            return null;
        }
        
        T current = foundAndParent[0];
        
        v = threadMap.get(current);
        
        // this is null for root
        T parent = foundAndParent[1];
            
        if (!hasARightChild(v)) {
            if (hasALeftChild(v)) {
                //http://adtinfo.org/libavl.html/Deleting-from-a-TBST.html#260
                // case 1 of Right thread and Left child
                
                // p is node to delete
                // q is parent
                T t = (T)current.getLeft();
                while ((t != null) && hasARightChild(t)) {
                    if (t.getRight() != null) {
                        t = (T)t.getRight();
                    } else {
                        break;
                    }
                }
                // current.right does not have a child link
                t.setRight(current.getRight());
                if (dirA[0] == 0) {
                    parent.setLeft(current.getLeft());
                } else {
                    parent.setRight(current.getLeft());
                }
            } else {
                //http://adtinfo.org/libavl.html/Deleting-from-a-TBST.html#260
                // case 2 of Right thread and Left thread,
                //   i.e. the node to delete is a leaf
                
                // p is node to delete
                // q is parent
                                
                if (dirA[0] == 0) {
                    parent.setLeft(current.getLeft());
                    //TODO: revisit this:
                    if (!parent.equals(root)) {
                        setLeftIsAChild(parent, false);
                    }
                } else {
                    parent.setRight(current.getRight());
                    //TODO: revisit this
                    if (!parent.equals(root)) {
                        setRightIsAChild(parent, false);
                    }
                }
            }
        } else {
            // current has a right child
            T r = (T) current.getRight();
            if (!hasALeftChild(r)) {
                //http://adtinfo.org/libavl.html/Deleting-from-a-TBST.html#260
                // case 3 current right child has a left thread
                r.setLeft(current.getLeft());
                boolean hasLeftChild = hasALeftChild(current);
                setLeftIsAChild(r, hasLeftChild);
                if (hasLeftChild) {
                    T t = (T) r.getLeft();
                    while (hasARightChild(t)) {
                        t = (T)t.getRight();
                    }
                    t.setRight(r);
                }
                
                if (parent == null) {
                    //NLK: added to handle root
                    root = r;
                } else {
                    if (dirA[0] == 0) {
                        parent.setLeft(r);
                    } else {
                        parent.setRight(r);
                    }
                }
            } else {
                //case 4:  current's right child has a left child
                T s = null;
                while (true) {
                    s = (T)r.getLeft();
                    if (!hasALeftChild(s)) {
                        break;
                    }
                    r = s;
                }
            }
        }
        
        threadMap.remove(current);
                
        n--;
        
        assert(threadMap.size() == n);
        
        return current;
    }
    
    //TODO:  create a balance method:
    //http://adtinfo.org/libavl.html/Balancing-a-TBST.html
    
    public T minimum() {
        return minimum(root);
    }
    
    @SuppressWarnings({"unchecked"})
    protected T minimum(T x) {
        if (x == null) {
            return null;
        }
        T nd = x;
        while (hasALeftChild(nd)) {
            nd = (T) nd.getLeft();
        }
        return nd;
    }
    
    public T maximum() {
        return maximum(root);
    }
    
    @SuppressWarnings({"unchecked"})
    protected T maximum(T x) {
        if (x == null) {
            return null;
        }
        T nd = x;
        while (hasARightChild(nd)) {
            nd = (T) nd.getRight();
        }
        return nd;
    }
    
    /**
     * smallest element in the tree with key greater
     * than node.key.
     * @param node
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    public T successor(T node) {
        
        if (hasARightChild(node)) {
            return minimum((T)node.getRight());
        }
        // node has a right thread
        T s = (node.getRight() != null) ?
            (T) node.getRight() : null;
        if (s == null) {
            // node is already the largest in the tree
            return node;
        }
        return s;
    }
    
    /*
     * the largest element in the tree with key smaller 
     * than node.key.
     * 
     * @param node
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    public T predecessor(T node) {

        if (hasALeftChild(node)) {
            return maximum((T)node.getLeft());
        }
        // node has a left thread
        T p = (node.getLeft() != null) ?
            (T) node.getLeft() : null;
        if (p == null) {
            // node is already the smallest in the tree
            return node;
        }
        return p;
    }
  
    public int getNumberOfNodes() {
        return n;
    }
    
    @SuppressWarnings({"unchecked"})
    public T first(T node) {
 
        if (node == null || root == null) {
            return null;
        }
        
        while (hasALeftChild(node)) {
            node = (T)node.getLeft();
        }
        
        return node;
    }
    
    @SuppressWarnings({"unchecked"})
    public T next (T node) {
  
        if (node == null) {
            return first(node);
        } else if (!hasARightChild(node)) {
            node = (T)node.getRight();
            return node;
        } else {
            node = (node.getRight() != null) ?
                (T)node.getRight() : null;
            while ((node != null) && hasALeftChild(node)) {
                node = (T)node.getLeft();
            }
            return node;
        }
    }
    
}
