package algorithms;

import algorithms.imageProcessing.HeapNode;
import java.lang.reflect.Array;
import java.util.Stack;

/**
   Binary Search Tree – has the property that for every node:
       The left subtree contains only nodes that have a key less
           than or equal to the node
       The right subtree contains only nodes that have a key
           greater than the node
   The left and right subtrees are also binary search trees

 * use a Node with left, right and parent.
 * use a binary search method impl as: while loop rot != null: if (root.value < value) root=root.getRight()...
 *
 * based upon pseudocode from Cormen et al. Introduction to Algorithms.  added in-order, pre-order and post-order traversal methods too.
 *
 * @author nichole
 */
public class BinarySearchTree<T extends HeapNode> {
    
    protected int n = 0;
    
    protected T root = null;
  
    public T minimum() {
        return minimum(root);
    }
    
    protected T minimum(T x) {
        if (x == null) {
            return null;
        }
        T nd = x;
        while (nd.getLeft() != null) {
            nd = (T) nd.getLeft();
        }
        return nd;
    }
    
    public T maximum() {
        return maximum(root);
    }
    
    protected T maximum(T x) {
        if (x == null) {
            return null;
        }
        T nd = x;
        while (nd.getRight() != null) {
            nd = (T) nd.getRight();
        }
        return nd;
    }
    
    /**
     * smallest element in the tree with key greater
     * than x.key.
     * @param x
     * @return 
     */
    public T successor(T x) {
        if (x.getRight() != null) {
            return minimum((T)x.getRight());
        }
        T y = (x.getParent() != null) ? (T) x.getParent() : null;
        while ((y != null) && 
            (y.getRight() != null && x.equals(y.getRight()))) {
            x = y;
            y = (y.getParent() != null) ? (T) y.getParent() : null;
        }
        return y;
    }
    
    /*
     * the largest element in the tree with key smaller 
     * than x.key.
     * 
     * a smaller number can either be a descendant to the left
     * or a parent or the left descendants of a parent.
     * 
     * This is not an efficient data structure for this operation
     * so will not implement it for this.
     * 
     * @param x
     * @return 
     *
    public T predecessor(T x) {
        the smallest smaller than x could be the same
        value key, so if found that, could stop the
        search early.
    }
    */
    
    public void insert(T z) {
        
        if (root == null) {
            root = z;
            n++;
            return;
        }
        
        T zParent = null;
        T x = root;
        
        while (x != null) {
            zParent = x;
            if (z.getKey() < x.getKey()) {
                x = (x.getLeft() != null) ? (T) x.getLeft() : null;
            } else if (z.getKey() == x.getKey()) {
                zParent = z;
                break;
            } else {
                x = (x.getRight() != null) ? (T) x.getRight() : null;
            }
        }
        z.setParent(zParent); 
        if (zParent == null) {
            root = z;
        } else if (z.getKey() < zParent.getKey()) {
            zParent.setLeft(z);           
        } else {
            zParent.setRight(z);
        }
        
        n++;
    }
    
    public void delete(T z) {
       
        T y = null;
        T x = null;
        if ((z.getLeft() == null) || (z.getRight() == null)) {
            y = z;
        } else {
            y = successor(z);
        }
        
        if (y.getLeft() != null) {
            x = (T) y.getLeft();
        } else {
            x = (y.getRight() != null) ? (T) y.getRight() : null;
        }
        
        if (x != null) {
            x.setParent(y.getParent());
        }
         
        //if (y.getParent() == null) {
        if (root.equals(z)) {
            root = x;
        } else if (y == y.getParent().getLeft()) {
            y.getParent().setLeft(x);
        } else {
            y.getParent().setRight(x);
        }
        
        if (y != z) {
            // copy y data into z
            z.setKey(y.getKey());
            z.setData(y.getData());
        }
        
        n--;
    }
    
    /**
     * runtime complexity is O(h) where h is height of
     * tree, which is usually lg_2(n)
     * @param z
     * @return 
     */
    public T search(T z) {
        if (root == null) {
            return null;
        }
        return search(root, z);
    }
    
    private T search(T tn, T z) {
        return search(tn, z.getKey());
    }
    
    private T search(T tn, long theKey) {
        while ((tn != null) && (tn.getKey() != theKey)) {
            if (theKey < tn.getKey()) {
                tn = (tn.getLeft() != null) ? (T) tn.getLeft() : null;
            } else if (theKey > tn.getKey()) {
                tn = (tn.getRight() != null) ? (T) tn.getRight() : null;
            }
        }
        return tn;
    }
    
    private T getRoot(T nd) {
        while ((nd != null) && (nd.getParent() != null)) {
            nd = (T) nd.getParent();
        }
        return nd;
    }
    
    public int getNumberOfNodes() {
        return n;
    }
    
    /**
     * visit each node using pattern left subtree, root, right subtree
     * in an iterative manner rather than invoking the method recursively.
     */
    protected T[] getInOrderTraversalIterative(T node) {
       
        Class cls = T.getType();
        T[] array = (T[]) Array.newInstance(cls, n);
        int count = 0;
        
        Stack<T> stack = new Stack<>();
               
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                
                stack.push(node);
                
                node = (node.getLeft() != null) ? (T) node.getLeft() : null;
            
            } else {
                
                node = stack.pop();
                
                array[count] = node;
                count++;
                
                //System.out.println(node.key);
                
                node = (node.getRight() != null) ? (T) node.getRight() : null;
            }
        }
        
        return array;
    }
    
    /**
     * visit each node using pattern: root, left subtree, right subtree
     * in an iterative manner rather than invoking the method recursively.
     */
    protected T[] printPreOrderTraversalIterative(T node) {
       
        Class cls = T.getType();
        T[] array = (T[]) Array.newInstance(cls, n);
        int count = 0;
        
        Stack<T> stack = new Stack<>();
        
        while (!stack.isEmpty() || (node != null)) {
            if (node != null) {
                
                array[count] = node;
                count++;
                //System.out.println(node.key);
                
                stack.push(node);
                
                node = (node.getLeft() != null) ? (T)node.getLeft() : null;
            
            } else {
                
                node = stack.pop();
                
                node = (node.getRight() != null) ? (T)node.getRight() : null;
            }
        }
        
        return array;
    }
    
    /**
     * visit each node using pattern: left subtree, right subtree, root subtree
     * in an iterative manner rather than invoking the method recursively.
     */
    protected T[] printPostOrderTraversalIterative(T node) {
        
        Class cls = T.getType();
        T[] array = (T[]) Array.newInstance(cls, n);
        int count = 0;
        
        if (node == null) {
            return array;
        }
        
        Stack<T> stack = new Stack<>();
        Stack<T> stack2 = new Stack<T>();
        stack.push(node);
        
        while (!stack.isEmpty()) {
            
            node = stack.pop();
            
            stack2.push(node);
            
            if (node.getLeft() != null) {
                stack.push((T)node.getLeft());
            }

            if (node.getRight() != null) {
                stack.push((T)node.getRight());
            }            
        }
        
        while (!stack2.isEmpty()) {
            
            node = stack2.pop();
            
            //process(node);
            array[count] = node;
            count++;
            //System.out.println(node.key);
        }
         
        return array;
    }
}
