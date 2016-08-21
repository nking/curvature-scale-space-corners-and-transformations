package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
 *  Compilation:  javac RangeSearch.java
 *  Execution:    java RangeSearch < words.txt
 *  
 *  Range search implemented using a randomized BST.
 *  
 *  % java RangeSearch < words.txt
 *  height:          33
 *  size:            20068
 *  min key:         a
 *  max key:         zygote
 *  integrity check: true
 *
 * [kevin, kfg]
 *  key
 *  keyboard
 *  keyed
 *  keyhole
 *  keynote
 *  keypunch
 *  keys
 *  keystone
 *  keyword
 *
 *  [paste, pasty]
 *  paste
 *  pasteboard
 *  pastel
 *  pasteup
 *  pastiche
 *  pastime
 *  pastor
 *  pastoral
 *  pastry
 *  pasture
 *  pasty

    class downloaded from http://algs4.cs.princeton.edu/92search/
    copyright for authors Robert Sedgewick and Kevin Wayne
    is GPLV3, http://algs4.cs.princeton.edu/faq/
 
 ******************************************************************************/

public class RangeSearch<Key extends Comparable<Key>, Value>  {

    private Node root;   // root of the BST

    // BST helper node data type
    private class Node {
        Key key;            // key
        Value val;          // associated data
        Node left, right;   // left and right subtrees
        int N;              // node count of descendents

        public Node(Key key, Value val) {
            this.key = key;
            this.val = val;
            this.N   = 1;
        }
    }

   /***************************************************************************
    *  BST search
    ***************************************************************************/

    public boolean contains(Key key) {
        return (get(key) != null);
    }

    // return value associated with the given key
    // if no such value, return null
    public Value get(Key key) {
        return get(root, key);
    }

    private Value get(Node x, Key key) {
        if (x == null) return null;
        int cmp = key.compareTo(x.key);
        if      (cmp == 0) return x.val;
        else if (cmp  < 0) return get(x.left,  key);
        else               return get(x.right, key);
    }

   /***************************************************************************
    *  randomized insertion
    ***************************************************************************/
    public void put(Key key, Value val) {
        root = put(root, key, val);
    }

    // make new node the root with uniform probability
    private Node put(Node x, Key key, Value val) {
        if (x == null) return new Node(key, val);
        int cmp = key.compareTo(x.key);
        if (cmp == 0) { x.val = val; return x; }
        if (StdRandom.bernoulli(1.0 / (size(x) + 1.0))) return putRoot(x, key, val);
        if (cmp < 0) x.left  = put(x.left,  key, val); 
        else         x.right = put(x.right, key, val); 
        // (x.N)++;
        fix(x);
        return x;
    }


    private Node putRoot(Node x, Key key, Value val) {
        if (x == null) return new Node(key, val);
        int cmp = key.compareTo(x.key);
        if      (cmp == 0) { x.val = val; return x; }
        else if (cmp  < 0) { x.left  = putRoot(x.left,  key, val); x = rotR(x); }
        else               { x.right = putRoot(x.right, key, val); x = rotL(x); }
        return x;
    }




   /***************************************************************************
    *  deletion
    ***************************************************************************/
    private Node joinLR(Node a, Node b) { 
        if (a == null) return b;
        if (b == null) return a;

        if (StdRandom.bernoulli((double) size(a) / (size(a) + size(b))))  {
            a.right = joinLR(a.right, b);
            fix(a);
            return a;
        }
        else {
            b.left = joinLR(a, b.left);
            fix(b);
            return b;
        }
    }

    private Node remove(Node x, Key key) {
        if (x == null) return null; 
        int cmp = key.compareTo(x.key);
        if      (cmp == 0) x = joinLR(x.left, x.right);
        else if (cmp  < 0) x.left  = remove(x.left,  key);
        else               x.right = remove(x.right, key);
        fix(x);
        return x;
    }

    // remove and return value associated with given key; if no such key, return null
    public Value remove(Key key) {
        Value val = get(key);
        root = remove(root, key);
        return val;
    }




   /***************************************************************************
    *  Range searching
    ***************************************************************************/


    // return all keys in given interval
    public Iterable<Key> range(Key min, Key max) {
        return range(new Interval<Key>(min, max));
    }
    public Iterable<Key> range(Interval<Key> interval) { 
        Queue<Key> list = new Queue<Key>();
        range(root, interval, list);
        return list;
    }
    private void range(Node x, Interval<Key> interval, Queue<Key> list) {
        if (x == null) return;
        if (!less(x.key, interval.min()))  range(x.left, interval, list);
        if (interval.contains(x.key))      list.enqueue(x.key);
        if (!less(interval.max(), x.key))  range(x.right, interval, list);
    }



   /***************************************************************************
    *  Utility functions
    ***************************************************************************/

    // return the smallest key
    public Key min() {
        Key key = null;
        for (Node x = root; x != null; x = x.left)
            key = x.key;
        return key;
    }
    
    // return the largest key
    public Key max() {
        Key key = null;
        for (Node x = root; x != null; x = x.right)
            key = x.key;
        return key;
    }


   /***************************************************************************
    *  useful binary tree functions
    ***************************************************************************/

    // return number of nodes in subtree rooted at x
    public int size() { return size(root); }
    private int size(Node x) { 
        if (x == null) return 0;
        else           return x.N;
    }

    // height of tree (empty tree height = 0)
    public int height() { return height(root); }
    private int height(Node x) {
        if (x == null) return 0;
        return 1 + Math.max(height(x.left), height(x.right));
    }


   /***************************************************************************
    *  helper BST functions
    ***************************************************************************/

    // fix subtree count field
    private void fix(Node x) {
        if (x == null) return;                 // check needed for remove
        x.N = 1 + size(x.left) + size(x.right);
    }

    // right rotate
    private Node rotR(Node h) {
        Node x = h.left;
        h.left = x.right;
        x.right = h;
        fix(h);
        fix(x);
        return x;
    }

    // left rotate
    private Node rotL(Node h) {
        Node x = h.right;
        h.right = x.left;
        x.left = h;
        fix(h);
        fix(x);
        return x;
    }


   /***************************************************************************
    *  Debugging functions that test the integrity of the tree
    ***************************************************************************/

    // check integrity of subtree count fields
    public boolean check() { return checkCount() && isBST(); }

    // check integrity of count fields
    private boolean checkCount() { return checkCount(root); }
    private boolean checkCount(Node x) {
        if (x == null) return true;
        return checkCount(x.left) && checkCount(x.right) && (x.N == 1 + size(x.left) + size(x.right));
    }


    // does this tree satisfy the BST property?
    private boolean isBST() { return isBST(root, min(), max()); }

    // are all the values in the BST rooted at x between min and max, and recursively?
    private boolean isBST(Node x, Key min, Key max) {
        if (x == null) return true;
        if (less(x.key, min) || less(max, x.key)) return false;
        return isBST(x.left, min, x.key) && isBST(x.right, x.key, max);
    } 



   /***************************************************************************
    *  helper comparison functions
    ***************************************************************************/

    private boolean less(Key k1, Key k2) {
        return k1.compareTo(k2) < 0;
    }


   /***************************************************************************
    *  test client
    ***************************************************************************/
    public static void main(String[] args) {
        /*
        int N = 0;
        RangeSearch<String, Integer> st = new RangeSearch<String, Integer>();
        while (!StdIn.isEmpty()) {
            String s = StdIn.readString();
            st.put(s, N++);
        }

        StdOut.println("height:          " + st.height());
        StdOut.println("size:            " + st.size());
        StdOut.println("min key:         " + st.min());
        StdOut.println("max key:         " + st.max());
        StdOut.println("integrity check: " + st.check());
        StdOut.println();

        StdOut.println(new Interval<String>("kevin", "kfg"));
        Iterable<String> list = st.range(new Interval<String>("kevin", "kfg"));
        for (String s : list)
            StdOut.println(s + " " + st.get(s));
        StdOut.println();

        StdOut.println(new Interval<String>("paste", "pasty"));
        list = st.range(new Interval<String>("paste", "pasty"));
        for (String s : list)
            StdOut.println(s);
        StdOut.println();
        */
    }

}
