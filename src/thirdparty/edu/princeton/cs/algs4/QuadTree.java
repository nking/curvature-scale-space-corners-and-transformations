package thirdparty.edu.princeton.cs.algs4;

/**
 * a QuadTree implementation from Sedgewick and Wayne
 * from 
 * http://algs4.cs.princeton.edu/92search/QuadTree.java.html
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
 */
public class QuadTree<T extends Comparable<T>, Value>  {
    
    private Node root;

    // helper node data type
    private class Node<S> {
        S x;
        S y; // x- and y- coordinates
        Node NW, NE, SE, SW;   // four subtrees
        Value value;           // associated data

        Node(S x, S y, Value value) {
            this.x = x;
            this.y = y;
            this.value = value;
        }
    }


  /***********************************************************************
    *  Insert (x, y) into appropriate quadrant
    ***************************************************************************/
    public void insert(T x, T y, Value value) {
        root = insert(root, x, y, value);
    }

    private Node<T> insert(Node<T> h, T x, T y, Value value) {
        if (h == null) return new Node<T>(x, y, value);
        //// if (eq(x, h.x) && eq(y, h.y)) h.value = value;  // duplicate
        else if ( less(x, h.x) &&  less(y, h.y)) h.SW = insert(h.SW, x, y, value);
        else if ( less(x, h.x) && !less(y, h.y)) h.NW = insert(h.NW, x, y, value);
        else if (!less(x, h.x) &&  less(y, h.y)) h.SE = insert(h.SE, x, y, value);
        else if (!less(x, h.x) && !less(y, h.y)) h.NE = insert(h.NE, x, y, value);
        return h;
    }


  /***********************************************************************
    *  Range search.
    ***************************************************************************/

    public void query2D(Interval2D<T> rect) {
        query2D(root, rect);
    }

    private void query2D(Node<T> h, Interval2D<T> rect) {
        if (h == null) return;
        T xmin = rect.intervalX.min();
        T ymin = rect.intervalY.min();
        T xmax = rect.intervalX.max();
        T ymax = rect.intervalY.max();
        if (rect.contains(h.x, h.y))
            System.out.println("    (" + h.x + ", " + h.y + ") " + h.value);
        if ( less(xmin, h.x) &&  less(ymin, h.y)) query2D(h.SW, rect);
        if ( less(xmin, h.x) && !less(ymax, h.y)) query2D(h.NW, rect);
        if (!less(xmax, h.x) &&  less(ymin, h.y)) query2D(h.SE, rect);
        if (!less(xmax, h.x) && !less(ymax, h.y)) query2D(h.NE, rect);
    }


   /***************************************************************************
    *  helper comparison functions
    ***************************************************************************/

    private boolean less(T k1, T k2) { return k1.compareTo(k2) <  0; }
    private boolean eq  (T k1, T k2) { return k1.compareTo(k2) == 0; }


   /***************************************************************************
    *  test client
    ***************************************************************************/
    public static void main(String[] args) {
        int M = Integer.parseInt(args[0]);   // queries
        int N = Integer.parseInt(args[1]);   // points

        QuadTree<Integer, String> st = new QuadTree<Integer, String>();

        // insert N random points in the unit square
        for (int i = 0; i < N; i++) {
            Integer x = (int) (100 * Math.random());
            Integer y = (int) (100 * Math.random());
            // StdOut.println("(" + x + ", " + y + ")");
            st.insert(x, y, "P" + i);
        }
        System.out.println("Done preprocessing " + N + " points");

        // do some range searches
        for (int i = 0; i < M; i++) {
            Integer xmin = (int) (100 * Math.random());
            Integer ymin = (int) (100 * Math.random());
            Integer xmax = xmin + (int) (10 * Math.random());
            Integer ymax = ymin + (int) (20 * Math.random());
            Interval<Integer> intX = new Interval<Integer>(xmin, xmax);
            Interval<Integer> intY = new Interval<Integer>(ymin, ymax);
            Interval2D<Integer> rect = new Interval2D<Integer>(intX, intY);
            System.out.println(rect + " : ");
            st.query2D(rect);
        }
        
    }

}
