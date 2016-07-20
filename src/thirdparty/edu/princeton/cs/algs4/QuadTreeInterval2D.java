package thirdparty.edu.princeton.cs.algs4;

import java.util.ArrayList;
import java.util.List;

/**
 * adapted from QuadTree implementation from Sedgewick and Wayne
 * from 
 * http://algs4.cs.princeton.edu/92search/QuadTree.java.html
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
 */
public class QuadTreeInterval2D<T extends Comparable<T>, Value>  {
                 //Interval2D type T
    
    private Node<T> root;

    // helper node data type
    private class Node<T extends Comparable<T>> {
        Interval2D<T> xy;
        Node<T> NW, NE, SE, SW;   // four subtrees
        Value value;           // associated data

        Node(Interval2D<T> box, Value value) {
            this.xy = box;
            this.value = value;
        }
    }


  /***********************************************************************
    *  Insert (x, y) into appropriate quadrant
    ***************************************************************************/
    public void insert(Interval2D<T> box, Value value) {
        root = insert(root, box, value);
    }

    private Node<T> insert(Node<T> h, Interval2D<T> box, Value value) {
        if (h == null) {
            return new Node<T>(box, value);
        }
        //// if (eq(x, h.x) && eq(y, h.y)) h.value = value;  // duplicate
        else if ( less(box.intervalX, h.xy.intervalX) 
            &&  less(box.intervalY, h.xy.intervalY)) 
            h.SW = insert(h.SW, box, value);
        else if ( less(box.intervalX, h.xy.intervalX) 
            && !less(box.intervalY, h.xy.intervalY)) 
            h.NW = insert(h.NW, box, value);
        else if (!less(box.intervalX, h.xy.intervalX) 
            &&  less(box.intervalY, h.xy.intervalY)) 
            h.SE = insert(h.SE, box, value);
        else if (!less(box.intervalX, h.xy.intervalX) 
            && !less(box.intervalY, h.xy.intervalY)) 
            h.NE = insert(h.NE, box, value);
        return h;
    }


  /***********************************************************************
    *  Range search.
    ***************************************************************************/

    public List<Interval2D<T>> query2D(Interval2D<T> rect) {
        
        List<Interval2D<T>> output = new ArrayList<Interval2D<T>>();
        
        query2D(root, rect, output);
        
        return output;
    }

    private void query2D(Node<T> h, Interval2D<T> rect,
        List<Interval2D<T>> output) {
        if (h == null) return;
        T xmin = rect.intervalX.min();
        T ymin = rect.intervalY.min();
        T xmax = rect.intervalX.max();
        T ymax = rect.intervalY.max();
        
        T xminH = h.xy.intervalX.min();
        T yminH = h.xy.intervalY.min();
        T xmaxH = h.xy.intervalX.max();
        T ymaxH = h.xy.intervalY.max();
        
        if (rect.compareTo(h.xy) == 0) {
            output.add(h.xy);
            System.out.println("    (" + h.xy + ") " + h.value);
        }
        if ( less(xmin, xminH) && less(ymin, yminH)) 
            query2D(h.SW, rect, output);
        if ( less(xmin, xminH) && !less(ymax, yminH)) 
            query2D(h.NW, rect, output);
        if (!less(xmax, xmaxH) &&  less(ymin, ymaxH)) 
            query2D(h.SE, rect, output);
        if (!less(xmax, xmaxH) && !less(ymax, ymaxH)) 
            query2D(h.NE, rect, output);
    }

   /***************************************************************************
    *  helper comparison functions
    ***************************************************************************/

    private boolean less(Interval<T> k1, Interval<T> k2) { 
        return k1.compareTo(k2) <  0; }
    private boolean eq  (Interval<T> k1, Interval<T> k2) { 
        return k1.compareTo(k2) == 0; }
    
    private boolean less(T k1, T k2) { 
        return k1.compareTo(k2) <  0; }
    private boolean eq  (T k1, T k2) { 
        return k1.compareTo(k2) == 0; }

}
