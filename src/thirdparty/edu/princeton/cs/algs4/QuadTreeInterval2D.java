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
        
        int cX = h.xy.intervalX.compareTo(box.intervalX);
        int cY = h.xy.intervalY.compareTo(box.intervalY);
        
        if ((cX < 0) && (cY < 0)) { 
            h.SW = insert(h.SW, box, value);
        } else if ((cX < 0) && (cY >= 0)) {
            h.NW = insert(h.NW, box, value);
        } else if ((cX >= 0) && (cY < 0)) {
            h.SE = insert(h.SE, box, value);
        } else if ((cX >= 0) && (cY >= 0)) { 
            h.NE = insert(h.NE, box, value);
        }
        
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

    private void query2D(Node<T> h, Interval2D<T> srch,
        List<Interval2D<T>> output) {
        
        if (h == null) return;
        
        int cX = srch.intervalX.compareTo(h.xy.intervalX);
        int cY = srch.intervalY.compareTo(h.xy.intervalY);

        if ((cX == 0) && (cY == 0)) {
            output.add(h.xy);
        }
        
        if (h.SW != null && (cX <= 0) && (cY <= 0)) 
            query2D(h.SW, srch, output);
        if (h.NW != null && (cX <= 0) && (cY >= 0)) 
            query2D(h.NW, srch, output);
        if (h.SE != null && (cX >= 0) && (cY <= 0)) 
            query2D(h.SE, srch, output);
        if (h.NE != null && (cX >= 0) && (cY >= 0)) 
            query2D(h.NE, srch, output);    
    }

}
