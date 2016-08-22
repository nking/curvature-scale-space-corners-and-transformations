package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
 Adapted from RangeSearch, specialized for Interval parameterization.
  
 * Compilation:  javac RangeSearch.java
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

public class IntervalRangeSearch<T extends Comparable<T>, Value> extends
    RangeSearch<Interval<T>, Value> {
    
    public Queue<Interval<T>> range0(Interval<T> interval) {
        Queue<Interval<T>> list = new Queue<Interval<T>>();
        System.out.println("root=" + root);
        System.out.println("srch interval=" + interval);
        range0(root, interval, list);
        return list;
    }
    
    private void range0(RangeSearchNode<Interval<T>, Value> x, 
        Interval<T> interval, Queue<Interval<T>> list) {
       
        if (x == null) return;
       
        boolean intersects = interval.intersects(x.key);
        if (intersects) {
            list.enqueue(x.key);
        }
        
        /*
        interval has min max for search.
        
        tree has left as the larger keys
        
                     x
              lft         rgt
           lft  rgt     lft  rgt
        
        
        or viewed by increasing values--->
               xmin--x--xmax smin
              rgt         lft
           rgt  lft     rgt  lft
        */
          
        // if x.max < interval.min
        //   search lft
        // if xmin > interval.max
        //   search rgt
        
        if ((x.left != null) && (intersects || 
            (x.key.max().compareTo(interval.min()) < 1)) ) {
            range0(x.left, interval, list);
        }
        
        if ((x.right != null) && (intersects || 
            (interval.max().compareTo(x.key.min()) < 1))) {
            range0(x.right, interval, list);
        }
    }

}
