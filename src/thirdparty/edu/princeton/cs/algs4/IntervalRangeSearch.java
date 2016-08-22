package thirdparty.edu.princeton.cs.algs4;

/******************************************************************************
 A specialization of RangeSearch for an Interval parameterized type.
  
 NOTE that the intervals cannot overlap.  If a put of an interval
 intersects with a key, the existing interval in the tree gets
 the value of the new interval, but the key range does not
 change.
 
 @author nichole
 
 * @param <T> the data type used in the Intervals
 * @param <Value> the data type of the key associated with each
 * tree interval.
 ******************************************************************************/
public class IntervalRangeSearch<T extends Comparable<T>, Value> extends
    RangeSearch<Interval<T>, Value> {
    
    public Queue<Interval<T>> range0(Interval<T> interval) {
        Queue<Interval<T>> list = new Queue<Interval<T>>();
        //System.out.println("root=" + root);
        //System.out.println("srch interval=" + interval);
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
