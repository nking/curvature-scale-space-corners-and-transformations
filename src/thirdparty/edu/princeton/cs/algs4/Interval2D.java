/*
 *  2-dimensional interval data type.
 * from Interval2D.java in algs4.jar
 * http://algs4.cs.princeton.edu/
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
   (see bottom of this file)
Note this version of the file was copied from
http://algs4.cs.princeton.edu/91primitives/Interval2D.java.html
as it has the parameterization that the jar file did not
have.
 */
package thirdparty.edu.princeton.cs.algs4;

/**
 *  The <tt>Interval2D</tt> class represents a closed two-dimensional interval,
 *  which represents all points (x, y) with both xmin <= x <= xmax and
 *  ymin <= y <= ymax.
 *  Two-dimensional intervals are immutable: their values cannot be changed
 *  after they are created.
 *  The class <code>Interval2D</code> includes methods for checking whether
 *  a two-dimensional interval contains a point and determining whether
 *  two two-dimensional intervals intersect.
 *  <p>
 *  For additional documentation, 
 *  see <a href="http://algs4.cs.princeton.edu/12oop">Section 1.2</a> of 
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne. 
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Interval2D<T extends Comparable<T>> 
    implements Comparable<Interval2D<T>> { 
    
    public final Interval<T> intervalX;   // x-interval
    public final Interval<T> intervalY;   // y-interval
   
    public Interval2D(Interval<T> intervalX, Interval<T> intervalY) {
        this.intervalX = intervalX;
        this.intervalY = intervalY;
    }

    // does this 2D interval a intersect b?
    public boolean intersects(Interval2D<T> b) {
        if (intervalX.intersects(b.intervalX)) return true;
        if (intervalY.intersects(b.intervalY)) return true;
        return false;
    }

    // does this 2D interval contain (x, y)?
    public boolean contains(T x, T y) {
        return intervalX.contains(x) && intervalY.contains(y);
    }

    @Override
    public int compareTo(Interval2D<T> other) {
       
        //other.max   min  max
        int c1 = intervalX.min().compareTo(other.intervalX.max());
        if (c1 > 0) {
            return -c1;
        }
        //    max  other.min
        c1 = intervalX.max().compareTo(other.intervalX.min());
        if (c1 < 0) {
            return -c1;
        }
        
        c1 = intervalY.min().compareTo(other.intervalY.max());
        if (c1 > 0) {
            return -c1;
        }
        //    max  other.min
        c1 = intervalY.max().compareTo(other.intervalY.min());
        if (c1 < 0) {
            return -c1;
        }
        
        return 0;
    }
    
    // return string representation
    public String toString() {
        return intervalX + " x " + intervalY;
    }

    // test client
    public static void main(String[] args) {
        Interval<Double> intervalX = new Interval<Double>(0.0, 1.0);
        Interval<Double> intervalY = new Interval<Double>(5.0, 6.0);
        Interval2D<Double> box1 = new Interval2D<Double>(intervalX, intervalY);
        intervalX = new Interval<Double>(-5.0, 5.0);
        intervalY = new Interval<Double>(3.0, 7.0);
        Interval2D<Double> box2 = new Interval2D<Double>(intervalX, intervalY);
        System.out.println("box1 = " + box1);
        System.out.println("box2 = " + box2);
        System.out.println(box1.contains(0.5, 5.5));
        System.out.println(!box1.contains(1.5, 5.5));
        System.out.println(box1.intersects(box2));
        System.out.println(box2.intersects(box1));
    }

}