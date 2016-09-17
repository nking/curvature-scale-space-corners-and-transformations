package thirdparty.edu.princeton.cs.algs4;

/**
 * from Interval.java in algs4.jar
 * http://algs4.cs.princeton.edu/
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
   (see bottom of this file)
 * http://algs4.cs.princeton.edu/92search/QuadTree.java.html
 * copyright for authors Robert Sedgewick and Kevin Wayne
 * is GPLV3, http://algs4.cs.princeton.edu/faq/
 */
public class Interval<T extends Comparable<T>> implements 
    Comparable<Interval<T>> {
    
    private final T min;    // min endpoint
    private final T max;    // max endpoint

    public Interval(T min, T max) {
        if (less(max, min)) {
            throw new RuntimeException(
            "Illegal argument: max=" + max + " min=" + min);
        }
        this.min = min;
        this.max = max;
    }

    // return min endpoint
    public T min() {
        return min;
    }

    // return max endpoint
    public T max() {
        return max;
    }

    // is x between min and max
    public boolean contains(T x) {
        return !less(x, min) && !less(max, x);
    }

    // does this interval a intersect interval b?
    public boolean intersects(Interval<T> b) {
        Interval<T> a  = this;
        if (less(a.max, b.min)) return false;
        if (less(b.max, a.min)) return false;
        return true;
    }

    // does this interval a equal interval b?
    public boolean equals(Interval<T> b) {
        Interval<T> a  = this;
        return a.min.equals(b.min) && a.max.equals(b.max);
    }


    // comparison helper functions
    private boolean less(T x, T y) {
        return x.compareTo(y) < 0;
    }

    @Override
    public int compareTo(Interval<T> other) {
        
        //other.max   min  max
        int c1 = min().compareTo(other.max());
        if (c1 > 0) {
            return -c1;
        }
        //    max  other.min
        c1 = max().compareTo(other.min());
        if (c1 < 0) {
            return -c1;
        }
        
        return 0;
    }
    
    // return string representation
    public String toString() {
        return "[" + min + ", " + max + "]";
    }

    // test client
    public static void main(String[] args) {
        int n = Integer.parseInt(args[0]);

        Interval<Integer> a = new Interval<Integer>(5, 17);
        Interval<Integer> b = new Interval<Integer>(5, 17);
        Interval<Integer> c = new Interval<Integer>(5, 18);
        System.out.println(a.equals(b));
        System.out.println(!a.equals(c));
        System.out.println(!b.equals(c));


        // generate n random points in [-1, 2] and compute
        // fraction that lies in [0, 1]
        Interval<Double> interval = new Interval<Double>(0.0, 1.0);
        int count = 0;
        for (int i = 0; i < n; i++) {
            double x = 3 * Math.random() - 1.0;
            if (interval.contains(x)) count++;
        }
        System.out.println("fraction = " + (1.0 * count / n));
    }

}
