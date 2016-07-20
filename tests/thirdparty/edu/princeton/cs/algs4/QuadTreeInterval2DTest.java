package thirdparty.edu.princeton.cs.algs4;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class QuadTreeInterval2DTest extends TestCase {
    
    public QuadTreeInterval2DTest() {
    }
    
    public void test0() {
        
        /*
          7
          6         & &
          5         & &
          4
          3
            -5      0 1    5
        */
        
        Interval<Double> intervalX1 = new Interval<Double>(0.0, 1.0);
        Interval<Double> intervalY1 = new Interval<Double>(5.0, 6.0);
        Interval2D<Double> box1 
            = new Interval2D<Double>(intervalX1, intervalY1);
        Interval<Double> intervalX2 = new Interval<Double>(-5.0, 5.0);
        Interval<Double> intervalY2 = new Interval<Double>(3.0, 7.0);
        Interval2D<Double> box2 = new Interval2D<Double>(
          intervalX2, intervalY2);
         
        QuadTreeInterval2D<Double, String> st2
            = new QuadTreeInterval2D<Double, String>();
         
        st2.insert(box1, "box1");
        
        List<Interval2D<Double>> results = st2.query2D(box2);
        
        assertEquals(1, results.size());
        
        assertTrue(results.get(0).equals(box1));
    }
    
    public void test1() {
        
        /*
          7               *
          6               *
          5               *
          4     |     |
          3     |     |
          2     |_____|      + + +
          1                  + + +
          0
            -5 -4 -3 -2 -1 0 1 2 3 4 5
        */
        Interval<Integer> x1 = new Interval<Integer>(-4, -2);
        Interval<Integer> y1 = new Interval<Integer>(2, 4);
        Interval2D<Integer> box1 = new Interval2D<Integer>(x1, y1);

        Interval<Integer> x2 = new Interval<Integer>(0, 0);
        Interval<Integer> y2 = new Interval<Integer>(5, 7);
        Interval2D<Integer> box2 = new Interval2D<Integer>(x2, y2);
        
        Interval<Integer> x3 = new Interval<Integer>(1, 3);
        Interval<Integer> y3 = new Interval<Integer>(1, 2);
        Interval2D<Integer> box3 = new Interval2D<Integer>(x3, y3);

        Interval<Integer> srchX = new Interval<Integer>(-2, 2);
        Interval<Integer> srchY = new Interval<Integer>(2, 5);
        Interval2D<Integer> srch = new Interval2D<Integer>(srchX, srchY);

        QuadTreeInterval2D<Integer, String> st2
            = new QuadTreeInterval2D<Integer, String>();
         
        st2.insert(box1, "box1");
        st2.insert(box2, "box2");
        st2.insert(box3, "box3");
        
        List<Interval2D<Integer>> results = st2.query2D(srch);
        
        assertEquals(3, results.size());
        
        Set<Interval2D<Integer>> expected = new
            HashSet<Interval2D<Integer>>();
        expected.add(box1);
        expected.add(box2);
        expected.add(box3);
        
        for (Interval2D<Integer> ii : results) {
            assertTrue(expected.remove(ii));
        }
        assertTrue(expected.isEmpty());
    }
}
