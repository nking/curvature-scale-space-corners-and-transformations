/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thirdparty.edu.princeton.cs.algs4;

import java.util.List;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

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
}
