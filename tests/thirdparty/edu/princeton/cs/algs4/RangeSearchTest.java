package thirdparty.edu.princeton.cs.algs4;

import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class RangeSearchTest extends TestCase {
    
    public RangeSearchTest() {
    }

    public void testRangeSearch() {
        
        IntervalRangeSearch<Integer, Integer> rt =
            new IntervalRangeSearch<Integer, Integer>();
        
        // the mmd2 ranges are non overlapping for inserts
        Interval<Integer> o1 = new Interval<Integer>(15, 20);
        Interval<Integer> o2 = new Interval<Integer>(21, 22);
        Interval<Integer> o3 = new Interval<Integer>(45, 60);
        rt.put(o1, Integer.valueOf(1));
        rt.put(o2, Integer.valueOf(2));
        rt.put(o3, Integer.valueOf(3));
                
        Set<Interval<Integer>> expected = new
            HashSet<Interval<Integer>>();
        expected.add(o1);
        expected.add(o2);
        
        Interval<Integer> s0 = new Interval<Integer>(19, 23);
        
        Queue<Interval<Integer>> queue = rt.range0(s0);
        int nq = queue.size();
        for (Interval<Integer> fnd : queue) {
            System.out.println("found=" + fnd);
            assertTrue(expected.remove(fnd));
        }
        assertEquals(2, nq);
        assertTrue(expected.isEmpty());
    }
    
    public void testPutIfLessThan() {
        
        IntervalRangeSearch<Integer, Integer> rt =
            new IntervalRangeSearch<Integer, Integer>();
        
        boolean didIns;
        Integer v;
        
        // the mmd2 ranges are non overlapping for inserts
        Interval<Integer> o1 = new Interval<Integer>(15, 20);
        Interval<Integer> o2 = new Interval<Integer>(21, 22);
        Interval<Integer> o3 = new Interval<Integer>(21, 60);
        v = Integer.valueOf(1);
        didIns = rt.putIfLessThan(o1, v, v);
        assertTrue(didIns);
        
        v = Integer.valueOf(2);
        didIns = rt.putIfLessThan(o2, v, v);
        assertTrue(didIns);
        
        v = Integer.valueOf(3);
        didIns = rt.putIfLessThan(o3, v, v);
        assertFalse(didIns);
        
        Set<Interval<Integer>> expected = new
            HashSet<Interval<Integer>>();
        expected.add(o1);
        expected.add(o2);
        
        Interval<Integer> s0 = new Interval<Integer>(19, 23);
        
        Queue<Interval<Integer>> queue = rt.range0(s0);
        int nq = queue.size();
        for (Interval<Integer> fnd : queue) {
            System.out.println("found=" + fnd);
            assertTrue(expected.remove(fnd));
        }
        assertEquals(2, nq);
        assertTrue(expected.isEmpty());
        
    }
    
}
