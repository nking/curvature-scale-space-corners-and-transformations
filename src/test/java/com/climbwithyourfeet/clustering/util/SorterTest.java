package com.climbwithyourfeet.clustering.util;

import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SorterTest extends TestCase {
    
    /**
     *
     * @param testName
     */
    public SorterTest(String testName) {
        super(testName);
    }

    /**
     *
     */
    public void testMergeSortByXThenY_List() {
        
        List<PairFloat> a = new ArrayList<PairFloat>();
        a.add(new PairFloat(100.f, 100.f));
        a.add(new PairFloat(100.f, 70.f));
        a.add(new PairFloat(10.f, 7.f));
        a.add(new PairFloat(1000.f, 7.f));
        
        Sorter.mergeSortByXThenY(a);
        
        assertTrue(a.size() == 4);
        PairFloat t = a.get(0);
        assertTrue(t.getX() == 10.f && t.getY() == 7.f);
        t = a.get(1);
        assertTrue(t.getX() == 100.f && t.getY() == 70.f);
        t = a.get(2);
        assertTrue(t.getX() == 100.f && t.getY() == 100.f);
        t = a.get(3);
        assertTrue(t.getX() == 1000.f && t.getY() == 7.f);
        
        Sorter.mergeSortByYThenX(a);
        assertTrue(a.size() == 4);
        t = a.get(0);
        assertTrue(t.getX() == 10.f && t.getY() == 7.f);
        t = a.get(1);
        assertTrue(t.getX() == 1000.f && t.getY() == 7.f);
        t = a.get(2);
        assertTrue(t.getX() == 100.f && t.getY() == 70.f);
        t = a.get(3);
        assertTrue(t.getX() == 100.f && t.getY() == 100.f);
    }
    
}
