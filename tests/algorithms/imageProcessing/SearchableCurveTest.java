package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * 
 * @author nichole
 */
public class SearchableCurveTest {
    
    public SearchableCurveTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSortByYThenX() {
        
        System.out.println("sortByYThenX");
        
        PairIntArray xy = new PairIntArray();
        xy.add(1, 1);
        xy.add(2, 2);
        xy.add(4, 3);
        xy.add(3, 3);
        xy.add(2, 4);
        xy.add(1, 5);
        SearchableCurve instance = new SearchableCurve(xy);
        /*                       \/ \/
        int[] x = new int[]{1, 2, 4, 3, 2, 1};
        int[] y = new int[]{1, 2, 3, 3, 4, 5};
        int[] indexes = new int[]{0, 1, 2, 3, 4, 5};
        */
        int[] x = instance.getX();
        int[] y = instance.getY();
        
        assertTrue((y[0] == 5) && (x[0] == 1));
        assertTrue((y[1] == 4) && (x[1] == 2));
        assertTrue((y[2] == 3) && (x[2] == 3));
        assertTrue((y[3] == 3) && (x[3] == 4));
        assertTrue((y[4] == 2) && (x[4] == 2));
        assertTrue((y[5] == 1) && (x[5] == 1));
    }

    @Test
    public void testFindClosestMatchBinarySearch() {
        
        System.out.println("findClosestMatchBinarySearch");
        
        PairIntArray xy = new PairIntArray();
        xy.add(1, 1);
        xy.add(2, 2);
        xy.add(4, 3);
        xy.add(3, 3);
        xy.add(2, 4);
        xy.add(1, 5);
        
        /*
        1, 5
        2, 4
        3, 3
        4, 3
        2, 2
        1, 1
        */
        
        SearchableCurve instance = new SearchableCurve(xy);
        
        int idx = instance.findClosestMatchBinarySearch(3, 3);
        assertTrue(2 == idx);
        int idx2 = instance.findClosestMatchForYTooHigh(3, 3);
        assertTrue(2 == idx2);
        
        idx = instance.findClosestMatchBinarySearch(4, 3);
        assertTrue(3 == idx);
        idx2 = instance.findClosestMatchForYTooHigh(4, 3);
        assertTrue(3 == idx2);
        
        idx = instance.findClosestMatchBinarySearch(2, 3);
        assertTrue((1 == idx) || (2 == idx) || (4 == idx));
        idx2 = instance.findClosestMatchForYTooHigh(2, 3);
        assertTrue((1 == idx2) || (2 == idx2) || (4 == idx2));
        
        idx = instance.findClosestMatchBinarySearch(2, 0);
        assertTrue(4 == idx);
        idx2 = instance.findClosestMatchForYTooHigh(2, 0);
        
        idx = instance.findClosestMatchBinarySearch(1, 1);
        assertTrue(5 == idx);
        idx2 = instance.findClosestMatchForYTooHigh(1, 1);
        assertTrue(5 == idx2);
    }

    @Test
    public void testFindClosestMatch() {
        
        System.out.println("findClosestMatch");
        
        PairIntArray xy = new PairIntArray();
        xy.add(1, 1);
        xy.add(2, 2);
        xy.add(4, 3);
        xy.add(3, 3);
        xy.add(2, 4);
        xy.add(1, 5);
        
        /*
        1, 5
        2, 4
        3, 3
        4, 3
        2, 2
        1, 1
        */
        
        SearchableCurve instance = new SearchableCurve(xy);
        
        PairInt result = instance.findClosestMatch(3, 3);
        assertTrue((result.getX() == 3) && (result.getY() == 3));
        
        result = instance.findClosestMatch(4, 3);
        assertTrue((result.getX() == 4) && (result.getY() == 3));
        
        result = instance.findClosestMatch(1, 1);
        assertTrue((result.getX() == 1) && (result.getY() == 1));
        
    }
   
}
