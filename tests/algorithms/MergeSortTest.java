package algorithms;

import algorithms.util.PairIntWithIndex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MergeSortTest extends TestCase {
    
    public MergeSortTest() {
    }
    
    public void testSortByDecr() throws Exception {
        
        int[] a = new int[]{1, 2, 3, 4, 5, 6};

    	MergeSort.sortByDecr(a);

    	int[] expectedA = new int[]{6, 5, 4, 3, 2, 1};
        
        assertTrue(Arrays.equals(expectedA, a));
    }
    
    public void testSortByXThenY() throws Exception {
        
        List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
        
        list.add(new PairIntWithIndex(3, 4, 0));
        list.add(new PairIntWithIndex(1, 5, 0));
        list.add(new PairIntWithIndex(3, 3, 1));
        
        MergeSort.sortByXThenY(list);
        
        assertEquals(1, list.get(0).getX());
        assertEquals(5, list.get(0).getY());
        assertEquals(0, list.get(0).getPixIndex());
        
        assertEquals(3, list.get(1).getX());
        assertEquals(3, list.get(1).getY());
        assertEquals(1, list.get(1).getPixIndex());
        
        assertEquals(3, list.get(2).getX());
        assertEquals(4, list.get(2).getY());
        assertEquals(0, list.get(2).getPixIndex());
    }
    
    public void testSortByYThenX() throws Exception {
        
        List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
        
        list.add(new PairIntWithIndex(3, 4, 0));
        list.add(new PairIntWithIndex(1, 4, 0));
        list.add(new PairIntWithIndex(3, 3, 1));
        
        MergeSort.sortByYThenX(list);
        
        assertEquals(3, list.get(0).getX());
        assertEquals(3, list.get(0).getY());
        assertEquals(1, list.get(0).getPixIndex());
        
        assertEquals(1, list.get(1).getX());
        assertEquals(4, list.get(1).getY());
        assertEquals(0, list.get(1).getPixIndex());
        
        assertEquals(3, list.get(2).getX());
        assertEquals(4, list.get(2).getY());
        assertEquals(0, list.get(2).getPixIndex());
    }
}
