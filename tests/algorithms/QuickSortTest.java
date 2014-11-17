package algorithms;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class QuickSortTest {
    
    public QuickSortTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testSort() {
                
                            //  0  1  2  3  4  5  6  7  8  9
                            //              |  |
        float[] a = new float[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};

        QuickSort.sort(a);

        for (int i = 0; i < a.length; i++) {            
            assertTrue(a[i] == i);
        }
    }
}
