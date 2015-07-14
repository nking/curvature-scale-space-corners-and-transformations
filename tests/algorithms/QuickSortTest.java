package algorithms;

import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class QuickSortTest extends TestCase {
    
    public QuickSortTest() {
    }
    
    public void testSort() {
                
                            //  0  1  2  3  4  5  6  7  8  9
                            //              |  |
        float[] a = new float[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};

        QuickSort.sort(a);

        for (int i = 0; i < a.length; i++) {            
            assertTrue(a[i] == i);
        }
    }
    
    public void testSort2() {
                
                            //  0  1  2  3  4  5  6  7  8  9
                            //              |  |
        float[] a = new float[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};
        int[] b = new int[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};
        int[] c = new int[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};

        QuickSort.sort(a, b, c, 0, a.length - 1);

        for (int i = 0; i < a.length; i++) {
            assertTrue(a[i] == i);
            assertTrue(b[i] == i);
            assertTrue(c[i] == i);
        }
    }
    
    public void testSortBy1stThen2ndThen3rd() throws Exception {
        
        /*  a  b  c
            4  5  6
            4  6  7
            4  3  1
            1  2  3
        */
        float[] a = new float[]{4, 4, 4, 1};
        float[] b = new float[]{5, 6, 3, 2};
        float[] c = new float[]{6, 7, 1, 3};
        
        QuickSort.sortBy1stThen2ndThen3rd(a, b, c, 0, a.length - 1);
        
        assertTrue(a[0] == 1); assertTrue(b[0] == 2); assertTrue(c[0] == 3);
        assertTrue(a[1] == 4); assertTrue(b[1] == 3); assertTrue(c[1] == 1);
        assertTrue(a[2] == 4); assertTrue(b[2] == 5); assertTrue(c[2] == 6);
        assertTrue(a[3] == 4); assertTrue(b[3] == 6); assertTrue(c[3] == 7);
    }
    
}
