package algorithms;

import java.security.SecureRandom;
import java.util.Arrays;
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
    
    public void testSortBy1stThen2nd() throws Exception {
        
        /*  a  b  c
            4  5  6
            4  6  7
            4  3  1
            1  2  3
        */
        double[] a = new double[]{4, 4, 4, 1};
        double[] b = new double[]{5, 6, 3, 2};
        int[] c = new int[]{6, 7, 1, 3};
        int[] d = new int[]{6, 7, 1, 3};
        
        QuickSort.sortBy1stThen2nd(a, b, c, d, 0, a.length - 1);
        
        assertTrue(a[0] == 1); assertTrue(b[0] == 2); assertTrue(c[0] == 3); assertTrue(d[0] == 3);
        assertTrue(a[1] == 4); assertTrue(b[1] == 3); assertTrue(c[1] == 1); assertTrue(d[1] == 1);
        assertTrue(a[2] == 4); assertTrue(b[2] == 5); assertTrue(c[2] == 6); assertTrue(d[2] == 6);
        assertTrue(a[3] == 4); assertTrue(b[3] == 6); assertTrue(c[3] == 7); assertTrue(d[3] == 7);
    }
    
    public void testSortByFirstArgument() {
        
        int[] a = new int[]{3, 7, 1, 5};
        Object[] b = new Object[]{"3", "7", "1", "5"};
        
        QuickSort.sortBy1stArg(a, b);
        
        assertTrue(a[0] == 1);
        assertEquals(b[0], "1");
        
        assertTrue(a[1] == 3);
        assertEquals(b[1], "3");
        
        assertTrue(a[2] == 5);
        assertEquals(b[2], "5");
        
        assertTrue(a[3] == 7);
        assertEquals(b[3], "7");
    }
    
    public void testSortByFirstArgument_2() {
        
        float[] a = new float[]{3, 7, 1, 5};
        int[] b = new int[]{0, 1, 2, 3};
        
        QuickSort.sortBy1stArg(a, b);
        
        assertTrue(a[0] == 1);
        assertTrue(b[0] == 2);
        
        assertTrue(a[1] == 3);
        assertTrue(b[1] == 0);
        
        assertTrue(a[2] == 5);
        assertTrue(b[2] == 3);
        
        assertTrue(a[3] == 7);
        assertTrue(b[3] == 1);
    }
    
    public void testSortByFirstArgument2() {
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int n = 100;
        int[] a = new int[n];
        String[] b = new String[n];
        
        for (int i = 0; i < n; ++i) {
            int v = sr.nextInt();
            String vStr = Integer.toString(v);
            a[i] = v;
            b[i] = vStr;
        }
        
        int[] expected = Arrays.copyOf(a, n);
        Arrays.sort(expected);
        
        QuickSort.sortBy1stArg(a, b);
        
        for (int i = 0; i < n; i++) {
            int v = a[i];
            String vStr = Integer.toString(v);
            assertTrue(expected[i] == v);
            assertEquals(b[i], vStr);
        }
    }
    
    public void testSortByFirstArgument3() {
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int n = 100;
        float[] a = new float[n];
        String[] b = new String[n];
        
        for (int i = 0; i < n; ++i) {
            float v = sr.nextFloat();
            String vStr = Float.toString(v);
            a[i] = v;
            b[i] = vStr;
        }
        
        float[] expected = Arrays.copyOf(a, n);
        Arrays.sort(expected);
        
        QuickSort.sortBy1stArg(a, b);
        
        for (int i = 0; i < n; i++) {
            float v = a[i];
            String vStr = Float.toString(v);
            assertTrue(expected[i] == v);
            assertEquals(b[i], vStr);
        }
    }
    
}
