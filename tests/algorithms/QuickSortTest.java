package algorithms;

import algorithms.imageProcessing.util.PairIntWithIndex0;
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
    
    public void testSortBy1stThen2nd_2() throws Exception {
        
        /*  a  b  c
            4  5  6
            4  6  7
            4  3  1
            1  2  3
        */
        float[] a = new float[]{4, 4, 4, 1};
        float[] b = new float[]{5, 6, 3, 2};
        
        QuickSort.sortBy1stThen2nd(a, b, 0, a.length - 1);
        
        assertTrue(a[0] == 1); assertTrue(b[0] == 2);
        assertTrue(a[1] == 4); assertTrue(b[1] == 3); 
        assertTrue(a[2] == 4); assertTrue(b[2] == 5);
        assertTrue(a[3] == 4); assertTrue(b[3] == 6);
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
    
    public void testSortByYThenX() throws Exception {
        
        PairIntWithIndex0[] points = new PairIntWithIndex0[4];
        points[0] = new PairIntWithIndex0(10, 5, 2);
        points[1] = new PairIntWithIndex0(10, 6, 3);
        points[2] = new PairIntWithIndex0(1, 1, 1);
        points[3] = new PairIntWithIndex0(1, 0, 0);
        
        QuickSort.sortByYThenX(points);
        
        for (int i = 0; i < points.length; ++i) {
            assertEquals(i, points[i].getPixIndex());
        }
        
    }
    
     public void testDescendingSort() {
                
                        //  0  1  2  3  4  5  6  7  8  9
                        //              |  |
        int[] a = new int[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};
        int[] b = new int[]{5, 4, 0, 9, 6, 8, 3, 7, 2, 1};

        QuickSort.descendingSort(a, b);

        int prev = a[0];
        for (int i = 0; i < a.length; i++) {
            assertEquals(i, b[i]);
            assertTrue(a[i] <= prev);
            prev = a[i];
        }
    }
     
     public void testDescendingSort2() {
                          
        double[] a = new double[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};

        QuickSort.descendingSort(a);

        double prev = a[0];
        for (int i = 0; i < a.length; i++) {
            assertEquals((double)i, a.length - a[i] - 1);
            assertTrue(a[i] <= prev);
            prev = a[i];
        }
    }
}
