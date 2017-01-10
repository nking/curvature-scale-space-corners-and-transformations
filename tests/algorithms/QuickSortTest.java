package algorithms;

import algorithms.imageProcessing.util.PairIntWithIndex0;
import algorithms.util.PairInt;
import algorithms.util.IntIntDouble;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
    
    public void testSortBy1stThen2nd_3() throws Exception {
        
        /*  a  b  c
            4  5  6
            4  6  7
            4  3  1
            1  2  3
        */
        TFloatList a = new TFloatArrayList(new float[]{4, 4, 4, 1});
        TFloatList b = new TFloatArrayList(new float[]{5, 6, 3, 2});
        TIntList c = new TIntArrayList(new int[]{2, 3, 1, 0});
        
        QuickSort.sortBy1stThen2nd(a, b, c);
        
        assertEquals(1.f, a.get(0));
        assertEquals(4.f, a.get(1));
        assertEquals(4.f, a.get(2));
        assertEquals(4.f, a.get(3));
        
        assertEquals(2.f, b.get(0));
        assertEquals(3.f, b.get(1));
        assertEquals(5.f, b.get(2));
        assertEquals(6.f, b.get(3));
        
        assertEquals(0, c.get(0));
        assertEquals(1, c.get(1));
        assertEquals(2, c.get(2));
        assertEquals(3, c.get(3));
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
        
        // ----
        TDoubleList c = new TDoubleArrayList();
        c.add(3);
        c.add(7);
        c.add(1);
        c.add(5);
        int[] d = new int[]{3, 7, 1, 5};
        
        QuickSort.sortBy1stArg(c, d);
        
        assertTrue(c.get(0) == 1);
        assertEquals(1, d[0]);
        
        assertTrue(c.get(1) == 3);
        assertEquals(3, d[1]);
        
        assertTrue(c.get(2) == 5);
        assertEquals(5, d[2]);
        
        assertTrue(c.get(3) == 7);
        assertEquals(7, d[3]);
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
     
     public void testDescendingSort1() {
         
         List<Float> a = new ArrayList<Float>();
         a.add(4f);
         a.add(5f);
         a.add(9f);
         a.add(0f);
         a.add(3f);
         a.add(1f);
         a.add(6f);
         a.add(2f);
         a.add(7f);
         a.add(8f);
         List<Float> b = new ArrayList<Float>();
         b.add(5f);
         b.add(4f);
         b.add(0f);
         b.add(9f);
         b.add(6f);
         b.add(8f);
         b.add(3f);
         b.add(7f);
         b.add(2f);
         b.add(1f);
         QuickSort.descendingSort(a, b);

         float prev = a.get(0).floatValue();
         for (int i = 0; i < a.size(); i++) {
             assertTrue(Math.abs(i - b.get(i).floatValue()) < 0.01);
             assertTrue(a.get(i).floatValue() <= prev);
             prev = a.get(i).floatValue();
         }
         
         // -----
         //descendingSort(TIntList a, List<? extends Object> b
         TIntList aa = new TIntArrayList(a.size());
         aa.add(4);
         aa.add(5);
         aa.add(9);
         aa.add(0);
         aa.add(3);
         aa.add(1);
         aa.add(6);
         aa.add(2);
         aa.add(7);
         aa.add(8);
         b = new ArrayList<Float>();
         b.add(5f);
         b.add(4f);
         b.add(0f);
         b.add(9f);
         b.add(6f);
         b.add(8f);
         b.add(3f);
         b.add(7f);
         b.add(2f);
         b.add(1f);
         QuickSort.descendingSort(aa, b);
         
         int prev2 = aa.get(0);
         for (int i = 0; i < a.size(); i++) {
             assertTrue(Math.abs(i - b.get(i)) < 0.01);
             assertTrue(aa.get(i) <= prev2);
             prev2 = aa.get(i);
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
 
    public void testSortByDimension1FirstSecond() {

        int[][] a = new int[2][4];
        a[0] = new int[4];
        a[1] = new int[4];
        
        a[0][0] = 4;
        a[0][1] = 3;
        a[0][2] = 1;
        a[0][3] = 1;
        
        a[1][0] = 4;
        a[1][1] = 3;
        a[1][2] = 2;
        a[1][3] = 1;        

        QuickSort.sortByDimension1FirstSecond(a);
        
        for (int i = 0; i < 4; ++i) {
            assertEquals(i + 1, a[1][i]);
        }
    }

    public void testSortByDimension1FirstSecondThird() {

        int[][] a = new int[3][4];
        a[0] = new int[4];
        a[1] = new int[4];
        a[2] = new int[4];
        
        a[0][0] = 4;
        a[0][1] = 3;
        a[0][2] = 1;
        a[0][3] = 1;
        
        a[1][0] = 4;
        a[1][1] = 3;
        a[1][2] = 1;
        a[1][3] = 1; 
        
        a[2][0] = 4;
        a[2][1] = 3;
        a[2][2] = 2;
        a[2][3] = 1;        

        QuickSort.sortByDimension1FirstSecondThird(a);
        
        for (int i = 0; i < 4; ++i) {
            assertEquals(i + 1, a[2][i]);
        }
    }

    public void sortByDecrYThenIncrX() {
    
        PairInt[] a = new PairInt[4];
        a[0] = new PairInt(7, 7);
        a[1] = new PairInt(0, 0);
        a[2] = new PairInt(4, 4);
        a[3] = new PairInt(2, 4);
        int[] indexes = new int[]{0, 1, 2, 3};

        QuickSort.sortByDecrYThenIncrX(a, indexes);
        
        PairInt[] aExp = new PairInt[4];
        aExp[0] = new PairInt(7, 7);
        aExp[3] = new PairInt(0, 0);
        aExp[1] = new PairInt(4, 4);
        aExp[2] = new PairInt(2, 4);
        int[] indexesExp = new int[]{0, 3, 1, 2};

        for (int i = 0; i < a.length; ++i) {
            assertEquals(indexesExp[i], indexes[i]);
            assertEquals(aExp[i], a[i]);
        }
    }

    public void testSortBy1stArg() {
        
        TIntList a = new TIntArrayList();
        TDoubleList b = new TDoubleArrayList();
        TIntList c = new TIntArrayList();
        
        for (int i = 0; i < 10; ++i) {
            a.add(9 - i);
            b.add(9 - i);
            c.add(9 - i);
        }
        
        QuickSort.sortBy1stArg(a, b, c);
        
        for (int i = 0; i < 10; ++i) {
            assertEquals(i, a.get(i));
            assertEquals((double)i, b.get(i));
            assertEquals(i, c.get(i));
        }
    }
    
    public void testSortBy1stArg1() {
        
        TFloatList a = new TFloatArrayList();
        TIntList b = new TIntArrayList();
        TIntList c = new TIntArrayList();
        
        for (int i = 0; i < 10; ++i) {
            a.add(9 - i);
            b.add(9 - i);
            c.add(9 - i);
        }
        
        QuickSort.sortBy1stArg(a, b, c);
        
        for (int i = 0; i < 10; ++i) {
            assertEquals((float)i, a.get(i));
            assertEquals(i, b.get(i));
            assertEquals(i, c.get(i));
        }
    }
    
    public void testSortByA() {
        
        int n = 10;
        
        IntIntDouble[] abc = new IntIntDouble[n];
    
        for (int i = 0; i < abc.length; ++i) {
            int v = n - 1 - i;
            abc[i] = new IntIntDouble(v, v, v);
        }
        
        QuickSort.sortByA(abc);
        
        for (int i = 0; i < abc.length; ++i) {
            assertEquals(i, abc[i].getA());
        }
    }
    
     public void testSortBy1stArg_2() {
        
        TIntList a = new TIntArrayList();
        TIntList c = new TIntArrayList();
        
        for (int i = 0; i < 10; ++i) {
            a.add(9 - i);
            c.add(9 - i);
        }
        
        QuickSort.sortBy1stArg(a, c);
        
        for (int i = 0; i < 10; ++i) {
            assertEquals(i, a.get(i));
            assertEquals(i, c.get(i));
        }
    }
    
    public void testSortBy1stThen2ndThen3rd_2() {

        TIntList a = new TIntArrayList();
        TIntList b = new TIntArrayList();
        TIntList c = new TIntArrayList();

        a.add(5);
        b.add(4);
        c.add(4);
        
        a.add(5);
        b.add(7);
        c.add(7);
        
        a.add(5);
        b.add(4);
        c.add(5);

        a.add(1);
        b.add(2);
        c.add(3);
        
        QuickSort.sortBy1stThen2ndThen3rd(a, b, c);
        
        assertEquals(1, a.get(0));
        assertEquals(2, b.get(0));
        assertEquals(3, c.get(0));
        
        assertEquals(5, a.get(1));
        assertEquals(4, b.get(1));
        assertEquals(4, c.get(1));
        
        assertEquals(5, a.get(2));
        assertEquals(4, b.get(2));
        assertEquals(5, c.get(2));
        
        assertEquals(5, a.get(3));
        assertEquals(7, b.get(3));
        assertEquals(7, c.get(3));
        
        /*
        a.add(1);
        b.add(2);
        c.add(3);
        
        a.add(5);
        b.add(4);
        c.add(4);
        
        a.add(5);
        b.add(4);
        c.add(5);

        a.add(5);
        b.add(7);
        c.add(7);        
        */
    }

    public void testSortByXAsc() {
        
        int[] x = new int[]{10, 1, 5, 30,  5};
        int[] y = new int[]{10, 1, 6, 30,  5};
        
        int[] expectedX = new int[]{1, 5, 5, 10, 30};
        int[] expectedY = new int[]{1, 5, 6, 10, 30};
        
        QuickSort.sortBy1stThen2nd(x, y);
        
        for (int i = 0; i < x.length; i++) {
            assertEquals(expectedX[i], x[i]);
            assertEquals(expectedY[i], y[i]);
        }
    }
    
    public static void testSortBy1stThen2ndThen3rd_5() {
        
        float[] a = new float[]{4, 1, 1, 1};
        float[] b = new float[]{4, 2, 2, 1};
        float[] c = new float[]{4, 3, 2, 1};
        int[] d = new int[]{3, 2, 1, 0};
        
        QuickSort.sortBy1stThen2ndThen3rd(a, b, c, d);
        
        assertTrue(Arrays.equals(new float[]{1, 1, 1, 4}, 
            a));
        assertTrue(Arrays.equals(new float[]{1, 2, 2, 4}, 
            b));
        assertTrue(Arrays.equals(new float[]{1, 2, 3, 4}, 
            c));
        assertTrue(Arrays.equals(new int[]{0, 1, 2, 3}, 
            d));
    }
    
    public void testSortByFirstArg_3() {
        
        TFloatList a = new TFloatArrayList();
        a.add(3); a.add(7); a.add(1); a.add(5);
        
        List b = new ArrayList();
        b.add("3"); b.add("7"); b.add("1"); b.add("5");
        
        QuickSort.sortBy1stArg(a, b);
        
        assertEquals(1.f, a.get(0));
        assertEquals("1", (String)b.get(0));
        
        assertEquals(3.f, a.get(1));
        assertEquals("3", (String)b.get(1));
        
        assertEquals(5.f, a.get(2));
        assertEquals("5", (String)b.get(2));
        
        assertEquals(7.f, a.get(3));
        assertEquals("7", (String)b.get(3));
        
    }
}
