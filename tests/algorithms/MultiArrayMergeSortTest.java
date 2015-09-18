package algorithms;

import java.util.Arrays;
import algorithms.util.PairIntArray;
import static junit.framework.Assert.assertTrue;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/sorting/
 * under MIT License (MIT), Nichole King 2013
 * and added to
 */
public class MultiArrayMergeSortTest extends TestCase {

    /**
     * @see TestCase#setUp()
     */
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     * @see TestCase#tearDown()
     */
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testSortBy1stArgThen2nd(){

        /*   10
         *
         *                *p12
         *    9
         *
         *
         *    8
         *
         *
         *    7
         *
         *
         *    6                             *p10 *p11
         *
         *                   *p9       *p8
         *    5         *p7
         *
         *                        *p6
         *    4                                            *p5        *p4
         *
         *     p3*
         *    3                                                 *p2
         *
         *
         *    2                                                   *p1
         *
         *
         *    1            *p0
         *
         *
         *    0    1    2    3    4    5    6    7    8    9   10   11   12
    	 */
    	float[] x = new float[13];
    	float[] y = new float[13];
        x[0] =  2.5f;
        y[0] =  1.0f;
        x[1] = 10.5f;
        y[1] =  2.0f;
        x[2] = 10.0f;
        y[2] =  3.0f;
        x[3] =  0.5f;
        y[3] =  3.2f;
        x[5] =  9.0f;
        y[5] =  4.0f;

        x[4] = 11.2f;
        y[4] =  4.0f;
        x[6] =  4.0f;
        y[6] =  4.2f;
        x[7] =  2.0f;
        y[7] =  5.0f;
        x[9] =  3.0f;
        y[9] =  5.3f;
        x[8] =  5.0f;
        y[8] =  5.3f;

        x[10] =  6.0f;
        y[10] =  6.0f;
        x[11] =  7.0f;
        y[11] =  6.0f;
        x[12] =  2.3f;
        y[12] =  9.3f;


        float[] ex = Arrays.copyOf(x, x.length);
    	float[] ey = Arrays.copyOf(y, y.length);
        ex[4] =  x[5];
        ey[4] =  y[5];
        ex[5] =  x[4];
        ey[5] =  y[4];
        
        ex[8] =  x[9];
        ey[8] =  y[9];
        ex[9] =  x[8];
        ey[9] =  y[8];

    	MultiArrayMergeSort.sortBy1stArgThen2nd(y, x);
    	assertTrue(x.length == ex.length);

    	for (int i=0; i < ex.length; i++) {
            assertTrue( Math.abs(ex[i] - x[i]) < 0.01);
            assertTrue( Math.abs(ey[i] - y[i]) < 0.01);
        }
    	
    	boolean caughtException = true;
        try {
            float[] x2 = null;
            float[] y2 = null;
            MultiArrayMergeSort.sortBy1stArgThen2nd(x2, y2);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            float[] y2 = null;
            MultiArrayMergeSort.sortBy1stArgThen2nd(new float[]{1f,2f}, y2);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MultiArrayMergeSort.sortBy1stArgThen2nd(new float[]{1f,2f}, new float[]{1f,2f,3f});
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
    }
    
    public void testSortBy1stArgThen2nd_int(){

        /*   10
         *
         *                *p12
         *    9
         *
         *
         *    8
         *
         *
         *    7
         *
         *
         *    6                             *p10 *p11
         *
         *                   *p9       *p8
         *    5         *p7
         *
         *                        *p6
         *    4                                            *p5        *p4
         *
         *     p3*
         *    3                                                 *p2
         *
         *
         *    2                                                   *p1
         *
         *
         *    1            *p0
         *
         *
         *    0    1    2    3    4    5    6    7    8    9   10   11   12
    	 */
    	int[] x = new int[13];
    	int[] y = new int[13];
        x[0] =  25;
        y[0] =  10;
        x[1] = 105;
        y[1] =  20;
        x[2] = 100;
        y[2] =  30;
        x[3] =  5;
        y[3] =  32;
        x[5] =  90;
        y[5] =  40;

        x[4] = 112;
        y[4] =  40;
        x[6] =  40;
        y[6] =  42;
        x[7] =  20;
        y[7] =  50;
        x[9] =  30;
        y[9] =  53;
        x[8] =  50;
        y[8] =  53;

        x[10] =  60;
        y[10] =  60;
        x[11] =  70;
        y[11] =  60;
        x[12] =  23;
        y[12] =  93;


        int[] ex = Arrays.copyOf(x, x.length);
    	int[] ey = Arrays.copyOf(y, y.length);
        ex[4] =  x[5];
        ey[4] =  y[5];
        ex[5] =  x[4];
        ey[5] =  y[4];
        
        ex[8] =  x[9];
        ey[8] =  y[9];
        ex[9] =  x[8];
        ey[9] =  y[8];

    	MultiArrayMergeSort.sortBy1stArgThen2nd(y, x);
    	assertTrue(x.length == ex.length);

    	for (int i=0; i < ex.length; i++) {
            assertTrue( Math.abs(ex[i] - x[i]) < 0.01);
            assertTrue( Math.abs(ey[i] - y[i]) < 0.01);
        }
    	
    	boolean caughtException = true;
        try {
            int[] a0 = null;
            int[] a1 = null;
            MultiArrayMergeSort.sortBy1stArgThen2nd(a0, a1);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            int[] a1 = null;
            MultiArrayMergeSort.sortBy1stArgThen2nd(new int[]{1, 2}, a1);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MultiArrayMergeSort.sortBy1stArgThen2nd(new int[]{1, 2}, new int[]{1, 2, 3});
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
    }
    
    public void testSortByYThenX(){

        /*   10
         *
         *                *p12
         *    9
         *
         *
         *    8
         *
         *
         *    7
         *
         *
         *    6                             *p10 *p11
         *
         *                   *p9       *p8
         *    5         *p7
         *
         *                        *p6
         *    4                                            *p5        *p4
         *
         *     p3*
         *    3                                                 *p2
         *
         *
         *    2                                                   *p1
         *
         *
         *    1            *p0
         *
         *
         *    0    1    2    3    4    5    6    7    8    9   10   11   12
    	 */
                
    	float[] x = new float[13];
    	float[] y = new float[13];
        x[0] =  2.5f;
        y[0] =  1.0f;
        x[1] = 10.5f;
        y[1] =  2.0f;
        x[2] = 10.0f;
        y[2] =  3.0f;
        x[3] =  0.5f;
        y[3] =  3.2f;
        x[5] =  9.0f;
        y[5] =  4.0f;

        x[4] = 11.2f;
        y[4] =  4.0f;
        x[6] =  4.0f;
        y[6] =  4.2f;
        x[7] =  2.0f;
        y[7] =  5.0f;
        x[9] =  3.0f;
        y[9] =  5.3f;
        x[8] =  5.0f;
        y[8] =  5.3f;

        x[10] =  6.0f;
        y[10] =  6.0f;
        x[11] =  7.0f;
        y[11] =  6.0f;
        x[12] =  2.3f;
        y[12] =  9.3f;
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < x.length; i++) {
            int xp = Math.round(x[i] * 10.f);
            int yp = Math.round(y[i] * 10.f);
            xy.add(xp, yp);
        }

        float[] ex = Arrays.copyOf(x, x.length);
    	float[] ey = Arrays.copyOf(y, y.length);
        ex[4] =  x[5];
        ey[4] =  y[5];
        ex[5] =  x[4];
        ey[5] =  y[4];
        
        ex[8] =  x[9];
        ey[8] =  y[9];
        ex[9] =  x[8];
        ey[9] =  y[8];
        
        PairIntArray exy = new PairIntArray();
        for (int i = 0; i < x.length; i++) {
            int xp = Math.round(ex[i] * 10.f);
            int yp = Math.round(ey[i] * 10.f);
            exy.add(xp, yp);
        }
        
    	MultiArrayMergeSort.sortByYThenX(xy);

    	for (int i=0; i < exy.getN(); i++) {
            assertTrue(exy.getX(i) == xy.getX(i));
            assertTrue(exy.getY(i) == xy.getY(i));
        }
    	
    }
    
    public void testSortByDecr() throws Exception {
        
        int[] a = new int[]{1, 2, 3, 4, 5, 6};
    	int[] b = new int[]{0, 1, 2, 3, 4, 5};

    	MultiArrayMergeSort.sortByDecr(a, b);
    	assertTrue(a.length == b.length);

    	int[] expectedA = new int[]{6, 5, 4, 3, 2, 1};
        int[] expectedB = new int[]{5, 4, 3, 2, 1, 0};
        
        assertTrue(Arrays.equals(expectedA, a));
        assertTrue(Arrays.equals(expectedB, b));
    }
    
    public void testSortByDecr2() throws Exception {
        
        double[] a = new double[]{1, 2, 3, 4, 5, 6};
    	int[] b = new int[]{0, 1, 2, 3, 4, 5};

    	MultiArrayMergeSort.sortByDecr(a, b);
    	assertTrue(a.length == b.length);

    	double[] expectedA = new double[]{6, 5, 4, 3, 2, 1};
        int[] expectedB = new int[]{5, 4, 3, 2, 1, 0};
        
        assertTrue(Arrays.equals(expectedA, a));
        assertTrue(Arrays.equals(expectedB, b));
    }
    
    public void testSortBy1stDescThen2ndAsc() throws Exception {
        
        int[] nSimilarSummary = new int[2];
        double[] costsSummary = new double[2];
        Integer[][] indexesSummary = new Integer[2][];
        int[] mainIndexSummary = new int[2];
        
        // ===================================
        nSimilarSummary[0] = 1;
        nSimilarSummary[1] = 4;
        costsSummary[0] = 3.;
        costsSummary[1] = 3;
        indexesSummary[0] = new Integer[]{Integer.valueOf(0)};
        indexesSummary[1] = new Integer[]{Integer.valueOf(1)};
        mainIndexSummary[0] = 0;
        mainIndexSummary[1] = 1;
        
        MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary,
            indexesSummary, mainIndexSummary);
        
        assertTrue(nSimilarSummary[0] == 4);
        assertTrue(costsSummary[0] == 3);
        assertTrue(mainIndexSummary[0] == 1);
        assertTrue(Arrays.equals(indexesSummary[0], new Integer[]{Integer.valueOf(1)}));
        assertTrue(nSimilarSummary[1] == 1);
        assertTrue(costsSummary[1] == 3);
        assertTrue(mainIndexSummary[1] == 0);
        assertTrue(Arrays.equals(indexesSummary[1], new Integer[]{Integer.valueOf(0)}));

        // ================================
        nSimilarSummary[1] = 1;
        nSimilarSummary[0] = 4;
        costsSummary[1] = 3.;
        costsSummary[0] = 3;
        indexesSummary[1] = new Integer[]{Integer.valueOf(0)};
        indexesSummary[0] = new Integer[]{Integer.valueOf(1)};
        mainIndexSummary[1] = 0;
        mainIndexSummary[0] = 1;
        
        MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary,
            indexesSummary, mainIndexSummary);
        
        assertTrue(nSimilarSummary[0] == 4);
        assertTrue(costsSummary[0] == 3);
        assertTrue(mainIndexSummary[0] == 1);
        assertTrue(Arrays.equals(indexesSummary[0], new Integer[]{Integer.valueOf(1)}));
        assertTrue(nSimilarSummary[1] == 1);
        assertTrue(costsSummary[1] == 3);
        assertTrue(mainIndexSummary[1] == 0);
        assertTrue(Arrays.equals(indexesSummary[1], new Integer[]{Integer.valueOf(0)}));
        
        //====================================================
        
        nSimilarSummary[0] = 1;
        nSimilarSummary[1] = 1;
        costsSummary[0] = 3.;
        costsSummary[1] = 400;
        indexesSummary[0] = new Integer[]{Integer.valueOf(0)};
        indexesSummary[1] = new Integer[]{Integer.valueOf(1)};
        mainIndexSummary[0] = 0;
        mainIndexSummary[1] = 1;
        
        MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary,
            indexesSummary, mainIndexSummary);
        
        assertTrue(nSimilarSummary[0] == 1);
        assertTrue(costsSummary[0] == 3);
        assertTrue(mainIndexSummary[0] == 0);
        assertTrue(Arrays.equals(indexesSummary[0], new Integer[]{Integer.valueOf(0)}));
        
        assertTrue(nSimilarSummary[1] == 1);
        assertTrue(costsSummary[1] == 400);
        assertTrue(mainIndexSummary[1] == 1);
        assertTrue(Arrays.equals(indexesSummary[1], new Integer[]{Integer.valueOf(1)}));
        
        // ==== reverse order of input and assert same results ===
        nSimilarSummary[1] = 1;
        nSimilarSummary[0] = 1;
        costsSummary[1] = 3.;
        costsSummary[0] = 400;
        indexesSummary[1] = new Integer[]{Integer.valueOf(0)};
        indexesSummary[0] = new Integer[]{Integer.valueOf(1)};
        mainIndexSummary[1] = 0;
        mainIndexSummary[0] = 1;
        
        MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary,
            indexesSummary, mainIndexSummary);
        
        assertTrue(nSimilarSummary[0] == 1);
        assertTrue(costsSummary[0] == 3);
        assertTrue(mainIndexSummary[0] == 0);
        assertTrue(Arrays.equals(indexesSummary[0], new Integer[]{Integer.valueOf(0)}));
        
        assertTrue(nSimilarSummary[1] == 1);
        assertTrue(costsSummary[1] == 400);
        assertTrue(mainIndexSummary[1] == 1);
        assertTrue(Arrays.equals(indexesSummary[1], new Integer[]{Integer.valueOf(1)}));
    }
    
    /*
    private static void mergeBy1stDescThen2ndAsc(double[] a1, int[] a2, 
        Integer[][] a3, int[] a4, int idxLo, int idxMid, int idxHi) {
    */
    public void testSortBy1stAscThen2ndDesc() throws Exception {
        
        double[] a1 = new double[]{ 5,  4,  3,  3,  2,  1};
        int[] a2    = new    int[]{10, 20, 40, 30, 50, 60};
        Integer[][] a3 = new Integer[6][];
        for (int i = 0; i < a3.length; ++i) {
            a3[i] = new Integer[]{i};
        }
        int[] a4 = new int[]{100, 200, 300, 400, 500, 550};
        MultiArrayMergeSort.sortBy1stAscThen2ndDesc(a1, a2, a3, a4, 0, a1.length - 1);
        
        assertTrue(a1[0] == 1); 
        assertTrue(a1[1] == 2); 
        assertTrue(a1[2] == 3);
        assertTrue(a1[3] == 3);
        assertTrue(a1[4] == 4);
        assertTrue(a1[5] == 5);
        
        assertTrue(a2[0] == 60); 
        assertTrue(a2[1] == 50); 
        assertTrue(a2[2] == 40);
        assertTrue(a2[3] == 30);
        assertTrue(a2[4] == 20);
        assertTrue(a2[5] == 10);
        
        assertTrue(a3[0][0] == 5); 
        assertTrue(a3[1][0] == 4); 
        assertTrue(a3[2][0] == 2);
        assertTrue(a3[3][0] == 3);
        assertTrue(a3[4][0] == 1);
        assertTrue(a3[5][0] == 0);
        
        assertTrue(a4[0] == 550); 
        assertTrue(a4[1] == 500); 
        assertTrue(a4[2] == 300);
        assertTrue(a4[3] == 400);
        assertTrue(a4[4] == 200);
        assertTrue(a4[5] == 100);
    }
    
    /**
     * Test suite
     * @return static Test
    */
    public static Test suite(){
        return new TestSuite(MultiArrayMergeSortTest.class);
    }

    /**
     * Set up a Junit test runner
     * @param args Not used.
    */
    public static void main(String[] args) {

        junit.textui.TestRunner.run(suite());
    }

}
