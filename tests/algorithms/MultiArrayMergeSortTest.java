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
            MultiArrayMergeSort.sortBy1stArgThen2nd(null, null);
        } catch (IllegalArgumentException t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MultiArrayMergeSort.sortBy1stArgThen2nd(new float[]{1f,2f}, null);
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
