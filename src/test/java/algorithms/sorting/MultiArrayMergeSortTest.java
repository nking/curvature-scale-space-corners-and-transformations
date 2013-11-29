package algorithms.sorting;

import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

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

    public void testSortByY() throws Exception {

        int nPoints = 1000;
    	float[] x = new float[nPoints];
    	float[] y = new float[nPoints];

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( System.currentTimeMillis() );

        for (int i = 0; i < nPoints; i++) {
            x[i] = sr.nextInt()*sr.nextFloat();
            y[i] = sr.nextInt()*sr.nextFloat();
        }

    	MultiArrayMergeSort.sortByY(x, y);
    	assertTrue(x.length == nPoints);

        float previousX = x[0];
        float previousY = y[0];
    	for (int i=1; i < x.length; i++) {
            assertTrue( y[i] >= previousY);
            previousX = x[i];
            previousY = y[i];
        }
    	
    	// ===   test   exceptions ====
    	boolean caughtException = true;
    	try {
    	    MultiArrayMergeSort.sortByY(null, y);
    	} catch (Throwable t) {
    	    caughtException = true;
    	}
    	assertTrue(caughtException);
    	
    	caughtException = true;
        try {
            MultiArrayMergeSort.sortByY(x, null);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MultiArrayMergeSort.sortByY(x, Arrays.copyOf(y, y.length - 3));
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        
        
        //===== test w/ int array too ====
        int[] sortedXIndexes = new int[x.length];
        for (int i = 0; i < nPoints; i++) {
            x[i] = sr.nextInt()*sr.nextFloat();
            y[i] = sr.nextInt()*sr.nextFloat();
            sortedXIndexes[i] = i;
        }
        
        float[] cpX = Arrays.copyOf(x, x.length);
        float[] cpY = Arrays.copyOf(y, y.length);        

        MultiArrayMergeSort.sortByY(x, y, sortedXIndexes, x.length);
        assertTrue(x.length == nPoints);

        previousX = x[0];
        previousY = y[0];
        for (int i=1; i < x.length; i++) {
            assertTrue( y[i] >= previousY);
            previousX = x[i];
            previousY = y[i];
            
            int idx = sortedXIndexes[i];
            assertTrue(cpX[idx] == x[i]);
            assertTrue(cpY[idx] == y[i]);
        }
        
        // ===   test   exceptions ====
        caughtException = true;
        try {
            MultiArrayMergeSort.sortByY(null, y, null, x.length);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MultiArrayMergeSort.sortByY(x, null, null, x.length);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = true;
        try {
            MultiArrayMergeSort.sortByY(x, y, null, x.length);
        } catch (Throwable t) {
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


        float[] ex = new float[13];
    	float[] ey = new float[13];
        ex[0] =  2.5f;
        ey[0] =  1.0f;
        ex[1] = 10.5f;
        ey[1] =  2.0f;
        ex[2] = 10.0f;
        ey[2] =  3.0f;
        ex[3] =  0.5f;
        ey[3] =  3.2f;
        ex[4] =  9.0f;
        ey[4] =  4.0f;

        ex[5] = 11.2f;
        ey[5] =  4.0f;
        ex[6] =  4.0f;
        ey[6] =  4.2f;
        ex[7] =  2.0f;
        ey[7] =  5.0f;
        ex[8] =  3.0f;
        ey[8] =  5.3f;
        ex[9] =  5.0f;
        ey[9] =  5.3f;

        ex[10] =  6.0f;
        ey[10] =  6.0f;
        ex[11] =  7.0f;
        ey[11] =  6.0f;
        ex[12] =  2.3f;
        ey[12] =  9.3f;


    	MultiArrayMergeSort.sortByYThenX(x, y);
    	assertTrue(x.length == ex.length);

    	for (int i=0; i < ex.length; i++) {
            assertTrue( Math.abs(ex[i] - x[i]) < 0.01);
            assertTrue( Math.abs(ey[i] - y[i]) < 0.01);
        }
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
