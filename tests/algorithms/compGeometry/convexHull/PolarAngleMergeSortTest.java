package algorithms.compGeometry.convexHull;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairIntArray;

import java.util.logging.Logger;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/convexHull/GrahamScan.java
 * under MIT License (MIT), Nichole King 2013
 */
public class PolarAngleMergeSortTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

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

    public void testReduceToUniquePolarAngles() throws Exception {

        /*            2,6
         *
    	 *     0,2    2,2
         *                   4,1
    	 *            2,0
    	 */
        // points should have been sorted by y and decr x already
    	float[] xx = new float[]{2,4,2,2,0};
    	float[] yy = new float[]{0,1,2,6,2};

        double[] ap = new double[xx.length];
        for (int i = 1; i < xx.length; i++) {
            ap[i] = AngleUtil.polarAngleCCW(
                (double)(xx[i] - xx[0]), 
                (double)(yy[i] - yy[0]));
        }

        int nUsable = PolarAngleMergeSort.reduceToUniquePolarAngles(xx[0], 
            yy[0], xx, yy, ap);

        double[] expectedxx = {2, 4, 2, 0};
    	double[] expectedyy = {0, 1, 6, 2};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        // =================================
        xx = new float[]{2,4,2,2,0};
        yy = new float[]{0,1,2,6,2};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }

    public void testReduceToUniquePolarAngles2() throws Exception {

        /*
         *    4                         p5
         *
         *
         *    3                    p4   <----- expect is 'removed' and n returned is reduced
         *          p2
         *
         *    2               p1 <----- expect is 'removed' and n returned is reduced
         *
         *
         *    1          p0
         *
         *    0
         *     0    1    2    3    4    5
    	 */

        float[] xx = new float[]{2, 3, 4, 5, 1};
    	float[] yy = new float[]{1, 2, 3, 4, 2.5f};

        double[] ap = new double[xx.length];
        for (int i = 1; i < xx.length; i++) {
            ap[i] = AngleUtil.polarAngleCCW(
                (double)(xx[i] - xx[0]), 
                (double)(yy[i] - yy[0]));
        }

        int nUsable = PolarAngleMergeSort.reduceToUniquePolarAngles(xx[0], yy[0], xx, yy, ap);

        double[] expectedxx = {2, 5, 1};
    	double[] expectedyy = {1, 4, 2.5};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        // =================================
        xx = new float[]{20, 30, 40, 50, 10};
        yy = new float[]{10, 20, 30, 40, 25};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }

    public void testReduceToUniquePolarAngles3() throws Exception {

        /*
         *    4
         *
         *
         *    3     p3             p2
         *
         *
         *    2          p1  <----- expect this is 'removed' and returned nItems truncated
         *
         *
         *    1               p0
         *
         *    0
         *     0    1    2    3    4
    	 */

        float[] xx = new float[]{3, 4, 2, 1};
    	float[] yy = new float[]{1, 3, 2, 3};

        double[] ap = new double[xx.length];
        for (int i = 1; i < xx.length; i++) {
            ap[i] = AngleUtil.polarAngleCCW(
                (double)(xx[i] - xx[0]), 
                (double)(yy[i] - yy[0]));
        }

        int nUsable = PolarAngleMergeSort.reduceToUniquePolarAngles(xx[0], yy[0], xx, yy, ap);

        double[] expectedxx = {3, 4, 1};
    	double[] expectedyy = {1, 3, 3};


    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        // =================================
        xx = new float[]{3, 4, 2, 1};
        yy = new float[]{1, 3, 2, 3};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }

    public void testSort1() {

        /*            2,6
         *
    	 *     0,2    2,2
         *                   7,1
    	 *            2,0
    	 */
        // points should have been sorted by y and decr x already
    	float[] xx = {2,7,0,2,2};
    	float[] yy = {0,1,2,2,6};

        int nUsable = PolarAngleMergeSort.sort(2, 0, xx, yy);

        float[] expectedxx = {2, 7, 2, 0};
        float[] expectedyy = {0, 1, 6, 2};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        //---------
        xx = new float[]{2,7,0,2,2};
    	yy = new float[]{0,1,2,2,6};
    	
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
        
    	// ======================================================================
    	/*            2,6
        *
        *
        *     0,2    2,2 3,2
        *                      7,1
        *            2,0
        *    7
        *    6   <>
        *    5
        *    4
        *    3
        *    <> <> <>
        *    1             <>
        *    0 1<> 3 4 5 6 7
        */
    	
    	xx = new float[] {2.0f, 7.0f, 0.0f, 2.0f, 3.0f, 2.0f};
    	yy = new float[] {0.0f, 1.0f, 2.0f, 2.0f, 2.0f, 6.0f};
    	
    	
        expectedxx = new float[]{2.0f, 7.0f, 3.0f, 2.0f, 2.0f, 0.0f};
        expectedyy = new float[]{0.0f, 1.0f, 2.0f, 2.0f, 6.0f, 2.0f};
        
        double[] polarAngle = new double[xx.length];
        for (int i = 1; i < xx.length; i++) {
            polarAngle[i] = AngleUtil.polarAngleCCW((double)(xx[i] - xx[0]), 
                (double)(yy[i] - yy[0]));
        }
    	PolarAngleMergeSort.sort(2.0f, 0.0f, xx, yy, 1, xx.length - 1, polarAngle);
        for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        //---------
        xx = new float[] {2.0f, 7.0f, 0.0f, 2.0f, 3.0f, 2.0f};
    	yy = new float[] {0.0f, 1.0f, 2.0f, 2.0f, 2.0f, 6.0f};
    	
        xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
        expectedxx = new float[]{2.0f, 7.0f, 3.0f, 2.0f, 2.0f, 0.0f};
        expectedyy = new float[]{0.0f, 1.0f, 2.0f, 2.0f, 6.0f, 2.0f};
        
        
        // ================================================
        
    	nUsable = PolarAngleMergeSort.sort(2, 0, xx, yy);

    	expectedxx = new float[]{2.0f, 7.0f, 3.0f, 2.0f, 0.0f};
        expectedyy = new float[]{0.0f, 1.0f, 2.0f, 6.0f, 2.0f};
        
        for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
    }

    public void testSort2() {

        // p0 has already been found to be the lowest y of the points and
        //   for that y, the lowest x value
        //
        // -- points should be sorted by counterclockwise polar angle from p0
        // -- more than one points w/ same angle are removed excepting the furthest from p0

    	/*
         *    4                         p5
         *
         *
         *    3                    p4   <----- expect is 'removed' and n returned is reduced
         *          p2
         *
         *    2               p1 <----- expect is 'removed' and n returned is reduced
         *
         *
         *    1          p0
         *
         *    0
         *     0    1    2    3    4    5
    	 */

        float[] xx = {2, 3, 1,   4, 5};
        float[] yy = {1, 2, 2.5f, 3, 4};

        int nUsable = PolarAngleMergeSort.sort(2, 1, xx, yy);

        float[] expectedxx = new float[]{2, 5, 1};
        float[] expectedyy = new float[]{1, 4, 2.5f};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        //------------
        xx = new float[]{20, 30, 10, 40, 50};
        yy = new float[]{10, 20, 25, 30, 40};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }

    public void testSort3() {

        // p0 has already been found to be the lowest y of the points and
        //   for that y, the lowest x value
        //
        // -- points should be sorted by counterclockwise polar angle from p0
        // -- more than one points w/ same angle are removed excepting the furthest from p0

    	/*
         *    4
         *
         *
         *    3     p3             p2
         *
         *
         *    2          p1  <----- expect this is 'removed' and returned nItems truncated
         *
         *
         *    1               p0
         *
         *    0
         *     0    1    2    3    4
    	 */

        float[] xx = {3, 2, 4, 1};
        float[] yy = {1, 2, 3, 3};

        int nUsable = PolarAngleMergeSort.sort(3, 1, xx, yy);

        float[] expectedxx = {3, 4, 1};
        float[] expectedyy = {1, 3, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        //------------
        xx = new float[]{3, 4, 1};
        yy = new float[]{1, 3, 3};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }

    public void testSort3e() {

        // p0 has already been found to be the lowest y of the points and
        //   for that y, the lowest x value
        //
        // -- points should be sorted by counterclockwise polar angle from p0
        // -- more than one points w/ same angle are removed excepting the furthest from p0

    	/*
         *    4
         *
         *
         *    3     p3             p2
         *
         *
         *    2          p1  <----- expect this is 'removed' and returned nItems truncated
         *
         *
         *    1               p0
         *
         *    0
         *     0    1    2    3    4
    	 */

        float[] xx = {3, 1, 4, 2};
        float[] yy = {1, 3, 3, 2};

        int nUsable = PolarAngleMergeSort.sort(3, 1, xx, yy);

        float[] expectedxx = {3, 4, 1};
        float[] expectedyy = {1, 3, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        //------------
        xx = new float[]{3, 1, 4, 2};
        yy = new float[]{1, 3, 3, 2};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }

    public void testSort() {

    	/*            2,5
    	 *     0,2    2,2    4,1
    	 *            2,0
    	 */

        float[] xx = {2,4,2,0,2};
    	float[] yy = {0,1,2,2,5};

        PolarAngleMergeSort.sort(2, 0, xx, yy);

        float[] expectedxx = {2, 4, 2, 0};
    	float[] expectedyy = {0, 1, 5, 2};

    	for (int i=0; i < expectedxx.length; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
        
        //------------
        xx = new float[]{2,4,2,0,2};
        yy = new float[]{0,1,2,2,5};
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < xx.length; i++) {
            xy.add(Math.round(xx[i]), Math.round(yy[i]));
        }
        
    }
    
    public void testSortExceptions() {

        boolean threwException;
        float[] xx = {2,4,2,0};
        float[] yy = {0,1,2,2,5};

        threwException = false;
        try {
            PolarAngleMergeSort.sort(2, 0, null, yy);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            PolarAngleMergeSort.sort(2, 0, xx, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            PolarAngleMergeSort.sort(2, 0, xx, yy);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        int n = 0;
        try {
            n = PolarAngleMergeSort.sort(2, 0, new float[]{1.0f}, new float[]{1.0f});
            assertTrue(n == 1);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(n == 1);
        assertFalse(threwException);

    }
    
    public void test0() throws Exception {
        // make coverage reports happy: work around for class built purely for static methods. constructor not directly used...
        PolarAngleMergeSort t = new PolarAngleMergeSort();
        assertNotNull(t);
    }

    /**
     * Test suite
     * @return static Test
    */
    public static Test suite(){
        return new TestSuite(PolarAngleMergeSortTest.class);
    }

    /**
     * Set up a Junit test runner
     * @param args Not used.
    */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }

}
