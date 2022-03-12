package algorithms.compGeometry.convexHull;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.sort.MiscSorter;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import java.util.logging.Logger;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

public class PolarAngleQuickSortTest extends TestCase {

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
    
    public void testSortCCWBy1stPoint() {
        /*            2,6
         *
         *
         *     0,2   2,2 3,2
         *                      7,1
         *            2,0
         *         
         */
        int n = 6;
        long[] x = new long[]{0, 2, 7, 2, 2, 3};
        long[] y = new long[]{2, 2, 1, 6, 0, 2};
        
        long[] ex = new long[]{2, 7, 0, 2, 3, 2};
        long[] ey = new long[]{0, 1, 2, 2, 2, 6};
        int[] origIdxs = MiscSorter.sortBy1stArg2(y, x);
                
        assertTrue(Arrays.equals(ex, x));
        assertTrue(Arrays.equals(ey, y));
        
        double[] outPolarAngle = new double[n];
        ex = new long[]{2, 7, 3, 2, 0};
        ey = new long[]{0, 1, 2, 6, 2};
        int nR = PolarAngleQuickSort.sortCCWBy1stPoint(x, y, outPolarAngle);
        assertEquals(5, nR);
        
        x = Arrays.copyOf(x, nR);
        y = Arrays.copyOf(y, nR);
        outPolarAngle = Arrays.copyOf(outPolarAngle, nR);
                
        assertTrue(Arrays.equals(ex, x));
        assertTrue(Arrays.equals(ey, y));
    }

    public void testReduceToUniquePolarAngles() throws Exception {

        /*            2,6
         *
    	 *     0,2    2,2
         *                   4,1
    	 *            2,0
    	 */
        // points should have been sorted by y and decr x already
        PairInt[] points = new PairInt[5];
        points[0] = new PairInt(2, 0);
        points[1] = new PairInt(4, 1);
        points[2] = new PairInt(2, 2);
        points[3] = new PairInt(2, 6);
        points[4] = new PairInt(0, 2);
        
    	//float[] xx = new float[]{2,4,2,2,0};
    	//float[] yy = new float[]{0,1,2,6,2};

        double[] ap = new double[points.length];
        for (int i = 1; i < points.length; i++) {
            ap[i] = AngleUtil.polarAngleCCW(
                (double)(points[i].getX() - points[0].getX()), 
                (double)(points[i].getY() - points[0].getY()));
        }

        int nUsable = PolarAngleQuickSort.reduceToUniquePolarAngles(
            points[0], points, ap);

        double[] expectedxx = {2, 4, 2, 0};
    	double[] expectedyy = {0, 1, 6, 2};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
        }
        
    }

    public void testReduceToUniquePolarAngles2() throws Exception {

        /*
         *    4                         p5
         *
         *
         *    3     p2             p4   <----- expect is 'removed' and n returned is reduced
         *          
         *
         *    2               p1 <----- expect is 'removed' and n returned is reduced
         *
         *
         *    1          p0
         *
         *    0
         *     0    1    2    3    4    5
    	 */

        PairInt[] points = new PairInt[5];
        points[0] = new PairInt(2, 1);
        points[1] = new PairInt(3, 2);
        points[2] = new PairInt(4, 3);
        points[3] = new PairInt(5, 4);
        points[4] = new PairInt(1, 3);
        
        //float[] xx = new float[]{2, 3, 4, 5, 1};
    	//float[] yy = new float[]{1, 2, 3, 4, 2.5f};

        double[] ap = new double[points.length];
        for (int i = 1; i < points.length; i++) {
            ap[i] = AngleUtil.polarAngleCCW(
                (double)(points[i].getX() - points[0].getX()), 
                (double)(points[i].getY() - points[0].getY()));
        }

        int nUsable = PolarAngleQuickSort.reduceToUniquePolarAngles(points[0],
            points, ap);

        double[] expectedxx = {2, 5, 1};
    	double[] expectedyy = {1, 4, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
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

        PairInt[] points = new PairInt[4];
        points[0] = new PairInt(3, 1);
        points[1] = new PairInt(4, 3);
        points[2] = new PairInt(2, 2);
        points[3] = new PairInt(1, 3);
        
        //float[] xx = new float[]{3, 4, 2, 1};
    	//float[] yy = new float[]{1, 3, 2, 3};

        double[] ap = new double[points.length];
        for (int i = 1; i < points.length; i++) {
            ap[i] = AngleUtil.polarAngleCCW(
                (double)(points[i].getX() - points[0].getX()), 
                (double)(points[i].getY() - points[0].getY()));
        }

        int nUsable = PolarAngleQuickSort.reduceToUniquePolarAngles(points[0], 
            points, ap);

        double[] expectedxx = {3, 4, 1};
    	double[] expectedyy = {1, 3, 3};


    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
        }
        
    }

    public void testSort1() {

        /*            2,6
         *
    	 *     0,2    2,2
         *                   7,1
    	 *            2,0
    	 */
        
        PairInt[] points = new PairInt[5];
        points[0] = new PairInt(2, 0);
        points[1] = new PairInt(7, 1);
        points[2] = new PairInt(0, 2);
        points[3] = new PairInt(2, 2);
        points[4] = new PairInt(2, 6);
        
        // points should have been sorted by y and decr x already
    	//float[] xx = {2,7,0,2,2};
    	//float[] yy = {0,1,2,2,6};

        int nUsable = PolarAngleQuickSort.sort(points[0], points);

        float[] expectedxx = {2, 7, 2, 0};
        float[] expectedyy = {0, 1, 6, 2};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
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
        
    	points = new PairInt[6];
        points[0] = new PairInt(2, 0);
        points[1] = new PairInt(7, 1);
        points[2] = new PairInt(0, 2);
        points[3] = new PairInt(2, 2);
        points[4] = new PairInt(3, 2);
        points[5] = new PairInt(2, 6);
        
        expectedxx = new float[]{2.0f, 7.0f, 3.0f, 2.0f, 2, 0.0f};
        expectedyy = new float[]{0.0f, 1.0f, 2.0f, 6.0f, 2, 2.0f};
        
        double[] polarAngle = new double[points.length];
        for (int i = 1; i < points.length; i++) {
            polarAngle[i] = AngleUtil.polarAngleCCW((double)(points[i].getX() - points[0].getX()), 
                (double)(points[i].getY() - points[0].getY()));
        }
        
    	PolarAngleQuickSort.sortByPolarAngle(points, 1, points.length - 1, polarAngle);
        for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
        }
        
        //---------
        points = new PairInt[6];
        points[0] = new PairInt(2, 0);
        points[1] = new PairInt(7, 1);
        points[2] = new PairInt(0, 2);
        points[3] = new PairInt(2, 2);
        points[4] = new PairInt(3, 2);
        points[5] = new PairInt(2, 6);
        
        // ================================================
        
    	nUsable = PolarAngleQuickSort.sort(points[0], points);

    	expectedxx = new float[]{2.0f, 7.0f, 3.0f, 2.0f, 0.0f};
        expectedyy = new float[]{0.0f, 1.0f, 2.0f, 6.0f, 2.0f};
        
        for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
        }
    }

    public void testSort1_2() {

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
        
        List<PairInt> points = new ArrayList<PairInt>();
        points.add(new PairInt(2, 0));
        points.add(new PairInt(7, 1));
        points.add(new PairInt(0, 2));
        points.add(new PairInt(2, 2));
        points.add(new PairInt(3, 2));
        points.add(new PairInt(2, 6));
             
        float[] expectedxx 
            = new float[]{2.0f, 7.0f, 3.0f, 2.0f, 2, 0.0f};
        float[] expectedyy 
            = new float[]{0.0f, 1.0f, 2.0f, 6.0f, 2, 2.0f};
        
        double[] polarAngle = new double[points.size()];
        for (int i = 1; i < points.size(); i++) {
            polarAngle[i] = 
                AngleUtil.polarAngleCCW((double)(
                points.get(i).getX() - points.get(0).getX()), 
                (double)(points.get(i).getY() - 
                points.get(0).getY()));
        }
        
    	PolarAngleQuickSort.sortByPolarAngle(
            points, 1, points.size() - 1, polarAngle);
        for (int i=0; i < points.size(); i++) {
            assertTrue( Math.abs(expectedxx[i] 
                - points.get(i).getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] 
                - points.get(i).getY()) < 0.01);
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
         *    3     p2             p4   <----- expect is 'removed' and n returned is reduced
         *          
         *
         *    2               p1 <----- expect is 'removed' and n returned is reduced
         *
         *
         *    1          p0
         *
         *    0
         *     0    1    2    3    4    5
    	 */

        PairInt[] points = new PairInt[5];
        points[0] = new PairInt(2, 1);
        points[1] = new PairInt(3, 2);
        points[2] = new PairInt(1, 3);
        points[3] = new PairInt(4, 3);
        points[4] = new PairInt(5, 4);
        //float[] xx = {2, 3, 1,   4, 5};
        //float[] yy = {1, 2, 2.5f, 3, 4};

        int nUsable = PolarAngleQuickSort.sort(points[0], points);

        float[] expectedxx = new float[]{2, 5, 1};
        float[] expectedyy = new float[]{1, 4, 3.f};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
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

        PairInt[] points = new PairInt[4];
        points[0] = new PairInt(3, 1);
        points[1] = new PairInt(2, 2);
        points[2] = new PairInt(4, 3);
        points[3] = new PairInt(1, 3);
        //float[] xx = {3, 2, 4, 1};
        //float[] yy = {1, 2, 3, 3};

        int nUsable = PolarAngleQuickSort.sort(points[0], points);

        float[] expectedxx = {3, 4, 1};
        float[] expectedyy = {1, 3, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
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
        PairInt[] points = new PairInt[4];
        points[0] = new PairInt(3, 1);
        points[1] = new PairInt(1, 3);
        points[2] = new PairInt(4, 3);
        points[3] = new PairInt(2, 2);
        //float[] xx = {3, 1, 4, 2};
        //float[] yy = {1, 3, 3, 2};

        int nUsable = PolarAngleQuickSort.sort(points[0], points);

        float[] expectedxx = {3, 4, 1};
        float[] expectedyy = {1, 3, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
        }
        
    }

    public void testSort() {

    	/*            2,5
    	 *     0,2    2,2    4,1
    	 *            2,0
    	 */

        PairInt[] points = new PairInt[5];
        points[0] = new PairInt(2, 0);
        points[1] = new PairInt(4, 1);
        points[2] = new PairInt(2, 2);
        points[3] = new PairInt(0, 2);
        points[4] = new PairInt(2, 5);
        
        //float[] xx = {2,4,2,0,2};
    	//float[] yy = {0,1,2,2,5};

        PolarAngleQuickSort.sort(points[0], points);

        float[] expectedxx = {2, 4, 2, 0};
    	float[] expectedyy = {0, 1, 5, 2};

    	for (int i=0; i < expectedxx.length; i++) {
            assertTrue( Math.abs(expectedxx[i] - points[i].getX()) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - points[i].getY()) < 0.01);
        }
        
    }
    
    public void testSortExceptions() {

        boolean threwException;
        
        PairInt[] points = new PairInt[4];
        points[0] = new PairInt(2, 0);
        points[1] = new PairInt(4, 1);
        points[2] = new PairInt(2, 2);
        points[3] = new PairInt(0, 2);
        //float[] xx = {2,4,2,0};
        //float[] yy = {0,1,2,2,5};

        threwException = false;
        try {
            PairInt[] t = null;
            PolarAngleQuickSort.sort(points[0], t);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        points = new PairInt[1];
        points[0] = new PairInt(1, 1);
        threwException = false;
        int n = 0;
        try {
            n = PolarAngleQuickSort.sort(new PairInt(2, 0), points);
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
        return new TestSuite(PolarAngleQuickSortTest.class);
    }

    /**
     * Set up a Junit test runner
     * @param args Not used.
    */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }

}
