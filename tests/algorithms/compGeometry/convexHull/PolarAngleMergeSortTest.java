package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.LinesAndAngles;
import java.util.logging.Logger;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

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

        /*            2,5
         *
    	 *     0,2    2,2
         *                   4,1
    	 *            2,0
    	 */
        // points should have been sorted by y and decr x already
    	double[] xx = {2,4,2,2,0};
    	double[] yy = {0,1,2,5,2};

        double[] ap = new double[xx.length];
        for (int i = 0; i < xx.length; i++) {
            ap[i] = LinesAndAngles.calculatePolarSineTheta(xx[0], yy[0], xx[i], yy[i]);
        }

        int nUsable = PolarAngleMergeSort.reduceToUniquePolarAngles(xx[0], yy[0], xx, yy, ap);

        double[] expectedxx = {2, 4, 2, 0};
    	double[] expectedyy = {0, 1, 5, 2};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
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

        double[] xx = {2, 3, 4, 5, 1};
    	double[] yy = {1, 2, 3, 4, 2.5};

        double[] ap = new double[xx.length];
        for (int i = 0; i < xx.length; i++) {
            ap[i] = LinesAndAngles.calculatePolarSineTheta(xx[0], yy[0], xx[i], yy[i]);
        }

        int nUsable = PolarAngleMergeSort.reduceToUniquePolarAngles(xx[0], yy[0], xx, yy, ap);

        double[] expectedxx = {2, 5, 1};
    	double[] expectedyy = {1, 4, 2.5};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
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

        double[] xx = {3, 4, 2, 1};
    	double[] yy = {1, 3, 2, 3};

        double[] ap = new double[xx.length];
        for (int i = 0; i < xx.length; i++) {
            ap[i] = LinesAndAngles.calculatePolarSineTheta(xx[0], yy[0], xx[i], yy[i]);
        }

        int nUsable = PolarAngleMergeSort.reduceToUniquePolarAngles(xx[0], yy[0], xx, yy, ap);

        double[] expectedxx = {3, 4, 1};
    	double[] expectedyy = {1, 3, 3};


    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
    }

    public void testSort1() {

        /*            2,5
         *
    	 *     0,2    2,2
         *                   4,1
    	 *            2,0
    	 */
        // points should have been sorted by y and decr x already
    	double[] xx = {2,4,0,2,2};
    	double[] yy = {0,1,2,2,5};

        int nUsable = PolarAngleMergeSort.sort(2, 0, xx, yy);

        double[] expectedxx = {2, 4, 2, 0};
    	double[] expectedyy = {0, 1, 5, 2};

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

        double[] xx = {2, 3, 1,   4, 5};
    	double[] yy = {1, 2, 2.5, 3, 4};

        int nUsable = PolarAngleMergeSort.sort(2, 1, xx, yy);

        double[] expectedxx = {2, 5, 1};
    	double[] expectedyy = {1, 4, 2.5};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
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

        double[] xx = {3, 2, 4, 1};
    	double[] yy = {1, 2, 3, 3};

        int nUsable = PolarAngleMergeSort.sort(3, 1, xx, yy);

        double[] expectedxx = {3, 4, 1};
    	double[] expectedyy = {1, 3, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
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

        double[] xx = {3, 1, 4, 2};
    	double[] yy = {1, 3, 3, 2};

        int nUsable = PolarAngleMergeSort.sort(3, 1, xx, yy);

        double[] expectedxx = {3, 4, 1};
    	double[] expectedyy = {1, 3, 3};

    	for (int i=0; i < nUsable; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
    }

    public void testSort() {

    	/*            2,5
    	 *     0,2    2,2    4,1
    	 *            2,0
    	 */

        double[] xx = {2,4,2,0,2};
    	double[] yy = {0,1,2,2,5};

        PolarAngleMergeSort.sort(2, 0, xx, yy);

        double[] expectedxx = {2, 4, 2, 0};
    	double[] expectedyy = {0, 1, 5, 2};

    	for (int i=0; i < expectedxx.length; i++) {
            assertTrue( Math.abs(expectedxx[i] - xx[i]) < 0.01);
            assertTrue( Math.abs(expectedyy[i] - yy[i]) < 0.01);
        }
    }

    /**
     * Test suite
     * @return static Test
    */
    public static Test suite(){
        System.out.println("Creating a TestSuite for PolarAngleMergeSort");
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
