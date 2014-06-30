package algorithms.compGeometry;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LinesAndAnglesTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public LinesAndAnglesTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of crossProduct method, of class LinesAndAngles.
     */
    public void testCrossProduct() {

        log.info("crossProduct");

        /*         |
         *       4 -
         *         |
         *       3 -  o2
         *         |
         * 1 o   2 -
         *         |
         *       1 -
         *         |
         *       0 o--|--|--|--|--|--|
         *         0  1  2  3  4
         *
         */
        int x1 = -3;
        int y1 = 2;
        int x2 = 1;
        int y2 = 3;
        int result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result < 0); // clockwise
        assertTrue(result == -11);

        x1 = 3;
        result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result >= 0); // counter clockwise
        assertTrue(result == 7);
    }

    /**
     * Test of linesIntersect method, of class LinesAndAngles.
     */
    public void testLinesIntersect() {
        log.info("linesIntersect");
        /*               @
         *       |
         *       |
         *          *   
         *       |
         *       |
         *       |
         *  -----@-------
         *       |
         *       *
         *       |
         */
        int x1 = 0;
        int y1 = -1;
        int x2 = 1;
        int y2 = 2;

        int x3 = 0;
        int y3 = 0;
        int x4 = 2;
        int y4 = 3;

        boolean result = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertTrue(result);
    }

    /**
     * Test of linesIntersect method, of class LinesAndAngles.
     */
    public void testLinesIntersect2() {
        log.info("result");
        /*    
         *       |
         *       |  @
         *       *   
         *       |
         *       |
         *       |
         *  -----*--@----
         *       |
         *       |
         *       |
         */
        int x1 = 0;
        int y1 = 0;
        int x2 = 0;
        int y2 = 2;

        int x3 = 1;
        int y3 = 0;
        int x4 = 1;
        int y4 = 3;

        boolean result = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertFalse(result);
    }
    
    public void testLinesIntersect3() {
        log.info("testLinesIntersect3");
        /*    
         *       |
         *       |  
         *       *   
         *       |
         *       |
         *       |
         *  -----*------
         *       |
         *       @
         *       |
         */
        int x1 = 0;
        int y1 = 0;
        int x2 = 0;
        int y2 = 2;

        int x3 = 0;
        int y3 = -1;
        int x4 = 0;
        int y4 = 3;

        boolean result = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertTrue(result);
    }
    public void testLinesIntersect4() {
        log.info("testLinesIntersect4");
        /*    
         *       |
         *       |
         *  ---@-*-*--@-
         *       |
         *       |
         */
        int x1 = 0;
        int y1 = 0;
        int x2 = 2;
        int y2 = 0;

        int x3 = -1;
        int y3 = 0;
        int x4 = 3;
        int y4 = 0;

        boolean result = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertTrue(result);


        int[] seed0 = new int[]{2, 1};
        int[] line0 = new int[]{3, 2, 4, 0};
        int[] seed1 = new int[]{4, 1};        
        result = LinesAndAngles.linesIntersect(line0[0], line0[1], line0[2], line0[3], seed0[0], seed0[1], seed1[0], seed1[1]);
        assertTrue(result);

    }
 
    public void testIntersectionOf2Lines() {

        /*
        0                 <p4>

        1

        2     <p1>

        3               P

        4                     <p2>

        5

        6          <p3>

        7

        8

        9

       10   1   2   3   4   5   6   7   8   9   10
         */

        float x1 = 2; float y1 = 2;
        float x2 = 6; float y2 = 4;

        float x3 = 3; float y3 = 6;
        float x4 = 5; float y4 = 0;

        double[] p = LinesAndAngles.intersectionOf2Lines(x1, y1, x2, y2, x3, y3, x4, y4);

        double[] p2 = LinesAndAngles.intersectionOf2Lines_2(x1, y1, x2, y2, x3, y3, x4, y4);

        assertNotNull(p);
        assertTrue(p.length == 2);
        assertTrue( Math.abs(p[0] - 4.) < 0.01);
        assertTrue( Math.abs(p[1] - 3.) < 0.01);

        assertNotNull(p2);
        assertTrue(p2.length == 2);
        assertTrue( Math.abs(p2[0] - 4.) < 0.01);
        assertTrue( Math.abs(p2[1] - 3.) < 0.01);


        p = LinesAndAngles.intersectionOf2Lines(x2, y2, x1, y1, x3, y3, x4, y4);
        p2 = LinesAndAngles.intersectionOf2Lines_2(x1, y1, x2, y2, x3, y3, x4, y4);

        assertNotNull(p);
        assertTrue(p.length == 2);
        assertTrue( Math.abs(p[0] - 4.) < 0.01);
        assertTrue( Math.abs(p[1] - 3.) < 0.01);

        assertNotNull(p2);
        assertTrue(p2.length == 2);
        assertTrue( Math.abs(p2[0] - 4.) < 0.01);
        assertTrue( Math.abs(p2[1] - 3.) < 0.01);


        p = LinesAndAngles.intersectionOf2Lines(x1, y1, x2, y2, x4, y4, x3, y3);
        p2 = LinesAndAngles.intersectionOf2Lines_2(x1, y1, x2, y2, x3, y3, x4, y4);

        assertNotNull(p);
        assertTrue(p.length == 2);
        assertTrue( Math.abs(p[0] - 4.) < 0.01);
        assertTrue( Math.abs(p[1] - 3.) < 0.01);

        assertNotNull(p2);
        assertTrue(p2.length == 2);
        assertTrue( Math.abs(p2[0] - 4.) < 0.01);
        assertTrue( Math.abs(p2[1] - 3.) < 0.01);

        p = LinesAndAngles.intersectionOf2Lines(x2, y2, x1, y1, x4, y4, x3, y3);
        p2 = LinesAndAngles.intersectionOf2Lines_2(x1, y1, x2, y2, x3, y3, x4, y4);

        assertNotNull(p);
        assertTrue(p.length == 2);
        assertTrue( Math.abs(p[0] - 4.) < 0.01);
        assertTrue( Math.abs(p[1] - 3.) < 0.01);

        assertNotNull(p2);
        assertTrue(p2.length == 2);
        assertTrue( Math.abs(p2[0] - 4.) < 0.01);
        assertTrue( Math.abs(p2[1] - 3.) < 0.01);
    }

    public void testIntersectionOf2Lines_2() {

        // results for 2 lines that do not intersect

        /*
        0      <p3>

        1

        2      <p1>            <p4>

        3               P

        4                      <p2>

        ...

       10   1   2   3   4   5   6   7   8   9   10
         */

        float x1 = 2; float y1 = 2;
        float x2 = 6; float y2 = 4;

        float x3 = 2; float y3 = 0;
        float x4 = 6; float y4 = 2;

        double[] p = LinesAndAngles.intersectionOf2Lines(x1, y1, x2, y2, x3, y3, x4, y4);

        assertNull(p);


        p = LinesAndAngles.intersectionOf2Lines(x2, y2, x1, y1, x3, y3, x4, y4);

        assertNull(p);


        p = LinesAndAngles.intersectionOf2Lines(x1, y1, x2, y2, x4, y4, x3, y3);

        assertNull(p);


        p = LinesAndAngles.intersectionOf2Lines(x2, y2, x1, y1, x4, y4, x3, y3);

        assertNull(p);
    }

    public void testCalculateSlope() {
        float x1 = 6.f;
        float y1 = 2.f;
        float x2 = 8.f;
        float y2 = 4.f;

        float slope = LinesAndAngles.calculateSlope(x1, y1, x2, y2);
        assertTrue(Math.abs(slope - 1.) < 0.01);
    }

    public void testCalculateSlope2() {

        float x1 = 8.f; float y1 = 4.f;
        float x2 = 2.f; float y2 = 6.f;

        float ymin = 0.f;
        float ymax = 10.f;

        float xmin = 0.f;
        float xmax = 10.f;

        float[] endpointsArray = LinesAndAngles.calculatePerpendicularBisectingSegment(x1, y1, x2, y2, xmin, xmax, ymin, ymax);

        assertTrue( Math.abs(endpointsArray[0] - 3.33) < 0.01);
        assertTrue( Math.abs(endpointsArray[1] - 0.) < 0.01);
        assertTrue( Math.abs(endpointsArray[2] - 6.66) < 0.01);
        assertTrue( Math.abs(endpointsArray[3] - 10.) < 0.01);
    }

    public void testDirection() throws Exception {

        /* 0   1   2   3   4   5   6   7   8   9   10  11

        0
   direction0=2.6         .     direction1=-2.4
        1                .
            seed0 <>    . <> seed1
        2              .

        3
                  <> seed2   direction2=2.25
        4   ..................
                  <> seed3   direction3=-2.25
        5
                   .
        6           .
  direction4=-1.5    .
        7    seed4  <>
                       .  <> seed5
        8               .      direction5=2.25
                         .
        9

        10
            0   1   2   3   4   5   6   7   8   9   10  11 */

        float[] seed0 = new float[]{2.0f, 1.5f};
        float[] line0 = new float[]{3.0f, 2.0f, 4.2f, 0.f};

        float[] seed1 = new float[]{4.5f, 1.5f};


        float[] seed2 = new float[]{2.0f, 3.5f};
        float[] line2 = new float[]{0.0f, 4.0f, 4.5f, 4.f};

        float[] seed3 = new float[]{2.0f, 4.5f};

        float[] seed4 = new float[]{2.0f, 7.0f};
        float[] line4 = new float[]{2.0f, 6.0f, 3.5f, 9.f};
        float[] seed5 = new float[]{3.5f, 7.5f};


        double direction0 = LinesAndAngles.direction(
            line0[0], line0[1], line0[2], line0[3], seed0[0], seed0[1]);

        double direction1 = LinesAndAngles.direction(
            line0[0], line0[1], line0[2], line0[3], seed1[0], seed1[1]);

        double direction2 = LinesAndAngles.direction(
            line2[0], line2[1], line2[2], line2[3], seed2[0], seed2[1]);

        double direction3 = LinesAndAngles.direction(
            line2[0], line2[1], line2[2], line2[3], seed3[0], seed3[1]);

        double direction4 = LinesAndAngles.direction(
            line4[0], line4[1], line4[2], line4[3], seed4[0], seed4[1]);

        double direction5 = LinesAndAngles.direction(
            line4[0], line4[1], line4[2], line4[3], seed5[0], seed5[1]);

        log.info("direction0=" + direction0 + "\ndirection1=" + direction1
            + "\ndirection2=" + direction2 + "\ndirection3=" + direction3
            + "\ndirection4=" + direction4 + "\ndirection5=" + direction5);
        
        
    }
    
    public void testSegmentIsWithinSegment() throws Exception {

        // ==== horizontal line
        /*   segment completely within other segment
         *    0   2   3   1
         */
        float x0, x1, x2, x3, y0, y1, y2, y3;
        boolean isWithin;
        float eps = 0.01f;

        x0 = 10.0f;
        y0 = 10.0f;
        x1 = 40.0f;
        y1 = y0;
        x2 = 20.f;
        y2 = y0;
        x3 = 30.f;
        y3 = y0;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertTrue(isWithin);

        // horizontal line where a segment is not completely contained within another
        /*
         *    0   2   1   3
         */
        x0 = 10.0f;
        y0 = 10.0f;
        x1 = 30.0f;
        y1 = y0;
        x2 = 20.f;
        y2 = y0;
        x3 = 40.f;
        y3 = y0;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        // horizontal line where a segment is completely outside the another
        /*
         *    0   1   2   3
         */
        x0 = 10.0f;
        y0 = 10.0f;
        x1 = 20.0f;
        y1 = y0;
        x2 = 30.f;
        y2 = y0;
        x3 = 40.f;
        y3 = y0;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);


        // ==== vertical line   segment completely within other segment
        /*    0
         *    2
         *    3
         *    1   */
        x0 = 40.0f;
        y0 = 40.0f;
        x1 = x0;
        y1 = 10.f;
        x2 = x0;
        y2 = 30.f;
        x3 = x0;
        y3 = 20.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertTrue(isWithin);


        /*    0   segment is not completely contained within another
         *    2
         *    1
         *    3   */
        x0 = 40.0f;
        y0 = 40.0f;
        x1 = x0;
        y1 = 20.f;
        x2 = x0;
        y2 = 30.f;
        x3 = x0;
        y3 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        /*    0     a segment is completely outside the another
         *    1
         *    2
         *    3   */
        x0 = 40.0f;
        y0 = 40.0f;
        x1 = x0;
        y1 = 30.f;
        x2 = x0;
        y2 = 20.f;
        x3 = x0;
        y3 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        // ==== diagonal line   segment completely within other segment
        /*           0
         *        2
         *     3
         *  1   */
        x0 = 40.0f;
        y0 = 40.0f;
        x1 = 10.f;
        y1 = 10.f;
        x2 = 30.f;
        y2 = 30.f;
        x3 = 20.f;
        y3 = 20.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertTrue(isWithin);


        /*          0   segment is not completely contained within another
         *        2
         *      1
         *    3   */
        x0 = 40.0f;
        y0 = 40.0f;
        x1 = 20.f;
        y1 = 20.f;
        x2 = 30.f;
        y2 = 30.f;
        x3 = 10.f;
        y3 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        /*          0     a segment is completely outside the another
         *        1
         *      2
         *    3   */
        x0 = 40.0f;
        y0 = 40.0f;
        x1 = 30.f;
        y1 = 30.f;
        x2 = 20.f;
        y2 = 20.f;
        x3 = 10.f;
        y3 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);


        // ==== diagonal line   segment completely within other segment
        /*  0
         *    2
         *      3
         *        1   */
        x0 = 10.0f;
        y0 = 40.0f;
        x2 = 20.f;
        y2 = 30.f;
        x3 = 30.f;
        y3 = 20.f;
        x1 = 40.f;
        y1 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertTrue(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertTrue(isWithin);


        /*  0   segment is not completely contained within another
         *    2
         *      1
         *        3   */
        x0 = 10.0f;
        y0 = 40.0f;
        x2 = 20.f;
        y2 = 30.f;
        x1 = 30.f;
        y1 = 20.f;
        x3 = 40.f;
        y3 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);


        /*  0   a segment is completely outside the another
         *    1
         *      2
         *        3   */
        x0 = 10.0f;
        y0 = 40.0f;
        x1 = 20.f;
        y1 = 30.f;
        x2 = 30.f;
        y2 = 20.f;
        x3 = 40.f;
        y3 = 10.f;
        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x2, y2, x3, y3, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x0, y0, x1, y1, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

        isWithin = LinesAndAngles.segmentIsWithinSegment(x1, y1, x0, y0, x3, y3, x2, y2, eps);
        assertFalse(isWithin);

    }

    /**
     * Test of crossProduct method, of class LinesAndAngles.
     */
    public void testDistSquared() {

        log.info("testDistSquared");

        /*         |
         *       4 -
         *         |
         *       3 -  o2
         *         |
         * 1 o   2 -
         *         |
         *       1 -
         *         |
         *       0 o--|--|--|--|--|--|
         *         0  1  2  3  4
         */
        int x1 = -3;
        int y1 = 2;
        int x2 = 1;
        int y2 = 3;
        double expected = (4.*4. + 1);
        double result = LinesAndAngles.distSquared(x1, y1, x2, y2);
        assertTrue( Math.abs(result - expected) < 0.01);

        x1 = 3;
        expected = (2.*2. + 1);
        result = LinesAndAngles.distSquared(x1, y1, x2, y2);
        assertTrue( Math.abs(result - expected) < 0.01);
    }

    /**
     *   *  |  *
     *      |
     *      |.....  <---- angle is w.r.t y=0, x=xc.  increases in CCW order
     *  *
     *         *
     *
     * @param xc
     * @param yc
     * @param radius
     * @param angleInDegreesFromYEQ0XGT0  angle in degrees, CCW from point y=0, x=xc
     * @return
     */
    static float[] calculateXAndYFromXcYcAndRadiusCCW(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        return calculateXAndYFromXcYcAndRadius(xc, yc, radius, 360 - angleInDegreesFromYEQ0XGT0);
    }
    public static float[] calculateXAndYFromXcYcAndRadius(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        double dx = radius * Math.cos(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));
        double dy = radius * Math.sin(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));

        float x = (float) (xc + dx);
        float y = (float) (yc - dy);
        return new float[]{x, y};
    }
}
