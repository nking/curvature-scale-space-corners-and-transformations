package algorithms.compGeometry;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PointInPolygonTest extends TestCase {

    public PointInPolygonTest(String testName) {
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

    public void getSquarePoly(float[] xPoly, float[] yPoly) {

        xPoly[0] = 0;
        yPoly[0] = 0;
        xPoly[1] = 10;
        yPoly[1] = 0;

        xPoly[2] = 10;
        yPoly[2] = 0;
        xPoly[3] = 10;
        yPoly[3] = 10;

        xPoly[4] = 10;
        yPoly[4] = 10;
        xPoly[5] = 0;
        yPoly[5] = 10;

        xPoly[6] = 0;
        yPoly[6] = 10;
        xPoly[7] = 0;
        yPoly[7] = 0;
    }

    public void getSquareHolePoly(float[] xPoly, float[] yPoly) {

        xPoly[0] = 0;
        yPoly[0] = 0;
        xPoly[1] = 10;
        yPoly[1] = 0;

        xPoly[2] = 10;
        yPoly[2] = 0;
        xPoly[3] = 10;
        yPoly[3] = 10;

        xPoly[4] = 10;
        yPoly[4] = 10;
        xPoly[5] = 0;
        yPoly[5] = 10;

        xPoly[6] = 0;
        yPoly[6] = 10;
        xPoly[7] = 0;
        yPoly[7] = 0;

        xPoly[8] = 2.5f;
        yPoly[8] = 2.5f;
        xPoly[9] = 7.5f;
        yPoly[9] = 2.5f;

        xPoly[10] = 7.5f;
        yPoly[10] = 2.5f;
        xPoly[11] = 7.5f;
        yPoly[11] = 7.5f;

        xPoly[12] = 7.5f;
        yPoly[12] = 7.5f;
        xPoly[13] = 2.5f;
        yPoly[13] = 7.5f;

        xPoly[14] = 2.5f;
        yPoly[14] = 7.5f;
        xPoly[15] = 2.5f;
        yPoly[15] = 2.5f;
    }

    public void getStrangePoly(float[] xPoly, float[] yPoly) {

        xPoly[0] = 0;
        yPoly[0] = 0;
        xPoly[1] = 2.5f;
        yPoly[1] = 2.5f;

        xPoly[2] = 2.5f;
        yPoly[2] = 2.5f;
        xPoly[3] = 0;
        yPoly[3] = 10;

        xPoly[4] = 2.5f;
        yPoly[4] = 7.5f;
        xPoly[5] = 7.5f;
        yPoly[5] = 7.5f;

        xPoly[6] = 7.5f;
        yPoly[6] = 7.5f;
        xPoly[7] = 10;
        yPoly[7] = 10;

        xPoly[8] = 10;
        yPoly[8] = 10;
        xPoly[9] = 10;
        yPoly[9] = 0;

        xPoly[10] = 10;
        yPoly[10] = 0;
        xPoly[11] = 2.5f;
        yPoly[11] = 2.5f;
    }

    public void getExagonPoly(float[] xPoly, float[] yPoly) {

        xPoly[0] = 3;
        yPoly[0] = 0;
        xPoly[1] = 7;
        yPoly[1] = 0;

        xPoly[2] = 7;
        yPoly[2] = 0;
        xPoly[3] = 10;
        yPoly[3] = 5;

        xPoly[4] = 10;
        yPoly[4] = 5;
        xPoly[5] = 7;
        yPoly[5] = 10;

        xPoly[6] = 7;
        yPoly[6] = 10;
        xPoly[7] = 3;
        yPoly[7] = 10;

        xPoly[8] = 3;
        yPoly[8] = 10;
        xPoly[9] = 0;
        yPoly[9] = 5;

        xPoly[10] = 0;
        yPoly[10] = 5;
        xPoly[11] = 3;
        yPoly[11] = 0;
    }

     public void getTestPts(float[] xPoly, float[] yPoly) {

         xPoly[0] = 5;
         yPoly[0] = 5;
         xPoly[1] = 5;
         yPoly[1] = 8;

         xPoly[2] = -10;
         yPoly[2] = 5;
         xPoly[3] = 0;
         yPoly[3] = 5;

         xPoly[4] = 10;
         yPoly[4] = 5;
         xPoly[5] = 8;
         yPoly[5] = 5;

         xPoly[6] = 10;
         yPoly[6] = 10;
     }

    /**
     * Test of isInSimpleCurve method, of class PointInPolygon.
     */
    public void testIsInSimpleCurve() {

        System.out.println("testIsInSimpleCurve");

        float[] squarePolyX = new float[8];
        float[] squarePolyY = new float[8];
        getSquarePoly(squarePolyX, squarePolyY);


        float[] squareHolePolyX = new float[16];
        float[] squareHolePolyY = new float[16];
        getSquareHolePoly(squareHolePolyX, squareHolePolyY);


        float[] strangePolyX = new float[12];
        float[] strangePolyY = new float[12];
        getStrangePoly(strangePolyX, strangePolyY);


        float[] exagonPolyX = new float[12];
        float[] exagonPolyY = new float[12];
        getExagonPoly(exagonPolyX, exagonPolyY);


        float[] testPointsX = new float[7];
        float[] testPointsY = new float[7];
        getTestPts(testPointsX, testPointsY);


        PointInPolygon instance = new PointInPolygon();

        // test square poly
        boolean[] expected = new boolean[] {
            true, true, false, false, true, true, false
        };
        int count=0;
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, squarePolyX, squarePolyY, squarePolyX.length);

            System.out.println("  result=" + result + " expect=" + expected[count]);
            assertTrue(result == expected[count]);
            count++;
        }

        // test square hole poly
        expected = new boolean[] {
            false, true, false, false, true, true, false
        };
        count=0;
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, squareHolePolyX, squareHolePolyY, squareHolePolyX.length);

            System.out.println("  result=" + result + " expect=" + expected[count]);
            assertTrue(result == expected[count]);
            count++;
        }

        // test strange poly
        expected = new boolean[] {
            true, false, false, false, true, true, false
        };
        count=0;
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, strangePolyX, strangePolyY, strangePolyX.length);

            System.out.println("  result=" + result + " expect=" + expected[count]);
            assertTrue(result == expected[count]);
            count++;
        }

        // test exagon poly
        expected = new boolean[] {
            true, true, false, false, false /*it's on exiting edge*/, false /*it's on exiting edge*/, false
        };
        count=0;
        /*for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, exagonPolyX, exagonPolyY, exagonPolyX.length);

            System.out.println("  result=" + result + " expect=" + expected[count]);
            assertTrue(result == expected[count]);
            count++;
        }*/
    }

    public void test10() throws Exception {

         /* 0   1   2   3   4   5   6   7   8   9   10  11
        0   -----------------------------------------
            |              .                        |
        1   |             .                         |
            | seed0 <>   . <> seed1                 |
        2   |           .                           |
            |          .                            |
        3   |         .                             |
                  <> seed2
        4   ..................
                  <> seed3
        5
                   .
        6           .
                     .
        7    seed4  <>
                       .  <> seed5
        8               .
                         .
        9

        10
            0   1   2   3   4   5   6   7   8   9   10  11 */

        /*
        float[] seed0 = new float[]{2.0f, 1.5f};
        float[] line0 = new float[]{3.0f, 2.0f, 4.2f, 0.f};

        float[] seed1 = new float[]{4.5f, 1.5f};


        float[] seed2 = new float[]{2.0f, 3.5f};
        float[] line2 = new float[]{0.0f, 4.0f, 4.5f, 4.f};

        float[] seed3 = new float[]{2.0f, 4.5f};

        float[] seed4 = new float[]{2.0f, 7.0f};
        float[] line4 = new float[]{2.0f, 6.0f, 3.5f, 9.f};
        float[] seed5 = new float[]{3.5f, 7.5f};


        PointInPolygon instance = new PointInPolygon();

        boolean result = instance.isInSimpleCurve(testPointX, testPointY, squarePolyX, squarePolyY, squarePolyX.length);
        */
    }

}
