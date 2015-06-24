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

        /*
        10********
          *      *
          *      *
        0 ********
          0     10
        */
        
        xPoly[0] = 0;
        yPoly[0] = 0;
        xPoly[1] = 10;
        yPoly[1] = 0;
        xPoly[2] = 10;
        yPoly[2] = 10;
        xPoly[3] = 0;
        yPoly[3] = 10;
        xPoly[4] = 0;
        yPoly[4] = 0;
        
    }

    /*
    note, the polygon is defined as a sequence of points and the first and
    last points are the same.  the definition is no-longer using redundant points
    such as pt0 to pt1 then pt1 to pt2 then pt2 to pt3.
    */
    public void getExagonPoly(float[] xPoly, float[] yPoly) {

        /*

       10

                *7,8--*5,6          10
              /          \
            /             \
        *9,10               *3,4    5
           \               /
             \           /
       0        *0,11-*1,2   
          0     3      7     10     0
        
        */
        xPoly[0] = 3;
        yPoly[0] = 0;
        xPoly[1] = 7;
        yPoly[1] = 0;
        xPoly[2] = 10;
        yPoly[2] = 5;
        xPoly[3] = 7;
        yPoly[3] = 10;
        xPoly[4] = 3;
        yPoly[4] = 10;
        xPoly[5] = 0;
        yPoly[5] = 5;
        xPoly[6] = 3;
        yPoly[6] = 0;
    }

     public void getTestPts(float[] xPoly, float[] yPoly) {

         /*
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (8, 5)
         (10, 10)
         */
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

        float[] squarePolyX = new float[5];
        float[] squarePolyY = new float[5];
        getSquarePoly(squarePolyX, squarePolyY);


        float[] exagonPolyX = new float[7];
        float[] exagonPolyY = new float[7];
        getExagonPoly(exagonPolyX, exagonPolyY);


        float[] testPointsX = new float[7];
        float[] testPointsY = new float[7];
        getTestPts(testPointsX, testPointsY);


        PointInPolygon instance = new PointInPolygon();

        // test square poly
        boolean[] expected = new boolean[] {
            true, true, false, true, true, true, true
        };
        /*
        10******** 2 2
          *      *
          *      *
        0 ********
          0     10
          0      1 1
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (8, 5)
         (10, 10)
        */
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, 
                squarePolyX, squarePolyY, squarePolyX.length);

            //System.out.println("  result=" + result + " expect=" + expected[i]);
            assertTrue(result == expected[i]);
        }
         /*

           10

                    *7,8--*5,6         10
                  /          \
                /             \
            *9,10               *3,4   5
               \               /
                 \           /
           0        *0,11-*1,2   
              0    3       7     10
        
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (8, 5)
         (10, 10)
        */
        // test exagon poly
        expected = new boolean[] {
            true, true, false, true, true, true, false
        };
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, 
                exagonPolyX, exagonPolyY, exagonPolyX.length);
            boolean expect = expected[i];
            System.out.println("  result=" + result + " expect=" + expect);
            assertTrue(result == expect);//i=5  test=(8,5)
        }
    }
  
    public void test10() throws Exception {

        float[] xPoly = new float[]{0.325f, 0.52f, 0.6f, 0.6f,  0.325f};
        float[] yPoly = new float[]{0.28f,  0.52f, 0.52f, 0.28f, 0.28f};
        
        float xPt = 0.3712f;
        float yPt = 0.4467f;
        
        PointInPolygon instance = new PointInPolygon();
        
        assertFalse(instance.isInSimpleCurve(xPt, yPt, xPoly, yPoly, xPoly.length));
       
        assertTrue(instance.isInSimpleCurve(0.45f, 0.40f, xPoly, yPoly, xPoly.length));
        
        assertFalse(instance.isInSimpleCurve(0.4f, 0.2f, xPoly, yPoly, xPoly.length));
        
        assertFalse(instance.isInSimpleCurve(0.65f, 0.2f, xPoly, yPoly, xPoly.length));
    }

}
