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

    public void testSquarePoly() throws Exception {

        /*
        10********
          *      *
          *      *
        0 ********
          0     10
        */
        float[] xPoly = new float[5];
        float[] yPoly = new float[5];
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
        
        float[] testPointsX = new float[9];
        float[] testPointsY = new float[9];
         testPointsX[0] = 5;
         testPointsY[0] = 5;
         testPointsX[1] = 5;
         testPointsY[1] = 8;

         testPointsX[2] = -10;
         testPointsY[2] = 5;
         testPointsX[3] = 0;
         testPointsY[3] = 5;

         testPointsX[4] = 10;
         testPointsY[4] = 5;
         testPointsX[5] = 8;
         testPointsY[5] = 5;

         testPointsX[6] = 10;
         testPointsY[6] = 10;
         
         testPointsX[7] = 15;
         testPointsY[7] = 5;
         testPointsX[8] = 5;
         testPointsY[8] = 15;
         
        // test square poly
        boolean[] expected = new boolean[] {
            true, true, false, true, true, true, true, false, false
        };
        
        /*
        10********
          *      *
          *      *
        0 ********
          0     10
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (8, 5)
         (10, 10)
         (15, 5)
         (5, 15)
        */
        PointInPolygon instance = new PointInPolygon();
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, 
                xPoly, yPoly, xPoly.length);

            //System.out.println("  result=" + result + " expect=" + expected[i]);
            assertTrue(result == expected[i]);
        }
    }

    /*
    note, the polygon is defined as a sequence of points and the first and
    last points are the same.  the definition is no-longer using redundant points
    such as pt0 to pt1 then pt1 to pt2 then pt2 to pt3.
    */
    public void testExagonPoly() throws Exception {

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
        float[] xPoly = new float[7];
        float[] yPoly = new float[7];
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
        
        float[] testPointsX = new float[9];
        float[] testPointsY = new float[9];
         testPointsX[0] = 5;
         testPointsY[0] = 5;
         testPointsX[1] = 5;
         testPointsY[1] = 8;

         testPointsX[2] = -10;
         testPointsY[2] = 5;
         testPointsX[3] = 0;
         testPointsY[3] = 5;

         testPointsX[4] = 10;
         testPointsY[4] = 5;
         testPointsX[5] = 5;
         testPointsY[5] = -5;

         testPointsX[6] = 10;
         testPointsY[6] = 10;
         
         testPointsX[7] = 15;
         testPointsY[7] = 5;
         testPointsX[8] = 5;
         testPointsY[8] = 15;
         
        boolean[] expected = new boolean[] {
            true, true, false, true, true, false, false, false, false
        };
        
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
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (5, -5)
         (10, 10)
         (15, 5)
         (5, 15)
        */
        PointInPolygon instance = new PointInPolygon();
        for (int i = 0; i < testPointsX.length; i++) {
            float testPointX = testPointsX[i];
            float testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, 
                xPoly, yPoly, xPoly.length);

            //System.out.println("[" + i + "]   result=" + result + " expect=" + expected[i]);
            
            assertTrue(result == expected[i]);
        }
    }
  
    public void testSquarePoly_int() throws Exception {

        /*
        10********
          *      *
          *      *
        0 ********
          0     10
        */
        int[] xPoly = new int[5];
        int[] yPoly = new int[5];
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
        
        int[] testPointsX = new int[9];
        int[] testPointsY = new int[9];
        testPointsX[0] = 5;
        testPointsY[0] = 5;
        testPointsX[1] = 5;
        testPointsY[1] = 8;

        testPointsX[2] = -10;
        testPointsY[2] = 5;
        testPointsX[3] = 0;
        testPointsY[3] = 5;

        testPointsX[4] = 10;
        testPointsY[4] = 5;
        testPointsX[5] = 8;
        testPointsY[5] = 5;

        testPointsX[6] = 10;
        testPointsY[6] = 10;

        testPointsX[7] = 15;
        testPointsY[7] = 5;
        testPointsX[8] = 5;
        testPointsY[8] = 15;
         
        // test square poly
        boolean[] expected = new boolean[] {
            true, true, false, true, true, true, true, false, false
        };
        
        /*
        10********
          *      *
          *      *
        0 ********
          0     10
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (8, 5)
         (10, 10)
         (15, 5)
         (5, 15)
        */
        
        PointInPolygon instance = new PointInPolygon();
        for (int i = 0; i < testPointsX.length; i++) {
            int testPointX = testPointsX[i];
            int testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, 
                xPoly, yPoly, xPoly.length);

            //System.out.println("[" + i + "]   result=" + result + " expect=" + expected[i]);
            
            assertEquals(expected[i], result);
        }
    }

    /*
    note, the polygon is defined as a sequence of points and the first and
    last points are the same.  the definition is no-longer using redundant points
    such as pt0 to pt1 then pt1 to pt2 then pt2 to pt3.
    */
    public void testExagonPoly_int() throws Exception {

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
        int[] xPoly = new int[7];
        int[] yPoly = new int[7];
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
        
        int[] testPointsX = new int[9];
        int[] testPointsY = new int[9];
         testPointsX[0] = 5;
         testPointsY[0] = 5;
         testPointsX[1] = 5;
         testPointsY[1] = 8;

         testPointsX[2] = -10;
         testPointsY[2] = 5;
         testPointsX[3] = 0;
         testPointsY[3] = 5;

         testPointsX[4] = 10;
         testPointsY[4] = 5;
         testPointsX[5] = 5;
         testPointsY[5] = -5;

         testPointsX[6] = 10;
         testPointsY[6] = 10;
         
         testPointsX[7] = 15;
         testPointsY[7] = 5;
         testPointsX[8] = 5;
         testPointsY[8] = 15;
         
        // test square poly
        boolean[] expected = new boolean[] {
            true, true, false, true, true, false, false, false, false
        };
        
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
         (5, 5)
         (5, 8)
         (-10, 5)
         (0, 5)
         (10, 5)
         (5, -5)
         (10, 10)
         (15, 5)
         (5, 15)
        */
        PointInPolygon instance = new PointInPolygon();
        //for (int i = 0; i < testPointsX.length; i++) {
        for (int i = 6; i < 7; i++) {
            int testPointX = testPointsX[i];
            int testPointY = testPointsY[i];

            boolean result = instance.isInSimpleCurve(testPointX, testPointY, 
                xPoly, yPoly, xPoly.length);

            //System.out.println("[" + i + "]   result=" + result + " expect=" + expected[i]);
            
            assertTrue(result == expected[i]);
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
