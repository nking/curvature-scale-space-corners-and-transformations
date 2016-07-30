package algorithms.compGeometry;

import java.util.logging.Logger;
import junit.framework.TestCase;

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
        double x1 = -3;
        double y1 = 2;
        double x2 = 1;
        double y2 = 3;
        double result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result < 0); // clockwise
        assertTrue(result == -11);

        x1 = 3;
        result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result >= 0); // counter clockwise
        assertTrue(result == 7);
    }
    
    public void testDirection() throws Exception {
        
        /*
         * 
         * 10                             p2
         *  9
         *  8                          /
         *  7     p3
         *  6      \
         *  5                     /
         *  4        \
         *  3
         *  2           \     /
         *  1              p1
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * direction of change from P1:P2 to P1:P3
         */
        
        float x2 = 10;
        float y2 = 10;
        float x1 = 5;
        float y1 = 1;
        float x3 = 2;
        float y3 = 7;
        
        double direction = LinesAndAngles.direction(x1, y1, x2, y2, x3, y3);
        
        assertTrue(direction > 0);
        
        // swap p1 and p3
        x3 = 10;
        y3 = 10;
        x2 = 2;
        y2 = 7;
        direction = LinesAndAngles.direction(x1, y1, x2, y2, x3, y3);
        
        assertTrue(direction < 0);

    }
    
    public void testDirection_int() throws Exception {
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                      
         *  4     p1     
         *  3
         *  2                 
         *  1              p2
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * direction of change from P1:P2 to P1:P3
         */
        
        // CCW sweep of P1:P2 to P1:P3 gives negative number
        int direction = LinesAndAngles.direction(
            2, 4, 
            5, 1, 
            3, 7);
        assertTrue(direction < 0);
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                      
         *  4     p1     
         *  3
         *  2                 
         *  1     p2       
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * direction of change from P1:P2 to P1:P3
         */
        
        // CCW sweep of P1:P2 to P1:P3 gives negative number
        direction = LinesAndAngles.direction(
            2, 4, 
            2, 1, 
            3, 7);
        assertTrue(direction < 0);
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                      
         *  4     p1     
         *  3
         *  2  p2             
         *  1         
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * direction of change from P1:P2 to P1:P3
         */
        
        // CCW sweep of P1:P2 to P1:P3 gives negative number
        direction = LinesAndAngles.direction(
            2, 4, 
            1, 2, 
            3, 7);
        assertTrue(direction > 0);
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                    p1  
         *  4            
         *  3
         *  2                 
         *  1              p2
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * direction of change from P1:P2 to P1:P3
         */
        // CCW sweep of P1:P2 to P1:P3 gives positive number
        direction = LinesAndAngles.direction(
            7, 5, 
            5, 1, 
            3, 7);
        assertTrue(direction > 0);
       
    }
    
    public void testPointIsInLine() throws Exception {
        
        boolean isInLine = LinesAndAngles.pointIsInLine(0, 5, 0, 0, 0, 10);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(1, 5, 0, 0, 0, 10);
        assertFalse(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(0.f, 5.f, 0.f, 0.f, 0.f, 10.f);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(2, 2, 0, 0, 4, 4);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(3, 0, 2, 0, 4, 0);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(3, 1, 2, 0, 4, 0);
        assertFalse(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(2, 0, 2, 0, 4, 0);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(4, 4, 0, 0, 4, 4);
        assertTrue(isInLine);
    }
    
    public void testCalcClockwiseAngle() {
        
        int x1, y1, x2, y2, x3, y3;
        double angle, expected; 
        double eps = 0.001;
        
        /*
           1     2
              3
        */
        x1 = 0;
        y1 = 10;
        x3 = 10;
        y3 = 0;
        x2 = 20;
        y2 = 10;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = Math.PI/2;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3   2
        */
        x1 = 0;
        y1 = 10;
        x3 = 10;
        y3 = 0;
        x2 = 20;
        y2 = 0;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = 3*Math.PI/4;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3 
                 2
        */
        x1 = 0;
        y1 = 20;
        x3 = 10;
        y3 = 10;
        x2 = 20;
        y2 = 0;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = Math.PI;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3 
              2
        */
        x1 = 0;
        y1 = 20;
        x3 = 10;
        y3 = 10;
        x2 = 10;
        y2 = 0;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = 5.*Math.PI/4.;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3 
           2
        */
        x1 = 0;
        y1 = 20;
        x3 = 10;
        y3 = 10;
        x2 = 0;
        y2 = 0;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = 6.*Math.PI/4.;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
              1     
           2  3 
        */
        x1 = 10;
        y1 = 10;
        x3 = 10;
        y3 = 0;
        x2 = 0;
        y2 = 0;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = 6.*Math.PI/4.;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
              3  2
              1
        */
        x1 = 0;
        y1 = 0;
        x3 = 0;
        y3 = 10;
        x2 = 10;
        y2 = 10;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        expected = 6.*Math.PI/4.;
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
    }
}
