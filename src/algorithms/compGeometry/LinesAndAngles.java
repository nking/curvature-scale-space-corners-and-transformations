package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import java.util.Arrays;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/LinesAndAngles.java
 * under MIT License (MIT), Nichole King 2013
 * 
 */
public class LinesAndAngles {
    
    public static double distSquared(double x1, double y1, double x2, double y2) {

        double dx2 = (x2 - x1);
        dx2 *= dx2;
        double dy2 = (y2 - y1);
        dy2 *= dy2;
        return dx2 + dy2;
    }
    
    /**
     * calculate the direction of change for the 2 vectors 
     * P1:P2 to P1:P3  returns negative when direction is clockwise,
     * else if zero the vectors are collinear, else if positive the
     * directions is counterclockwise.
     * 
     *           p2
     *  p3      / 
     *   \    /
     *     p1
     * 
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return 
     */
    public static double direction(float x1, float y1, float x2, float y2, 
        float x3, float y3) {

        double d = ((x2 - x1)*(y3 - y1)) - ((y2 - y1)*(x3 - x1));
        
        return d;
    }
    
    /**
     <pre>
      determine the cross product and return negative number if the 2nd
      set of values, p1, are clockwise from the first, p1, w.r.t. origin (0,0).
      
          * P2
          .
          .
          .   * P1      &lt;--- P2 is counterclockwise from P1 w.r.t. origin o
          .
          o
      
                * P2
                .
        P1 *    .       &lt;--- P2 is clockwise from P1 w.r.t. origin o
                .
                .
                o
     </pre>
     
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static double crossProduct(double x1, double y1, double x2, 
        double y2) {
        
        return ((x1*y2) - (x2*y1));
    }
    
}