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
    
    public static double lengthOfLine(double x1, double y1, double x2, double y2) {
        double d = Math.sqrt( distSquared( x1,  y1,  x2,  y2) );
        return d;
    }
    
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
    
    /**
     * create an interpolated polygon of area around the peak for the top full width topFraction of max
     * immediately surrounding the peak.
     * For example, this is used to find the centroid of the peak of a Generalized Extreme Value distribution.
     *
     * @param x
     * @param y
     * @param xe errors for array x, can be null
     * @param ye errors for array y, can be null
     * @param topFraction
     * @return
     */
    public static ArrayPair createPolygonOfTopFWFractionMax(float[] x, 
        float[] y, float[] xe, float[] ye, float topFraction) {

        float[] xTopPolygon = new float[x.length + 3];
        float[] yTopPolygon = new float[y.length + 3];
        
        float[] xTopPolygonError = (xe == null) ? null : new float[x.length + 3];
        float[] yTopPolygonError = (ye == null) ? null : new float[x.length + 3];
        
        int yPeakIndex = MiscMath.findYMaxIndex(y);

        float yFractionLimit = y[yPeakIndex]*topFraction;

        int polygonCount = 1; // leave space for the interpolated 1st point

        // these are w.r.t. x and y
        int firstIndex = -1;
        int lastIndex = -1;

        for (int i = 0; i < x.length; i++) {

            if (i <= yPeakIndex) {

                if (y[i] >= yFractionLimit) {

                    if (firstIndex == -1) {
                        firstIndex = i;
                    }

                    xTopPolygon[polygonCount] = x[i];
                    yTopPolygon[polygonCount] = y[i];
                    if (xTopPolygonError != null) {
                        xTopPolygonError[polygonCount] = xe[i];
                        yTopPolygonError[polygonCount] = ye[i];
                    }
                    polygonCount++;
                }
            } else if (y[i] >= yFractionLimit) {

                xTopPolygon[polygonCount] = x[i];
                yTopPolygon[polygonCount] = y[i];
                if (xTopPolygonError != null) {
                    xTopPolygonError[polygonCount] = xe[i];
                    yTopPolygonError[polygonCount] = ye[i];
                }
                polygonCount++;

            } else {
                // we are below the yFractionLimit
                lastIndex = i - 1;
                break;
            }
        }
        if (lastIndex == -1) {
            lastIndex = x.length - 1;
        }

        // compress the array.  we will add two items beyond polygonCount.
        //    one for the interpolated yFractionLimit, and the other for a
        //    copy of the first point to close the polygon
        if ((polygonCount + 1) < xTopPolygon.length) {
            xTopPolygon = Arrays.copyOf(xTopPolygon, (polygonCount + 2));
            yTopPolygon = Arrays.copyOf(yTopPolygon, (polygonCount + 2));
            
            if (xTopPolygonError != null) {
                xTopPolygonError = Arrays.copyOf(xTopPolygonError, (polygonCount + 2));
                yTopPolygonError = Arrays.copyOf(yTopPolygonError, (polygonCount + 2));
            }
        }

        /* <pre>
         * (0)                     | (1)                     | (2)                   | (3)
         *           peak          |           peak          |         peak          |            peak
         *                         |                         |                       |
         *                         |     pt_i+1    pt_j-1    |                       |
         *     pt_i        pt_j    |   pt_i           pt_j   |                       |
         *    ------------------   |  ---------------------  |  -------------------  |  ----------------------
         *                         |                         |                       |
         *   pt_0             pt_k | pt_0              pt_k  |   pt_0      pt_k      |  &lt;no point&gt;      &lt;no point&gt;
         *                         |                         |                       |
         *
         *                                                                           | (4)     peak
         *                                                                           |
         *                                                                           |      pt_i       pt_j
         *                                                                           |  ----------------------
         *                                                                           |
         *                                                                           | &lt;no point&gt;      &lt;no point&gt;
         * </pre>
         */

        // interpolate first point
        yTopPolygon[0] = yFractionLimit;
        if (xTopPolygonError != null) {
            // TODO: correct this
            yTopPolygonError[0] = ye[0];
        }
        
        if (firstIndex > 0) {
            /* case (0), case (1), and case (2)
             *  y[pt_i] - y[pt_0]    y[pt_i] - y[limit]                               y[pt_i] - y[limit]
             *  ----------------- =  ------------------  ====>  x[pt_i] - x[limit] =  ------------------ * x[pt_i] - x[pt_0]
             *  x[pt_i] - x[pt_0]    x[pt_i] - x[limit]                               y[pt_i] - y[pt_0]
             */
            float tmp = (y[firstIndex] - yFractionLimit)/(y[firstIndex] - y[firstIndex - 1]);
            xTopPolygon[0] = x[firstIndex] - (tmp * (x[firstIndex] - x[firstIndex - 1]));
            
            if (xTopPolygonError != null) {
                xTopPolygonError[0] = (float)Math.sqrt(Math.pow(xe[firstIndex], 2) + 
                    Math.pow(xe[firstIndex - 1], 2));
            }
            
        } else if ((firstIndex == 0) && (yPeakIndex == 0)) {
            /* case (3)
             * there is no point before the peak, so we create a point right below as a straight line
             */
             xTopPolygon[0] = x[yPeakIndex];

             if (xTopPolygonError != null) {
                 xTopPolygonError[0] = xe[yPeakIndex];
             }
             
        } else if ((firstIndex == 0) && (yPeakIndex > 0)) {
            /* case (4)
             * There is a point before the peak, so we extend line.
             *   y[firstIndex+1] - y[firstIndex]     y[firstIndex] - y[limit]
             *   -------------------------------  = -------------------------
             *   x[firstIndex+1] - x[firstIndex]     x[firstIndex] - x[limit]
             *
             *                                 y[firstIndex] - y[limit]
             * x[firstIndex] - x[limit] =  ------------------------------- * x[firstIndex+1] - x[firstIndex]
             *                             y[firstIndex+1] - y[firstIndex]
             *
             *                                  y[firstIndex] - y[limit]
             * x[limit] = x[firstIndex] -   ------------------------------ * x[firstIndex+1] - x[firstIndex]
             *                              y[firstIndex+1] - y[firstIndex]
             */
            float tmp = (y[firstIndex] - yFractionLimit)/(y[firstIndex + 1] - y[firstIndex]);
            xTopPolygon[0] = x[firstIndex] - (tmp * (x[firstIndex + 1] - x[firstIndex]));

            if (xTopPolygonError != null) {
                xTopPolygonError[0] = (float)Math.sqrt(Math.pow(xe[firstIndex], 2) + 
                    Math.pow(xe[firstIndex + 1], 2));
            }
            
        } else {
            throw new IllegalStateException("have not solved for this first point! ");
        }

        // interpolate the point after lastIndex

        yTopPolygon[polygonCount] = yFractionLimit;
        if (xTopPolygonError != null) {
            // TODO: correct this
            yTopPolygonError[polygonCount] = ye[0];
        }
        
        if (lastIndex < (x.length - 1)) {
            /* case (0), case (1), and case (2)
             *  y[pt_j] - y[pt_k]    y[pt_j] - y[limit]                               y[pt_j] - y[limit]
             *  ----------------- =  ------------------  ====>  x[pt_j] - x[limit] =  ------------------ * x[pt_j] - x[pt_k]
             *  x[pt_j] - x[pt_k]    x[pt_j] - x[limit]                               y[pt_j] - y[pt_k]
             */
            float slope = (y[lastIndex] - yFractionLimit)/(y[lastIndex] - y[lastIndex + 1]);
            xTopPolygon[polygonCount] = x[lastIndex] - (slope * (x[lastIndex] - x[lastIndex + 1]));

            if (xTopPolygonError != null) {
                xTopPolygonError[0] = (float)Math.sqrt(Math.pow(xe[lastIndex], 2) + 
                    Math.pow(xe[lastIndex + 1], 2));
            }
            
        } else if ((lastIndex == (x.length - 1)) && (lastIndex == yPeakIndex)) {
            /* case (3)
             * there is no point after the peak, so we create a point right below as a straight line
             */
             xTopPolygon[polygonCount] = x[yPeakIndex];

             if (xTopPolygonError != null) {
                 xTopPolygonError[polygonCount] = xe[yPeakIndex];
             }
            
        } else if ((lastIndex == (x.length - 1)) && (yPeakIndex < lastIndex)) {
            /* case (4)
             * There is a point after the peak, so we extend line.
             *   y[lastIndex-1] - y[lastIndex]     y[lastIndex] - y[limit]
             *   -------------------------------  = -------------------------
             *   x[lastIndex-1] - x[lastIndex]     x[lastIndex] - x[limit]
             *
             *                                y[lastIndex] - y[limit]
             * x[lastIndex] - x[limit] =  ------------------------------- * x[lastIndex-1] - x[lastIndex]
             *                             y[lastIndex-1] - y[lastIndex]
             *
             *                                  y[lastIndex] - y[limit]
             * x[limit] = x[lastIndex] -   ------------------------------ * x[lastIndex-1] - x[lastIndex]
             *                              y[lastIndex-1] - y[lastIndex]
             */
            float tmp = (y[lastIndex] - yFractionLimit)/(y[lastIndex - 1] - y[lastIndex]);
            xTopPolygon[polygonCount] = x[lastIndex] - (tmp * (x[lastIndex - 1] - x[lastIndex]));

            if (xTopPolygonError != null) {
                xTopPolygonError[polygonCount] = (float)Math.sqrt(Math.pow(xe[lastIndex], 2) + 
                    Math.pow(xe[lastIndex - 1], 2));
            }
            
        } else {
            throw new IllegalStateException("have not solved for this later point! ");
        }

        xTopPolygon[xTopPolygon.length - 1] = xTopPolygon[0];
        yTopPolygon[xTopPolygon.length - 1] = yTopPolygon[0];

        if (xTopPolygonError != null) {
            xTopPolygonError[xTopPolygon.length - 1] = xTopPolygonError[0];
            yTopPolygonError[xTopPolygon.length - 1] = yTopPolygonError[0];
        }
        
        ArrayPair xy = null;
        if (xTopPolygonError == null) {
            xy = new ArrayPair(xTopPolygon, yTopPolygon);
        } else {
            xy = new ArrayPair(xTopPolygon, yTopPolygon, xTopPolygonError, yTopPolygonError);
        }

        return xy;
    }
    
    public static float[] calcAreaAndCentroidOfSimplePolygon(float[] xPolygon, 
        float[] yPolygon) {

        if (xPolygon.length < 3) {
            return null;
        }

        float area = 0;

        boolean allPointsArePositive = true;

        // from wikipedia
        for (int i = 0; i < (xPolygon.length - 1); i++) {
            float a = xPolygon[i] * yPolygon[i+1];
            float b = xPolygon[i+1] * yPolygon[i];
            area += (a-b);

            if ((xPolygon[i] < 0) || (yPolygon[i] < 0)) {
                allPointsArePositive = false;
            }
        }
        if ((xPolygon[xPolygon.length - 1] < 0) || (yPolygon[xPolygon.length - 1] < 0)) {
            allPointsArePositive = false;
        }
        if ((area < 0) && allPointsArePositive) {
            area *= -1;
        }
        area *= 0.5f;


        float xc = 0;
        for (int i = 0; i < (xPolygon.length - 1); i++) {
            float a = xPolygon[i] * yPolygon[i+1];
            float b = xPolygon[i+1] * yPolygon[i];

            xc += ( (xPolygon[i] + xPolygon[i+1]) * (a - b) );
        }
        xc *= (1./(6.*area));
        if ((xc < 0) && allPointsArePositive) {
            xc *= -1;
        }


        float yc = 0;
        for (int i = 0; i < (xPolygon.length - 1); i++) {
            float a = xPolygon[i] * yPolygon[i+1];
            float b = xPolygon[i+1] * yPolygon[i];

            yc += ( (yPolygon[i] + yPolygon[i+1]) * (a - b) );
        }
        //yc = yPolygon[0] - yc; // make it relative to zero point
        yc *= (1./(6.*area));
        if ((yc < 0) && allPointsArePositive) {
            yc *= -1;
        }

        return new float[]{area, xc, yc};
    }

}