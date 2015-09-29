package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class LinesAndAngles {

    /**
     <pre>
      determine the cross product and return negative number if the 2nd
      set of values, p1, are clockwise from the first, p1, w.r.t. origin (0,0).
      
          * P2
          .
          .
          .   * P1     &lt;--- P2 is counterclockwise from P1 w.r.t. origin o
          .
          o
      
                * P2
                .
        P1 *    .      &lt;--- P2 is clockwise from P1 w.r.t. origin o
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
    public static int crossProduct(int x1, int y1, int x2, int y2) {
        return ((x1*y2) - (x2*y1));
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static double crossProduct(double x1, double y1, double x2, double y2) {
        return ((x1*y2) - (x2*y1));
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @param x4
     * @param y4
     * @return
     */
    public static boolean linesIntersect(int x1, int y1,
        int x2, int y2, int x3, int y3, int x4, int y4) {

        int d1 = direction(x3, y3, x4, y4, x1, y1);
        int d2 = direction(x3, y3, x4, y4, x2, y2);
        int d3 = direction(x1, y1, x2, y2, x3, y3);
        int d4 = direction(x1, y1, x2, y2, x4, y4);

        if (
            (((d1 > 0) && (d2 < 0)) || ((d1 < 0) && (d2 > 0)))
            &&
            (((d3 > 0) && (d4 < 0)) || ((d3 < 0) && (d4 > 0)))
        ) {
            return true;
        } else if ( ((d1 == 0) && onSegment(x3, y3, x4, y4, x1, y1)) ) {
            return true;
        } else if ( ((d2 == 0) && onSegment(x3, y3, x4, y4, x2, y2)) ) {
            return true;
        } else if ( ((d3 == 0) && onSegment(x1, y1, x2, y2, x3, y3)) ) {
            return true;
        } else if ( ((d4 == 0) && onSegment(x1, y1, x2, y2, x4, y4)) ) {
            return true;
        }
        return false;
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return
     */
    public static int direction(int x1, int y1, int x2, int y2, int x3, int y3) {

        int x31 = x3 - x1;
        int y31 = y3 - y1;

        int x21 = x2 - x1;
        int y21 = y2 - y1;

        return crossProduct(x31, y31, x21, y21);
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return
     */
    public static double direction(float x1, float y1, float x2, float y2, float x3, float y3) {

        double x31 = x3 - x1;
        double y31 = y3 - y1;

        double x21 = x2 - x1;
        double y21 = y2 - y1;

        return crossProduct(x31, y31, x21, y21);
    }

    static boolean onSegment(int x1, int y1,
        int x2, int y2, int x3, int y3) {

        int minx12 = (x1 < x2) ? x1 : x2;
        int miny12 = (y1 < y2) ? y1 : y2;

        int maxx12 = (x1 > x2) ? x1 : x2;
        int maxy12 = (y1 > y2) ? y1 : y2;

        if ( (minx12 <= x3) && (x3 <= maxx12) && (miny12 <= y3) && (y3 <= maxy12) ) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * calculate the angle from the "horizon" at point p1 to the line p1 to point p2 
     * 
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static double calculatePolarSineTheta(double x1, double y1, double x2, double y2) {
        // log.info("(" + x1 + "," + y1 + ")," + "(" + x2 + "," + y2 + ")");
        // determine quadrant of points

        if (x2 == x1) {
            if (y2 == y1) {
                return 0.0f;
            }
    		return 1.0;//0.0f??
    	} else if (x2 > x1) {
             //
             //       p2 o
             //         /
             //      c /       sin(t) = b/c
             //       /   b
             //      /
             //  p1 o....
             //
             //   p1 o....
             //       \   c
             //      b \       sin(t) = b/c
             //         \
             //       p2 o
             //
            double c = lengthOfLine(x1, y1, x2, y2);
            double b = lengthOfLine(x2, y1, x2, y2);
            return (b/c);
        } else {
             //                                                       .
             //   p2 o              p2 o
             //       \                |\c        sin(t) = b/c + sin(pi)
             //        \              b| \
             //      p1 o------        ...o------
             //                           p1
             //            p1
             //         .. o-----
             //        b| /          sin(t) = b/c + sin(pi)
             //         |/ c
             //      p2 o
            double c = lengthOfLine(x1, y1, x2, y2);
            double b = lengthOfLine(x2, y1, x1, y1);

            // the angle is actually (1-b/c) + sin90 = (1-b/c) + 1 = (b/c)

            return (b/c) + 1.0;
        }
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static double calculatePolarSineTheta(float x1, float y1, float x2, float y2) {
        // log.info("(" + x1 + "," + y1 + ")," + "(" + x2 + "," + y2 + ")");
        // determine quadrant of points
        if (x2 == x1) {
            if (y2 == y1) {
                return 0.0f;
            }
                return 1.0;//0.0f??
        } else if (x2 > x1) {
             //
             //       p2 o
             //         /
             //      c /       sin(t) = b/c
             //       /   b
             //      /
             //  p1 o....
             //
             //   p1 o....
             //       \   c
             //      b \       sin(t) = b/c
             //         \
             //       p2 o
             //
            double c = lengthOfLine(x1, y1, x2, y2);
            double b = lengthOfLine(x2, y1, x2, y2);
            return (b/c);
        } else {
            //                                                       .
            //   p2 o              p2 o
            //       \                |\c        sin(t) = b/c + sin(pi)
            //        \              b| \
            //      p1 o------        ...o------
            //                           p1
            //            p1
            //         .. o-----
            //        b| /          sin(t) = b/c + sin(pi)
            //         |/ c
            //      p2 o
           double c = lengthOfLine(x1, y1, x2, y2);
           double b = lengthOfLine(x2, y1, x1, y1);

           // the angle is actually (1-b/c) + sin90 = (1-b/c) + 1 = (b/c)

           return (b/c) + 1.0;
       }
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static double lengthOfLine(double x1, double y1, double x2, double y2) {
        double d = Math.sqrt( distSquared( x1,  y1,  x2,  y2) );
        return d;
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static double distSquared(double x1, double y1, double x2, double y2) {

        double dx2 = (x2 - x1);
        dx2 *= dx2;
        double dy2 = (y2 - y1);
        dy2 *= dy2;
        return dx2 + dy2;
    }

    /** <pre>from http://en.wikipedia.org/wiki/Line-line_intersection
          
              (x1, y1)
                   *
                    \    (x4, y4)
                     \      *
                      \    /
                       \  /
                        \/
                      P /\
                       /  \
                      /    \
                     /      \
                    *        \
              (x3, y3)        \
                               \
                                *
                            (x2, y2)
           </pre>
         * 
         * @param x1
         * @param y1
         * @param x2
         * @param y2
         * @param x3
         * @param y3
         * @param x4
         * @param y4
         * @return the intersection of the 2 lines as a 2 item double array holding xTopPolygon then yTopPolygon, else
         * if lines do not intersect, will return null;
         */
    public static double[] intersectionOf2Lines(float x1, float y1, float x2, float y2,
        float x3, float y3, float x4, float y4) {

        /*
         *        | x1  y1 |   | x1  1 |
         *        | x2  y2 |   | x2  1 |
         *        |                    |
         *        | x3  y3 |   | x3  1 |
         *        | x4  y4 |   | x4  1 |
         * P_x =  ----------------------------
         *        | x1   1 |   | y1  1 |
         *        | x2   1 |   | y2  1 |
         *        |                    |
         *        | x3   1 |   | y3  1 |
         *        | x4   1 |   | y4  1 |
         *
         *
         *        | x1  y1 |   | y1  1 |
         *        | x2  y2 |   | y2  1 |
         *        |                    |
         *        | x3  y3 |   | y3  1 |
         *        | x4  y4 |   | y4  1 |
         * P_y =  ----------------------------
         *        | x1   1 |   | y1  1 |
         *        | x2   1 |   | y2  1 |
         *        |                    |
         *        | x3   1 |   | y3  1 |
         *        | x4   1 |   | y4  1 |
         *
         */

        double x1y2minusy1x2 = (x1*y2) - (y1*x2);

        double x3y4minusy3x4 = (x3*y4) - (y3*x4);

        double x3minusx4 = x3 - x4;

        double y3minusy4 = y3 - y4;

        double y1minusy2 = y1 - y2;

        double x1minusx2 = x1 - x2;

        double p_denom = (x1minusx2*y3minusy4) - (y1minusy2*x3minusx4);

        if (p_denom == 0.0) {
            return null;
        }

        double p_x = ((x1y2minusy1x2 * x3minusx4) - (x3y4minusy3x4 * x1minusx2))/p_denom;

        double p_y = ((x1y2minusy1x2 * y3minusy4) - (x3y4minusy3x4 * y1minusy2))/p_denom;

        return new double[]{p_x, p_y};
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @param x4
     * @param y4
     * @return
     */
    public static double[] intersectionOf2Lines_2(float x1, float y1, float x2, float y2,
        float x3, float y3, float x4, float y4) {

        double tc1=(x2-x1);

        double tc2=(y2-y1);

        double sc1=(x4-x3);

        double sc2=(y4-y3);

        double con1=(x3-x1);

        double con2=(y3-y1);

        double p_denom = (tc2*sc1-tc1*sc2);

        if (p_denom == 0.0) {
            return null;
        }

        double con=tc2*con1-tc1*con2;
        double s = con/p_denom;

        double ax = (s*sc1);
        double ay = (s*sc2);
        double p_x = x3 - ax;
        double p_y = y3 - ay;

        return new double[]{p_x, p_y};
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static float calculateSlope(float x1, float y1, float x2, float y2) {
        return ( y2 - y1 ) / (x2 - x1);
    }

    /**
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static float calculateTangentSlope(float x1, float y1, float x2, float y2) {
        return -1.f*( x2 - x1 ) / (y2 - y1);
    }

    /**
     *
     * @param x0
     * @param y0
     * @param y
     * @param slope
     * @return
     */
    public static float calculateX(float x0, float y0, float y, float slope) {
        return ( x0 + ((y - y0)/slope) );
    }

    //  (yTopPolygon-y0)/(xTopPolygon-x0) = slope

    /**
     *
     * @param x0
     * @param y0
     * @param x
     * @param slope
     * @return
     */
        public static float calculateY(float x0, float y0, float x, float slope) {
        return ( y0 + slope*(x - x0) );
    }

    /**
     *
     * @param xi
     * @param yi
     * @param xj
     * @param yj
     * @param xMin
     * @param xMax
     * @param yMin
     * @param yMax
     * @return
     */
    public static float[] calculatePerpendicularBisectingSegment(float xi, 
        float yi, float xj, float yj, float xMin, float xMax, float yMin, 
        float yMax) {

        float xmidpoint = (xi + xj)/2.f;
        float ymidpoint = (yi + yj)/2.f;

        float yLeft = yMax;
        float yRight = yMin;
        float xLeft, xRight;

        if (xi == xj) {
            yLeft = ymidpoint;
            yRight = ymidpoint;

            xLeft = xMin;

            xRight = xMax;

        } else if (yi == yj) {
            xLeft = xmidpoint;
            xRight = xmidpoint;

            yLeft = yMin;
            yRight = yMax;

        } else {

            float slope = calculateTangentSlope(xi, yi, xj, yj);

            xLeft = calculateX(xmidpoint, ymidpoint, yLeft, slope);
            xRight = calculateX(xmidpoint, ymidpoint, yRight, slope);

            if (xLeft > xRight) {
                float tmp = xLeft;
                xLeft = xRight;
                xRight = tmp;
                tmp = yLeft;
                yLeft = yRight;
                yRight = tmp;
            }

            if (xLeft > xMax) {
                xLeft = xMax;
                yLeft = calculateY(xmidpoint, ymidpoint, xLeft, slope);
            }
            if (xLeft < xMin) {
                xLeft = xMin;
                yLeft = calculateY(xmidpoint, ymidpoint, xLeft, slope);
            }
            if (xRight > xMax) {
                xRight = xMax;
                yRight = calculateY(xmidpoint, ymidpoint, xRight, slope);
            }
            if (xRight < xMin) {
                xRight = xMin;
                yRight = calculateY(xmidpoint, ymidpoint, xRight, slope);
            }
            if (yRight > yMax) {
                yRight = yMax;
                xRight = calculateX(xmidpoint, ymidpoint, yRight, slope);
            }
            if (yRight < yMin) {
                yRight = yMin;
                xRight = calculateX(xmidpoint, ymidpoint, yRight, slope);
            }
            if (yLeft > yMax) {
                yLeft = yMax;
                xLeft = calculateX(xmidpoint, ymidpoint, yLeft, slope);
            }
            if (yLeft < yMin) {
                yLeft = yMin;
                xLeft = calculateX(xmidpoint, ymidpoint, yLeft, slope);
            }
        }

        return new float[]{xLeft, yLeft, xRight, yRight};
    }
   
    /**
     *
     * @param x0
     * @param y0
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @param generalEps
     * @return
     */
    public static boolean segmentIsWithinSegment(float x0, float y0, float x1, 
        float y1, float x2, float y2, float x3, float y3, float generalEps) {

        // check that slopes are the same.
        // then range checks
        float segment01Slope = (y1 - y0)/(x1 - x0);
        float segment23Slope = (y3 - y2)/(x3 - x2);

        float t01 = (segment01Slope < 0) ? -1.f*segment01Slope : segment01Slope;
        float t23 = (segment23Slope < 0) ? -1.f*segment23Slope : segment23Slope;
        if (Math.abs(t01 - t23) > generalEps) {
            return false;
        }

        //re-order by incr xTopPolygon and yTopPolygon
        if ( (x0 > x1) || ((x0 == x1) && (y0 > y1)) ) {
            float tmp = x0;
            x0 = x1;
            x1 = tmp;
            tmp = y0;
            y0 = y1;
            y1 = tmp;
        }
        if ( (x2 > x3) || ((x2 == x3) && (y2 > y3)) ) {
            float tmp = x2;
            x2 = x3;
            x3 = tmp;
            tmp = y2;
            y2 = y3;
            y3 = tmp;
        }
        // Now we have x0 < x1   and   x2 < x3


        //  look for  0  2  3  1  in xTopPolygon  and matching pattern in yTopPolygon

        if ( (x0 <= x2) && (x3 <= x1) ) {
            if ( (y0 <= y1) && (y2 <= y3) && (y0 <= y2) && (y3 <= y1) ) {
                return true;
            } else if ((y0 >= y1) && (y2 >= y3) && (y0 >= y2) && (y3 >= y1)) {
                return true;
            }
        } else {
            //  look for  2  0  1  3  in xTopPolygon  and matching pattern in yTopPolygon
            if (segment01Slope > 0) {
                if ( (x2 <= x0) && (x1 <= x3) ) {
                    if ( (y1 <= y0) && (y3 <= y2) && (y2 <= y0) && (y1 <= y3) ) {
                        return true;
                    }
                }
            } else {
                if ( (x2 <= x0) && (x1 <= x3) ) {
                    if ( (y0 >= y1) && (y2 >= y3) && (y2 >= y0) && (y1 >= y3) ) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    /**
     *
     * @param xPolygon
     * @param yPolygon
     * @return
     */
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
         *   pt_0             pt_k | pt_0              pt_k  |   pt_0      pt_k      |  <no point>      <no point>
         *                         |                         |                       |
         *
         *                                                                           | (4)     peak
         *                                                                           |
         *                                                                           |      pt_i       pt_j
         *                                                                           |  ----------------------
         *                                                                           |
         *                                                                           | <no point>      <no point>
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

}
