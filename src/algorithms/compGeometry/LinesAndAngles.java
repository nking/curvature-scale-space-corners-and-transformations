package algorithms.compGeometry;

/**
 *
 * @author nichole
 */
public class LinesAndAngles {

    /**
     * determine the cross product and return negative number if the 2nd
     * set of values, p1, are clockwise from the first, p1, w.r.t. origin (0,0).
     *
     *    o P2
     *    .
     *    .
     *    .   o P1      <--- P2 is counterclockwise from P1
     *    .
     *    0
     *
     *          o P2
     *          .
     *  P1 o    .       <--- P2 is clockwise from P1
     *          .
     *          .
     *          0
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static int crossProduct(int x1, int y1, int x2, int y2) {
        return ((x1*y2) - (x2*y1));
    }
    public static double crossProduct(double x1, double y1, double x2, double y2) {
        return ((x1*y2) - (x2*y1));
    }
    public static double crossProduct(float x1, float y1, float x2, float y2) {
        return ((x1*y2) - (x2*y1));
    }

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

    public static int direction(int x1, int y1, int x2, int y2, int x3, int y3) {

        int x31 = x3 - x1;
        int y31 = y3 - y1;

        int x21 = x2 - x1;
        int y21 = y2 - y1;

        return crossProduct(x31, y31, x21, y21);
    }

    public static double direction(double x1, double y1, double x2, double y2, double x3, double y3) {

        double x31 = x3 - x1;
        double y31 = y3 - y1;

        double x21 = x2 - x1;
        double y21 = y2 - y1;

        return crossProduct(x31, y31, x21, y21);
    }

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

    public static double calculatePolarSineTheta(double x1, double y1, double x2, double y2) {
        //System.out.println("(" + x1 + "," + y1 + ")," + "(" + x2 + "," + y2 + ")");
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

    public static double calculatePolarSineTheta(float x1, float y1, float x2, float y2) {
        //System.out.println("(" + x1 + "," + y1 + ")," + "(" + x2 + "," + y2 + ")");
        // determine quadrant of points
    	if (x2 == x1) {
            if (y2 == y1) {
                return 0.0f;
            }
    		return 0.0f;//1.0;
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

    /** from http://en.wikipedia.org/wiki/Line-line_intersection
         *
         *    (x1, y1)
         *         *
         *          \    (x4, y4)
         *           \      *
         *            \    /
         *             \  /
         *              \/
         *            P /\
         *             /  \
         *            /    \
         *           /      \
         *          *        \
         *    (x3, y3)        \
         *                     \
         *                      *
         *                  (x2, y2)
         * @param x1
         * @param y1
         * @param x2
         * @param y2
         * @param x3
         * @param y3
         * @param x4
         * @param y4
         * @return the intersection of the 2 lines as a 2 item double array holding x then y, else
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

    public static double[] intersectionOf2Segments(float x1, float y1, float x2, float y2,
        float x3, float y3, float x4, float y4, float eps) {

        double[] p = intersectionOf2Lines(x1, y1, x2, y2, x3, y3, x4, y4);

        if (p == null) {
            return null;
        }

        /*         \
         *           \
         *            *\        the intersection is not within bounds of both lines.
         *         /     \
         *       /
         *     /
         *
         *
         */
        if (x3 > x4) {
            float tmp = x3;
            x3 = x4;
            x4 = tmp;
            tmp = y3;
            y3 = y4;
            y4 = tmp;
        }
        if (x1 > x2) {
            float tmp = x1;
            x1 = x2;
            x2 = tmp;
            tmp = y1;
            y1 = y2;
            y2 = tmp;
        }

        if (x1 > x2) {
            //  x2|-------|x1
            if ((p[0] < (x2 - eps)) || (p[0] > (x1 + eps))) {
                return null;
            }
        } else {
            //  x1|-------|x2
            if ((p[0] < (x1 - eps)) || (p[0] > (x2 + eps))) {
                return null;
            }
        }
        if (x3 > x4) {
            //  x4|-------|x3
            if ((p[0] < (x4 - eps)) || (p[0] > (x3 + eps))) {
                return null;
            }
        } else {
            //  x3|-------|x4
            if ((p[0] < (x3 - eps)) || (p[0] > (x4 + eps))) {
                return null;
            }
        }

        if (y1 > y2) {
            //  y2|-------|y1
            if ((p[1] < (y2 - eps)) || (p[1] > (y1 + eps))) {
                return null;
            }
        } else {
            //  y1|-------|y2
            if ((p[1] < (y1 - eps)) || (p[1] > (y2 + eps))) {
                return null;
            }
        }
        if (y3 > y4) {
            //  y4|-------|y3
            if ((p[1] < (y4 - eps)) || (p[1] > (y3 + eps))) {
                return null;
            }
        } else {
            //  y3|-------|y4
            if ((p[1] < (y3 - eps)) || (p[1] > (y4 + eps))) {
                return null;
            }
        }

        return p;
    }

    public static float calculateSlope(float x1, float y1, float x2, float y2) {
        return ( y2 - y1 ) / (x2 - x1);
    }
    public static float calculateTangentSlope(float x1, float y1, float x2, float y2) {
        return -1.f*( x2 - x1 ) / (y2 - y1);
    }

    public static float calculateX(float x0, float y0, float y, float slope) {
        return ( x0 + ((y - y0)/slope) );
    }

    //  (y-y0)/(x-x0) = slope
    public static float calculateY(float x0, float y0, float x, float slope) {
        return ( y0 + slope*(x - x0) );
    }

    public static float[] calculatePerpendicularBisectingSegment(float xi, float yi, float xj, float yj,
        float xMin, float xMax, float yMin, float yMax) {

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

    public static double distToLine(float x0, float y0, float x1, float y1,
        float xPoint, float yPoint) {

        float xLineDiff = (x1 - x0);
        float xPointDiff = (xPoint - x0);
        float yLineDiff = (y1 - y0);
        float yPointDiff = (yPoint - y0);

        double prDist = Math.sqrt(  xPointDiff*xPointDiff + yPointDiff*yPointDiff );

        /*
         *     (xPoint, yPoint)
         *        *
         *       /
         *   PR /                        perpendicular distance = PR sin theta
         *     /theta
         *    *-------------------*
         *  (x0,y0)            (x1,y1)
         *
         *  |A||B|cos(theta) = A.B = (x1-x0)(xPoint-x0) + (y1-y0)(yPoint-y0)
         *  |A|= sqrt( (x1-x0)^2 + (y1-y0)^2)
         *  |B|= sqrt( (xPoint-x0)^2 + (yPoint-y0)^2)
         *
         *   cos(theta) = (x1-x0)(xPoint-x0) + (y1-y0)(yPoint-y0) / (sqrt( (x1-x0)^2 + (y1-y0)^2) * sqrt( (xPoint-x0)^2 + (yPoint-y0)^2))
         */

        double a = (xLineDiff * xPointDiff + yLineDiff * yPointDiff);
        double b = Math.sqrt(xLineDiff*xLineDiff + yLineDiff*yLineDiff)
            * Math.sqrt(xPointDiff*xPointDiff + yPointDiff*yPointDiff);
        double theta = Math.acos( a/b );

        return prDist*Math.sin(theta);
    }

    /*
     *
     * determine angle subtended by line and point where we are returning the angle
     * (NOTE:  set the first point as the one near the ngle of interest)
     *
     *      (x1,y1)
     *         /
     *        /      *(xPoint,yPoint)
     *       /
     *      /theta
     *   (x0,y0)
     *
     *
     */
    public static double angleBetweenPointAndLine(float x0, float y0, float x1, float y1,
        float xPoint, float yPoint) {

        /*
         * Use Law of Cosines
         *
         *                   (x1, y1)
         *                      *
         *                    .
         *        side c    .  \/.
         *                .    B               c^2 = a^2 + b^2 - 2ab cos C
         *              .         .            or  a^2 = b^2 + c^2 - 2bc cos A
         *            .
         *          .              . side a
         * (x0,y0)*  _\A
         *            .             .
         *                .
         *                    .      .
         *            side b      . C/
         *                             *(xp,yp)
         *
         *  cosine A = (b^2 + c^2 - a^2)/2bc
         */

        double bsq = distSquared(x0, y0, xPoint, yPoint);
        double csq = distSquared(x0, y0, x1, y1);
        double asq = distSquared(xPoint, yPoint, x1, y1);

        double cosA = ( bsq + csq - asq)/( 2.*Math.sqrt(bsq) * Math.sqrt(csq));
        double theta = Math.acos( cosA );// value is in range the 0.0 through pi

        // use a convention for the sign
        theta = ((xPoint - x0) > 0) ? theta : -1.*theta;

        return theta;
    }

    /**
     * For Clockwise, solving theta = C in diagram below.
     * For Counter Clockwise, solving theta = A in diagram below.
     *
     *                   (xp, yp)
     *                      *
     *                    .
     *        side c    .  \/.             Law of cosines
     *                .    B               c^2 = a^2 + b^2 - 2ab cos C
     *              .         .            or  a^2 = b^2 + c^2 - 2bc cos A
     *            .
     *          .              . side a
     * (x0,y0)*  _\A
     *            .             .
     *                .
     *                    .      .
     *            side b      . C/
     *                             *(x1,y1)
     *
     *   Right triangle cosine and sine
     *     r /
     *      /
     *     /theta    x = r * cos(theta)
     *    -------
     *       x
     *
     * @param x0
     * @param y0
     * @param x1
     * @param y1
     * @param xPoint
     * @param yPoint
     * @param calculateClockWise
     * @return the angle between the edge 0:1 and edge 1:Point which is C for Clockwise, else A for Counter Clockwise
     */
    public static double angleBetweenPointAndLine(float x0, float y0, float x1, float y1,
        float xPoint, float yPoint, boolean calculateClockWise) {

        if (x1 < x0) {
            throw new IllegalArgumentException("Expecting x0 <= x1.");
        }

        float slope01 = (y1 - y0)/(x1 - x0);
        float slope1p = (yPoint - y1)/(xPoint - x1);

        // FIRST handling the cases where edge 0:1 and/or 1:point are horizontal or vertical
        //   to use right triangle forumulas

        double angle = angleBetweenPointAndLineWithRightTriangle(x0, y0, x1, y1, xPoint, yPoint, calculateClockWise);

        if (angle > -1) {
        	return angle;
        }

        // ELSE transform points to reference frame where edge 0:1 is horizontal
        // And calculate angles using Law Of Cosines.
        // Then see if 'CW' or 'CCW' w.r.t 0:1 means an inside or outside angle should be
        //    returned.

        // Transform points 1 and P so that 0:1 is horizontal.  Using the origin as center of rotation.

        double thetaTransformCW = 0.;
       /*   Right triangle cosine and sine
        *     r /
        *      /
        *     /theta    x = r * cos(theta)
        *    -------
        *       x
        */
        if (x0 == x1) {
        	thetaTransformCW = -90.;
        } else if (y0 == y1) {
        	thetaTransformCW = 0.;
        } else {
        	thetaTransformCW = -1*Math.atan2((y1-y0), (x1-x0));
        }

        // P'_x = origin_x + ((x-origin_x)*math.cos(theta) - (y-origin_y)*math.sin(theta))
        // P'_y = origin_y + ((y-origin_y)*math.cos(theta) + (x-origin_x)*math.sin(theta))
        float x00 = x0;
        float y00 = y0;

        float x11 = x0 + (float)((x1 - x0)*Math.cos(thetaTransformCW) - (y1 - y0)*Math.sin(thetaTransformCW));
        float y11 = y0 + (float)((y1 - y0)*Math.cos(thetaTransformCW) + (x1 - x0)*Math.sin(thetaTransformCW));

        float xpp = x0 + (float)((xPoint - x0)*Math.cos(thetaTransformCW) - (yPoint - y0)*Math.sin(thetaTransformCW));
        float ypp = y0 + (float)((yPoint - y0)*Math.cos(thetaTransformCW) + (xPoint - x0)*Math.sin(thetaTransformCW));

        /*
         *   right triangles have already been handled, so we can limit cases by (xpp, ypp) position
         *      relative to 0:1
         * _________________   _________________   _________________
         * |case 0:         |  |case 1:         |  |case 2:         |
         * |                |  |                |  |                |
         * | P              |  |       P        |  |             P  |   <=== inner angles in polygon
         * |    *     *     |  |    *     *     |  |    *     *     |
         * |    0     1     |  |    0     1     |  |    0     1     |
         * |                |  |                |  |                |
         * ------------------  ------------------  ------------------
         *
         * _________________   _________________   _________________
         * |case 3:         |  |case 4:         |  |case 5:         |
         * |                |  |                |  |                |
         * |                |  |                |  |                |  <=== the inner angle is measured
         * |    *     *     |  |    *     *     |  |    *     *     |       but we want the outer angle
         * |    0     1     |  |    0     1     |  |    0     1     |       for the answer
         * | P              |  |       P        |  |             P  |
         * ------------------  ------------------  ------------------
         *
         */

        double bsq = distSquared(x00, y00, x11, y11);
        double csq = distSquared(x00, y00, xpp, ypp);
        double asq = distSquared(xpp, ypp, x11, y11);

        if (calculateClockWise) {
            // c^2 = a^2 + b^2 - 2ab cos C
            // cos C = (a^2 + b^2 - c^2)/2ab
            double cosC = (asq + bsq - csq)/(2.* Math.sqrt(bsq)*Math.sqrt(asq));
            double c = Math.acos(cosC);
            double theta = c;

            if (ypp > y00) {
                // case 0, 1, 2
                return theta;
            } else {
                return 2. * Math.PI - theta;
            }
        } else {
            // a^2 = c^2 + b^2 - 2bc cos A
            // cos A = (c^2 + b^2 - a^2)/2bc
            double cosA = (csq + bsq - asq)/(2.* Math.sqrt(bsq)*Math.sqrt(csq));
            double a = Math.acos(cosA);
            double theta = a;

            if (ypp > y00) {
                // case 0, 1, 2
                return theta;
            } else {
                return 2. * Math.PI - theta;
            }
        }
    }

    /**
     * use right triangles when possible to measure
     *   either the angle between edge 0:1 and 1:point when calculateClockWise=true
     *   OR the angle between edge 0:1 and 0:point when calculateClockWise=false.
     * It includes calculations to learn whether the angle is outside or inside the
     * triangle etc.
     *
     * @param x0
     * @param y0
     * @param x1
     * @param y1
     * @param xPoint
     * @param yPoint
     * @param calculateClockWise
     * @return
     */
    private static double angleBetweenPointAndLineWithRightTriangle(float x0,
		float y0, float x1, float y1, float xPoint, float yPoint, boolean calculateClockWise) {

    	//TODO:  simplify this one day...

    	if (x1 < x0) {
            throw new IllegalArgumentException("Expecting x0 < x1 or (x0 == x1 and y1 > y1).");
        }

        if (calculateClockWise) {
        	// measuring angle between 0:1 and 1:Point
            float slope01 = (y1 - y0)/(x1 - x0);
            float slope1p = (yPoint - y1)/(xPoint - x1);

            if (slope1p == 0) {
                if (slope01 == 0) {
                    if (xPoint > x1) {
                        return Math.PI;
                    } else {
                        return 2.*Math.PI;
                    }
                }
                // can use  sine(theta) = y/r and complement or cosine(theta) = x/r
                double r = Math.sqrt(distSquared(x0, y0, x1, y1));
                double y =  (y0 - y1);

                if ((xPoint > x1) && (y0 > y1)) {
                	double angle = Math.asin(y/r);
                    return Math.PI - angle;
                } else if ((xPoint > x1) && (yPoint <= y1) && (x1 > x0)) {
                	y *= -1;
                    double angle = Math.asin(y/r);
                    return Math.PI + angle;
                } else if ((xPoint > x1) && (yPoint <= y1)) {
                	 y *= -1;
                     double angle = Math.asin(y/r);
                     return (2.*Math.PI) - angle;

                } else if ((xPoint < x1) && (y0 > y1)) {
                    double angle = Math.asin(y/r);
                    return 2.*Math.PI - angle;

                } else if ((xPoint <= x1) && (yPoint >= y1)) {
                    y *= -1;
                    double angle = Math.asin(y/r);
                    return angle;
                }
            } else if (slope01 == 0) {
                // can use cosine(theta) = x/r and complement
                double r = Math.sqrt(distSquared(xPoint, yPoint, x1, y1));
                double x = (xPoint - x1);
                if ((xPoint < x1) && (yPoint > y1)) {// COVERAGE HERE?
                	double angle = Math.acos(x / r);
                    return Math.PI - angle;
                } else if ((xPoint < x1) && (yPoint < y1)){
                	double angle = Math.acos(x / r);
                    return Math.PI + angle;
                } else if (xPoint == x1) {
                    double angle = Math.acos(x / r);
                    if (yPoint < y0) {
                        return Math.PI + angle;
                    } else {
                        return angle;
                    }
                } else {
                    // xPoint > x1
                    if (yPoint > y1) {
                        double angle = Math.acos(x / r);
                        return Math.PI - angle;
                    }
                }
            } else if ((x1 == x0) && (xPoint == x1)) {
                if (yPoint < y1) {
                    return 2.*Math.PI;
                } else {
                    return Math.PI;
                }
            } else if (x1 == x0) {

                double r = Math.sqrt(distSquared(xPoint, yPoint, x1, y1));

                if (xPoint < x1) {
                    double x = (x1 - xPoint);
                    double angle = Math.acos(x / r);
                    if (yPoint <= y1) {
                        return Math.PI/2. - angle;
                    } else {
                        return Math.PI/2. + angle;
                    }
                } else {
                    // xPoint > x1
                    double x = (xPoint - x1);
                    double angle = Math.acos(x / r);
                    if (yPoint <= y1) {
                        return (3.*Math.PI/2.) + angle;
                    } else {
                        return (3.*Math.PI/2.) - angle;
                    }
                }
            }

            if (Math.abs(slope1p - slope01) < 0.01) { //TODO:  calc an eps instead of 0.01
                if (xPoint > x1) {
                    return Math.PI;
                } else {
                    return 0;
                }
            }

        } else {

        	// =============  Counter Clockwise ===============

        	// measuring angle between 0:1 and 1:Point
            float slope01 = (y1 - y0)/(x1 - x0);
            float slope0p = (yPoint - y0)/(xPoint - x0);

            if (slope0p == 0) {
                if (slope01 == 0) {
                    if (xPoint < x0) {
                        return Math.PI;
                    } else {
                        return 2.*Math.PI;
                    }
                }

                // can use  sine(theta) = y/r and complement or cosine(theta) = x/r
                double r = Math.sqrt(distSquared(x0, y0, x1, y1));
                double y =  (y0 - y1);

                /*     1                        1
                 *          1                      1
                 *     0:P        OR          P:0
                 *          1                      1
                 *      1                        1
                 */
                if ((xPoint >= x0) && (yPoint <= y1)) {
                    y *= -1;
                    double angle = Math.asin(y/r);
                    return 2.*Math.PI - angle;
                } else if ((xPoint >= x0) && (yPoint > y1)) {
                	double angle = Math.asin(y/r);
                    return angle;
                } else if ((xPoint < x0) && (yPoint <= y1)) {
                	y *= -1;
                    double angle = Math.asin(y/r);
                    return Math.PI - angle;
                } else if ((xPoint > x0) && (yPoint > y1)) {
                    double angle = Math.asin(y/r);
                    return (2.*Math.PI) - angle;
                }
            } else if (slope01 == 0) {
                /*     P  P  P
                 *       0:1
                 *     P  P  P */
                double r = Math.sqrt(distSquared(xPoint, yPoint, x0, y0));
                double y =  (yPoint - y0);
                if (yPoint > y0) {
                    double angle = Math.asin(y / r);
                    if (xPoint >= x0) {
                        return angle;
                    } else {
                        return Math.PI - angle;
                    }
                } else {
                    y *= -1;
                    double angle = Math.asin(y / r);
                    if (xPoint <= x0) {
                        return Math.PI + angle;
                    } else {
                        return 2.*Math.PI - angle;
                    }
                }
            } else if ((x1 == x0) && (xPoint == x1)) {
                if (yPoint < y0) {
                    return Math.PI;
                } else {
                    return 0;
                }
            } else if (x1 == x0) {
                /*     P
                 *
                 *     1
                 *  P     P
                 *     0
                 *  P     P
                 *     P
                 */
                double r = Math.sqrt(distSquared(xPoint, yPoint, x0, y0));
                double y = (yPoint - y0);

                if ((xPoint <= x0) && (yPoint >= y0)) {
                    double angle = Math.asin(y / r);
                    return Math.PI/2. - angle;
                } else if ((xPoint <= x0) && (yPoint < y0)) {
                    y *= -1;
                    double angle = Math.asin(y / r);
                    return Math.PI/2. + angle;
                } else if ((xPoint > x0) && (yPoint >= y0)) {
                    double angle = Math.asin(y / r);
                    return (3.*Math.PI/2.) + angle;
                } else {
                    y *= -1;
                    double angle = Math.asin(y / r);
                    return (3.*Math.PI/2.) - angle;
                }
            }

            if (Math.abs(slope0p - slope01) < 0.01) { //TODO:  calc an eps instead of 0.01
                if (xPoint < x0) {
                    return Math.PI;
                } else {
                    return 0;
                }
            }
        }
        return -1;
	}

	public static boolean pointIsInLine(float xPoint, float yPoint,
        float xL, float yL, float xR, float yR, double epsilon) {

        // determine slope for line and plugin xPoint.
        //    does the resulting y = yPoint within epsilon?

        /*
         *   yL - yR            yL - yPoint
         *   ------- = slope =  -----------
         *   xL - xR            xL - xPoint
         *
         *   yL - yPoint = slope*(xL - xPoint)
         *
         */
        float slope = ( yL - yR ) / ( xL - xR );

        float yPointPredicted = yL - (slope * (xL - xPoint));

        if (Math.abs( yPointPredicted - yPoint) < epsilon) {
            return true;
        }
        return false;
    }

    public static boolean segmentIsWithinSegment(float x0, float y0, float x1, float y1,
        float x2, float y2, float x3, float y3, float generalEps) {

        // check that slopes are the same.
        // then range checks
        float segment01Slope = (y1 - y0)/(x1 - x0);
        float segment23Slope = (y3 - y2)/(x3 - x2);

        float t01 = (segment01Slope < 0) ? -1.f*segment01Slope : segment01Slope;
        float t23 = (segment23Slope < 0) ? -1.f*segment23Slope : segment23Slope;
        if (Math.abs(t01 - t23) > generalEps) {
            return false;
        }

        //re-order by incr x and y
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


        //  look for  0  2  3  1  in x  and matching pattern in y

        if ( (x0 <= x2) && (x3 <= x1) ) {
            if ( (y0 <= y1) && (y2 <= y3) && (y0 <= y2) && (y3 <= y1) ) {
                return true;
            } else if ((y0 >= y1) && (y2 >= y3) && (y0 >= y2) && (y3 >= y1)) {
                return true;
            }
        } else {
            //  look for  2  0  1  3  in x  and matching pattern in y
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


    public static float[] calcAreaAndCentroidOfSimplePolygon(float[] xHull, float[] yHull) {

        if (xHull.length < 3) {
            return null;
        }

        float area = 0;

        boolean allPointsArePositive = true;

        // from wikipedia
        for (int i = 0; i < (xHull.length - 1); i++) {
            float a = xHull[i] * yHull[i+1];
            float b = xHull[i+1] * yHull[i];
            area += (a-b);

            if ((xHull[i] < 0) || (yHull[i] < 0)) {
                allPointsArePositive = false;
            }
        }
        if ((xHull[xHull.length - 1] < 0) || (yHull[xHull.length - 1] < 0)) {
            allPointsArePositive = false;
        }
        if ((area < 0) && allPointsArePositive) {
            area *= -1;
        }
        area *= 0.5f;


        float xc = 0;
        for (int i = 0; i < (xHull.length - 1); i++) {
            float a = xHull[i] * yHull[i+1];
            float b = xHull[i+1] * yHull[i];

            xc += ( (xHull[i] + xHull[i+1]) * (a - b) );
        }
        xc *= (1./(6.*area));
        if ((xc < 0) && allPointsArePositive) {
            xc *= -1;
        }


        float yc = 0;
        for (int i = 0; i < (xHull.length - 1); i++) {
            float a = xHull[i] * yHull[i+1];
            float b = xHull[i+1] * yHull[i];

            yc += ( (yHull[i] + yHull[i+1]) * (a - b) );
        }
        //yc = yHull[0] - yc; // make it relative to zero point
        yc *= (1./(6.*area));
        if ((yc < 0) && allPointsArePositive) {
            yc *= -1;
        }

        return new float[]{area, xc, yc};
    }
}
