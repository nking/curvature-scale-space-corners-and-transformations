package algorithms.compGeometry;

import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.util.PairInt;

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
     * direction is counterclockwise.
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
     * @param p1
     * @param p2
     * @param p3
     * @return 
     */
    public static <T extends PairInt>  double direction(T p1, T p2, T p3) {
        
        int x1 = p1.getX();
        int y1 = p1.getY();
        int x2 = p2.getX();
        int y2 = p2.getY();
        int x3 = p3.getX();
        int y3 = p3.getY();
        
        double d = ((x2 - x1)*(y3 - y1)) - ((y2 - y1)*(x3 - x1));
        
        return d;
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
     * @param p1
     * @param p2
     * @param p3
     * @return 
     */
    public static double direction(PairIntWithIndex p1, PairIntWithIndex p2, 
        PairIntWithIndex p3) {

        double d = ((p2.getX() - p1.getX())*(p3.getY() - p1.getY())) -
            ((p2.getY() - p1.getY())*(p3.getX() - p1.getX()));
                
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
     * from pseudocode in "Introduction to Algorithms" by Cormen et al.
     * 
     * runtime complexity is O(N lg N) where N is the number of line segments.
     * worse case runtime complexity is O(N^2).
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

        /*
        direction is positive when direction is clockwise,
        else if zero the vectors are collinear, else if negative if the
        direction is counterclockwise.
         * 
         *           p3
         *  p2      / 
         *   \    /
         *     p1
        */
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
        } else if ( (d1 == 0) && onSegment(x3, y3, x4, y4, x1, y1) ) {
            return true;
        } else if ( (d2 == 0) && onSegment(x3, y3, x4, y4, x2, y2) ) {
            return true;
        } else if ( (d3 == 0) && onSegment(x1, y1, x2, y2, x3, y3) ) {
            return true;
        } else if ( (d4 == 0) && onSegment(x1, y1, x2, y2, x4, y4) ) {
            return true;
        }
        return false;
    }
   
    /**
     * calculate the direction of change for the 2 vectors 
     * P1:P2 to P1:P3  returns positive when direction is clockwise,
     * else if zero the vectors are collinear, else if negative the
     * direction is counterclockwise.
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
    public static int direction(int x1, int y1, int x2, int y2, int x3, int y3) {

        int x31 = x3 - x1;
        int y31 = y3 - y1;

        int x21 = x2 - x1;
        int y21 = y2 - y1;

        int cp = ((x31*y21) - (x21*y31));
        
        return cp;
    }
    
    /**
     * determine the cross product and return negative number if the 2nd
     * set of values, p1, are clockwise from the first, p1, w.r.t. origin (0,0).
     <pre>
          o P2
          .
          .
          .   o P1       P2 is counterclockwise from P1
          .
          0
     ---------------------------------------------------------
                o P2
                .
        P1 o    .        P2 is clockwise from P1
                .
                .
                0
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
    test whether (x3, y3) is between (x1,y1) and (x2,y2) on same segment.
    */
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
     * from pseudocode in "Introduction to Algorithms" by Cormen et al.
     * 
     * runtime complexity is O(N lg N) where N is the number of line segments.
     * worse case runtime complexity is O(N^2).
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
    public static boolean linesIntersect(float x1, float y1,
        float x2, float y2, float x3, float y3, float x4, float y4) {

        double d1 = -1*direction(x3, y3, x4, y4, x1, y1);
        double d2 = -1*direction(x3, y3, x4, y4, x2, y2);
        double d3 = -1*direction(x1, y1, x2, y2, x3, y3);
        double d4 = -1*direction(x1, y1, x2, y2, x4, y4);

        if (
            (((d1 > 0) && (d2 < 0)) || ((d1 < 0) && (d2 > 0)))
            &&
            (((d3 > 0) && (d4 < 0)) || ((d3 < 0) && (d4 > 0)))
        ) {
            return true;
        } else if ( (d1 == 0) && onSegment(x3, y3, x4, y4, x1, y1) ) {
            return true;
        } else if ( (d2 == 0) && onSegment(x3, y3, x4, y4, x2, y2) ) {
            return true;
        } else if ( (d3 == 0) && onSegment(x1, y1, x2, y2, x3, y3) ) {
            return true;
        } else if ( (d4 == 0) && onSegment(x1, y1, x2, y2, x4, y4) ) {
            return true;
        }
        return false;
    }
    
    /**
     * determine the cross product and return negative number if the 2nd
     * set of values, p1, are clockwise from the first, p1, w.r.t. origin (0,0).
     <pre>
          o P2
          .
          .
          .   o P1       P2 is counterclockwise from P1
          .
          0
     ---------------------------------------------------------
                o P2
                .
        P1 o    .        P2 is clockwise from P1
                .
                .
                0
      </pre>
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
     */
    public static float crossProduct(float x1, float y1, float x2, float y2) {
        return ((x1*y2) - (x2*y1));
    }
    
    static boolean onSegment(float x1, float y1,
        float x2, float y2, float x3, float y3) {

        float minx12 = (x1 < x2) ? x1 : x2;
        float miny12 = (y1 < y2) ? y1 : y2;

        float maxx12 = (x1 > x2) ? x1 : x2;
        float maxy12 = (y1 > y2) ? y1 : y2;

        if ( (minx12 <= x3) && (x3 <= maxx12) && (miny12 <= y3) && (y3 <= maxy12) ) {
            return true;
        } else {
            return false;
        }
    }

    public static boolean pointIsInLine(int xPt, int yPt, int x2, int y2, int x3, int y3) {
        
        int diffY = y3 - y2;
        int diffX = x3 - x2;
        
        int numer = (diffY*xPt) - (diffX*yPt) + (x3*y2) - (y3*x2);
        
        //double denom = Math.sqrt((diffX*diffX) + (diffY*diffY));
        //double dist = Math.abs(numer)/denom;
        
        if (numer != 0) {
            return false;
        }
        
        // check that pt is between bounds
        if (x2 < x3) {
            if (xPt < x2 || xPt > x3) {
                return false;
            }
        } else if (x2 > x3) {
            if (xPt > x2 || xPt < x3) {
                return false;
            }
        }
        if (y2 < y3) {
            if (yPt < y2 || yPt > y3) {
                return false;
            }
        } else if (y2 > y3) {
            if (yPt > y2 || yPt < y3) {
                return false;
            }
        }
        return true;
    }
    
    public static boolean pointIsInLine(float xPt, float yPt, float x2, float y2, 
        float x3, float y3) {
        
        //TODO: this one needs a tolerance
        
        float diffY = y3 - y2;
        float diffX = x3 - x2;
        
        float numer = (diffY*xPt) - (diffX*yPt) + (x3*y2) - (y3*x2);
        
        //double denom = Math.sqrt((diffX*diffX) + (diffY*diffY));
        //double dist = Math.abs(numer)/denom;
        
        if (numer != 0) {
            return false;
        }
        
        // check that pt is between bounds
        if (x2 < x3) {
            if (xPt < x2 || xPt > x3) {
                return false;
            }
        } else if (x2 > x3) {
            if (xPt > x2 || xPt < x3) {
                return false;
            }
        }
        if (y2 < y3) {
            if (yPt < y2 || yPt > y3) {
                return false;
            }
        } else if (y2 > y3) {
            if (yPt > y2 || yPt < y3) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Calculate the angle of segment P3:P1 sweeping clockwise
     * to segment P3:P2.
     * 
     * Internally, the method uses the law of cosines and
     * the direction method.
     * <pre>
     *   P1     P2
     * 
     *      P3
     * </pre>
     * NOTE that if given a single point, it will return NaN
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return 
     */
    public static double calcClockwiseAngle(int x1, int y1, 
        int x2, int y2, int x3, int y3) {

        double d = direction(x3, y3, x2, y2, x1, y1);

        double angleA = calcAngle(x1, y1, x2, y2, x3, y3);
       
        if (d > 0) {
            angleA = 2.*Math.PI - angleA;
        }
        
        return angleA;
    }
 
    /**
     * Calculate the angle of segment P3:P1 sweeping clockwise
     * to segment P3:P2.
     * 
     * Internally, the method uses the law of cosines and
     * the direction method.
     * <pre>
     *   P1     P2
     * 
     *      P3
     * </pre>
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return 
     */
    public static double calcAngle(int x1, int y1, 
        int x2, int y2, int x3, int y3) {

        // c
        int d1d3 = distanceSqEucl(x1, y1, x3, y3);

        // b
        int d2d3 = distanceSqEucl(x2, y2, x3, y3);

        // a
        int d1d2 = distanceSqEucl(x1, y1, x2, y2);

        //cosine A = (b^2 + c^2 - a^2)/2bc
        double cosA = (d1d3 + d2d3 - d1d2)
            / (2 * Math.sqrt((double)d1d3 * (double)d2d3));

        double angleA = Math.acos(cosA);
        
        return angleA;
    }
 
    protected static int distanceSqEucl(int x1, int y1, 
        int x2, int y2) {
        
        int diffX = x1 - x2;
        int diffY = y1 - y2;
        return (diffX * diffX + diffY * diffY);
    }
    
    /**
     * NOT YET TESTED
     * given a polar angle and radius, calculate the line perpendicular to that
     * point where it intersects with the image boundary.
     * @param thetaDegress
     * @param radius
     * @param imageWidth
     * @param imageHeight
     * @return new int[] {x1, y1, x2, y2}
     */
    public static int[] calcPolarLineEndPoints(int thetaDegrees, int radius,
        int imageWidth, int imageHeight) {
        
        if (thetaDegrees == 0) {
            return new int[]{radius, 0, radius, imageHeight - 1};
        } else if (thetaDegrees == 180) {
            return new int[]{-radius, 0, -radius, imageHeight - 1};
        } else if (thetaDegrees == 90) {
            return new int[]{0, radius, imageWidth - 1, radius};
        } else if (thetaDegrees == 270) {
            return new int[]{0, -radius, imageWidth - 1, -radius};
        } else if (thetaDegrees < 90) {
            return calcPolarLineEndPointsForThetaLT90(thetaDegrees, radius,
                imageWidth, imageHeight);
        } else if (thetaDegrees < 180) {
            int[] ep = calcPolarLineEndPointsForThetaLT90(
                180 - thetaDegrees, radius,
                imageWidth, imageHeight);
            ep[0] *= -1;
            ep[2] *= -1;
            return ep;
        } else if (thetaDegrees < 270) {
            int[] ep = calcPolarLineEndPointsForThetaLT90(
                thetaDegrees - 180, radius,
                imageWidth, imageHeight);
            for (int i = 0; i < ep.length; ++i) {
                ep[i] *= -1;
            }
            return ep;
        } else {
            int[] ep = calcPolarLineEndPointsForThetaLT90(
                360 - thetaDegrees, radius,
                imageWidth, imageHeight);
            ep[1] *= -1;
            ep[3] *= -1;
            return ep;
        }
    } 
    
    private static int[] calcPolarLineEndPointsForThetaLT90(int thetaDegrees, 
        int radius, int imageWidth, int imageHeight) {

        int lx = imageWidth - 1;
        int ly = imageHeight - 1;
        
        double a = (double) thetaDegrees * Math.PI / 180.;
        double b = (Math.PI / 2.) - a;
        
        double xInter1 = (double) radius / Math.cos(a);
        double yInter1 = xInter1 * Math.tan(b);
        if (xInter1 < imageWidth) {
            if (yInter1 < imageHeight) {
                return new int[]{0, (int) Math.round(yInter1), 
                    (int) Math.round(xInter1), 0};
            }
            double y2 = yInter1 - ly;
            // yInter1/y2 = xInter1/x2 ==> x2 = xInter1 * y2 / yInter1
            double x2 = (y2 * xInter1) / (yInter1);
            return new int[]{(int) Math.round(x2), ly, 
                (int) Math.round(xInter1), 0};
        }
        // xInt > imageWidth - 1
        double x2 = xInter1 - lx;
        double y2 = x2 * Math.tan(b);
        if (yInter1 < imageHeight) {
            return new int[]{0, (int) Math.round(yInter1), 
                lx, (int) Math.round(y2)};
        }
        // yInt > imageHeight - 1
        double y3 = yInter1 - ly;
        double x3 = (y3 * xInter1) / (yInter1);
        
        return new int[]{(int) Math.round(x3), ly, 
            lx, (int) Math.round(y2)};
    }
}
