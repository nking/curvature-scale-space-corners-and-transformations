package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;

/**
 * class to hold a few details from the creation of corners that
 * help later to match corners.
 * 
 * @author nichole
 */
public class CornerRegion {
    
    protected final int edgeIdx;
    
    protected final int kMaxIdx;
    
    protected final float[] k;
    
    protected final int[] x;
    
    protected final int[] y;
    
    protected double orientation = Double.MIN_VALUE;
    
    /**
     * constructor with edge index, the index with the maximum curvature in it
     * with respect to nPoints.
     * @param theEdgeIndex
     * @param maxCurvatureIndex
     * @param nPoints 
     */
    public CornerRegion(int theEdgeIndex, int nPoints, int maxCurvatureIndex) {
        
        this.edgeIdx = theEdgeIndex;
        
        if (nPoints < 0) {
            throw new IllegalArgumentException("nPoints must be 0 or larger");
        }
        
        if (nPoints > 0) {
            if (nPoints < 3) {
                throw new IllegalArgumentException(
                "nPoints, if larger than 0, must be at least 3");
            }
            if ((maxCurvatureIndex == 0) || (maxCurvatureIndex > (nPoints - 2))) {
                throw new IllegalArgumentException(
                    "maxCurvatureIndex must be > 0 and less than nPoints-1");
            }
        }
        
        this.kMaxIdx = maxCurvatureIndex;
        
        k = new float[nPoints];
        x = new int[nPoints];
        y = new int[nPoints];
    }
    
    public void set(int index, float kValue, int xCoordinate, int yCoordinate) {
        
        if (index < 0) {
            throw new IllegalArgumentException("index cannot be < 0");
        }
        if (index > (x.length - 1)) {
            throw new IllegalArgumentException("index cannot be larger than " +
                x.length);
        }
        
        k[index] = kValue;
        x[index] = xCoordinate;
        y[index] = yCoordinate;
    }
    
    public float[] getK() {
        return k;
    }
    public int[] getX() {
        return x;
    }
    public int[] getY() {
        return y;
    }
    
    public int getXBeforeMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return x[kMaxIdx - 1];
    }
    
    public int getYBeforeMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return y[kMaxIdx - 1];
    }
    
    public int getXAfterMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return x[kMaxIdx + 1];
    }
    
    public int getYAfterMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return y[kMaxIdx + 1];
    }
    
    /**
     * get the relative angle of this corner (a.k.a. dominant orientation)
     * that represents the angle perpendicular to the corner along the edge it 
     * was extracted from (units are radians)
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of degrees.
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException
     * an exception is thrown if dx=0 or dy=0 for all points owned by this 
     * instance.
     */
    public float getRelativeOrientationInDegrees() throws CornerRegionDegneracyException {
        
        double rotRadians = getRelativeOrientation();
        
        return (float)(rotRadians * 180./Math.PI);
    }
    
    /**
     * get the relative angle of this corner (a.k.a. dominant orientation)
     * that represents the angle perpendicular to the corner along the edge it 
     * was extracted from (units are radians)
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of radians.
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException
     * an exception is thrown if dx=0 or dy=0 for all points owned by this 
     * instance.
     */
    public double getRelativeOrientation() throws CornerRegionDegneracyException {
        
        if (orientation == Double.MIN_VALUE) {
            orientation = calculateOrientation();
        }
        
        return orientation;
    }
    
    /**
     * NOT READY FOR USE YET.  NEEDs MORE TESTING
     * calculate the angle perpendicular to the maximum of curvature.
     * The curvature has to be large enough so that a change in the neighboring
     * points is present (dx and dy cannot both be zero for all the points
     * given) else an exception is thrown.
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of radians.
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException
     * an exception is thrown if dx=0 or dy=0 for all points owned by this 
     * instance.
     */
    protected double calculateOrientation() throws 
        CornerRegionDegneracyException {
        
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        int dx0 = x[kMaxIdx] - x[kMaxIdx - 1];
        int dy0 = y[kMaxIdx] - y[kMaxIdx - 1];
        
        double theta0 = AngleUtil.polarAngleCCW(dx0, dy0);
        if (theta0 != 0) {
            theta0 = (2.*Math.PI) - theta0;
        }
        
        /* determine wheter to add or subtract 90 for vector by determining
        which direction proceeding from point in question is further from
        the point on the other side of the maximum.
        the further is the direction to choose for the resulting change by 90
        degrees.
        
        using quadrants to understand which direction the change is:
        QIII | QIV
        -----|-----
        QII  | QI
        */
        
        // a point extending the line from kMaxIdx to kMaxIdx - 1 is in this quadrant
        int q0 = calculateQuadrant(x[kMaxIdx - 1] - dx0, y[kMaxIdx - 1] - dy0,
            x[kMaxIdx - 1], y[kMaxIdx - 1]);
        
        int prevX0 = x[kMaxIdx - 1] - dx0;
        int prevY0 = y[kMaxIdx - 1] + dy0;
        
        int nextX0 = x[kMaxIdx - 1] + dx0;
        int nextY0 = y[kMaxIdx - 1] - dy0;
        
        int ref0Idx = kMaxIdx+1;
        
        double distSqP0 = (prevX0 - x[ref0Idx])*(prevX0 - x[ref0Idx])
            + (prevY0 - y[ref0Idx])*(prevY0 - y[ref0Idx]);
        
        double distSqN0 = (nextX0 - x[ref0Idx])*(nextX0 - x[ref0Idx])
            + (nextY0 - y[ref0Idx])*(nextY0 - y[ref0Idx]);
        
        if (distSqP0 == distSqN0) {
            
            while (distSqP0 == distSqN0) {
                ref0Idx++;
                if (ref0Idx > (x.length - 1)) {
                    throw new CornerRegionDegneracyException(
                    "need more neighboring points because the slopes are all vertical");
                }
                distSqP0 = (prevX0 - x[ref0Idx])*(prevX0 - x[ref0Idx])
                    + (prevY0 - y[ref0Idx])*(prevY0 - y[ref0Idx]);
                distSqN0 = (nextX0 - x[ref0Idx])*(nextX0 - x[ref0Idx])
                    + (nextY0 - y[ref0Idx])*(nextY0 - y[ref0Idx]);
            }
        }
        
        int q0Furthest = 0;
        
        boolean furthestIsPrev = true;
        
        if (distSqP0 < distSqN0) {
            furthestIsPrev = false;
            q0Furthest = calculateQuadrant(nextX0, nextY0,
                x[kMaxIdx - 1], y[kMaxIdx - 1]);
        } else if (distSqP0 > distSqN0) {
            q0Furthest = calculateQuadrant(prevX0, prevY0,
                x[kMaxIdx - 1], y[kMaxIdx - 1]);
        }
        
        double perp0 = theta0;
        if (q0 < q0Furthest) {
            perp0 -= (Math.PI/2.);
        } else if (q0 > q0Furthest) {
            perp0 += (Math.PI/2.);
        } else {
            double ang0Ext = AngleUtil.polarAngleCCW(
                x[kMaxIdx - 1] - dx0 - x[kMaxIdx - 1], 
                y[kMaxIdx - 1] - dy0 - y[kMaxIdx - 1]);
            if (furthestIsPrev) {
                double ang0P = AngleUtil.polarAngleCCW(
                    prevX0 - x[kMaxIdx - 1], prevY0 - y[kMaxIdx - 1]);
                if (ang0Ext > ang0P) {
                    perp0 -= (Math.PI/2.);
                } else {
                    perp0 += (Math.PI/2.);
                }
            } else {
                double ang0N = AngleUtil.polarAngleCCW(
                    nextX0 - x[kMaxIdx - 1], nextY0 - y[kMaxIdx - 1]);
                if (ang0Ext > ang0N) {
                    perp0 -= (Math.PI/2.);
                } else {
                    perp0 += (Math.PI/2.);
                }
            }
        }
        
        if (perp0 >= 2*Math.PI) {
            perp0 = perp0 - 2*Math.PI;
        } else if (perp0 < 0) {
            perp0 += 2*Math.PI;
        }
        
        //---------------------------------
        int dx1 = x[kMaxIdx + 1] - x[kMaxIdx];
        int dy1 = y[kMaxIdx + 1] - y[kMaxIdx];
        
        double theta1 = AngleUtil.polarAngleCCW(dx1, dy1);
        if (theta1 != 0) {
            theta1 = (2.*Math.PI) - theta1;
        }
        
        // a point extending the line from kMaxIdx to kMaxIdx + 1 is in this quadrant
        int q1 = calculateQuadrant(x[kMaxIdx + 1] + dx1, y[kMaxIdx + 1] + dy1,
            x[kMaxIdx + 1], y[kMaxIdx + 1]);
        
        int prevX1 = x[kMaxIdx + 1] - dx1;
        int prevY1 = y[kMaxIdx + 1] + dy1;
        
        int nextX1 = x[kMaxIdx + 1] + dx1;
        int nextY1 = y[kMaxIdx + 1] - dy1;
        
        int ref1Idx = kMaxIdx - 1;
        
        double distSqP1 = (prevX1 - x[ref1Idx])*(prevX1 - x[ref1Idx])
            + (prevY1 - y[ref1Idx])*(prevY1 - y[ref1Idx]);
        
        double distSqN1 = (nextX1 - x[ref1Idx])*(nextX1 - x[ref1Idx])
            + (nextY1 - y[ref1Idx])*(nextY1 - y[ref1Idx]);
        
        if (distSqP1 == distSqN1) {
            
            while (distSqP1 == distSqN1) {
                ref1Idx++;
                if (ref1Idx > (x.length - 1)) {
                    throw new CornerRegionDegneracyException(
                    "need more neighboring points because the slopes are all vertical");
                }
                distSqP1 = (prevX1 - x[ref1Idx])*(prevX1 - x[ref1Idx])
                    + (prevY1 - y[ref1Idx])*(prevY1 - y[ref1Idx]);
                distSqN1 = (nextX1 - x[ref1Idx])*(nextX1 - x[ref1Idx])
                    + (nextY1 - y[ref1Idx])*(nextY1 - y[ref1Idx]);
            }
        }
        
        int q1Furthest = 0;
        
        furthestIsPrev = true;
        
        if (distSqP1 < distSqN1) {
            furthestIsPrev = false;
            q1Furthest = calculateQuadrant(nextX1, nextY1,
                x[kMaxIdx + 1], y[kMaxIdx + 1]);
        } else if (distSqP1 > distSqN1) {
            q1Furthest = calculateQuadrant(prevX1, prevY1,
                x[kMaxIdx + 1], y[kMaxIdx + 1]);
        }
        
        double perp1 = theta1;
        if (q1 < q1Furthest) {
            perp1 -= (Math.PI/2.);
        } else if (q1 > q1Furthest) {
            perp1 += (Math.PI/2.);
        } else {
            double ang1Ext = AngleUtil.polarAngleCCW(
                x[kMaxIdx + 1] - dx1 - x[kMaxIdx + 1], 
                y[kMaxIdx + 1] - dy1 - y[kMaxIdx + 1]);
            if (furthestIsPrev) {
                double ang1P = AngleUtil.polarAngleCCW(
                    x[kMaxIdx + 1] - prevX1, y[kMaxIdx + 1] - prevY1);
                if (ang1Ext < ang1P) {
                    perp1 -= (Math.PI/2.);
                } else {
                    perp1 += (Math.PI/2.);
                }
            } else {
                double ang1N = AngleUtil.polarAngleCCW(
                    x[kMaxIdx + 1] - nextX1, y[kMaxIdx + 1] - nextY1);
                if (ang1Ext < ang1N) {
                    perp1 -= (Math.PI/2.);
                } else {
                    perp1 += (Math.PI/2.);
                }
            }
        }
        
        if (perp1 >= 2*Math.PI) {
            perp1 = perp1 - 2*Math.PI;
        } else if (perp1 < 0) {
            perp1 += 2*Math.PI;
        }
      
        // weight by the respective adjacent curvature strengths
        double k0 = k[kMaxIdx - 1];
        double k1 = k[kMaxIdx + 1];
        double kTot = k0 + k1;
        
        double weighted = (float)(((k0/kTot)*perp0) + ((k1/kTot)*perp1));
        
        return weighted;
    }

    /**
     *  QIII | QIV
        -----|-----
        QII  | QI
     * @param xp
     * @param yp
     * @param xc
     * @param yc
     * @return 
     */
    private int calculateQuadrant(int xp, int yp, int xc, int yc) {
        
        int q = 4;
        if ((xp < xc) && (yp >= yc)) {
            q = 3;
        } else if ((xp < xc) && (yp < yc)) {
            q = 2;
        } else if ((xp >= xc) && (yp < yc)) {
            q = 1;
        }
        return q;
    }
    
    public static class CornerRegionDegneracyException extends Exception {

        public CornerRegionDegneracyException(String message) {
            super(message);
        }
        
        public CornerRegionDegneracyException(String message, Throwable cause) {
            super(message, cause);
        }
        
        public CornerRegionDegneracyException(Throwable cause) {
            super(cause);
        }
    }
}
