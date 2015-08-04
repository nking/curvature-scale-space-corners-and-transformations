package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairIntArray;

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
     * with respect to nPoints.  A minimum of 5 points is recommended and
     * large maximum curvature, k[kMaxIdx].
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
    public int getEdgeIdx() {
        return edgeIdx;
    }
    public int getKMaxIdx() {
        return kMaxIdx;
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
     * points is present (slopes cannot both be the same on both sides of the
     * maximum for all the points given) else an exception is thrown.
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of radians.
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException
     * an exception is thrown if the slopes are the same for all points owned 
     * by this instance, that is the points represent a line.  This can happen
     * when the radius of curvature is very large, but that produces a small
     * curvature value (k=1/R).
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

        int dx1 = x[kMaxIdx + 1] - x[kMaxIdx];
        int dy1 = y[kMaxIdx + 1] - y[kMaxIdx];

        double theta1 = AngleUtil.polarAngleCCW(dx1, dy1);
        if (theta1 != 0) {
            theta1 = (2.*Math.PI) - theta1;
        }

        /* determine whether to add or subtract 90 for each vector:

        the centroid of the middle and neighboring points determines the
        opposite side of the vector pointing outwards from the edge.
        */

        PairIntArray xy = new PairIntArray(3);
        xy.add(x[kMaxIdx - 1], y[kMaxIdx - 1]);
        xy.add(x[kMaxIdx], y[kMaxIdx]);
        xy.add(x[kMaxIdx + 1], y[kMaxIdx + 1]);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[] centroidXY = curveHelper.calculateXYCentroids(xy);

        if ((dx1 == dx0) && (dy1 == dy0)) {

            // this is a straight line so far between points at kMaxIdx and
            // either side of it

            int ref0Idx = kMaxIdx - 1;
            int ref1Idx = kMaxIdx + 1;
            ref0Idx--;
            ref1Idx++;

            if ((ref0Idx < 0) && (ref1Idx > (x.length - 1))) {
                throw new CornerRegionDegneracyException(
                "need more neighboring points because the slopes are all the same");
            }

            while ((ref0Idx > -1) || (ref1Idx < x.length)) {

                if (ref0Idx > -1) {
                    int dxRef0 = x[kMaxIdx - 1] - x[ref0Idx];
                    int dyRef0 = y[kMaxIdx - 1] - y[ref0Idx];
                    if ((dxRef0 != dx0) || (dyRef0 != dy0)) {
                        xy.add(x[ref0Idx], y[ref0Idx]);
                        centroidXY = curveHelper.calculateXYCentroids(xy);
                        break;
                    }
                    ref0Idx--;
                }
                if (ref1Idx < x.length) {
                    int dxRef1 = x[ref1Idx] - x[kMaxIdx + 1];
                    int dyRef1 = y[ref1Idx] - y[kMaxIdx + 1];
                    if ((dxRef1 != dx1) || (dyRef1 != dy1)) {
                        xy.add(x[ref1Idx], y[ref1Idx]);
                        centroidXY = curveHelper.calculateXYCentroids(xy);
                        break;
                    }
                    ref1Idx++;
                }

                if ((ref0Idx < 0) && (ref1Idx > (x.length - 1))) {
                    throw new CornerRegionDegneracyException(
                    "need more neighboring points because the slopes are all the same");
                }
            }
        }

        /*
        The angle vectors for theta0 and theta1 should point away from the
        centroid.
        */
        
        //TODO: not sure will keep this weighting.  needs testing
        // weight by the respective adjacent curvature strengths
        double k0 = k[kMaxIdx - 1];
        double k1 = k[kMaxIdx + 1];
        double kTot = k0 + k1;

        /*
        // determine both angles separately, then calc weighted average:
        double perp0 = calculatePerpendicularAngleAwayFromCentroid(theta0,
            x[kMaxIdx - 1], y[kMaxIdx - 1], x[kMaxIdx], y[kMaxIdx], centroidXY);

        double perp1 = calculatePerpendicularAngleAwayFromCentroid(theta1,
            x[kMaxIdx], y[kMaxIdx], x[kMaxIdx + 1], y[kMaxIdx + 1], centroidXY);
       
        double weighted = (float)(((k0/kTot)*perp0) + ((k1/kTot)*perp1));
        */
        
        // determine weighted average theta first, then calculate the perpendicular
        // angle pointing out from the edge at the maximum of curvature:
        double weightedTheta = (float)(((k0/kTot)*theta0) + ((k1/kTot)*theta1));
        
        double weighted = calculatePerpendicularAngleAwayFromCentroid(
            weightedTheta, x[kMaxIdx - 1], y[kMaxIdx - 1], x[kMaxIdx], 
            y[kMaxIdx], centroidXY);
        
        return weighted;
    }
    
    /**
     * given theta and the point (xp, yp), determine which direction and hence
     * polar angle (clockwise) is perpendicular away from the centroid.
     * The reference point (xm, ym) is the point from which theta was also
     * calculated, which is probably the point for kMaxIdx.
     * @param theta
     * @param xp
     * @param yp
     * @param xm
     * @param ym
     * @param centroidXY
     * @return
     */
    protected double calculatePerpendicularAngleAwayFromCentroid(
        double theta, int xp, int yp, int xm, int ym, double[] centroidXY) {

        /*
        rotate the point (xm, ym) around (xp, yp) 90 degrees and -90 degrees.
        The rotated point which is furthest from the centroid is the
        direction of the vector pointing away from the centroid.
        */

        /*
        math.cos(math.pi/2) = 0
        math.sin(math.pi/2) = 1
        math.sin(-math.pi/2) = -1

        double xr = centroidX + ((y - centroidY) * sine(angle)));
        double yr = centroidY + ((-(x - centroidX) * sine(angle)))
        */

        int xmRot90 = xp + (ym - yp);
        int ymRot90 = yp + (-(xm - xp));

        int xmRotNegative90 = xp  - (ym - yp);
        int ymRotNegative90 = yp + (xm - xp);

        double distSqRot90 = (xmRot90 - centroidXY[0]) * (xmRot90 - centroidXY[0])
            + (ymRot90 - centroidXY[1]) * (ymRot90 - centroidXY[1]);

        double distSqRotNegative90 =
            (xmRotNegative90 - centroidXY[0]) * (xmRotNegative90 - centroidXY[0])
            + (ymRotNegative90 - centroidXY[1]) * (ymRotNegative90 - centroidXY[1]);

        double perp = theta;

        if (distSqRot90 > distSqRotNegative90) {
            perp += Math.PI/2.;
        } else {
            perp -= Math.PI/2.;
        }

        if (perp >= 2*Math.PI) {
            perp = perp - 2*Math.PI;
        } else if (perp < 0) {
            perp += 2*Math.PI;
        }

        return perp;
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
