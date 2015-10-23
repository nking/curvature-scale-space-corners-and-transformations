package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairIntArray;
import java.util.Arrays;

/**
 * class to hold a few details from the creation of corners that
 * help later to match corners.  Note that the equals compares the
 * contents of x and y only, so this corner can be used in Collections
 * to establish equals or not, but does not have any other comparison 
 * attributes.  Note also that the equals assumption means that the user
 * has to manage possible conflicts if they place in the same 
 * Set instances from more than one image.
 *
 * @author nichole
 */
public class CornerRegion {

    protected final int edgeListIdx;

    private int idxWithinCurve = -1;
    
    protected final int kMaxIdx;

    protected final float[] k;

    protected final int[] x;

    protected final int[] y;

    protected double orientation = Double.MIN_VALUE;
    
    //TODO: a temporary work around until the classes for curvature are refactored
    // is that this instance may hold fake neighboring curvature values
    // if set by the junction finding method.
    protected boolean dummyValuesInKNeighbors = false;

    /**
     * constructor with edge index, the index with the maximum curvature in it
     * with respect to nPoints.  A minimum of 5 points is recommended and
     * for nPoints=5, a minimum k of 0.2 is needed as the maximum curvature.
     * <pre>
     *   solid angle where r = radius of curvature.  k=1/r.
     *             .
     *            /|\
     *           / | \
     *          / r-h \r
     *         /   |   \
     *        .----|----.
     *             -     bottom portion is a triangle           w
     *                                                  .----.-----.
     *                                                    .  |h .
     *                                                       .
     *   the curvature is too small to determine slopes from neighboring
     *   points when h is less than 1 pixel and w is 3 or more.
     * 
     *   limit to k for h=1.0 and w=3:
     * 
     *        r^2 = (r-h)^2 + w^2
     * 
     *        r^2 = r^2 - 2*r*h + h^2 + w^2
     *        2*r*h = h^2 + w^2
     *          r = (h^2 + w^2)/(2*h)
     *          r = 5  which is k = 0.2  
     * 
     * Therefore, for k smaller than 0.2 won't see changes in slope in the
     * neighboring 2 points on either side.
     * 
     * </pre>
     * @param theEdgeIndex
     * @param maxCurvatureIndex
     * @param nPoints
     */
    public CornerRegion(int theEdgeIndex, int nPoints, int maxCurvatureIndex) {

        this.edgeListIdx = theEdgeIndex;

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
        return edgeListIdx;
    }
    public int getKMaxIdx() {
        return kMaxIdx;
    }
    
    public void setFlagThatNeighborsHoldDummyValues() {
        dummyValuesInKNeighbors = true;
    }
    public boolean kValuesAreDummy() {
        return dummyValuesInKNeighbors;
    }

    /**
     * get the relative angle of this corner (a.k.a. dominant orientation)
     * that represents the angle perpendicular to the corner along the edge it
     * was extracted from (units are radians)
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of degrees.
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException
     * an exception is thrown if the slopes are the same for all points owned
     * by this instance, that is the points represent a line.  This can happen
     * when the radius of curvature is very large (== a small k) and any
     * change is over a larger number of pixels than present here.
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
     * an exception is thrown if the slopes are the same for all points owned
     * by this instance, that is the points represent a line.  This can happen
     * when the radius of curvature is very large (== a small k) and any
     * change is over a larger number of pixels than present here.
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
     * when the radius of curvature is very large (== a small k) and any
     * change is over a larger number of pixels than present here.
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

        //TODO: not sure will keep this weighting.  needs testing.
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

        /*
        // determine weighted average theta first, then calculate the perpendicular
        // angle pointing out from the edge at the maximum of curvature:
        
        /*needs to account for averaging when one angle is near 360 and the
        other is 0 or greater
        */

        double theta = AngleUtil.getAngleAverageInRadians(theta0, theta1);
                
        double perp = curveHelper.calculatePerpendicularAngleAwayFromCentroid(
            theta, x[kMaxIdx - 1], y[kMaxIdx - 1], x[kMaxIdx],
            y[kMaxIdx], centroidXY);

        return perp;
    }
    
    public void setIndexWithinCurve(int theIndex) {
        idxWithinCurve = theIndex;
    }
    
    public int getIndexWithinCurve() {
        return idxWithinCurve;
    }

    public CornerRegion copy() {
        
        CornerRegion cr = new CornerRegion(edgeListIdx, x.length, kMaxIdx);
        cr.orientation = orientation;
        cr.idxWithinCurve = idxWithinCurve;
        System.arraycopy(k, 0, cr.getK(), 0, k.length);
        System.arraycopy(x, 0, cr.getX(), 0, x.length);
        System.arraycopy(y, 0, cr.getY(), 0, y.length);
        
        return cr;
    }
    
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        for (int i = 0; i < x.length; ++i) {
            sb.append(String.format("k[%d]=%.2f  x,y=(%d, %d)  idx=%d\n", i, k[i],
                x[i], y[i], idxWithinCurve));
        }
        
        return sb.toString();
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof CornerRegion)) {
            return false;
        }
            
        CornerRegion other = (CornerRegion)obj;
            
        if (x.length != other.getX().length) {
            return false;
        }
        
        if (!Arrays.equals(x, other.getX())) {
            return false;
        }
     
        return Arrays.equals(y, other.getY());
    }
    
    @Override
    public int hashCode() {

        int hash = fnvHashCode();

        return hash;
    }

    int fnv321aInit = 0x811c9dc5;
    int fnv32Prime = 0x01000193;

    /**
     * hash = offset_basis
     * for each octet_of_data to be hashed
     *     hash = hash xor octet_of_data
     *     hash = hash * FNV_prime
     * return hash
     *
     * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
     */
    protected int fnvHashCode() {

        int sum = fnv321aInit;

        sum = includeInHashSum(sum, edgeListIdx);
        
        sum = includeInHashSum(sum, kMaxIdx);

        for (int xCoord : x) {
            sum = includeInHashSum(sum, xCoord);
        }
        for (int yCoord : y) {
            sum = includeInHashSum(sum, yCoord);
        }
        for (float kValue : k) {
            int kBits = Float.floatToIntBits(kValue);
            sum = includeInHashSum(sum, kBits);
        }
        
        // not including orientation which may or may not be set because it's
        // calculated only upon need
        
        return sum;
    }

    private int includeInHashSum(int sum, int variable) {
        
        // xor the bottom with the current octet.
        sum ^= variable;

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;
        
        return sum;
    }

    public static class CornerRegionDegneracyException extends Exception {

        protected static final long serialVersionUID = 456789;
        
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
