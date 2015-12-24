package algorithms.imageProcessing;

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

        /* orientation calculation is now handled in IntensityFeatures so
        this rule no longer applies.
        TODO: refactor to current usage regarding orientation.
        if (nPoints > 0) {
            if (nPoints < 3) {
                throw new IllegalArgumentException(
                "nPoints, if larger than 0, must be at least 3");
            }
            if ((maxCurvatureIndex == 0) || (maxCurvatureIndex > (nPoints - 2))) {
                throw new IllegalArgumentException(
                    "maxCurvatureIndex must be > 0 and less than nPoints-1");
            }
        }*/

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
        
        int ref0Idx = kMaxIdx - 1;
        int ref1Idx = kMaxIdx + 1;

        int dx0 = x[kMaxIdx] - x[ref0Idx];
        int dy0 = y[kMaxIdx] - y[ref0Idx];

        int dx1 = x[ref1Idx] - x[kMaxIdx];
        int dy1 = y[ref1Idx] - y[kMaxIdx];

        // extending reference points out further if gradients are same
        if ((dx1 == dx0) && (dy1 == dy0)) {

            // this is a straight line so far between points at kMaxIdx and
            // either side of it

            while ((ref0Idx > 0) || (ref1Idx < (x.length - 1))) {

                if (ref0Idx > 0) {
                    ref0Idx--;
                    dx0 = x[kMaxIdx] - x[ref0Idx];
                    dy0 = y[kMaxIdx] - y[ref0Idx];
                    dx1 = x[ref1Idx] - x[kMaxIdx];
                    dy1 = y[ref1Idx] - y[kMaxIdx];
                    if ((dx1 != dx0) || (dy1 != dy0)) {
                        break;
                    }
                }
                if (ref1Idx < (x.length - 1)) {
                    ref1Idx++;
                    dx0 = x[kMaxIdx] - x[ref0Idx];
                    dy0 = y[kMaxIdx] - y[ref0Idx];
                    dx1 = x[ref1Idx] - x[kMaxIdx];
                    dy1 = y[ref1Idx] - y[kMaxIdx];
                    if ((dx1 != dx0) || (dy1 != dy0)) {
                        break;
                    }
                }
            }            
            
            if ((dx1 == dx0) && (dy1 == dy0)) {
                throw new CornerRegionDegneracyException(
                "need more neighboring points because the slopes are all the same");
            }
        }

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double perp = curveHelper.calculateAngleTangentToMidpoint(
            x[ref0Idx], y[ref0Idx], x[kMaxIdx], y[kMaxIdx], x[ref1Idx], y[ref1Idx]);

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
    
    @Override
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
