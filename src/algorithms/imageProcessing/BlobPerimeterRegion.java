package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class BlobPerimeterRegion {
    
    protected final int edgeIdx;

    private final int prevX;
    
    private final int x;
    
    private final int nextX;

    private final int prevY;
    
    private final int y;

    private final int nextY;
    
    protected final Set<PairInt> blob;
    
    protected double orientation = Double.MIN_VALUE;
    
    public BlobPerimeterRegion(final int theEdgeIndex, final int prevXCoord, 
        final int prevYCoord, final int xCoord, final int yCoord, 
        final int nextXCoord, final int nextYCoord,
        Set<PairInt> theBlob) {
        
        edgeIdx = theEdgeIndex;
        
        prevX = prevXCoord;
        
        prevY = prevYCoord;
        
        x = xCoord;
        
        y = yCoord;
        
        nextX = nextXCoord;
        
        nextY = nextYCoord;
        
        blob = theBlob;
    }
    
    /**
     * get the angle tangent to the point (a.k.a. dominant orientation)
     * (units are radians).
     * @return the angle perpendicular to the point location along the perimeter
     * of the blob in units of degrees.
     */
    public float getRelativeOrientationInDegrees() {

        double rotRadians = getRelativeOrientation();

        return (float)(rotRadians * 180./Math.PI);
    }
    
    /**
     * get the relative angle of this corner (a.k.a. dominant orientation)
     * that represents the angle perpendicular to the corner along the edge it
     * was extracted from (units are radians)
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of radians.
     */
    public double getRelativeOrientation() {

        if (orientation == Double.MIN_VALUE) {
            orientation = calculateOrientation();
        }

        return orientation;
    }
    
    public void overrideRelativeOrientation(final double thetaInRadians) {
        
        this.orientation = thetaInRadians;
    }
    
    /**
     * NOT READY FOR USE YET.  NEEDs MORE TESTING
     * calculate the angle perpendicular to the maximum of curvature.
     * The curvature has to be large enough so that a change in the neighboring
     * points is present (slopes cannot both be the same on both sides of the
     * maximum for all the points given) else an exception is thrown.
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of radians.
     */
    protected double calculateOrientation() {
        
        int dx0 = x - prevX;
        int dy0 = y - prevY;

        double theta0 = AngleUtil.polarAngleCCW(dx0, dy0);
        if (theta0 != 0) {
            theta0 = (2.*Math.PI) - theta0;
        }

        int dx1 = nextX - x;
        int dy1 = nextY - y;

        double theta1 = AngleUtil.polarAngleCCW(dx1, dy1);
        if (theta1 != 0) {
            theta1 = (2.*Math.PI) - theta1;
        }

        /* determine whether to add or subtract 90 for each vector:

        the centroid of the middle and neighboring points determines the
        opposite side of the vector pointing outwards from the edge.
        */

        PairIntArray xy = new PairIntArray(3);
        xy.add(prevX, prevY);
        xy.add(x, y);
        xy.add(nextX, nextY);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[] centroidXY = curveHelper.calculateXYCentroids(xy);

        // may need special handling for ((dx1 == dx0) && (dy1 == dy0))

        // determine weighted average theta first, then calculate the perpendicular
        // angle pointing out from the edge at the maximum of curvature:
        
        double theta = AngleUtil.getAngleAverageInRadians(theta0, theta1);
                
        double perp = curveHelper.calculatePerpendicularAngleAwayFromCentroid(
            theta, prevX, prevY, x, y, centroidXY, blob);

        return perp;
    }

    /**
     * @return the prevX
     */
    public int getPrevX() {
        return prevX;
    }

    /**
     * @return the x
     */
    public int getX() {
        return x;
    }

    /**
     * @return the nextX
     */
    public int getNextX() {
        return nextX;
    }

    /**
     * @return the prevY
     */
    public int getPrevY() {
        return prevY;
    }

    /**
     * @return the y
     */
    public int getY() {
        return y;
    }

    /**
     * @return the nextY
     */
    public int getNextY() {
        return nextY;
    }

}
