package algorithms.imageProcessing.features;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.Arrays;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class BlobPerimeterRegion extends CornerRegion {
    
    /**
     * this is released after the orientation is calculated.
    */
    protected Set<PairInt> blob = null;
        
    public BlobPerimeterRegion(final int theEdgeIndex, final int prevXCoord, 
        final int prevYCoord, final int xCoord, final int yCoord, 
        final int nextXCoord, final int nextYCoord,
        Set<PairInt> theBlob, final float sigma) {
        
        super(theEdgeIndex, 3, 1);

        dummyValuesInKNeighbors = true;
        
        float dummyKValue = -99;
        
        set(0, dummyKValue, prevXCoord, prevYCoord);
        
        set(1, sigma, xCoord, yCoord);
        
        set(2, dummyKValue, nextXCoord, nextYCoord);
        
        blob = theBlob;
    }
    
    public BlobPerimeterRegion(final int theEdgeIndex, final int prevXCoord, 
        final int prevYCoord, final int xCoord, final int yCoord, 
        final int nextXCoord, final int nextYCoord,
        double theOrientation, final float sigma) {
        
        super(theEdgeIndex, 3, 1);

        dummyValuesInKNeighbors = true;
        
        float dummyKValue = -99;
        
        set(0, dummyKValue, prevXCoord, prevYCoord);
        
        set(1, sigma, xCoord, yCoord);
        
        set(2, dummyKValue, nextXCoord, nextYCoord);
        
        orientation = theOrientation;
    }
    
    public void overrideRelativeOrientation(final double thetaInRadians) {
        
        this.orientation = thetaInRadians;
        
        this.blob = null;
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
    @Override
    protected double calculateOrientation() {
        
        int dx0 = x[1] - x[0];
        int dy0 = y[1] - y[0];

        double theta0 = AngleUtil.polarAngleCCW(dx0, dy0);
        if (theta0 != 0) {
            theta0 = (2.*Math.PI) - theta0;
        }

        int dx1 = x[2] - x[1];
        int dy1 = y[2] - y[1];

        double theta1 = AngleUtil.polarAngleCCW(dx1, dy1);
        if (theta1 != 0) {
            theta1 = (2.*Math.PI) - theta1;
        }

        /* determine whether to add or subtract 90 for each vector:

        the centroid of the middle and neighboring points determines the
        opposite side of the vector pointing outwards from the edge.
        */

        PairIntArray xy = new PairIntArray(3);
        xy.add(x[0], y[0]);
        xy.add(x[1], y[1]);
        xy.add(x[2], y[2]);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[] centroidXY = curveHelper.calculateXYCentroids(xy);

        // may need special handling for ((dx1 == dx0) && (dy1 == dy0))

        // determine weighted average theta first, then calculate the perpendicular
        // angle pointing out from the edge at the maximum of curvature:
        
        double theta = AngleUtil.getAngleAverageInRadians(theta0, theta1);
                
        double perp = curveHelper.calculatePerpendicularAngleAwayFromCentroid(
            theta, x[0], y[0], x[1], y[1], centroidXY, blob);

        return perp;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof BlobPerimeterRegion)) {
            return false;
        }
            
        BlobPerimeterRegion other = (BlobPerimeterRegion)obj;
            
        if (x.length != other.getX().length) {
            return false;
        }
        
        if (!Arrays.equals(x, other.getX())) {
            return false;
        }
     
        return Arrays.equals(y, other.getY());
    }
    
    public float getSigma() {
        return k[1];
    }
    
    @Override
    public BlobPerimeterRegion copy() {
        
        if (blob == null) {
            
            BlobPerimeterRegion c = new BlobPerimeterRegion(edgeListIdx, 
                x[0], y[0], x[1], y[1], x[2], y[2], orientation, k[1]);
            c.setIndexWithinCurve(this.getIndexWithinCurve());
            if (dummyValuesInKNeighbors) {
                c.setFlagThatNeighborsHoldDummyValues();
            }
            return c;
        }
        
        BlobPerimeterRegion c = new BlobPerimeterRegion(edgeListIdx, 
            x[0], y[0], x[1], y[1], x[2], y[2], blob, k[1]);
        c.setIndexWithinCurve(this.getIndexWithinCurve());
        if (dummyValuesInKNeighbors) {
            c.setFlagThatNeighborsHoldDummyValues();
        }
        
        return c;
    }
}
