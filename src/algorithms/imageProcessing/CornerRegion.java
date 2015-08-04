package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;

/**
 * class to hold a few details from the creation of corners that
 * help later to match corners.
 * 
 * @author nichole
 */
public class CornerRegion {
    
    protected final int edgeIndex;
    
    protected final int kMaxIndex;
    
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
        
        this.edgeIndex = theEdgeIndex;
        
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
        
        this.kMaxIndex = maxCurvatureIndex;
        
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
        
        return x[kMaxIndex - 1];
    }
    
    public int getYBeforeMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return y[kMaxIndex - 1];
    }
    
    public int getXAfterMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return x[kMaxIndex + 1];
    }
    
    public int getYAfterMax() {
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        return y[kMaxIndex + 1];
    }
    
    /**
     * get the relative angle of this corner (a.k.a. dominant orientation)
     * that represents the angle perpendicular to the corner along the edge it 
     * was extracted from (units are radians)
     * @return orientation in radians
     */
    public float getRelativeOrientationInDegrees() {
        
        double rotRadians = getRelativeOrientation();
        
        return (float)(rotRadians * 180./Math.PI);
    }
    
    /**
     * get the relative angle of this corner (a.k.a. dominant orientation)
     * that represents the angle perpendicular to the corner along the edge it 
     * was extracted from (units are radians)
     * @return orientation in radians
     */
    public double getRelativeOrientation() {
        
        if (orientation == Double.MIN_VALUE) {
            orientation = calculateOrientation();
        }
        
        return orientation;
    }
    
    /**
     * calculate the angle perpendicular to the maximum of curvature.
     * @return the angle perpendicular to the maximum of curvature at the
     * location along the edge in units of radians.
     */
    protected double calculateOrientation() {
        
        if (x.length == 0) {
            throw new IllegalStateException("this is an empty instance");
        }
        
        float dx0 = x[kMaxIndex] - x[kMaxIndex - 1];
        float dy0 = y[kMaxIndex] - y[kMaxIndex - 1];
        
        double theta0 = AngleUtil.polarAngleCCW(dx0, dy0);
        if (theta0 != 0) {
            theta0 = (2.*Math.PI) - theta0;
        }
        
        float dx1 = x[kMaxIndex + 1] - x[kMaxIndex];
        float dy1 = y[kMaxIndex + 1] - y[kMaxIndex];
        
        double theta1 = AngleUtil.polarAngleCCW(dx1, dy1);
        if (theta1 != 0) {
            theta1 = (2.*Math.PI) - theta1;
        }
        
        double perp0 = (theta0 <= Math.PI) ? (theta0 - (Math.PI/2.)) :
            (theta0 + (Math.PI/2.));
        
        if (perp0 >= 2*Math.PI) {
            perp0 = perp0 - 2*Math.PI;
        }
        
        double perp1 = (theta1 <= Math.PI) ? (theta1 - (Math.PI/2.)) :
            (theta1 + (Math.PI/2.));
        
        if (perp1 >= 2*Math.PI) {
            perp1 = perp1 - 2*Math.PI;
        }
        
        // weight by the respective adjacent curvature strengths
        double k0 = k[kMaxIndex - 1];
        double k1 = k[kMaxIndex + 1];
        double kTot = k0 + k1;
        
        return (float)(((k0/kTot)*perp0) + ((k1/kTot)*perp1));
    }
}
