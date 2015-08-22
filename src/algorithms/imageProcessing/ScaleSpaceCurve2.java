package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.Arrays;

/**
 * Holds X(t, sigma), Y(t, sigma), k(t, sigma) and the array t where t is 
 * the range of indexes for the curve in edges, normalized to values that 
 * range from 0 to 1.   Also holds information on inflection points (where
 * the curvature is zero).
 * 
 * @author nichole
 */
public class ScaleSpaceCurve {
    
    /**
     * sigma is the scale parameter used in a Gaussian kernel to convolve with
     * the curves (which is edges[edgesIdx]).
     */
    private final float sigma;

    /*
    x and y points of the curve.  the color field in xy can be used to indicate
    a closed curve.
    */
    private final PairIntArrayWithColor xy;

    /**
     * the computed curvature for each point in the curve.
     */
    private final float[] k;

    /**
     * the values of the point indexes rescaled to have values between 0 and 1.
     * For example, the point at x[10] has a t value from t[10]
     */
    private final float[] t;
    
    private final int size;

    /**
     * Holds the indexes for the points where the curvature is 0 for an edge.
     * The enclosed indexes are used with x, y, k, or t.
     */
    private int[] kIsZeroIdx;
    private int[] kIsZeroX;
    private int[] kIsZeroY;

    /**
     * the number of usable points in kIsZeroIdx. the array may be longer that
     * this number, but those values are not valid.
     */
    private int kIsZeroIdxSize;
    
    public ScaleSpaceCurve(float theSigma, PairIntArray curve, 
        boolean curveIsClosed) {

        sigma = theSigma;

        k = new float[curve.getN()];
        
        xy = (curve instanceof PairIntArrayWithColor) ? 
            (PairIntArrayWithColor) curve : new PairIntArrayWithColor(curve);
        
        if (curveIsClosed) {
            xy.setColor(1);
        }
        
        size = curve.getN();
        
        t = new float[curve.getN()];
        for (int i = 0; i < curve.getN(); i++) {
            t[i] = i / ((float) curve.getN());
        }

        // this will be reduced in size later:
        kIsZeroIdx = new int[curve.getN()];
        kIsZeroX = new int[kIsZeroIdx.length];
        kIsZeroY = new int[kIsZeroIdx.length];

        kIsZeroIdxSize = 0;
    }
    
    /**
     * @return the sigma
     */
    public float getSigma() {
        return sigma;
    }

    public float getX(int idx) {
        if (idx < 0 || idx > (xy.getN() - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of xy");
        }
        return xy.getX(idx);
    }
    
    public float getY(int idx) {
        if (idx < 0 || idx > (xy.getN() - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of xy");
        }
        return xy.getY(idx);
    }
    
    public void setXY(int idx, int xValue, int yValue) {
        if (idx < 0 || idx > (xy.getN() - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of xy");
        }
        xy.set(idx, xValue, yValue);
    }

    public float getK(int idx) {
        if (idx < 0 || idx > (k.length - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of array k");
        }
        return k[idx];
    }
    
    public void setK(int idx, float value) {
        if (idx < 0 || idx > (k.length - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of array k");
        }
        k[idx] = value;
    }
    
    /**
     * @return the t
     */
    public float[] getT() {
        return t;
    }

    /**
     * @return the kIsZeroIdx
     */
    public int[] getKIsZeroIdx() {
        return kIsZeroIdx;
    }
    
    /**
     * return the x coordinates where the curvature is zero.
     * @return 
     */
    public int[] getKIsZeroX() {
        return kIsZeroX;
    }
    
    /**
     * return the y coordinates where the curvature is zero.
     * @return 
     */
    public int[] getKIsZeroY() {
        return kIsZeroY;
    }

    /**
     * @return the kIsZeroIdxSize
     */
    public int getKIsZeroIdxSize() {
        return kIsZeroIdxSize;
    }
    
    public void addKIsZeroIdx(int idxForAZeroValue, int xCoord, int yCoord) {
        
        kIsZeroIdx[kIsZeroIdxSize] = idxForAZeroValue;
        kIsZeroX[kIsZeroIdxSize] = xCoord;
        kIsZeroY[kIsZeroIdxSize] = yCoord;
        
        kIsZeroIdxSize++;
    }
    
    public void compressKIsZeroIdx() {
        kIsZeroIdx = Arrays.copyOf(kIsZeroIdx, kIsZeroIdxSize);
        kIsZeroX = Arrays.copyOf(kIsZeroX, kIsZeroIdxSize);
        kIsZeroY = Arrays.copyOf(kIsZeroY, kIsZeroIdxSize);
    }
    
    public int getSize() {
        return size;
    }
    
    public float[] getK() {
        return k;
    }
    
    public int[] getX() {
        return xy.getX();
    }
    
    public int[] getY() {
        return xy.getY();
    }
    
    public boolean curveIsClosed() {
        return (xy.getColor() == 1);
    }
    
    public PairIntArray getXYCurve() {
        return xy;
    }
}
