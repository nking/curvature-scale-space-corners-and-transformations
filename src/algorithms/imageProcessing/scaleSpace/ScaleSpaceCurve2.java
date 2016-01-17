package algorithms.imageProcessing.scaleSpace;

import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.Arrays;

/**
 * extends ScaleSpaceCurve to hold the gaussian 2nd derivatives in x, and y
 * 
 * @author nichole
 */
public class ScaleSpaceCurve2 extends ScaleSpaceCurve {
        
    private final float[] x2ndDeriv;
    private final float[] y2ndDeriv;

    /**
     * the number of usable points in kIsZeroIdx. the array may be longer that
     * this number, but those values are not valid.
     */
    private int kIsZeroIdxSize;
    
    public ScaleSpaceCurve2(float theSigma, PairIntArray curve, 
        boolean curveIsClosed) {

        super(theSigma, curve, curveIsClosed);
        
        x2ndDeriv = new float[curve.getN()];
        y2ndDeriv = new float[curve.getN()];
    }
 
    public void setX2ndDeriv(int idx, float value) {
        if (idx < 0 || idx > (x2ndDeriv.length - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of array x2ndDeriv");
        }
        x2ndDeriv[idx] = value;
    }
    public void setY2ndDeriv(int idx, float value) {
        if (idx < 0 || idx > (y2ndDeriv.length - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of array y2ndDeriv");
        }
        y2ndDeriv[idx] = value;
    }
    public float getX2ndDeriv(int idx) {
        if (idx < 0 || idx > (x2ndDeriv.length - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of array x2ndDeriv");
        }
        return x2ndDeriv[idx];
    }
    public float getY2ndDeriv(int idx) {
        if (idx < 0 || idx > (y2ndDeriv.length - 1)) {
            throw new IllegalArgumentException("idx is out of bounds of array y2ndDeriv");
        }
        return y2ndDeriv[idx];
    }
    
}
