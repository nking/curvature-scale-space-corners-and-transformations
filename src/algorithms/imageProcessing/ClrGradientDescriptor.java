package algorithms.imageProcessing;

import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class ClrGradientDescriptor implements GradientDescriptor {
    
    //TODO: use more compact data structures after the general logic
    // is working and tested
    
    protected static int sentinel = Integer.MIN_VALUE;
    
    protected final int[] red;
    protected final int[] green;
    protected final int[] blue;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    public ClrGradientDescriptor(int[] r, int[] g, int[] b) {
        if (r == null) {
            throw new IllegalArgumentException("r cannot be null");
        }
        if (g == null) {
            throw new IllegalArgumentException("g cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (r.length != g.length || r.length != b.length) {
            throw new IllegalArgumentException("r, g, and b must be same length");
        }
        this.red = r;
        this.green = g;
        this.blue = b;
    }
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof ClrGradientDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ClrGradientDescriptor");
        }
        
        ClrGradientDescriptor other = (ClrGradientDescriptor)otherDesc;
        
        if (this.red.length != other.red.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
                
        float ssdR = MiscMath.calculateSSD(red, other.red, sentinel);
        float ssdG = MiscMath.calculateSSD(green, other.green, sentinel);
        float ssdB = MiscMath.calculateSSD(blue, other.blue, sentinel);
           
        float avg = (float)(ssdR + ssdG + ssdB)/3.f;
        
        return avg;
    }
    
    /**
     * Determine the sum squared error within this descriptor using 
     * auto-correlation and the assumption that the value at the middle index 
     * is the value from the original central pixel.
     * Note that the value is persisted after one calculation.  Any use
     * of normalization should happen before this is first invoked.
     * (A function with a force calculation argument can be made if necessary though).
     * @return 
     */
    @Override
    public float sumSquaredError() {
        
        if (!Float.isNaN(sumSquaredError)) {
            return sumSquaredError;
        }
        
        int n = red.length;
        int midIdx = n >> 1;
        
        int rVC = red[midIdx];
        int gVC = green[midIdx];
        int bVC = blue[midIdx];
        
        if (rVC == sentinel || gVC == sentinel || bVC == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central values for the array are somehow sentinels");
        }
        
        float sqErrR = MiscMath.sumSquaredError(red, sentinel);
        float sqErrG = MiscMath.sumSquaredError(green, sentinel);
        float sqErrB = MiscMath.sumSquaredError(blue, sentinel);
        
        float avg = (float)(sqErrR + sqErrG + sqErrB)/3.f;
        
        this.sumSquaredError = avg;
        
        return sumSquaredError;
    }
}
