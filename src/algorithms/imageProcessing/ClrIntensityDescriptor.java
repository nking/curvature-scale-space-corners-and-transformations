package algorithms.imageProcessing;

import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class ClrIntensityDescriptor implements IntensityDescriptor {
    
    //TODO: use more compact data structures after the general logic
    // is working and tested
    
    protected static int sentinel = Integer.MIN_VALUE;
    
    protected final int[] red;
    protected final int[] green;
    protected final int[] blue;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    public ClrIntensityDescriptor(int[] r, int[] g, int[] b) {
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
    
    /**
     * NOT IMPLEMENTED
     * apply a normalization to pixel values such that 
     * I[pixel] = (I[pixel] - mean(all I))/standardDeviation(all I).
     * The method invoked a second time does not change the internal values.
     */
    @Override
    public void applyNormalization() {
        
        if (hasBeenNormalized) {
            return;
        }
        
        /*
        TODO: implement this...
        
        histogram equalization at the pre-processing stage of the entire image
        can stretch the range of values over the available range, but for
        images containing a small intersection of content that might not be
        a helpful operation.
        
        corrections at the block level for illumination probably need to 
        be derived at a larger level with knowledge of the illumination
        source...
        */
        
        /*
        float[] meanAndStDevR = MiscMath.getAvgAndStDevIgnoreForSentinel(red, 
            red.length, sentinel);
        float[] meanAndStDevG = MiscMath.getAvgAndStDevIgnoreForSentinel(green, 
            green.length, sentinel);
        float[] meanAndStDevB = MiscMath.getAvgAndStDevIgnoreForSentinel(blue, 
            blue.length, sentinel);        
        */
        
        hasBeenNormalized = true;
    }

    @Override
    public boolean isNormalized() {
        return false;
        //return hasBeenNormalized;
    }
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof ClrIntensityDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ClrIntensityDescriptor");
        }
        
        ClrIntensityDescriptor other = (ClrIntensityDescriptor)otherDesc;
        
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
