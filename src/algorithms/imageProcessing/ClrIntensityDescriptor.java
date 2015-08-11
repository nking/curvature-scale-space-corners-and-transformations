package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class ClrIntensityDescriptor implements IntensityDescriptor {
    
    //TODO: use more compact data structures after the general logic
    // is working and tested
    
    public static int sentinel = Integer.MIN_VALUE;
    
    protected final int[] red;
    protected final int[] green;
    protected final int[] blue;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    /**
     * the index within arrays red, green and blue that the central pixel
     * value is stored in.
     */
    protected final int centralIndex;
    
    public ClrIntensityDescriptor(int[] r, int[] g, int[] b, int centralPixelIndex) {
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
        this.centralIndex = centralPixelIndex;
    }
    
    /**
     * NOT IMPLEMENTED
     * apply a normalization to pixel values 
     */
    @Override
    public void applyNormalization() {
        
        if (hasBeenNormalized) {
            return;
        }
        
        /*        
        histogram equalization at the pre-processing stage of the entire image
        can stretch the range of values over the available range, but for
        images containing a small intersection of content that might not be
        a helpful operation.
        
        corrections at the block level for illumination probably need to 
        be derived at a larger level with knowledge of the illumination
        source...
        */
       
        throw new UnsupportedOperationException("not yet implemented");
       
        //hasBeenNormalized = true;
    }

    @Override
    public boolean isNormalized() {
        return hasBeenNormalized;
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
           
        float avg = (ssdR + ssdG + ssdB)/3.f;
        
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
                
        int rVC = red[centralIndex];
        int gVC = green[centralIndex];
        int bVC = blue[centralIndex];
        
        if (rVC == sentinel || gVC == sentinel || bVC == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central values for the array are somehow sentinels");
        }
        
        float sqErrR = MiscMath.sumSquaredError(red, sentinel, centralIndex);
        float sqErrG = MiscMath.sumSquaredError(green, sentinel, centralIndex);
        float sqErrB = MiscMath.sumSquaredError(blue, sentinel, centralIndex);
        
        float avg = (sqErrR + sqErrG + sqErrB)/3.f;
        
        this.sumSquaredError = avg;
        
        return sumSquaredError;
    }

    @Override
    public int getCentralIndex() {
        return centralIndex;
    }

    protected int[] getInternalRedArrayCopy() {
        return Arrays.copyOf(red, red.length);
    }

    @Override
    public String toString() {
        String str1 = "red=" + Arrays.toString(red);
        String str2 = "\ngreen=" + Arrays.toString(green);
        String str3 = "\nblue=" + Arrays.toString(blue);
        return str1 + str2 + str3;
    }
    
}
