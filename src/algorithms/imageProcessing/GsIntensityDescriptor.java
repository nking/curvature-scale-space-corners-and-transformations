package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class GsIntensityDescriptor implements IntensityDescriptor {

    public static float sentinel = Float.MIN_VALUE;
    
    protected final float[] grey;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    /**
     * the index within array grey that the central pixel
     * value is stored in.
     */
    protected final int centralIndex;
    
    public GsIntensityDescriptor(float[] intensities, int centralPixelIndex) {
        this.grey = intensities;
        this.centralIndex = centralPixelIndex;
    }
    
    /**
     * NOT YET IMPLEMENTED
     * apply a normalization to pixel values to try to reduce the differences 
     * due to images of the same region due to lighting or perspective 
     * for example.
     * The method invoked a second time does not change the internal values.
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
        
        if (!(otherDesc instanceof GsIntensityDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type GsIntensityDescriptor");
        }
        
        GsIntensityDescriptor other = (GsIntensityDescriptor)otherDesc;
        
        if (this.grey.length != other.grey.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        float ssd = MiscMath.calculateSSD(grey, other.grey, sentinel);
                
        return ssd;
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
        
        float vc = grey[centralIndex];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        float sqErr = MiscMath.sumSquaredError(grey, sentinel, centralIndex);
            
        this.sumSquaredError = sqErr;
        
        return sumSquaredError;
    }

    protected float[] getInternalArrayCopy() {
        return Arrays.copyOf(grey, grey.length);
    }

    @Override
    public int getCentralIndex() {
        return centralIndex;
    }

    @Override
    public String toString() {
        return Arrays.toString(grey);
    }
    
}
