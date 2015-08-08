package algorithms.imageProcessing;

import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class GsGradientDescriptor implements GradientDescriptor {
    
    protected static int sentinel = Integer.MIN_VALUE;
    
    protected final int[] grey;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    public GsGradientDescriptor(int[] binnedCells) {
        this.grey = binnedCells;
    }
    
    // TODO: may change to use histograms or another technique
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof GsGradientDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type GsGradientDescriptor");
        }
        
        GsGradientDescriptor other = (GsGradientDescriptor)otherDesc;
        
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
        
        int n = grey.length;
        int midIdx = n >> 1;
        
        int vc = grey[midIdx];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        float sqErr = MiscMath.sumSquaredError(grey, sentinel);
            
        this.sumSquaredError = sqErr;
        
        return sumSquaredError;
    }

}
