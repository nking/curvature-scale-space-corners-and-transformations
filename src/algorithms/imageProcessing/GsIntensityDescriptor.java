package algorithms.imageProcessing;

import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class GsIntensityDescriptor implements IntensityDescriptor {

    protected static int sentinel = Integer.MIN_VALUE;
    
    protected final int[] grey;
    
    protected float sumSquaredError = Float.NaN;
    
    protected boolean hasBeenNormalized = false;
    
    public GsIntensityDescriptor(int[] intensities) {
        this.grey = intensities;
    }
    
    /**
     * NOT YET TESTED
     * apply a normalization to pixel values such that 
     * I[pixel] = (I[pixel] - mean(all I))/standardDeviation(all I).
     * The method invoked a second time does not change the internal values.
     */
    @Override
    public void applyNormalization() {
        
        if (hasBeenNormalized) {
            return;
        }
        
        float[] meanAndStDev = MiscMath.getAvgAndStDevIgnoreForSentinel(grey, 
            grey.length, sentinel);
        
        for (int i = 0; i < grey.length; ++i) {
            
            if (grey[i] == sentinel) {
                continue;
            }
            
            grey[i] -= meanAndStDev[0];
            grey[i] /= meanAndStDev[1];
        }
        
        hasBeenNormalized = true;
    }

    @Override
    public boolean isNormalized() {
        return hasBeenNormalized;
    }
    
    @Override
    public float calculateSSD(IntensityDescriptor otherDesc) {
        
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
                
        int count = 0;
        double sum = 0;
        
        for (int i = 0; i < this.grey.length; ++i) {
            
            int v1 = grey[i];
            if (v1 == sentinel) {
                continue;
            }
            int v2 = other.grey[i];
            if (v2 == sentinel) {
                continue;
            }
            
            sum += (v1 - v2)*(v1 - v2);
            
            count++;
        }
        sum /= (double)count;
                
        return (float)sum;
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
        
        int count = 0;
        
        double sum = 0;
        
        for (int i = 0; i < this.grey.length; ++i) {
            
            int v1 = grey[i];
            if (v1 == sentinel) {
                continue;
            }
            
            int diff = grey[i] - vc;
        
            sum += (diff * diff);
            count++;
        }
        sum /= (double)count;
               
        this.sumSquaredError = (float)sum;
        
        return sumSquaredError;
    }

}
