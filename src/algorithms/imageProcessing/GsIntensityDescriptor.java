package algorithms.imageProcessing;

import static algorithms.imageProcessing.ClrIntensityDescriptor.sentinel;
import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class GsIntensityDescriptor implements IntensityDescriptor {

    protected static int sentinel = Integer.MIN_VALUE;
    
    protected final int[] grey;
    
    protected boolean hasBeenNormalized = false;
    
    public GsIntensityDescriptor(int[] intensities) {
        this.grey = intensities;
    }
    
    /**
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
    
}
