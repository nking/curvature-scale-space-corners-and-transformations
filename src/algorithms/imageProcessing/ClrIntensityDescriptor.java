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
     * apply a normalization to pixel values such that 
     * I[pixel] = (I[pixel] - mean(all I))/standardDeviation(all I).
     * The method invoked a second time does not change the internal values.
     */
    @Override
    public void applyNormalization() {
        
        if (hasBeenNormalized) {
            return;
        }
        
        float[] meanAndStDevR = MiscMath.getAvgAndStDevIgnoreForSentinel(red, 
            red.length, sentinel);
        float[] meanAndStDevG = MiscMath.getAvgAndStDevIgnoreForSentinel(green, 
            green.length, sentinel);
        float[] meanAndStDevB = MiscMath.getAvgAndStDevIgnoreForSentinel(blue, 
            blue.length, sentinel);
        
        for (int i = 0; i < red.length; ++i) {
            
            if (red[i] == sentinel) {
                continue;
            }
            red[i] -= meanAndStDevR[0];
            red[i] /= meanAndStDevR[1];
            
            green[i] -= meanAndStDevG[0];
            green[i] /= meanAndStDevG[1];
            
            blue[i] -= meanAndStDevB[0];
            blue[i] /= meanAndStDevB[1];
            
        }
        
        hasBeenNormalized = true;
    }

    @Override
    public boolean isNormalized() {
        return hasBeenNormalized;
    }
}
