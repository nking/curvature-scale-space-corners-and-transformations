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
    
    @Override
    public float calculateSSD(IntensityDescriptor otherDesc) {
        
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
                
        int count = 0;
        double sumR = 0;
        double sumG = 0;
        double sumB = 0;
        
        for (int i = 0; i < this.red.length; ++i) {
            
            int r1 = red[i];
            if (r1 == sentinel) {
                continue;
            }
            int r2 = other.red[i];
            if (r2 == sentinel) {
                continue;
            }
            
            int g1 = green[i];
            int b1 = blue[i];
            
            int g2 = other.green[i];
            int b2 = other.blue[i];
            
            sumR += (r1 - r2)*(r1 - r2);
            sumG += (g1 - g2)*(g1 - g2);
            sumB += (b1 - b2)*(b1 - b2);
            
            count++;
        }
        sumR /= (double)count;
        sumG /= (double)count;
        sumB /= (double)count;
        
        float avg = (float)(sumR + sumG + sumB)/3.f;
        
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
        
        int count = 0;
        
        double sumR = 0;
        double sumG = 0;
        double sumB = 0;
        
        for (int i = 0; i < this.red.length; ++i) {
            
            int r1 = red[i];
            if (r1 == sentinel) {
                continue;
            }
            
            int diffR = red[i] - rVC;
            int diffG = green[i] - gVC;
            int diffB = blue[i] - bVC;
        
            sumR += (diffR * diffR);
            sumG += (diffG * diffG);
            sumB += (diffB * diffB);
            count++;
        }
        sumR /= (double)count;
        sumG /= (double)count;
        sumB /= (double)count;
        
        float avg = (float)(sumR + sumG + sumB)/3.f;
        
        this.sumSquaredError = avg;
        
        return sumSquaredError;
    }

}
