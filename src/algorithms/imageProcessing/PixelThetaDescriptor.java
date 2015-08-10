package algorithms.imageProcessing;

import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class PixelThetaDescriptor extends ThetaDescriptor {
    
    protected final int[] a;
    
    protected float sumSquaredError = Float.NaN;
    
    /**
     * the index within array a that the central pixel
     * value is stored in.
     */
    protected final int centralIndex;
    
    public PixelThetaDescriptor(int[] intensities, int centralPixelIndex) {
        this.a = intensities;
        this.centralIndex = centralPixelIndex;
    }
    
    @Override
    public float calculateDifference(ThetaDescriptor otherDesc) {
        return calculateSSD(otherDesc);
    }
    
    @Override
    public float calculateError() {
        return sumSquaredError();
    }
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof PixelThetaDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ThetaDescriptor");
        }
        
        PixelThetaDescriptor other = (PixelThetaDescriptor)otherDesc;
        
        if (this.a.length != other.a.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        float ssd = MiscMath.calculateAngular360SSD(a, other.a, sentinel);
                
        return ssd;
    }

    @Override
    public float sumSquaredError() {
        
        if (!Float.isNaN(sumSquaredError)) {
            return sumSquaredError;
        }
                
        int vc = a[centralIndex];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        float sqErr = MiscMath.sumSquaredAngular360Error(a, sentinel, centralIndex);
            
        this.sumSquaredError = sqErr;
        
        return sumSquaredError;
    }

    @Override
    public int getCentralIndex() {
        return centralIndex;
    }
    
}
