package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class PixelThetaDescriptor extends ThetaDescriptor {
    
    protected final int[] a;
    
    protected float sumSquaredError = Float.NaN;
    
    public PixelThetaDescriptor(int[] intensities) {
        this.a = intensities;
    }
}
