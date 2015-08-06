package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public interface IntensityDescriptor {
    
    public void applyNormalization();
    
    public boolean isNormalized();
    
    public float calculateSSD(IntensityDescriptor otherDesc);
    
    public float sumSquaredError();
    
}
