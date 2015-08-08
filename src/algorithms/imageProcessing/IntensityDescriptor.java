package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public interface IntensityDescriptor extends IDescriptor {
    
    public void applyNormalization();
    
    public boolean isNormalized();
    
}
