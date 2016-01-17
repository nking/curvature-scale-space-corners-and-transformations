package algorithms.imageProcessing.features;

/**
 *
 * @author nichole
 */
public interface IntensityDescriptor extends IDescriptor {
    
    public void applyNormalization();
    
    public boolean isNormalized();
    
}
