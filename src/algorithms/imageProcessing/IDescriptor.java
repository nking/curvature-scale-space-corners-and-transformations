package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public interface IDescriptor {

    public float calculateSSD(IDescriptor otherDesc);

    public float sumSquaredError();
    
}
