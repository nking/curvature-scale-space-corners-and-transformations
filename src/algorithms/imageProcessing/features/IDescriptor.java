package algorithms.imageProcessing.features;

/**
 *
 * @author nichole
 */
public interface IDescriptor {

    public int getCentralIndex();
    
    public float calculateSSD(IDescriptor otherDesc);

    public float sumSquaredError();
    
    public float calculateCosineSimilarity(IDescriptor otherDesc);
}
