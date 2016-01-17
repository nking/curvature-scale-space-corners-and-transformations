package algorithms.imageProcessing.features;

/**
 *
 * @author nichole
 */
public abstract class ThetaDescriptor implements IDescriptor {

    public static final int sentinel = Integer.MIN_VALUE;
    
    public abstract float calculateDifference(ThetaDescriptor otherDesc);
    
    public abstract float calculateError();

    public abstract float[] calculateMeanAndStDev(ThetaDescriptor descTheta2);
}
