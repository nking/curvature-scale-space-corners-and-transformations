package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class HistogramThetaDescriptor extends ThetaDescriptor {
    
    protected final float[][] cellHistograms;
    
    protected float sumSquaredError = Float.NaN;
    
    /**
     * the index within arrays cellHistograms that the central pixel
     * value is stored in.
     */
    protected final int centralIndex;
    
    public HistogramThetaDescriptor(float[][] theCellHistograms, 
        int centralPixelIndex) {
        
        this.cellHistograms = theCellHistograms;
        this.centralIndex = centralPixelIndex;
    }
    
    @Override
    public float calculateDifference(ThetaDescriptor otherDesc) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    @Override
    public float calculateError() {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof HistogramThetaDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ThetaDescriptor");
        }
        
        HistogramThetaDescriptor other = (HistogramThetaDescriptor)otherDesc;
        
        if (this.cellHistograms.length != other.cellHistograms.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        throw new UnsupportedOperationException("not yet implemented");
    }

    @Override
    public float sumSquaredError() {
        
        if (!Float.isNaN(sumSquaredError)) {
            return sumSquaredError;
        }
        
        throw new UnsupportedOperationException("not yet implemented");
    }

    @Override
    public int getCentralIndex() {
        return centralIndex;
    }
    
}
