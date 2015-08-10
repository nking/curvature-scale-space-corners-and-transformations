package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class HistogramThetaDescriptor extends ThetaDescriptor {
    
    protected final float[][] cellHistograms;
    
    protected float sumSquaredError = Float.NaN;
    
    public HistogramThetaDescriptor(float[][] theCellHistograms) {
        this.cellHistograms = theCellHistograms;
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
    
}
