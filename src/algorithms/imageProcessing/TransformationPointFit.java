package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPointFit {
    
    private TransformationParameters parameters = null;
    
    private int nMatchedPoints = 0;
    
    private double meanDistFromModel = Double.MAX_VALUE;
    
    private double stDevFromMean = Double.MAX_VALUE;
    
    private double tolerance = Double.MAX_VALUE;
    
    public TransformationPointFit(TransformationParameters theParameters, 
        int numberOfMatchedPoints, double theMeanDistFromModel, 
        double theStDevFromMean, double theTolerance) {
        
        parameters = theParameters;
        
        nMatchedPoints = numberOfMatchedPoints;
        
        meanDistFromModel = theMeanDistFromModel;
        
        stDevFromMean = theStDevFromMean;
        
        tolerance = theTolerance;
    }

    /**
     * @return the parameters
     */
    public TransformationParameters getParameters() {
        return parameters;
    }

    
    public int getNumberOfMatchedPoints() {
        return nMatchedPoints;
    }
    
    public double getRotationInRadians() {
        if (parameters == null) {
            return 0;
        }
        
        return (parameters.getRotationInDegrees() * Math.PI/180.);
    }
    
    public double getScale() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getScale();
    }
    
    public double getTranslationX() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getTranslationX();
    }
    
    public double getTranslationY() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getTranslationY();
    }
    
    /**
     * tolerance used when including only residuals whose absolute value
     * is less than tolerance.
     * @return 
     */
    public double getTolerance() {
        return tolerance;
    }

    /**
     * @return the meanDistFromModel
     */
    public double getMeanDistFromModel() {
        return meanDistFromModel;
    }

    /**
     * @return the stDevFromMean
     */
    public double getStDevFromMean() {
        return stDevFromMean;
    }
   
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("nMatchedPoints=").append(Integer.toString(nMatchedPoints))
            .append(" meanDistFromModel=")
            .append(Double.toString(meanDistFromModel))
            .append(" stDevFromMean=")
            .append(Double.toString(stDevFromMean))
            .append(" tolerance=")
            .append(Double.toString(tolerance))
            .append(" ")
            .append(parameters.toString());
        
        return sb.toString();
    }

}
