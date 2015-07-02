package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPointFit {
    
    private final TransformationParameters parameters;
    
    private final int nMatchedPoints;
    
    private int nMaxMatchable = 0;
    
    private final double meanDistFromModel;
    
    private final double stDevFromMean;
    
    private float transTolX = Float.MAX_VALUE;
    
    private float transTolY = Float.MAX_VALUE;
    
    public TransformationPointFit(TransformationParameters theParameters, 
        int numberOfMatchedPoints, double theMeanDistFromModel, 
        double theStDevFromMean, float theTranslationXTolerance,
        float theTranslationYTolerance) {
        
        parameters = theParameters;
        
        nMatchedPoints = numberOfMatchedPoints;
        
        meanDistFromModel = theMeanDistFromModel;
        
        stDevFromMean = theStDevFromMean;
        
        transTolX = theTranslationXTolerance;
        
        transTolY = theTranslationYTolerance;
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
    
    public float getRotationInRadians() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getRotationInRadians();
    }
    
    public float getScale() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getScale();
    }
    
    public float getTranslationX() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getTranslationX();
    }
    
    public float getTranslationY() {
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
    public float getTranslationXTolerance() {
        return transTolX;
    }
    
    /**
     * tolerance used when including only residuals whose absolute value
     * is less than tolerance.
     * @return 
     */
    public float getTranslationYTolerance() {
        return transTolY;
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
    
    /**
     * @return the nMaxMatchable
     */
    public int getNMaxMatchable() {
        return nMaxMatchable;
    }

    public void setTranslationXTolerance(float theTolerance) {
        transTolX = theTolerance;
    }
    
    public void setTranslationYTolerance(float theTolerance) {
        transTolY = theTolerance;
    }
    
    /**
     * @param maximumNumberMatchable the nMaxMatchable to set
     */
    public void setMaximumNumberMatchable(int maximumNumberMatchable) {
        this.nMaxMatchable = maximumNumberMatchable;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("nMatchedPoints=").append(Integer.toString(nMatchedPoints))
            .append(" nMaxMatchable=")
            .append(Double.toString(nMaxMatchable))
            .append(" meanDistFromModel=")
            .append(Double.toString(meanDistFromModel))
            .append(" stDevFromMean=")
            .append(Double.toString(stDevFromMean))
            .append(" translationXTolerance=")
            .append(Float.toString(transTolX))
            .append(" translationYTolerance=")
            .append(Float.toString(transTolY))
            .append(" ")
            .append(parameters.toString());
        
        return sb.toString();
    }

}
