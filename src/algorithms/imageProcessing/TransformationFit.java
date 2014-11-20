package algorithms.imageProcessing;

import algorithms.util.PairIntArray;

/**
 *
 * @author nichole
 */
public class TransformationFit {
    
    private TransformationParameters parameters = null;
    
    private PairIntArray[] transformedEdges1 = null;
    
    private double chiSqSum = Double.MAX_VALUE;
        
    public TransformationFit(TransformationParameters theParameters, 
        PairIntArray[] theTransformedEdges1) {
        
        parameters = theParameters;
        
        transformedEdges1 = theTransformedEdges1;
    }

    /**
     * @return the parameters
     */
    public TransformationParameters getParameters() {
        return parameters;
    }

    /**
     * @return the transformedEdges1
     */
    public PairIntArray[] getTransformedEdges1() {
        return transformedEdges1;
    }

    /**
     * @return the chiSqSum
     */
    public double getChiSqSum() {
        return chiSqSum;
    }

    /**
     * @param theChiSqSum the chi squared sum to set
     */
    public void setChiSqSum(double theChiSqSum) {
        this.chiSqSum = theChiSqSum;
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
   
}
