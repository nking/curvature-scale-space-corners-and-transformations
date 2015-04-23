package algorithms.imageProcessing;

import java.util.Map;

/**
 * a custom coefficient built specifically for a statement within
 * an ANDedClauses instance.
 * 
 * @author nichole
 */
public interface CustomCoeff {
    
    public double evaluate(ColorData data, float[] coefficients,
        Map<Integer, Float> customCoefficients);
  
}
