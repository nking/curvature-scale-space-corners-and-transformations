package algorithms.util;

/**
 *
 * @author nichole
 */
public class CoefficientWrapper {
    
    private final float[] coefficients;
    
    public CoefficientWrapper(float[] coeff) {
        this.coefficients = coeff;
    }
    
    public float[] getCoefficients() {
        return coefficients;
    }
}
