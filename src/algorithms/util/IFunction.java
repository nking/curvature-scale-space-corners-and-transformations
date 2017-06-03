package algorithms.util;

/**
 * interface for classes providing an implementation
 * of a function and its derivative.
 * 
 * @author nichole
 */
public interface IFunction {
    
    /** evaluates the objective at this set of coefficients.
    */
    double f (double[] coeffs);
    
    /** estimates the gradient
    */
    double[] der(double[] coeffs);
}
