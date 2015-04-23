package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class CustomCoeff00 implements CustomCoeff {
    
    /**
     * <pre>
     * evaluates 
     *     coefficients[4] + (ABSOLUTE_CONTRAST - coefficients[5]) * (coefficients[6]);
     * 
     * where it's expected that the equation is near
     *     (1.5 + ((Math.abs(contrastV) - 0.5) * (-2.0)) )
     * </pre>
     * 
     * @param data
     * @param coefficients
     * @return 
     */
    public double evaluate(ColorData data, float[] coefficients) {
        
        double absContr = data.getParameter(PARAM.ABSOLUTE_CONTRAST);
        
        double r = coefficients[4] + (absContr - coefficients[5]) * (coefficients[6]);
        
        // 1.5 = cofficients[4]
        // 0.5 = cofficients[5]
        //-2.0 = cofficients[6]
        //
        //(1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
        
        return r;
    }
    
}
