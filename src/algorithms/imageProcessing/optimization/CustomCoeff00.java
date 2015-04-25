package algorithms.imageProcessing.optimization;

import java.util.Map;

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
     * @param customCoefficients
     * @return 
     */
    @Override
    public double evaluate(ColorData data, float[] coefficients,
        Map<Integer, Float> customCoefficients) {
        
        double absContr = data.getParameter(PARAM.ABSOLUTE_CONTRAST);
        
        float coeff4 = customCoefficients.get(Integer.valueOf(4));
        float coeff5 = customCoefficients.get(Integer.valueOf(5));
        float coeff6 = customCoefficients.get(Integer.valueOf(6));
        
        double r = coeff4 + ((absContr - coeff5) * coeff6);
        
        // 1.5 = cofficients[4]
        // 0.5 = cofficients[5]
        //-2.0 = cofficients[6]
        //
        //(1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
        
        return r;
    }
    
}
