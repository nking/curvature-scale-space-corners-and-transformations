package algorithms.imageProcessing;

import java.util.Map;

/**
 *
 * @author nichole
 */
public class CustomCoeff01 implements CustomCoeff {
    
    /**
     * <pre>
     * evaluates 
     *     
     *     coefficients[1] + diffCIEX;
     * 
     * where it's expected that the equation is near
     *     15 * diffCIEX
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
        
        double diffCIEX = data.getParameter(PARAM.DIFF_CIEX);
        
        float coeff = customCoefficients.get(Integer.valueOf(1));
       
        return coeff * diffCIEX;
    }
    
}
