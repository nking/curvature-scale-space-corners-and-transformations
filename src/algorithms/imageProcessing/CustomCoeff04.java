package algorithms.imageProcessing;

import java.util.Map;

/**
 *
 * @author nichole
 */
public class CustomCoeff04 implements CustomCoeff {
    
    /**
     * <pre>
     * evaluates 
     *     
     *     coefficients[7] + diffCIEY;
     * 
     * where it's expected that the equation is near
     *     15 * diffCIEY
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
        
        double diffCIEY = data.getParameter(PARAM.DIFF_CIEY);
        
        float coeff = customCoefficients.get(Integer.valueOf(7));
       
        return coeff * diffCIEY;
    }
    
}
