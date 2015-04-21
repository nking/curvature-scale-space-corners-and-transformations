package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ANDedClauses {

    protected final PARAM[] params1;
    protected final PARAM[] params2;
    protected final COMPARISON[] gtOrLT;
    protected final float[] coefficients;
    protected final SKYCONDITIONAL[] skyConditional;
    protected final int n;
    
    public ANDedClauses(int nClauses) {
        
        this.n = nClauses;
        params1 = new PARAM[n];
        params2 = new PARAM[n];
        gtOrLT = new COMPARISON[n];
        coefficients = new float[n];
        skyConditional = new SKYCONDITIONAL[n];
    }
    
    public ANDedClauses(SKYCONDITIONAL[] skyC, PARAM[] p1, PARAM[] p2, 
        COMPARISON[] comp, float[] c) {
        
        this.n = skyC.length;
        skyConditional = skyC;
        params1 = p1;
        params2 = p2;
        gtOrLT = comp;
        coefficients = c;
    }
    
    /**
     * evaluate the clause based upon the given data and return true if
     * the clauses are all true, else return false as quickly as possible.
     * @param data
     * @return 
     */
    public boolean evaluate(ColorData data) {
              
        for (int i = 0; i < n; i++) {
            
            double param1 = data.getParameter(params1[i]);
            double param2 = data.getParameter(params2[i]);
            double coeff = coefficients[i];
            
            if (skyConditional[i].ordinal() == SKYCONDITIONAL.ALL.ordinal()) {
                
                if (gtOrLT[i].ordinal() == COMPARISON.GREATER_THAN.ordinal()) {
                    if (!((param1 / param2) > coeff)) {
                        return false;
                    }
                } else {
                    if (!((param1 / param2) < coeff)) {
                        return false;
                    }
                }
                
            } else if ((skyConditional[i].ordinal() == 
                SKYCONDITIONAL.RED.ordinal()) && data.skyIsRed()) {
                
                if (gtOrLT[i].ordinal() == COMPARISON.GREATER_THAN.ordinal()) {
                    if (!((param1 / param2) > coeff)) {
                        return false;
                    }
                } else {
                    if (!((param1 / param2) < coeff)) {
                        return false;
                    }
                }
                
            } else if ((skyConditional[i].ordinal() == 
                SKYCONDITIONAL.BLUE.ordinal()) && !data.skyIsRed()) {
                
                if (gtOrLT[i].ordinal() == COMPARISON.GREATER_THAN.ordinal()) {
                    if (!((param1 / param2) > coeff)) {
                        return false;
                    }
                } else {
                    if (!((param1 / param2) < coeff)) {
                        return false;
                    }
                }
            }
        }
        
        return true;
    }
}
