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
        
        if (skyC == null) {
            throw new IllegalArgumentException("skyC cannot be null");
        }
        if (p1 == null) {
            throw new IllegalArgumentException("p1 cannot be null");
        }
        if (p2 == null) {
            throw new IllegalArgumentException("p2 cannot be null");
        }
        if (comp == null) {
            throw new IllegalArgumentException("comp cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        
        this.n = skyC.length;
        
        if (p1.length != n) {
            throw new IllegalArgumentException(
                "p1 is not the same length as skyC");
        }
        
        if (p2.length != n) {
            throw new IllegalArgumentException(
                "p2 is not the same length as skyC");
        }
        
        if (comp.length != n) {
            throw new IllegalArgumentException(
                "comp is not the same length as skyC");
        }
        
        if (c.length != n) {
            throw new IllegalArgumentException(
                "c is not the same length as skyC");
        }
        
        skyConditional = skyC;
        params1 = p1;
        params2 = p2;
        gtOrLT = comp;
        coefficients = c;
    }
    
    public void set(int index, SKYCONDITIONAL skyC, PARAM p1, PARAM p2, 
        COMPARISON comp, float c) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        skyConditional[index] = skyC;
        params1[index] = p1;
        params2[index] = p2;
        gtOrLT[index] = comp;
        coefficients[index] = c;
    }
    
    public float getCoefficients(int index) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        return coefficients[index];
    }
    
    public COMPARISON getGtOrLT(int index) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        return gtOrLT[index];
    }
    
    public PARAM getParams1(int index) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        return params1[index];
    }
    
    public PARAM getParams2(int index) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        return params2[index];
    }
    
    public SKYCONDITIONAL getSKYCONDITIONAL(int index) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        return skyConditional[index];
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
            
            if ((skyConditional[i].ordinal() == 
                SKYCONDITIONAL.RED.ordinal()) && !data.skyIsRed()) {
                
                return false;
             
            } else if ((skyConditional[i].ordinal() == 
                SKYCONDITIONAL.BLUE.ordinal()) && data.skyIsRed()) {
                
                return false;
                
            } if (skyConditional[i].ordinal() == SKYCONDITIONAL.ALL.ordinal()) {
                
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
