package algorithms.imageProcessing.optimization;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

/**
 *
 * @author nichole
 */
public class ANDedClauses implements Cloneable {

    protected final PARAM[] params1;
    protected final PARAM[] params2;
    protected final COMPARISON[] gtOrLT;
    protected final float[] coefficients;
    protected final SKYCONDITIONAL skyConditional;
    protected final int n;
    
    protected Map<Integer, CustomCoeff> customCoefficients = 
        new HashMap<Integer, CustomCoeff>();
    
    /**
     * these are coefficients not present in float[] coefficients
     * that are needed for customCoefficients.  
     * The invoker of the constructor is responsible for keeping track of 
     * filling both.
     */
    protected final Map<Integer, Float> customCoefficientVariables = 
        new HashMap<Integer, Float>();
    
    public ANDedClauses(int nClauses, SKYCONDITIONAL skyC) {
        
        this.n = nClauses;
        skyConditional = skyC;
        params1 = new PARAM[n];
        params2 = new PARAM[n];
        gtOrLT = new COMPARISON[n];
        coefficients = new float[n];
    }
    
    public ANDedClauses(SKYCONDITIONAL skyC, PARAM[] p1, PARAM[] p2, 
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
        
        this.n = p1.length;
        
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
    
    public void set(int index, PARAM p1, PARAM p2, COMPARISON comp, float c) {
        
        if ((index < 0) || (index > (n - 1))) {
            throw new IllegalArgumentException(
            "index is out of bounds of arrays size n");
        }
        
        params1[index] = p1;
        params2[index] = p2;
        gtOrLT[index] = comp;
        coefficients[index] = c;
    }
    
    public void setACustomCoefficient(int index, CustomCoeff cCoeff) {
        
        customCoefficients.put(Integer.valueOf(index), cCoeff);
    }
    
    /**
     * Set a coefficient needed by an instance of CustomCoeff.
     * 
     * @param coeffNumber
     * @param coefficient 
     */
    public void setCustomCoefficientVariable(int coeffNumber, Float coefficient) {
        
        customCoefficientVariables.put(Integer.valueOf(coeffNumber), coefficient);
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
    
    public SKYCONDITIONAL getSKYCONDITIONAL() {
        
        return skyConditional;
    }
    
    /**
     * evaluate the clause based upon the given data and return true if
     * the clauses are all true, else return false as quickly as possible.
     * @param data
     * @return 
     */
    public boolean evaluate(ColorData data) {
              
        if ((skyConditional.ordinal() == SKYCONDITIONAL.RED.ordinal()) 
            && !data.skyIsRed()) {        
            return false;
        } else if ((skyConditional.ordinal() == SKYCONDITIONAL.BLUE.ordinal()) 
            && data.skyIsRed()) {
            return false;
        }
        
        // iterate over each clause
        for (int i = 0; i < n; i++) {
            
            double param1 = data.getParameter(params1[i]);
            double param2 = data.getParameter(params2[i]);
            double coeff = coefficients[i];
            
            CustomCoeff cCoeff = customCoefficients.get(Integer.valueOf(i));
            if (cCoeff != null) {
                coeff = cCoeff.evaluate(data, coefficients, 
                    customCoefficientVariables);
            }
            
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
        
        return true;
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        return copy();
    }
        
    public ANDedClauses copy() {
        
        ANDedClauses c = new ANDedClauses(skyConditional, 
            Arrays.copyOf(params1, params1.length), 
            Arrays.copyOf(params2, params2.length),
            Arrays.copyOf(gtOrLT, gtOrLT.length), 
            Arrays.copyOf(coefficients, coefficients.length));
        
        if (!customCoefficients.isEmpty()) {
            Iterator<Entry<Integer, CustomCoeff>> iter = 
                customCoefficients.entrySet().iterator();
            while (iter.hasNext()) {
                Entry<Integer, CustomCoeff> entry = iter.next();
                c.setACustomCoefficient(entry.getKey(), entry.getValue());
            }
        }
        
        if (!customCoefficientVariables.isEmpty()) {
            Iterator<Entry<Integer, Float>> iter = 
                customCoefficientVariables.entrySet().iterator();
            while (iter.hasNext()) {
                Entry<Integer, Float> entry = iter.next();
                c.setCustomCoefficientVariable(entry.getKey(), entry.getValue());
            }
        }
        
        return c;
    }
}
