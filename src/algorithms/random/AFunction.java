package algorithms.random;

import algorithms.util.IFunction;

/**
 *
 * @author nichole
 */
public abstract class AFunction implements IFunction {
    
    private final int nCoeffs1;
    private final int nCoeffs2;
    
    public AFunction(int nParams1, int nParams2) {
        this.nCoeffs1 = nParams1;
        this.nCoeffs2 = nParams2;
    }
    
    public int getNumberOfCoeffs1() {
        return nCoeffs1;
    }
    
    public int getNumberOfCoeffs2() {
        return nCoeffs2;
    }
    
}
