package algorithms.random;

import thirdparty.fsu.random.QMCHaltonAdvanced;

/**
 *
 * @author nichole
 */
public class QuasiMonteCarlo {
    
    public double haltonAdvanced(AFunction f, int nIterations) {
        
        int nDim = f.getNumberOfCoeffs1();
        int[] base = new int[nDim];
        int seed[] = new int[nDim];
        //TODO: edit the step size for the problem or allow configuration
        int[] step_vec = new int[] { 0, 5, 1000, 1000000 };
        int step;
        double[] r = new double[nDim];
  
        QMCHaltonAdvanced qmcA = new QMCHaltonAdvanced();
        qmcA.halton_dim_num_set(nDim);
        step = step_vec[0];
        qmcA.halton_step_set(step);
        
        for (int i = 0; i < nDim; i++) {
            base[i] = qmcA.prime(i + 1);
        }
        qmcA.halton_base_set(base);
        
        double sum = 0;
        
        for (int ii = 0; ii < nIterations; ii++) {
            
            qmcA.halton(r);
            
            //System.out.println("R=" + Arrays.toString(r));
            
            sum += f.f(r);
        }
        sum /= (double)nIterations;
        
        return sum;
    }
    
    public double[] haltonAdvanced(int nDimensions) {
        
        double[] r = new double[nDimensions];
        
        haltonAdvanced(nDimensions, r);
        
        return r;
    }
    
    /**
     * initialize an instance of the qmc and return it.
     * 
     * To use it:
     *    initialize a double array of size nPoints and feed it to
     *    the halton instance as many times as needed for new points.
     * 
     *    double[] output = new double[nDimension];
     *    qmcA.halton(output);     
     * 
     * @param nDimension the number of points to generate for each iteration.
     * @return 
     */
    public QMCHaltonAdvanced initalize(int nDimension) {
        
        int nDim = nDimension;
        int[] base = new int[nDim];
        int seed[] = new int[nDim];
        //TODO: edit the step size for the problem or allow configuration
        int[] step_vec = new int[] { 0, 5, 1000, 1000000 };
        int step;
  
        QMCHaltonAdvanced qmcA = new QMCHaltonAdvanced();
        qmcA.halton_dim_num_set(nDim);
        step = step_vec[0];
        qmcA.halton_step_set(step);
        
        for (int i = 0; i < nDim; i++) {
            base[i] = qmcA.prime(i + 1);
        }
        qmcA.halton_base_set(base);
        
        return qmcA;
    }
    
    public void haltonAdvanced(int nPoints, double[] output) {

        if (output.length != nPoints) {
            throw new IllegalArgumentException("nPoints and output.length must "
                + "be the same");
        }
                
        QMCHaltonAdvanced qmcA = initalize(nPoints);
        
        qmcA.halton(output);        
    }
    
    public double lds(AFunction function, int nPoints) {
        
        int nParams = function.getNumberOfCoeffs1();
        
        double sum = 0;
        
        LowDiscrepancySequences lds = new LowDiscrepancySequences();
        
        double[] pt = new double[nParams];
        
        for (int i = 0; i < nPoints; ++i) {
            
            if (nParams == 2) {
                lds.haltonPoint(i, pt);
            } else {
                for (int j = 0; j < nParams; ++j) {
                    pt[j] = lds.halton(j, j + 1);
                }
            }
            
            sum += function.f(pt);
        }
        
        sum /= (double)nPoints;
        
        return sum;
    }
}
