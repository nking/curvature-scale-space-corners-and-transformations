package algorithms.search.global;

import algorithms.random.AFunction;
import algorithms.random.QuasiMonteCarlo;
import java.util.Arrays;
import thirdparty.dlib.optimization.LBFGSOptimization;
import thirdparty.dlib.optimization.LBFGSSearchStrategy;
import thirdparty.dlib.optimization.ObjectiveDeltaStopStrategy;
import thirdparty.fsu.random.QMCHaltonAdvanced;

/**
 * a global search to find the minimum of a given function within
 * the bounds given.
 * NRegions are sampled using Quasi Monte Carlo
 * and each region is searched using local search LBFGs.
 * 
 * LBFGS is a fast local search - it uses the
 * gradient provide by the function, which can be the first derivative or
 * a calculation using finite difference method.
 * 
 * Note that if one knows the function is convex, one can use LBFGS alone
 * to find the global minimum.

     The algorithm uses low discripeancy sequences, a.k.a. quasi-random sequences.
   from "Low Discrepancy Sequences for Monte Carlo Simulations on Reconfigurable Platforms"
   2008, Dalal, Stefan, and Harwayne-Gidansky
     The Cooper Union for the Advancement of Science and Art
     51 Astor Place, New York, NY 10003
       ...Low-discrepancy sequences, also known as “quasirandom” sequences, are numbers 
       that are better equidistributed in a given volume than pseudo-random
       numbers. Evaluation of high-dimensional integrals is commonly required in 
       scientific fields as well as other areas (such as finance), and is 
       performed by stochastic Monte Carlo simulations. Simulations which use 
       quasirandom numbers can achieve faster convergence and better accuracy 
       than simulations using conventional pseudo-random numbers. Such 
       simulations are called Quasi-Monte Carlo.
 * 
 * @author nichole
 */
public class QMCLBFGS {

    /**
     * search for the minimum of function in nRegions within the bounds 
     * using LBFGs and return the best
     * result as coordinates and the y value.
     * 
     * @param function function and the derivative or gradient of the function.
     * @param nRegions number of regions to sample over the boundary region
     * and then search using LBFGS
     * @param startBounds start of region to search in each dimension
     * @param stopBounds stop of region to search in in each dimension
     * @return an array of the coefficients for the minimum followed by
     * the minimum.
     */
    public double[] search(AFunction function, int nRegions,
        double[] startBounds, double[] stopBounds) {
        
        int nDim = function.getNumberOfCoeffs1();
        
        if (nDim < 1) {
            throw new IllegalArgumentException("nDim in function must be 1 or"
                + " larger");
        } 
        if (startBounds.length != stopBounds.length) {
            throw new IllegalArgumentException("bounds arrays must have same "
                + " length");
        }
        
        if (startBounds.length != nDim) {
            throw new IllegalArgumentException("bounds arrays must have same "
                + " length as function number of dimensions");
        }
        
        //output array for samples
        double[] cache = new double[nDim];
        
        double[] range = new double[nDim];
        for (int i = 0; i < nDim; ++i) {
            double b0 = startBounds[i];
            double b1 = stopBounds[i];
            if (b1 <= b0) {
                throw new IllegalArgumentException("stop bounds must be larger "
                    + " than start bounds");
            }
            range[i] = b1 - b0;
        }
        
        QuasiMonteCarlo qmc = new QuasiMonteCarlo();
        
        double[] bestResults = new double[nDim];
        Arrays.fill(bestResults, Double.POSITIVE_INFINITY);
        double bestMin = Double.MAX_VALUE;
        
        QMCHaltonAdvanced qmcA = qmc.initalize(nDim);
    
        for (int i = 0; i < nRegions; ++i) {
                    
            qmcA.halton(cache);
            
            // scale the sample to the bounds
            for (int j = 0; j < startBounds.length; ++j) {
                cache[j] = cache[j] + (cache[j] * range[j]);
            }
            
            System.out.println("input coeff=" + Arrays.toString(cache));
            
            LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
            ObjectiveDeltaStopStrategy stopStrategy
                = new ObjectiveDeltaStopStrategy(1.e-5, 100);

            double fLower = -10;

            LBFGSOptimization opt = new LBFGSOptimization();
            double min = opt.findMin(searchStrategy,
                stopStrategy, function, cache, fLower);

            System.out.format("  i=%d min=%f r=%s\n", i, min,
                Arrays.toString(cache));

            if (min < bestMin) {
                bestMin = min;
                System.arraycopy(cache, 0, bestResults, 0, cache.length);
            }
        }
        
        bestResults = Arrays.copyOf(bestResults, bestResults.length + 1);
        bestResults[bestResults.length - 1] = bestMin;
        
        return bestResults;
    }
   
    
}
