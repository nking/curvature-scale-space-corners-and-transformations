package algorithms.curves;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;

import algorithms.curves.GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM;
import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;

/**
   Class to find the minimum difference between a GEV model function(k, sigma, mu)
   and normalized data with errors.
   
   To make a solution faster than the very robust but slow downhill simplex method,
   would like to use the partial derivatives of the GEV 
   to make an iterative solution for chi-square minimization 
   of a non-linear, non-quadratic GEV model's difference from the data.
   
   This solution uses a preconditioner matrix with ICU0 function to help
   determine the next stop and direction. Note that the suggested changes 
   might be applied by the NonQuadraticConguteSolver as
   a fraction of the suggested change returned by the lineSearch method
 
  Useful for implementing the code below was reading:
   
   
   http://en.wikipedia.org/wiki/Fletcher-Reeves#Nonlinear_conjugate_gradient
   
   and
   
   http://netlib.org/linalg/html_templates/node64.html#figdilu
   
   and browsing various books on optimization
   
    =============================================================================================
    
    Minimize f(vars)
         
     The function f is called the objective function or cost function.
     
      The vector x is an n-vector of independent variables: 
         vars = [var_1, var_2, …, var_n]^T is a member of set Real numbers. 
         The variables var_1, …, var_n are often referred to as decision variables. 
     
      The optimization problem above can be viewed as a decision problem that involves 
      finding the “best” vector var of the decision variables over all possible vectors in Ω. 
      By the “best” vector we mean the one that results in the-smallest value of 
      the objective function. 
      
      This vector is called the minimizer of f over Ω. 
      It is possible that there may be many minimizers. In this case, finding 
      any of the minimizers will suffice.
          
       Df is the first derivative of f(vars) and is 
           [partial deriv f/partial deriv var_1, partial deriv f/partial deriv var_1, ...]
       
       ∇f = the gradient of f. 
       ∇f = (Df)^T
       
       F(vars) is the 2nd derivative of f and is sometimes called the Hessian.
                                 | ∂^2f/∂^2_var_1       ...    ∂^2f/∂_var_n d_var_1
          F(vars) = D^2f(vars) = |         ...          ...         ...
                                 | ∂^2f/∂_var_1 d_var_n ...    ∂^2f/∂^2_var_n
       
       Example:  Let f(x1, x2) = 5(x_1) + 8(x_2) + (x_1)(x_2) − (x_1)^2 − 2(x_2)^2
             Df(x) = (∇f(x))^T = [∂f/dx_1, ∂f/dx2] = [5 + x_2 - 2x_1,  8 + x_1 - 4x_2]
             
                                 | -2  1 |
             F(x) = D^2f(x)    = |  1 -4 |
             
       --------------------------------------------------------------------------------
       
       f(k, sigma, mu) = the GEV function
       
       Df = [∂f/∂k, ∂f/∂sigma, ∂f/∂mu]
       
                     |   ∂f/∂k   |
       ∇f = (Df)^T = | ∂f/∂sigma |
                     |   ∂f/∂mu  |
            
                                              | ∂^2f/∂k∂k       ∂^2f/∂sigma∂k        ∂^2f/∂mu∂k     |
       F(k, sigma, mu) = D^2f(k, sigma, mu) = | ∂^2f/∂k∂sigma   ∂^2f/∂sigma∂sigma    ∂^2f/∂mu∂sigma |
                                              | ∂^2f/∂k∂mu      ∂^2f/∂sigma∂mu       ∂^2f/∂mu∂mu    |
            
       M = D^2f(k, sigma, mu)
       
                                | ∂^2f/∂k∂k        ∂^2f/∂k∂sigma       ∂^2f/∂k∂mu     |     |   ∂f/∂k   |
       M^(-1) * ∇f = M^T * ∇f = | ∂^2f/∂sigma∂k    ∂^2f/∂sigma∂sigma   ∂^2f/∂sigma∂mu |  *  | ∂f/∂sigma |
                                | ∂^2f/∂mu∂k       ∂^2f/∂mu∂sigma      ∂^2f/∂mu∂mu    |     |   ∂f/∂mu  |
                                
                                | (∂^2f/∂k∂k) * (∂f/∂k) +  (∂^2f/∂k∂sigma) * (∂f/∂sigma) + (∂^2f/∂k∂mu) * (∂f/∂mu)          |
                              = | (∂^2f/∂sigma∂k) * (∂f/∂k) + (∂^2f/∂sigma∂sigma) * (∂f/∂sigma) + (∂^2f/∂sigma∂mu)*(∂f/∂mu) |
                                | (∂^2f/∂mu∂k)*(∂f/∂k) + (∂^2f/∂mu∂sigma)*(∂f/∂sigma) + (∂^2f/∂mu∂mu)*(∂f/∂mu)              |
                                
       Can use Incomplete Cholesky factorization with fill 0 on M^(-1) to create the pre-conditioning matrix.
       
           http://netlib.org/linalg/html_templates/node64.html#figdilu
           
           Let S be the non-zero set ({i,j} : a_i_j != 0 )
           
            for i = 1, 2, ...
                set d_i_i = a_i_i
            for i = 1, 2, ...
                set d_i_i = 1/d_i_i
                for j = i + 1, i + 2, ...                       
                if (i, j) in set S and (j, i) in set S then
                    set d_j_j = d_j_j - a_j_i * d_i_i * a_i_j
            
            set d(1,1) = (∂^2f/∂k∂k)
            set d(2,2) = (∂^2f/∂sigma∂sigma)
            set d(3,3) = (∂^2f/∂mu∂mu)
            
            i = 1:
                    set d(1,1) = 1./d(1,1) = 1./(∂^2f/∂k∂k)
                
                    | (1,1) (1,2) (1,3) | = | 1/(∂^2f/∂k∂k)   
                d = | (2,1) (2,2) (2,3) | = |                    (∂^2f/∂sigma∂sigma)  
                    | (3,1) (3,2) (3,3) | = |                                             (∂^2f/∂mu∂mu) |
                  
                j=2:
                    set d(2,2) = (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma)
                j=3:
                    set d(3,3) = (∂^2f/∂mu∂mu) - (∂^2f/∂mu∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂mu)
                    
            i = 2:
                    set d(2,2) = 1./d(2,2) = ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
                    
                    | (1,1) (1,2) (1,3) | = | 1/(∂^2f/∂k∂k)   
                d = | (2,1) (2,2) (2,3) | = |                latest d(2,2)  
                    | (3,1) (3,2) (3,3) | = |                                             (∂^2f/∂mu∂mu) |
                    
                j=3:
                    set d(3,3) = d(3,3) - a(3,2) * d(2,2) * a(2,3)
                               =
                                 ( (∂^2f/∂mu∂mu) - (∂^2f/∂mu∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂mu) )
                                 -
                                 (
                                    (∂^2f/∂mu∂sigma)
                                    *
                                    ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
                                    *
                                    ∂^2f/∂sigma∂mu
                                 )
                                  
            i = 3:
                    set d(3,3) = 1./d(3,3)
                    
                        = 1./(
                                 ( (∂^2f/∂mu∂mu) - (∂^2f/∂mu∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂mu) )
                                 -
                                 (
                                    (∂^2f/∂mu∂sigma)
                                    *
                                    ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
                                    *
                                    ∂^2f/∂sigma∂mu
                                 )
                             ) 
                             
                             
    Then using the ICU0 matrix as preconditioner:
           
                                        | d(1,1)   0        0      |     |   ∂f/∂k   |
        (M_icuo)^(-1) * ∇f = M^T * ∇f = | 0        d(2,2)   0      |  *  | ∂f/∂sigma |
                                        | 0        0        d(3,3) |     |   ∂f/∂mu  |
                                    
                                        | d(1,1) * (∂f/∂k)     |
                                      = | d(2,2) * (∂f/∂sigma) |
                                        | d(3,3) * (∂f/∂mu)    |

             where d(1,1) is 1./(∂^2f/∂k∂k)
                   d(2,2) is ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
                   d(3,3) is 1./(
                                     ( (∂^2f/∂mu∂mu) - (∂^2f/∂mu∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂mu) )
                                     -
                                     (
                                        (∂^2f/∂mu∂sigma)
                                        *
                                        ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
                                        *
                                        ∂^2f/∂sigma∂mu
                                     )
                                 ) 
                                 
       Note that below in the code, r is ∇f.
       
  @author nichole
 */
public class NonQuadraticConjugateGradientSolver extends AbstractCurveFitter {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected int maxIterations = 200;
    
    protected PolygonAndPointPlotter plotter;
    
    protected static float eps = 1e-8f;
    
    protected static float convergedEps = 0.00001f;
    
    public NonQuadraticConjugateGradientSolver(float[] xPoints, float[] yPoints,
        float[] xErrPoints, float[] yErrPoints) {

        super(xPoints, yPoints, xErrPoints, yErrPoints);
        
        try {
            plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
        } catch (IOException e) {
            log.severe(e.getMessage());
        }
    }    

    public void setMaximumNumberOfIterations(int maxNumber) {
        this.maxIterations = maxNumber;
    }
    
    /**
     * fit the x and y data using built in default ranges for the parameters k, sigma, and mu.
     * 
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {
        
        float kMin = 0.001f;
        float kMax = 2.0f;
        float sigmaMin = 0.025f;
        float sigmaMax = 0.5f;
        float muMin = 0.0001f;
        float muMax = 0.3f;
        
        // if yScale is small, muMax should be as large as 0.6 roughly
        if (yScale < 100) {
            muMax = 0.5f;
            sigmaMax = 0.5f;
        } else {
            int yMaxIndex = MiscMath.findYMaxIndex(y);
            float xAtYMax = x[yMaxIndex];
            muMin = xAtYMax/2.f;
            muMax = xAtYMax * 2.0f;
            kMax = 1.0f;
            if (yScale > 200) {
                kMin = 0.2f;
            }
        }
        
        return fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
    }

    /**
     * find the best fitting GEV by solving for each parameter in set {k, sigma, mu} separately
     * rather then minimizing the function for suggested changes by all derivatives at once.
     * 
     * So far, this is resulting in the best fits, but is sensitive to the starting point
     * and has not been tested over a wide range of data distributions.
     * 
     * The range of values given to this method by TwoPointVoidStats are those found to be most useful
     * for representing the range of normalized GEV curves that match the datasets given to it.
     * k < 0 are not fit because the distributions are not physical for the expected datasets,
     * though that can be changed if needed.
     * 
     * @param kMin  minimum range of value of k, the shape parameter
     * @param kMax  maximum range of value of k, the shape parameter
     * @param sigmaMin  minimum range of value of sigma, the scale parameter
     * @param sigmaMax  maximum range of value of sigma, the scale parameter
     * @param muMin  minimum range of value of mu, the location parameter
     * @param muMax  maximum range of value of mu, the location parameter
     * @return
     * @throws FailedToConvergeException
     */
    public GEVYFit fitCurveParametersSeparately(float kMin, float kMax, float sigmaMin, float sigmaMax,
        float muMin, float muMax) throws FailedToConvergeException {
        
        if (kMin < 0) {
            throw new IllegalArgumentException("kMin must be larger than zero");
        }
        if (kMin > kMax) {
            throw new IllegalArgumentException("kMin must be less than kMax");
        }
        if (muMin < 0) {
            throw new IllegalArgumentException("muMin must be larger than zero. mu is usually near the peak of the normalized histogram's x value.");
        }
        if (muMin > muMax) {
            throw new IllegalArgumentException("muMin must be less than muMax");
        }
        /*if (muMax > xmax) {
            throw new IllegalArgumentException("muMax must be less than the maximum value of x in the histogram (" + xmax + ")");
        }*/
        if (sigmaMin < 0) {
            throw new IllegalArgumentException("sigmaMin must be larger than zero");
        }
        if (sigmaMin > sigmaMax) {
            throw new IllegalArgumentException("sigmaMin must be less than sigmaMax");
        }
                
        /*
        float kVar = kMin;
        float sigmaVar = sigmaMin;
        float muVar = muMin;
        */
        /*
        float kVar = kMax;
        float sigmaVar = sigmaMax;
        float muVar = muMax;
        */
        
        float kVar = (kMax + kMin)/2.f;
        float sigmaVar = (sigmaMax + sigmaMin)/2.f;
        float muVar = (muMax + muMin)/2.f;
        
        // the variables k, sigma, and mu
        float[] vars    = new float[]{kVar, sigmaVar, muVar};
        float[] varsMin = new float[]{kMin, sigmaMin, muMin};
        float[] varsMax = new float[]{kMax, sigmaMax, muMax};
     
        int varStopIdx = vars.length - 1;
        
        // r is current residual.  it holds deltaK, deltaSigma, and deltaMu
        float[] r = new float[3];
        DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, 0, varStopIdx);
        
        // chiSqSumForLineSearch[0] holds current best chiSqSum for the last change in vars
        // chiSqSumForLineSearch[1] holds the return value from lineSearch
        float[] chiSqSumForLineSearch = new float[2];
        
        chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);
        
        int nSameSequentially = 0;
        float epsChiSame = 1e-5f;
        float lastChiSqSum = Float.MAX_VALUE;
        int nAltSolutionCount = 0;
        
        int nIter = 0;
        while ((nIter < maxIterations) && (chiSqSumForLineSearch[0] > convergedEps)) {
            
            /*if ((nIter > 1) && residualsAreSame(rPrev, r)) {
                break;
            }*/
            
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
            float chiSqSum = DerivGEV.chiSqSum(yGEV, y, ye);
            
            if (debug) {
                try {
                    String label = String.format(
                       "k=%4.4f <1.8>  s=%4.4f <0.85>  m=%4.4f <0.441>  n=%d  chi=%4.6f",
                        vars[0], vars[1], vars[2], nIter, chiSqSum);
                    plotFit(yGEV, label);
                } catch (IOException e) {
                    System.err.println(e.getMessage());
                }
            }
            
            // if solution has stalled, suggest deltas and apply the accepted
            if ((nIter > 3) && (nSameSequentially > 2)) {
                
                int idx = (nAltSolutionCount % 3);

                DerivGEV.exploreChangeInVars(vars, varsMin, varsMax, x, y, ye, r, chiSqSumForLineSearch, idx);
                
                // OR temporarily allow changes that increase chiSqSum
                
                nAltSolutionCount++;
                
                // apply changes
                if (chiSqSumForLineSearch[1] < chiSqSumForLineSearch[0]) {
                    for (int k = 0; k <= varStopIdx; k++) {
                        vars[k] = vars[k] + r[k];
                        chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                    }
                    nSameSequentially = 0;
                    lastChiSqSum = chiSqSumForLineSearch[0];
                    nIter++;
                    continue;
                }
            }
            
            for (int k = 0; k <= varStopIdx; k++) {
            
                if (nIter > 0) {
                    
                    // populate r with the best fitting derivatives for vars[]
                    DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, k, k);
                    
                    log.finest("   ->r[" + k + "]=" + r[k]  + "  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                }
                
                // line search finds the fraction of the derivatives in p to apply to the GEV to reduce the chi sq sum
                float alpha = lineSearch(r, vars, varsMin, varsMax, chiSqSumForLineSearch, k, k);
                if (alpha <= eps) {
                    log.finest("       r[" + k + "]=" + r[k] + "  last chiSqSum=" + chiSqSumForLineSearch[1]);
                    continue;
                }
                float ap = alpha*r[k];
                
                float tmpVar = vars[k] + ap;                
                
                if (!chiSqSumIsNotAcceptable(chiSqSumForLineSearch[0], chiSqSumForLineSearch[1])) {
                    vars[k] = tmpVar;
                    chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                }
                
                log.finest("    alpha=" + alpha + "  -> vars[" + k + "]=" + vars[k] + "  chiSqSum=" + chiSqSumForLineSearch[0] + " nIter=" + nIter);
            }
            
            if (Math.abs(lastChiSqSum - chiSqSumForLineSearch[0]) < epsChiSame) {
                nSameSequentially++;
            } else {
                nSameSequentially = 0;
            }
            lastChiSqSum = chiSqSumForLineSearch[0];
            
            nIter++;
        }
        
        float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);

        float chisqsum = calculateChiSquareSum(yGEV, WEIGHTS_DURING_CHISQSUM.ERRORS);

        GEVYFit yfit = new GEVYFit();
        yfit.setChiSqSum(chisqsum);
        yfit.setK(vars[0]);
        yfit.setSigma(vars[1]);
        yfit.setMu(vars[2]);
        yfit.setYFit(yGEV);
        yfit.setX(x);
        yfit.setXScale(xScale);
        yfit.setYScale(yScale);
        yfit.setYDataErrSq( calcYErrSquareSum() ); 
        
        log.info("number of times the alt solution was needed = " + nAltSolutionCount);

        return yfit;
    }
  
    /**
     * fit the x, y data with a GEV whose parameters are within the given ranges for k,
     * sigma, and mu.  The method attempts to fit for changes in k, sigma, and mu all
     * at once for each iteration.
     * 
     * @param kMin  minimum range of value of k, the shape parameter
     * @param kMax  maximum range of value of k, the shape parameter
     * @param sigmaMin  minimum range of value of sigma, the scale parameter
     * @param sigmaMax  maximum range of value of sigma, the scale parameter
     * @param muMin  minimum range of value of mu, the location parameter
     * @param muMax  maximum range of value of mu, the location parameter
     * @return
     * @throws FailedToConvergeException
     */
    public GEVYFit fitCurve(float kMin, float kMax, float sigmaMin, float sigmaMax,
        float muMin, float muMax) throws FailedToConvergeException {

        if (kMin < 0) {
            throw new IllegalArgumentException("kMin must be larger than zero");
        }
        if (kMin > kMax) {
            throw new IllegalArgumentException("kMin must be less than kMax");
        }
        if (muMin < 0) {
            throw new IllegalArgumentException("muMin must be larger than zero. mu is usually near the peak of the normalized histogram's x value.");
        }
        if (muMin > muMax) {
            throw new IllegalArgumentException("muMin must be less than muMax");
        }
        /*
        if (muMax > xmax) {
            throw new IllegalArgumentException("muMax must be less than the maximum value of x in the histogram (" + xmax + ")");
        }*/
        if (sigmaMin < 0) {
            throw new IllegalArgumentException("sigmaMin must be larger than zero");
        }
        if (sigmaMin > sigmaMax) {
            throw new IllegalArgumentException("sigmaMin must be less than sigmaMax");
        }
        
        float kVar = kMin;
        float sigmaVar = sigmaMin;
        float muVar = muMin;
        
        /*
        float kVar = kMax;
        float sigmaVar = sigmaMax;
        float muVar = muMax;
        */
        /*
        float kVar = (kMax + kMin)/2.f;
        float sigmaVar = (sigmaMax + sigmaMin)/2.f;
        float muVar = (muMax + muMin)/2.f;
        */
        
        // the variables k, sigma, and mu
        float[] vars = new float[]{kVar, sigmaVar, muVar};
        float[] varsMin = new float[]{kMin, sigmaMin, muMin};
        float[] varsMax = new float[]{kMax, sigmaMax, muMax};
     
        int varStopIdx = vars.length - 1;
        
        // r is current residual.  it holds deltaK, deltaSigma, and deltaMu
        float[] r = new float[3];
        DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, 0, varStopIdx);
        
        // chiSqSumForLineSearch[0] holds current best chiSqSum for the last change in vars
        // chiSqSumForLineSearch[1] holds the return value from lineSearch
        float[] chiSqSumForLineSearch = new float[2];
        
        chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);
        
        int nSameSequentially = 0;
        float epsChiSame = 1e-5f;
        float lastChiSqSum = Float.MAX_VALUE;
        int nAltSolutionCount = 0;
        
        int nIter = 0;
        while ((nIter < maxIterations) && (chiSqSumForLineSearch[0] > convergedEps)) {
            
            /*if ((nIter > 1) && residualsAreSame(rPrev, r)) {
                break;
            }*/
           
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
            float chiSqSum = DerivGEV.chiSqSum(yGEV, y, ye);
                       
            if (debug) {
                try {
                    String label = String.format(
                       "k=%4.4f <1.8>  s=%4.4f <0.85>  m=%4.4f <0.441>  n=%d  chi=%4.6f",
                        vars[0], vars[1], vars[2], nIter, chiSqSum);
                    plotFit(yGEV, label);
                } catch (IOException e) {
                    System.err.println(e.getMessage());
                }
            }
            
            // if solution has stalled, apply changes that improve the solution in small steps
            if ((nIter > 3) && (nSameSequentially > 2)) {
                
                int idx = (nAltSolutionCount % 3);
                
                DerivGEV.exploreChangeInVars(vars, varsMin, varsMax, x, y, ye, r, chiSqSumForLineSearch, idx);
                // OR temporarily allow changes that increase chiSqSum
                
                nAltSolutionCount++;
                
                // apply changes
                if (chiSqSumForLineSearch[1] < chiSqSumForLineSearch[0]) {
                    for (int k = 0; k <= varStopIdx; k++) {
                        vars[k] = vars[k] + r[k];
                        chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                    }
                    nSameSequentially = 0;
                    lastChiSqSum = chiSqSumForLineSearch[0];
                    nIter++;
                    continue;
                }
            }         
            
            if (nIter > 0) {
                // populate r with the best fitting derivatives for vars[]
                DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, 0, varStopIdx);      
            }
            
            // line search finds the fraction of the derivatives in p to apply to the GEV to reduce the chi sq sum
            float alpha = lineSearch(r, vars, varsMin, varsMax, chiSqSumForLineSearch, 0, varStopIdx);
 
            if (alpha <= eps) {
                // need 2nd deriv pre-conditioning
                break;
            }
                        
            if (!chiSqSumIsNotAcceptable(chiSqSumForLineSearch[0], chiSqSumForLineSearch[1])) {
                for (int k = 0; k <= varStopIdx; k++) {
                    float ap = alpha*r[k];
                    vars[k] = vars[k] + ap;
                    
                    log.finest("  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                }
                chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
            }
            
            if (Math.abs(lastChiSqSum - chiSqSumForLineSearch[0]) < epsChiSame) {
                nSameSequentially++;
            } else {
                nSameSequentially = 0;
            }
            lastChiSqSum = chiSqSumForLineSearch[0];
            
            nIter++;
        }
        
        float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);

        float chisqsum = calculateChiSquareSum(yGEV, WEIGHTS_DURING_CHISQSUM.ERRORS);

        GEVYFit yfit = new GEVYFit();
        yfit.setChiSqSum(chisqsum);
        yfit.setK(vars[0]);
        yfit.setSigma(vars[1]);
        yfit.setMu(vars[2]);
        yfit.setYFit(yGEV);
        yfit.setX(x);
        yfit.setXScale(xScale);
        yfit.setYScale(yScale);
        yfit.setYDataErrSq( calcYErrSquareSum() ); 
        
        return yfit;
    }
    
    /**
     * f(x_k + α*p_k) ≤ f(x_k) + c_1*α*((∇f_k)^T)*p_k  where c_1 is 0 or 1
     *     result will be applied as vars[k] = vars[k] + alpha*p[k]
     *
     * 
     * @param r suggested changes to apply to vars
     * @param vars the array holding current values for k, sigma, and mu
     * @param varsMin the minimum acceptable value of items in vars array
     * @param varsMax the maximum acceptable value of items in vars array
     * @param chiSqSum an array whose first item is the instance's current best chiSqSum corresponding to values in vars
     *         and whose second item is space to return the last value for chiSqSum computed here for the successful
     *         alpha.  an array with mixed uses is an attempt to optimize for runtime. this may change to return
     *         the alpha and chiSqSum in an object in the future.
     * @param idx0 start of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @param idx1 stop of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @return
     */
    protected float lineSearch(float[] r, float[] vars, float[] varsMin, float[] varsMax,
        float[] chiSqSum, int idx0, int idx1) throws FailedToConvergeException {
                
        float low = 0;
        float alpha = 1;
        float high = 1;
        
        float[] tmpVars = Arrays.copyOf(vars, vars.length);
                
        float[] yGEV = null;
        
        chiSqSum[1] = Float.MAX_VALUE;
                
        boolean notFound = true;
        while (notFound && (Math.abs(high - low) > eps)) {
            
            boolean failed = false;
            
            for (int i = idx0; i <= idx1; i++) {
                // solve alpha by trial and error starting with max value then use bisect:
                //     chiSqSum(... + alpha*p[i])  <=  chiSqSum + (c_1 * alpha * r[i]*p[i]).  
                //     where c_1 is 0 or 1.  c_1 being 0 doesn't make sense...
                float ap = alpha*r[i];                
                tmpVars[i] = vars[i] + ap;
                
                if ((tmpVars[i] < varsMin[i]) || (tmpVars[i] > varsMax[i])) {
                    failed = true;
                    break;
                }
                
                switch(i) {
                    case 0:
                        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, tmpVars[0], tmpVars[1], tmpVars[2]);
                        break;
                    case 1:
                        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, tmpVars[0], tmpVars[1], tmpVars[2]);
                        break;
                    case 2:
                        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, tmpVars[0], tmpVars[1], tmpVars[2]);                        
                        break;
                }
                
                float rght = chiSqSum[0] + alpha * r[i]*r[i];
                
                chiSqSum[1] = DerivGEV.chiSqSum(yGEV, y, ye);
                                
                if ((chiSqSum[1] > rght) || chiSqSumIsNotAcceptable(chiSqSum[0], chiSqSum[1])) {
                    failed = true;
                    break;
                } else {
                    log.finest("   alpha=" + alpha + "  lft=" + chiSqSum[1] + " rght=" + rght + " i=" + i + " prev chiSqSum=" + chiSqSum[0]);
                }
            }
            
            if (!failed) {
                // if in high side, let solution climb?
                return alpha;
            }
            
            // use bisect with pattern to always try high side first
            
            float bisector = (high + low)/2.f;
            if (alpha == 1) {
                alpha = 0.75f;
            } else if (alpha > bisector) {
                // try the low side
                alpha = 0.25f*(high + low)/2.f;
            } else {
                // reduce the high side since neither high nor low partition was successful
                high = bisector;
                alpha = 0.75f*(high + low)/2.f;
            }
                        
        }
        return alpha;
    }

    private boolean chiSqSumIsNotAcceptable(float bestChSqSum, float compareChSqSum) {
        
        boolean t1 = ((compareChSqSum/bestChSqSum) > 1.001);
        
        /*if (compareChSqSum < 0.002 && bestChSqSum < 0.002) {
            t1 = false;
        }*/
        
        return t1 ;
    }

    protected boolean isAcceptableMin(float[] vars, float eps, float kMin,
        float kMax, float sigmaMin, float sigmaMax, float muMin, float muMax) {
        
        if (vars[0] >= kMin && vars[0] <= kMax) {
            if (vars[1] >= sigmaMin && vars[1] <= sigmaMax) {
                if (vars[2] >= muMin && vars[2] <= muMax) {
                    return true;
                }
            }
        }
        
        // return false
        // temporarily allow it to explore all space to see if it runs away
        return true;
    }
    
    protected boolean isAcceptableMin(float[] vars, float eps, float kMin,
        float kMax, float sigmaMin, float sigmaMax) {
        
        if (vars[0] >= kMin && vars[0] <= kMax) {
            if (vars[1] >= sigmaMin && vars[1] <= sigmaMax) {
                return true;
            }
        }
        
        // return false
        // temporarily allow it to explore all space to see if it runs away
        return true;
    }

    protected void plotFit(float[] yGEV, String label) throws IOException {

        try {
            plotter.addPlot(x, y, xe, ye, x, yGEV, label);
            String filePath = plotter.writeFile2();
            System.out.println("*** filePath=" + filePath);
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }
    }
      
    protected boolean residualsAreSame(float[] rPrev, float[] r) {
        float limit = eps;//Float.MIN_VALUE;
        if ( (Math.abs(rPrev[0] - r[0]) < limit) && (Math.abs(rPrev[1] - r[1]) < limit)
            && (Math.abs(rPrev[2] - r[2]) < limit) ) {
            return true;
        }
        return false;
    }
}
