package algorithms.curves;

import java.io.IOException;
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
   determine the next stop and direction. 
 
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
     * get an array of kMin, kMax, sigmaMin, sigmaMax, muMin, muMax default ranges that
     * are tuned for distributions in which the peak of the GEV is usually < 0.4 in a normalized
     * range of x values.  If the fits are poor and the peak of the data is located near
     * 0.5 or later, you'll probably need to widen the shape parameter sigma range and use
     * the method fitCurveParametersSeparately directly.
     * 
     * @return
     */
    public float[] getMinMaxRanges() {
        
        int yMaxIndex = MiscMath.findYMaxIndex(y);
        float xAtYMax = x[yMaxIndex];
        
        float kMin = 0.0001f;
        float kMax = 2.0f;
        float sigmaMin = 0.025f;
        float sigmaMax = 3.0f;
        float muMin = xAtYMax/2.5f;
        float muMax = xAtYMax * 3.0f;
        if (muMax > 1.0) {
            muMax = xAtYMax * 1.5f;
        }
        
        if (yScale > 100) {
            //kMax = 1.0f;
            /*if (yScale > 200) {
                kMin = 0.2f;
            }*/
        }
        
        return new float[] {kMin, kMax, sigmaMin, sigmaMax, muMin, muMax};
    }
    
    /**
     * fit the x and y data using built in default ranges for the parameters k, sigma, and mu.
     * 
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {
        
        float[] minMaxes = getMinMaxRanges();
        
        return fitCurveParametersSeparately(minMaxes[0], minMaxes[1], minMaxes[2], minMaxes[3], minMaxes[4], minMaxes[5]);
    }
    
    /**
     * fit the x and y data using built in default ranges for the parameters k, sigma, and mu.
     * 
     */
    public GEVYFit fitCurveKGreaterThanZeroAllAtOnce(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {
        
        float[] minMaxes = getMinMaxRanges();
        
        return fitCurveParametersAllAtOnce(minMaxes[0], minMaxes[1], minMaxes[2], minMaxes[3], minMaxes[4], minMaxes[5]);
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
        
        // for most solutions, choosing middle of range is a good starting point,
        //   but for those that do not improved after several iterations,
        //   they should be restarted with the min variables.
        
        boolean hasTriedMinStarts = false;
        
        /*
        float kVar = kMin;
        float sigmaVar = sigmaMin;
        float muVar = muMin;
        */
        
        boolean hasTriedMaxStarts = false;
        boolean hasTriedAltSteps = false;
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
        
        // chiSqSumForLineSearch[0] holds current best chiSqSum for the last change in vars
        // chiSqSumForLineSearch[1] holds the return value from DerivGEV if step > 0 was admissable
        float[] chiSqSumForLineSearch = new float[2];
        
        chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);
        
        // r is current residual.  it holds deltaK, deltaSigma, and deltaMu
        float[] r = new float[3];
        
        int nSameSequentially = 0;
        float epsChiSame = 1e-5f;
        float lastChiSqSum = Float.MAX_VALUE;
        
        GEVYFit bestYFit = new GEVYFit();
        
        int nIter = 0;
                
        while ((nIter < maxIterations) /*&& (chiSqSumForLineSearch[0] > convergedEps)*/) {
            
            /*if ((nIter > 1) && residualsAreSame(rPrev, r)) {
                break;
            }*/
            
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
            
            float chiSqSum = DerivGEV.chiSqSum(yGEV, y, ye);
            
            if (debug) {
                try {
                    String label = String.format("k=%4.4f <*>  s=%4.4f <*>  m=%4.4f <*>  n=%d  chi=%4.8f yscl=%.0f",
                        vars[0], vars[1], vars[2], nIter, chiSqSum, this.yScale);
                    plotFit(yGEV, label);
                } catch (IOException e) {
                    System.err.println(e.getMessage());
                }
            }
            
            // if solution has stalled, start with the min start point, else the max start points
            if ((nIter > 3) && (nSameSequentially > 2)) {
                
                if (!hasTriedAltSteps) {
                    
                    bestYFit = compareFits(bestYFit, vars, chiSqSumForLineSearch);
                
                    for (int k = 0; k < vars.length; k++) {
                        DerivGEV.exploreChangeInVars(vars, varsMin, varsMax, x, y, ye, r, chiSqSumForLineSearch, k);
                        if (r[k] != 0) {
                            vars[k] += r[k];
                            chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                        }
                    }
                    
                    bestYFit = compareFits(bestYFit, vars, chiSqSumForLineSearch);
                    
                    hasTriedAltSteps = true;
                    nSameSequentially = 0;
                    lastChiSqSum = Float.MAX_VALUE;
                    
                } else if (true && !hasTriedMinStarts && (chiSqSumForLineSearch[0] > 0.001f)) {
                    
                    bestYFit = compareFits(bestYFit, vars, chiSqSumForLineSearch);
                    
                    kVar = kMin;
                    sigmaVar = sigmaMin;
                    muVar = muMin;
                    
                    vars[0] = kVar;
                    vars[1] = sigmaVar;
                    vars[2] = muVar;
                    nSameSequentially = 0;
                    lastChiSqSum = Float.MAX_VALUE;

                    chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);

                    hasTriedMinStarts = true;
                    hasTriedAltSteps = false;
                    
                    nIter = 0;

                } else {
                    // let end or let continue
                    nIter = maxIterations;
                }
            }
            
            for (int k = 0; k <= varStopIdx; k++) {
            
                if (nIter > 0) {
                    
                    // populate r with the best fitting derivatives for vars[]
                    DerivGEV.derivsThatMinimizeChiSqSum(vars, varsMin, varsMax, chiSqSumForLineSearch,
                        x, y, ye, r, k, k);
                    
                    if (r[k] != 0) {
                        vars[k] += r[k];
                        chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                    }
                    
                    log.finest("   ->r[" + k + "]=" + r[k]  + "  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                }
            }
            
            if (Math.abs(lastChiSqSum - chiSqSumForLineSearch[0]) < epsChiSame) {
                nSameSequentially++;
            } else {
                nSameSequentially = 0;
            }
            lastChiSqSum = chiSqSumForLineSearch[0];
            
            nIter++;
        }
        
        bestYFit = compareFits(bestYFit, vars, chiSqSumForLineSearch);
        
        if (debug) {
            try {
                String label = String.format("k=%4.4f <*>  s=%4.4f <*>  m=%4.4f <*>  n=%d  chi=%4.8f  yscl=%.0f",
                    vars[0], vars[1], vars[2], nIter, bestYFit.getChiSqSum(), this.yScale);
                plotFit(bestYFit.getYFit(), label);
            } catch (IOException e) {
                System.err.println(e.getMessage());
            }
        }
        
        return bestYFit;
    }
  
    protected GEVYFit compareFits(GEVYFit bestYFit, float[] vars, float[] chiSqSumForLineSearch) {
        
        if ((bestYFit == null) || (chiSqSumForLineSearch[0] < bestYFit.getChiSqSum())) {
            
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
    
            bestYFit = new GEVYFit();
            bestYFit.setChiSqSum(chiSqSumForLineSearch[0]);
            bestYFit.setK(vars[0]);
            bestYFit.setSigma(vars[1]);
            bestYFit.setMu(vars[2]);
            bestYFit.setYFit(yGEV);
            bestYFit.setX(x);
            bestYFit.setXScale(xScale);
            bestYFit.setYScale(yScale);
        }
        
        return bestYFit;
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
    public GEVYFit fitCurveParametersAllAtOnce(float kMin, float kMax, float sigmaMin, float sigmaMax,
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
        
        boolean hasTriedMinStarts = false;
        /*
        float kVar = kMin;
        float sigmaVar = sigmaMin;
        float muVar = muMin;
        */
        
        boolean hasTriedMaxStarts = false;
        /*
        float kVar = kMax;
        float sigmaVar = sigmaMax;
        float muVar = muMax;
        */
        float kVar = (kMax + kMin)/2.f;
        float sigmaVar = (sigmaMax + sigmaMin)/2.f;
        float muVar = (muMax + muMin)/2.f;
        
        // the variables k, sigma, and mu
        float[] vars = new float[]{kVar, sigmaVar, muVar};
        float[] varsMin = new float[]{kMin, sigmaMin, muMin};
        float[] varsMax = new float[]{kMax, sigmaMax, muMax};
     
        int varStopIdx = vars.length - 1;
        
        // r is current residual.  it holds deltaK, deltaSigma, and deltaMu
        float[] r = new float[3];
        
        // chiSqSumForLineSearch[0] holds current best chiSqSum for the last change in vars
        // chiSqSumForLineSearch[1] holds the return value from DerivGEV if step > 0 was admissable
        float[] chiSqSumForLineSearch = new float[2];
        
        chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);
        
        int nSameSequentially = 0;
        float epsChiSame = 1e-5f;
        float lastChiSqSum = Float.MAX_VALUE;
        
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
                       "k=%4.4f <1.8>  s=%4.4f <0.85>  m=%4.4f <0.441>  n=%d  chi=%4.8f  yscl=%.0f",
                        vars[0], vars[1], vars[2], nIter, chiSqSum, this.yScale);
                    plotFit(yGEV, label);
                } catch (IOException e) {
                    System.err.println(e.getMessage());
                }
            }
            
            // if solution has stalled, restart with min values, else max values, else exit loop
            if ((nIter > 3) && (nSameSequentially > 2)) {
                
                if (!hasTriedMinStarts && (chiSqSum > 0.01f)) {
                    kVar = kMin;
                    sigmaVar = sigmaMin;
                    muVar = muMin;
                    hasTriedMinStarts = true;
                    vars[0] = kVar;
                    vars[1] = sigmaVar;
                    vars[2] = muVar;
                    nSameSequentially = 0;
                    lastChiSqSum = Float.MAX_VALUE;
                    
                    chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);
                    
                    nIter = 0;
                
                } else if (!hasTriedMaxStarts && (chiSqSum > 0.01f)) {
                    kVar = kMax;
                    sigmaVar = sigmaMax;
                    muVar = muMax;
                    hasTriedMaxStarts = true;
                    vars[0] = kVar;
                    vars[1] = sigmaVar;
                    vars[2] = muVar;
                    nSameSequentially = 0;
                    lastChiSqSum = Float.MAX_VALUE;

                    chiSqSumForLineSearch[0] = DerivGEV.chiSqSum(vars[0], vars[1], vars[2], x, y, ye);
                    
                    nIter = 0;
                
                } else {
                    // let end or let continue
                    nIter = maxIterations;
                }
            }         
            
            if (nIter > 0) {
                
                // populate r with the best fitting derivatives for vars[]
                DerivGEV.derivsThatMinimizeChiSqSum(vars, varsMin, varsMax, chiSqSumForLineSearch,
                    x, y, ye, r, 0, varStopIdx); 
                
                for (int k = 0; k <= varStopIdx; k++) {
                    if (r[k] != 0) {
                        vars[k] += r[k];
                        chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                    }
                }
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
            log.finest("*** filePath=" + filePath);
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
