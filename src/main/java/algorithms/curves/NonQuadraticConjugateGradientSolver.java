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
   
   This solution uses a conjugate gradient method with a Polak-Robiere function to help
   determine the next stop and direction.
   
   NOTE:
   
      The preconditioning matrix (2nd derivatives) is in progress, but not yet present.
      The solution does not yet converge.
 
 
   Useful for implementing the code below was reading:
   
   
   http://en.wikipedia.org/wiki/Fletcher-Reeves#Nonlinear_conjugate_gradient
   
       1) Calculate the steepest direction: deltax_n = -d/dx f(x_n)
       2) Compute beta_n according to one of the formulas below (Polak-Robiere here)
       3) Update the conjugate direction: s_n = deltax_n + (beta_n)*(s_n-1)
       4) Perform a line search: alpha_n = arg min_{\alpha} f(x_n + alpha s_n)
       5) Update the position: x_{n+1} = x_{n} + (alpha_{n})*(s_{n})
   
   
   Polak-Robiere function:
   
               ∇(x_n)^T * (∇x_n - ∇x_n-1)
    beta_n =  ----------------------------
                 ∇(x_n-1)^T * ∇(x_n-1)   
    
    In the equation above, use beta = max (0, beta_n)
    
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
           
           The preconditioner has not yet been applied below.
           
           
  @author nichole
 */
public class NonQuadraticConjugateGradientSolver {

    protected final float[] x;
    
    // original y scaled to a max value of 1
    protected final float[] y;

    protected final float[] xe;
    protected final float[] ye;
    
    protected final float xmin; 
    protected final float xmax;
    protected final float ymin;
    protected final float ymax;

    protected float xScale = -1;
    
    // the factor by which the y axis can be multiplied to return it to the original values
    protected float yScale = -1;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;
    
    protected final GeneralizedExtremeValue gev;
    
    protected PolygonAndPointPlotter plotter;
    
    protected static float eps = 1e-8f;
    
    public NonQuadraticConjugateGradientSolver(float[] xPoints, float[] yPoints,
        float[] xErrPoints, float[] yErrPoints) {

        if (xPoints == null) {
            throw new IllegalArgumentException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalArgumentException("yPoints cannot be null");
        }
        if (xErrPoints == null) {
            throw new IllegalArgumentException("xErrPoints cannot be null");
        }
        if (yErrPoints == null) {
            throw new IllegalArgumentException("yErrPoints cannot be null");
        }

        float[] tmp = Arrays.copyOf(xPoints, xPoints.length);
        this.xScale = scaleDataTo1(tmp);
        this.x = tmp;

        tmp = Arrays.copyOf(yPoints, yPoints.length);
        this.yScale = scaleDataTo1(tmp);
        this.y = tmp;

        tmp = Arrays.copyOf(xErrPoints, xErrPoints.length);
        scaleDataTo1(tmp, xScale);
        this.xe = tmp;

        tmp = Arrays.copyOf(yErrPoints, yErrPoints.length);
        scaleDataTo1(tmp, yScale);
        this.ye = tmp;

        this.gev = new GeneralizedExtremeValue(x, y, xe, ye);
        
        xmin = x[0];
        xmax = x[x.length - 1];
        ymin = 0.0f;
        ymax = MiscMath.findMax(y);
        
        try {
            plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
        } catch (IOException e) {
            log.severe(e.getMessage());
        }
    }
    
    protected final float scaleDataTo1(float[] a) {
        float max = MiscMath.findMax(a);
        scaleDataTo1(a, max);
        return max;
    }
    protected final void scaleDataTo1(float[] a, float scale) {
        if (a == null) {
            return;
        }
        for (int i = 0; i < a.length; i++) {
            a[i] = a[i]/scale;
        }
    }

    //TODO:  when this class and GEVChiSquareMinimization are abstracted, this method should be abstract in the base class and implemented here
    
    public GEVYFit fitCurve(float kMin, float kMax, float sigmaMin, float sigmaMax,
        float muMin, float muMax) throws FailedToConvergeException {

        //TODO:  check that muMin and muMax are within bounds of x
        
        float kVar = kMin;
        float sigmaVar = sigmaMin;
        float muVar = muMin;
        
        /*
        float kVar = kMax;
        float sigmaVar = sigmaMax;
        float muVar = muMax;
        
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
        
        //TODO: derive M, the preconditioner
        
        float[] rPrev = new float[r.length];
        
        // p is search direction
        float[] p = Arrays.copyOf(r, r.length);    
        //TODO:  assign to p the inverse of M preconditioner * r
        
        float bestChiSqSum = calculateChiSquareSum(
            GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]), 
            WEIGHTS_DURING_CHISQSUM.ERRORS);
        
        int maxIterations = 200;
        int nIter = 0;
        while ( nIter < maxIterations) {
            
            if ((nIter > 0) && residualsAreSame(rPrev, r)) {
                System.out.println("residuals are same nIter=" + nIter);
                /*if (varStopIdx == (vars.length - 1)) {
                    // this holds mu fixed to current value in vars[2]
                    varStopIdx--;
                } */
            }
            
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
            float chiSqSum = calculateChiSquareSum(yGEV, WEIGHTS_DURING_CHISQSUM.ERRORS);
            
            try {
            String label = String.format("k=%4.4f <1.8>  s=%4.4f <0.85>  m=%4.4f <0.441>  n=%d  chiSqSum=%3.4f", vars[0], vars[1], vars[2], nIter, chiSqSum);
            plotFit(yGEV, label);
            } catch (IOException e) {
                System.err.println(e.getMessage());
            }
            
            if (nIter > 0) {
                for (int k = 0; k < r.length; k++) {
                    rPrev[k] = r[k];
                }
                // populate r with the best fitting derivatives for vars[]
                DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, 0, varStopIdx);
                for (int k = 0; k <= varStopIdx; k++) {  
                    
                    // Polak=Ribiere function
                    float beta =  (float) ((r[k] * (r[k] - rPrev[k])) / Math.pow(rPrev[k], 2));
                    if (beta < 0 || Float.isInfinite(beta) || Float.isNaN(beta)) {
                        beta = 0;
                    }
                    //p_k = r_k +β_k*(p_(k−1)).
                    p[k] = r[k] + beta*p[k];
                    
                    System.out.println("r[" + k + "]=" + r[k] + "  p[" + k + "]=" + p[k] + "  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                }                
            }
            // line search finds the fraction of the derivatives in p to apply to the GEV to reduce the chi sq sum
            float alpha = lineSearch(r, p, vars, varsMin, varsMax, bestChiSqSum, 0, varStopIdx);
 
            if (alpha <= eps) {
                // need 2nd deriv pre-conditioning
                break;
            }
            for (int k = 0; k <= varStopIdx; k++) {
                float ap = alpha*p[k];
                vars[k] = vars[k] + ap;
 System.out.println("  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
            }
            
            if (chiSqSum < bestChiSqSum) {
                bestChiSqSum = chiSqSum;
            }
            
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
     * find the best fitting GEV by solving for each parameter in set {k, sigma, mu} separately
     * rather then minimizing the function for suggested changes by all derivatives at once.
     * 
     * So far, this is resulting in the best fits, but is sensitive to the starting point
     * and has not been tested over a wide range of data distributions.
     * 
     * TODO: A pre-conditioner is being added to use the 2nd derivatives for better solutions or
     * faster solutions.
     * 
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param muMin
     * @param muMax
     * @return
     * @throws FailedToConvergeException
     */
    public GEVYFit fitCurveParametersSeparately(float kMin, float kMax, float sigmaMin, float sigmaMax,
        float muMin, float muMax) throws FailedToConvergeException {
        
        /*
        float kVar = kMin;
        float sigmaVar = sigmaMin;
        float muVar = muMin;
        
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
        DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, 0, varStopIdx);

        //TODO: derive M, the preconditioner
        
        float[] rPrev = new float[r.length];
        
        // p is search direction
        float[] p = Arrays.copyOf(r, r.length);    
        //TODO:  assign to p the inverse of M preconditioner * r
        
        float bestChiSqSum = calculateChiSquareSum(
            GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]), 
            WEIGHTS_DURING_CHISQSUM.ERRORS);
        
        int maxIterations = 200;
        int nIter = 0;
        while ( nIter < maxIterations) {
            
            if ((nIter > 0) && residualsAreSame(rPrev, r)) {
                System.out.println("residuals are same nIter=" + nIter);
                /*if (varStopIdx == (vars.length - 1)) {
                    // this holds mu fixed to current value in vars[2]
                    varStopIdx--;
                }*/
            }
            
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
            float chiSqSum = calculateChiSquareSum(yGEV, WEIGHTS_DURING_CHISQSUM.ERRORS);
            
            try {
            String label = String.format("k=%4.4f <1.8>  s=%4.4f <0.85>  m=%4.4f <0.441>  n=%d  chiSqSum=%3.4f", vars[0], vars[1], vars[2], nIter, chiSqSum);
            plotFit(yGEV, label);
            } catch (IOException e) {
                System.err.println(e.getMessage());
            }
            
            for (int k = 0; k < r.length; k++) {
                rPrev[k] = r[k];
            }
            
            for (int k = 0; k <= varStopIdx; k++) {
            
                if (nIter > 0) {
                    
                    // populate r with the best fitting derivatives for vars[]
                    DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r, k, k);
                    
                    // Polak=Ribiere function
                    float beta =  (float) ((r[k] * (r[k] - rPrev[k])) / Math.pow(rPrev[k], 2));
                    if (beta < 0 || Float.isInfinite(beta) || Float.isNaN(beta)) {
                        beta = 0;
                    }
                    //p_k = r_k +β_k*(p_(k−1)).
                    p[k] = r[k] + beta*p[k];
                        
                    System.out.println("r[" + k + "]=" + r[k] + "  p[" + k + "]=" + p[k] + "  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                }
                
             // line search finds the fraction of the derivatives in p to apply to the GEV to reduce the chi sq sum
                float alpha = lineSearch(r, p, vars, varsMin, varsMax, bestChiSqSum, k, k);
                if (alpha <= eps) {
                    // need 2nd deriv pre-conditioning
                    break;
                }
                float ap = alpha*p[k];
                vars[k] = vars[k] + ap;
     System.out.println("  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
            }
            
            if (chiSqSum < bestChiSqSum) {
                bestChiSqSum = chiSqSum;
            }
            
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
    
    public float calcYErrSquareSum() {
        float sum = 0.f;
        for (int i = 0; i < y.length; i++) {
            float z = yScale*ye[i];
            z *= z;
            sum += z;
        }
        return sum;
    }
  
    /**
     * f(x_k + α*p_k) ≤ f(x_k) + c_1*α*((∇f_k)^T)*p_k  where c_1 is 0 or 1
     *     result will be applied as vars[k] = vars[k] + alpha*p[k]
     *
     * 
     * @param r
     * @param p
     * @param idx0 start of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @param idx1 stop of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @return
     */
    protected float lineSearch(float[] r, float[] p, float[] vars, float[] varsMin, float[] varsMax,
        float chiSqSum, int idx0, int idx1) throws FailedToConvergeException {
                
        float low = 0;
        float alpha = 1;
        float high = 1;
        
        float[] tmpVars = new float[vars.length];
                
        boolean notFound = true;
        while (notFound && (Math.abs(high - low) > eps)) {
            
            boolean failed = false;
            
            for (int i = idx0; i <= idx1; i++) {
                // solve alpha by trial and error starting with max value then use bisect:
                //     chiSqSum(... + alpha*p[i])  <=  chiSqSum + (c_1 * alpha * r[i]*p[i]).  
                //     where c_1 is 0 or 1.  c_1 being 0 doesn't make sense...
                float ap = alpha*p[i];                
                tmpVars[i] = vars[i] + ap;
                
                if ((tmpVars[i] < varsMin[i]) || (tmpVars[i] > varsMax[i])) {
                    failed = true;
                    break;
                }
                
                float rght = chiSqSum + alpha * r[i]*p[i];
                
                float[] yGEV = null;
                switch(i) {
                    case 0:
                        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, tmpVars[0], vars[1], vars[2]);
                        break;
                    case 1:
                        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, tmpVars[0], tmpVars[1], vars[2]);
                        break;
                    case 2:
                        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, tmpVars[0], tmpVars[1], tmpVars[2]);                        
                        break;
                }
                
                float lft = calculateChiSquareSum(yGEV, WEIGHTS_DURING_CHISQSUM.ERRORS);
                
                boolean t2 = (lft > (chiSqSum + 0.01f)) && ((0.01f/chiSqSum) > 0.02f);
                
                if ((lft > rght) || ((lft/chiSqSum) > 1.5) || t2 ) {
                    failed = true;
                    break;
                } else {
                    System.out.println("alpha=" + alpha + "  lft=" + lft + " rght=" + rght + " i=" + i + " failed=" + failed);
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
            
            // TODO: could adapt the pattern to find the max of high partition when successful by raising low
            
        }
        return alpha;
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
    
    // TODO:  this needs to be in abstract base class  ******************
    
    /**
     * calculate the chi square of the yModel using the method for weights.
     * The best results for the fit are usually from using errors for the
     * fit.
     *
     * @param yModel
     * @param wdc
     * @return
     */
    public float calculateChiSquareSum(float[] yNormalizedModel, WEIGHTS_DURING_CHISQSUM wdc) {

        if (yNormalizedModel == null) {
            return Float.POSITIVE_INFINITY;
        }

        float[] w = calcWeights(wdc, yNormalizedModel);

        float chiSum = 0.f;

        for (int i = 0; i < yNormalizedModel.length; i++) {

            float z = yScale*(yNormalizedModel[i] - y[i]);
            z *= z*w[i];

            chiSum += z;
        }

        return chiSum;
    }

    /**
     * create the weight array used by the chi square sum method for the given
     * weight method.
     *
     * @param wdc
     * @param yModel
     * @return
     */
    float[] calcWeights(WEIGHTS_DURING_CHISQSUM wdc, float[] yModel) {

        float[] w = new float[yModel.length];

        if (wdc != null) {
            if (wdc.ordinal() == WEIGHTS_DURING_CHISQSUM.INVERSE_Y.ordinal()) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = 1.0f/(yScale*y[i]);
                }
            } else if (wdc.ordinal() == WEIGHTS_DURING_CHISQSUM.MODEL_Y.ordinal()) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = yModel[i]*yScale;
                }
            } else {
                boolean hasValidValues = hasValidValues(ye);
                if (hasValidValues) {
                    for (int i = 0; i < w.length; i++) {
                        w[i] = ye[i]*yScale;
                    }
                } else {
                    throw new IllegalStateException("dy has invalid values");
                }
            }
        } else {
            // defaults to using errors
            boolean hasValidValues = hasValidValues(ye);
            if (hasValidValues) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = ye[i]*yScale;
                }
            } else {
                throw new IllegalStateException("dy has invalid values");
            }
        }

        return w;
    }

    protected boolean hasValidValues(float[] a) {
        if (a == null) {
            return false;
        }
        for (int i = 0; i < a.length; i++) {
            if (Float.isNaN(a[i])) {
                return false;
            }
        }
        return true;
    }
    
    protected boolean residualsAreSame(float[] rPrev, float[] r) {
        float limit = 0.001f;
        if ( (Math.abs(rPrev[0] - r[0]) < limit) && (Math.abs(rPrev[1] - r[1]) < limit)
            && (Math.abs(rPrev[2] - r[2]) < limit) ) {
            return true;
        }
        return false;
    }
}
