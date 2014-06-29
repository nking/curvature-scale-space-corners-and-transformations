package algorithms.curves;

import java.io.IOException;
import java.util.logging.Logger;

import algorithms.curves.GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM;
import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <pre>
   Class to find the minimum difference between a GEV model function(k, sigma, mu)
   and normalized data with errors.
   
   To make a solution faster than the very robust but slow downhill simplex method,
   would like to use the partial derivatives of the GEV 
   to make an iterative solution for chi-square minimization 
   of a non-linear, non-quadratic GEV model's difference from the data.
   
   This solution uses a pre-conditioner matrix with ICU0 function to help
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
         vars = [var<sub>1</sub>, var<sub>2</sub>, &hellip;, var<sub>n</sub>]<sup>T</sup> is a member of set Real numbers. 
         The variables var<sub>1</sub>, &hellip;, var<sub>n</sub> are often referred to as decision variables. 
     
      The optimization problem above can be viewed as a decision problem that involves 
      finding the best vector var of the decision variables over all possible vectors in &#937;. 
      By the best vector we mean the one that results in the-smallest value of 
      the objective function. 
      
      This vector is called the minimizer of f over &#937;. 
      It is possible that there may be many minimizers. In this case, finding 
      any of the minimizers will suffice.
          
       Df is the first derivative of f(vars) and is 
           [partial deriv f/partial deriv var<sub>1</sub>, partial deriv f/partial deriv var<sub>1</sub>, &hellip;]
       
       &nabla;f = the gradient of f. 
       &nabla;f = (Df)<sup>T</sup>
       
       F(vars) is the 2nd derivative of f and is sometimes called the Hessian.
                                 | &#8706;<sup>2</sup>f/&#8706;<sup>2</sup>var<sub>1</sub>       &#8230;    &#8706;<sup>2</sup>f/&#8706;var<sub>n</sub>&#8706;var<sub>1</sub>
          F(vars) = D<sup>2</sup>f(vars) =   |    &#8230;           &#8230;     &#8230;
                                 | &#8706;<sup>2</sup>f/&#8706;var<sub>1</sub>&#8706;var<sub>n</sub>   &#8230;    &#8706;<sup>2</sup>f/&#8706;<sup>2</sup>var<sub>n</sub>
       
       Example:  Let f(x1, x2) = 5(x<sub>1</sub>) + 8(x<sub>2</sub>) + (x<sub>1</sub>)(x<sub>2</sub>) &#8722; (x<sub>1</sub>)<sup>2</sup> &#8722; 2(x<sub>2</sub>)<sup>2</sup>
             Df(x) = (&nabla;f(x))<sup>T</sup> = [&#8706;f/&#8706;x<sub>1</sub>, &#8706;f/&#8706;x<sup>2</sup>] = [5 + x<sub>2</sub> &#8722; 2x<sub>1</sub>,  8 + x<sub>1</sub> &#8722; 4x<sub>2</sub>]
             
                                | -2  1 |
             F(x) = D<sup>2</sup>f(x)    = |  1 -4 |
             
       --------------------------------------------------------------------------------
       
       f(k, sigma, mu) = the GEV function
       
       Df = [&#8706;/&#8706;k, &#8706;/&#8706;sigma, &#8706;/&#8706;mu]
       
                     |   &#8706;/&#8706;k   |
       &nabla;f = (Df)<sup>T</sup> =   | &#8706;/&#8706;sigma |
                     |   &#8706;/&#8706;mu  |
            
                                              | &#8706;<sup>2</sup>/&#8706;k&#8706;k     &#8706;<sup>2</sup>/&#8706;sigma&#8706;k        &#8706;<sup>2</sup>/&#8706;mu&#8706;k     |
       F(k, sigma, mu) = D<sup>2</sup>f(k, sigma, mu) =   | &#8706;<sup>2</sup>/&#8706;k&#8706;sigma   &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma    &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma |
                                              | &#8706;<sup>2</sup>/&#8706;k&#8706;mu      &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu       &#8706;<sup>2</sup>/&#8706;mu&#8706;mu    |
            
       M = D<sup>2</sup>f(k, sigma, mu)
       
                                | &#8706;<sup>2</sup>/&#8706;k&#8706;k        &#8706;<sup>2</sup>/&#8706;k&#8706;sigma       &#8706;<sup>2</sup>/&#8706;k&#8706;mu     |     |   &#8706;/&#8706;k   |
       M^(-1) * &nabla;f = M<sup>T</sup> * &nabla;f =  | &#8706;<sup>2</sup>/&#8706;sigma&#8706;k    &#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma   &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu |  *  | &#8706;/&#8706;sigma |
                                | &#8706;<sup>2</sup>/&#8706;mu&#8706;k       &#8706;<sup>2</sup>/&#8706;mu&#8706;sigma      &#8706;<sup>2</sup>/&#8706;mu&#8706;mu    |     |   &#8706;/&#8706;mu  |
                                
                                | (&#8706;<sup>2</sup>/&#8706;k&#8706;k * (&#8706;/&#8706;k) +  (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) * (&#8706;/&#8706;sigma) + (&#8706;<sup>2</sup>/&#8706;k&#8706;mu) * (&#8706;/&#8706;mu)          |
                              = | (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k) * (&#8706;/&#8706;k) + (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) * (&#8706;/&#8706;sigma) + (&#8706;<sup>2</sup>/&#8706;sigma&#8706;mu)*(&#8706;/&#8706;mu) |
                                | (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*(&#8706;/&#8706;k) + (&#8706;<sup>2</sup>/&#8706;mu&#8706;sigma)*(&#8706;/&#8706;sigma) + (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu)*(&#8706;/&#8706;mu)              |
                                
       Can use Incomplete Cholesky factorization with fill 0 on M<sup>(-1)</sup> to create the pre-conditioning matrix.
       
           http://netlib.org/linalg/html_templates/node64.html#figdilu
           
           Let S be the non-zero set ({i,j} : a<sub>ij</sub> != 0 )
           
            for i = 1, 2, &#8230;
                set d<sub>ii</sub> = a<sub>ii</sub>
            for i = 1, 2, &#8230;
                set d<sub>ii</sub> = 1/d<sub>ii</sub>
                for j = i + 1, i + 2, &#8230;                       
                if (i, j) in set S and (j, i) in set S then
                    set d<sub>jj</sub> = d<sub>jj</sub> - a<sub>ji</sub> * d<sub>ii</sub> * a<sub>ij</sub>
            
            set d(1,1) = (&#8706;<sup>2</sup>/&#8706;k&#8706;k)
            set d(2,2) = (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma)
            set d(3,3) = (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu)
            
            i = 1:
                    set d(1,1) = 1./d(1,1) = 1./(&#8706;<sup>2</sup>/&#8706;k&#8706;k)
                
                    | (1,1) (1,2) (1,3) | = | 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k)   
                d = | (2,1) (2,2) (2,3) | = |                     (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma)  
                    | (3,1) (3,2) (3,3) | = |                                             (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) |
                  
                j=2:
                    set d(2,2) = (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma)
                j=3:
                    set d(3,3) = (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) - (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;mu)
                    
            i = 2:
                    set d(2,2) = 1./d(2,2) = ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
                    
                    | (1,1) (1,2) (1,3) | = | 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k)   
                d = | (2,1) (2,2) (2,3) | = |                latest d(2,2)  
                    | (3,1) (3,2) (3,3) | = |                                             (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) |
                    
                j=3:
                    set d(3,3) = d(3,3) - a(3,2) * d(2,2) * a(2,3)
                               =
                                 ( (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) - (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;mu) )
                                 -
                                 (
                                    (&#8706;<sup>2</sup>/&#8706;mu&#8706;sigma)
                                    *
                                    ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
                                    *
                                    &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu
                                 )
                                  
            i = 3:
                    set d(3,3) = 1./d(3,3)
                    
                        = 1./(
                                 ( (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) - (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;mu) )
                                 -
                                 (
                                    (&#8706;<sup>2</sup>/&#8706;mu&#8706;sigma)
                                    *
                                    ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
                                    *
                                    &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu
                                 )
                             ) 
                             
                             
    Then using the ICU0 matrix as preconditioner:
           
                                        | d(1,1)   0        0      |     |   &#8706;/&#8706;k   |
        (M_icuo)<sup>(-1)</sup> * &nabla;f = M<sup>T</sup> * &nabla;f =     | 0        d(2,2)   0      |  *  | &#8706;/&#8706;sigma |
                                        | 0        0        d(3,3) |     |   &#8706;/&#8706;mu  |
                                    
                                        | d(1,1) * (&#8706;/&#8706;k)     |
                                      = | d(2,2) * (&#8706;/&#8706;sigma) |
                                        | d(3,3) * (&#8706;/&#8706;mu)    |

             where d(1,1) is 1./(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k)
                   d(2,2) is ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
                   d(3,3) is 1./(
                                     ( (&#8706;<sup>2</sup>/&#8706;mu&#8706;mu) - (&#8706;<sup>2</sup>/&#8706;mu&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;mu) )
                                     -
                                     (
                                        (&#8706;<sup>2</sup>/&#8706;mu&#8706;sigma)
                                        *
                                        ( 1./ ( (&#8706;<sup>2</sup>/&#8706;sigma&#8706;sigma) - (&#8706;<sup>2</sup>/&#8706;sigma&#8706;k)*( 1/(&#8706;<sup>2</sup>/&#8706;k&#8706;k, that is &#8706;<sup>2</sup>/&#8706;k&#8706;k) ) * (&#8706;<sup>2</sup>/&#8706;k&#8706;sigma) ) )
                                        *
                                        &#8706;<sup>2</sup>/&#8706;sigma&#8706;mu
                                     )
                                 ) 
                                 
       Note that below in the code, r is &nabla;f.
  </pre>
  @author nichole
 */
public class NonQuadraticConjugateGradientSolver extends AbstractCurveFitter {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected int maxIterations = 200;
    
    protected PolygonAndPointPlotter plotter;
    
    protected static float eps = 1e-8f;
    
    protected static float convergedEps = 0.00001f;
    
    protected boolean calcStepParamsSeparately = true;
        
    public NonQuadraticConjugateGradientSolver(float[] xPoints, float[] yPoints,
        float[] xErrPoints, float[] yErrPoints) {

        super(xPoints, yPoints, xErrPoints, yErrPoints);
    }    

    public void setMaximumNumberOfIterations(int maxNumber) {
        this.maxIterations = maxNumber;
    }
    
    @Override
    public void setDebug(boolean doUseDebug) {
        
        super.setDebug(doUseDebug);
        
        if (debug) {
            try {
                plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
            } catch (IOException e) {
                log.severe(e.getMessage());
            }
        }
    }
    
    /**
     * fit the x and y data using default ranges for the parameters k, sigma, and mu.
     * The generated curves use parameter step sizes determined from using
     * the k, sigma, and mu derivatives separately and keeping the best result.
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) 
        throws FailedToConvergeException, IOException {
        
        calcStepParamsSeparately = true;
        
        return fitCurveKGreaterThanZero();
    }
    
    /**
     * fit the x and y data using default ranges for the parameters k, sigma, and mu.
     * The generated curves use parameter step sizes determined from using all
     * of the k, sigma, and mu derivatives at once.
     */
    public GEVYFit fitCurveKGreaterThanZeroAllAtOnce(WEIGHTS_DURING_CHISQSUM weightMethod) 
        throws FailedToConvergeException, IOException {
        
        calcStepParamsSeparately = false;
        
        return fitCurveKGreaterThanZero();
    }
    
    protected GEVYFit fitCurveKGreaterThanZero() throws FailedToConvergeException, IOException {
        
        String filePath = ResourceFinder.findFileInTestResources("sim_curve_params_01.txt");
        
        File f = new File(filePath);
        
        BufferedReader reader = null;
        FileReader fr = null;
        
        List<Float> kParams = new ArrayList<Float>();
        List<Float> sParams = new ArrayList<Float>();
        List<Float> mParams = new ArrayList<Float>();
        
        try {
            fr = new FileReader(f);
            reader = new BufferedReader(fr);
            String line = reader.readLine();
            line = reader.readLine();
            while (line != null) {
                
                String[] params = line.split("\\s+");
                
                kParams.add(Float.valueOf(params[0]));
                sParams.add(Float.valueOf(params[1]));
                mParams.add(Float.valueOf(params[2]));
                                
                line = reader.readLine();
            }
            
        } finally {
            if (fr != null) {
                fr.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
        
        int n = kParams.size();
        int bestFitIndex = -1;
        GEVYFit bestFit = null;
        
        for (int i = 0; i < n - 1; i++) {
            
            float kMin = kParams.get(i).floatValue();
            float kMax = kParams.get(i + 1).floatValue();
            float sMin = sParams.get(i).floatValue();
            float sMax = sParams.get(i + 1).floatValue();
            float mMin = mParams.get(i).floatValue();
            float mMax = mParams.get(i + 1).floatValue();
            
            if (sMax < sMin) {
                sMax = sMin;
            }
            if (mMax < mMin) {
                mMax = mMin;
            }
            
            GEVYFit yFit = fitCurveParameters(kMin, kMax, sMin, sMax, mMin, mMax);
            
            if (yFit != null) {
                if (bestFit == null) {
                    bestFit = yFit;
                    bestFitIndex = i;
                } else if (yFit.chiSqStatistic < bestFit.chiSqStatistic) {
                    bestFit = yFit;
                    bestFitIndex = i;
                }
            }
        }
        
        return bestFit;
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
    public GEVYFit fitCurveParametersSeparately(float kMin, float kMax, 
        float sigmaMin, float sigmaMax, float muMin, float muMax) throws 
        FailedToConvergeException {
        
        calcStepParamsSeparately = true;
       
        return fitCurveParameters(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
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
    public GEVYFit fitCurveParametersAllAtOnce(float kMin, float kMax, 
        float sigmaMin, float sigmaMax, float muMin, float muMax) throws 
        FailedToConvergeException {

        calcStepParamsSeparately = false;
       
        return fitCurveParameters(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
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
    public GEVYFit fitCurveParameters(float kMin, float kMax, float sigmaMin, float sigmaMax,
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
        
        float kVar = (kMax + kMin)/2.f;
        float sigmaVar = (sigmaMax + sigmaMin)/2.f;
        float muVar = (muMax + muMin)/2.f;
        
        // the variables k, sigma, and mu
        float[] vars    = new float[]{kVar, sigmaVar, muVar};
        float[] varsMin = new float[]{kMin, sigmaMin, muMin};
        float[] varsMax = new float[]{kMax, sigmaMax, muMax};
     
        float[] prevVars = new float[]{kVar, sigmaVar, muVar};
        
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
            
            if (nSameSequentially > 0) {                            
                break;
            }
         
            if (nIter > 0) {
                
                if (calcStepParamsSeparately) {
                    
                    for (int k = 0; k <= varStopIdx; k++) {
                        // populate r with the best fitting derivatives for vars[]
                        DerivGEV.derivsThatMinimizeChiSqSum(vars, varsMin, varsMax, chiSqSumForLineSearch,
                            x, y, ye, r, k, k);

                        if (r[k] != 0) {
                            vars[k] += r[k];
                            chiSqSumForLineSearch[0] = chiSqSumForLineSearch[1];
                        }

                        log.finest("   ->r[" + k + "]=" + r[k]  + "  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                    }
                    
                } else {
                
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
            }
            
            if (debug) {
                float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
                float chiSqSum = DerivGEV.chiSqSum(yGEV, y, ye);
                try {
                    String label = String.format("k=%4.4f <*>  s=%4.4f <*>  m=%4.4f <*>  n=%d  chi=%4.8f yscl=%.0f",
                        vars[0], vars[1], vars[2], nIter, chiSqSum, this.yScale);
                    plotFit(yGEV, label);
                } catch (IOException e) {
                    System.err.println(e.getMessage());
                }
            }
            
            if ((nIter > 0) && 
                ((Math.abs(lastChiSqSum - chiSqSumForLineSearch[0]) < epsChiSame) 
                || (Arrays.equals(prevVars, vars)))
                ) {
                nSameSequentially++;
            } else {
                nSameSequentially = 0;
            }
            lastChiSqSum = chiSqSumForLineSearch[0];
            
            System.arraycopy(vars, 0, prevVars, 0, vars.length);
            
            nIter++;
        }
        
        bestYFit = compareFits(bestYFit, vars, chiSqSumForLineSearch);
        
        return bestYFit;
    }
  
    protected GEVYFit compareFits(GEVYFit bestYFit, float[] vars, float[] chiSqSumForLineSearch) {
        
        if ((bestYFit == null) || (chiSqSumForLineSearch[0] < bestYFit.getChiSqSum())) {
            
            float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, vars[0], vars[1], vars[2]);
    
            float degreesOfFreedom = yGEV.length - 3 - 1;

            float chiSqStatistic = chiSqSumForLineSearch[0]/degreesOfFreedom;
            
            bestYFit = new GEVYFit();
            bestYFit.setChiSqSum(chiSqSumForLineSearch[0]);
            bestYFit.setChiSqStatistic(chiSqStatistic);
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

    protected void plotFit(float[] yGEV, String label) throws IOException {

        if (plotter == null) {
            throw new IllegalStateException("set debug in order to use plotFit");
        }
        
        try {
            plotter.addPlot(x, y, xe, ye, x, yGEV, label);
            String filePath = plotter.writeFile2();
            log.finest("*** filePath=" + filePath);
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }
    }
    
}
