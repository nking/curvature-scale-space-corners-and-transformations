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
   of a non-linear GEV model's difference from the data.
   
   This solution uses a conjugate gradient method with a Polak-Robiere to help
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
   
              delta (x_n)^T (delta x_n - delta x_n-1)
    beta_n =  ---------------------------------------
                  delta (x_n-1)^T delta (x_n-1)
    
    In the equation above, use beta = max (0, beta_n)
    
    
  @author nichole
 */
public class PolakRibiere {

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
    
    public PolakRibiere(float[] xPoints, float[] yPoints,
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
        
        /* Minimize f(vars)
         *
         * The function f is called the objective function or cost function.
         *
         * The vector x is an n-vector of independent variables: 
         *    vars = [var_1, var_2, …, var_n]^T is a member of set Real numbers. 
         *    The variables var_1, …, var_n are often referred to as decision variables. 
         *
         * The optimization problem above can be viewed as a decision problem that involves 
         * finding the “best” vector var of the decision variables over all possible vectors in Ω. 
         * By the “best” vector we mean the one that results in the-smallest value of 
         * the objective function. 
         * 
         * This vector is called the minimizer of f over Ω. 
         * It is possible that there may be many minimizers. In this case, finding 
         * any of the minimizers will suffice.
         *     
         *  Df is the first derivative of f(vars) and is 
         *      [partial deriv f/partial deriv var_1, partial deriv f/partial deriv var_1, ...]
         *  
         *  ∇f = the gradient of f. 
         *  ∇f = (Df)^T
         *  
         *  F(vars) is the 2nd derivative of f and is sometimes called the Hessian.
         *                            | d^2f/d^2_var_1       ...    d^2f/d_var_n d_var_1
         *     F(vars) = D^2f(vars) = |         ...          ...         ...
         *                            | d^2f/d_var_1 d_var_n ...    d^2f/d^2_var_n
         *  
         *  Example:  Let f(x1, x2) = 5(x_1) + 8(x_2) + (x_1)(x_2) − (x_1)^2 − 2(x_2)^2
         *        Df(x) = (∇f(x))^T = [df/dx_1, df/dx2] = [5 + x_2 - 2x_1,  8 + x_1 - 4x_2]
         *        
         *                            | -2  1 |
         *        F(x) = D^2f(x)    = |  1 -4 |
         *
         */
        
        float kVar = (kMax + kMin)/2.f;
        float sigmaVar = (sigmaMax + sigmaMin)/2.f;
        float muVar = (muMax + muMin)/2.f;
        
        // the variables k, sigma, and mu
        float[] vars = new float[]{kVar, sigmaVar, muVar};
     
        // r is current residual.  it holds deltaK, deltaSigma, and deltaMu
        float[] r = new float[3];
        DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r);
        
        float[] rPrev = new float[r.length];
        
        // p is search direction
        float[] p = Arrays.copyOf(r, r.length);        
        
        int maxIterations = 100;
        int nIter = 0;
        float eps = 0.000001f;
        while ( (nIter < maxIterations) && isAcceptableMin(vars, eps, kMin, kMax, sigmaMin, sigmaMax, muMin, muMax)) {
            
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
                DerivGEV.derivsThatMinimizeChiSqSum(vars[2], vars[0], vars[1], x, y, ye, r); //TODO: could improve runtime by giving chiSqSum? same scale?
                for (int k = 0; k < r.length; k++) {  
                    
                    // Polak=Ribiere function
                    float beta = (float) Math.max((r[k] * (r[k] - rPrev[k]) / Math.pow(rPrev[k], 2)), 0);
                    //p_k = r_k +β_k*(p_(k−1)).
                    p[k] = r[k] + beta*p[k];
                    
                    System.out.println("r[" + k + "]=" + r[k] + "  p[" + k + "]=" + p[k] + "  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
                }                
            }
            // line search finds the fraction of the derivatives in p to apply to the GEV to reduce the chi sq sum
            float alpha = lineSearch(r, p, vars, chiSqSum);
            if (alpha == 0) {
                // need 2nd deriv pre-conditioning!!
                break;
            }
            for (int k = 0; k < r.length; k++) {
                float ap = alpha*p[k];
                vars[k] = vars[k] + ap;
 System.out.println("  vars[" + k + "]=" + vars[k] + " nIter=" + nIter);
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
     * @return
     */
    protected float lineSearch(float[] r, float[] p, float[] vars, float chiSqSum) throws FailedToConvergeException {
                
        float low = 0;
        float alpha = 1;
        float high = 1;
        
        float[] tmpVars = new float[vars.length];
                
        boolean notFound = true;
        while (notFound && (Math.abs(high - low) > 0.001)) {
            
            boolean failed = false;
            
            for (int i = 0; i < r.length; i++) {
                // solve alpha by trial and error starting with max value then use bisect:
                //     chiSqSum(... + alpha*p[i])  <=  chiSqSum + (c_1 * alpha * r[i]*p[i]).  
                //     where c_1 is 0 or 1.  c_1 being 0 doesn't make sense...
                float ap = alpha*p[i];                
                tmpVars[i] = vars[i] + ap;
                
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
                
                if ((lft > rght) || (lft > chiSqSum)) {
                    failed = true;
                    break;
                } else {
                    System.out.println("alpha=" + alpha + "  lft=" + lft + " rght=" + rght + " i=" + i + " failed=" + failed);
                }
            }
            
            if (!failed) {
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
}
