package algorithms.curves;

import java.util.Arrays;

import algorithms.misc.MiscMath;

public class DerivGEV {

    /*
     *                          (   (      ( x-mu))-(1/k))
     *                          (-1*(1 + k*(-----))      )
     *                 1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
     * y = y_const * ----- * exp                           * (1 + k*(-----))
     *               sigma                                   (      (sigma))
     *
     * mu is  the location parameter
     * sigma is the scale parameter and is > 0
     * k is the shape parameter
     * 
     * 
     * Components needed in the derivatives:
     * 
     *   Let z = (1 + k*( (x-mu)/sigma )
     *   
     *   then deriv of z w.r.t x is
     *   dzdx = k/sigma
     *   
     *   dzdk = (x-mu)/sigma
     *   
     *   dzdsigma =  -1 * k * (x-mu) * (sigma^-2) 
     *   
     *   dzdmu = -1*k/sigma
     *   
     *   
     *   deriv of -z^(-1/k) w.r.t. x is
     *      (1/k) * z^(-1 - (1/k)) * dzdx
     *   
     *   deriv of -z^(-1/k) w.r.t. k 
     *       use pattern: f(x) = u(x) ^(v(x))
                 ln ( f(x) ) = ln ( u(x) ^(v(x)) )
              
                 df(x)dx            du(x)dx
                 -------  =  v(x) * -------  +  dv(x)dx * ln(u(x))
                  f(x)               u(x)
                 
                                           du(x)dx
                 df(x)dx = f(x) * ( v(x) * -------  +  dv(x)dx * ln(u(x)) )
                                            u(x)
                
     *           df(k)dk = -z^(-1/k) * ( (-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
     *           
     *       *Because the ln( negative number) is not defined (though can be approximated with taylor series)
     *        need a different method of creating the derivative of a number with a power that is a function of x.
     *        This derivative is only needed with the derivative with respect to k.
     *        
              One solution is to calculate the GEV with slightly different values of k near the given k to get delta GEV/deltaK.
     *
     *   
     *   deriv of -z^(-1/k) w.r.t. sigma is
     *       (1/k) * z^(-1 - (1/k)) * dzdsigma 
     *   
     *   deriv of -z^(-1/k) w.r.t. mu is 
     *       (1/k) * z^(-1 - (1/k)) * dzdmu 
     *   
     *   
     *   Let f1 = the first exponential in y
     *          = exp(-1*(z^(-1/k)))
     *   
     *   df1dx     = f1 * (1/k) * z^(-1 - (1/k)) * dzdx
     *   
     *                                      
     *   df1dk     = f1 * -z^(-1/k) * ( (-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
     *                                         
     *                                  
     *   df1dsigma = f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
     *                                             
     *   
     *   df1dmu    = f1 * (1/k) * z^(-1 - (1/k)) * dzdmu
     *   
     *   
     *   
     *   Let f2 = z^(-1-(1/k))
     *   
     *   df2dx = (-1-(1/k)) * z^(-2-(1/k)) * dzdx 
     *   
     *                                  
     *   df2dk: u(k) = z and v(k) = (-1-(1/k))
     *      using 
                                        du(k)dk
                 df2dk = f2 * ( v(k) * -------  +  dv(k)dk * ln(u(k)) )
                                         u(k)
                       = f2 * ( (-1-(1/k)) * dzdk/z + (1/k^2) * ln(z) )
     *       
     *   df2dsigma = (-1-(1/k)) * z^(-2-(1/k)) * dzdsigma
     *   
     *   df2dmu = (-1-(1/k)) * z^(-2-(1/k)) * dzdmu
     *   
     *   
     *   Then putting it all together:
     *   
     *          yconst
     *   yfit = ------ * f1 * f2
     *           sigma
     *           
     *           yconst
     *   dydx =  ------ * ( f1 * df2dx + f2 * df1dx )
     *           sigma
     *          
     *           yconst
     *   dydk =  ------ * ( f1 * df2dk + f2 * df1dk )
     *           sigma
     *           
     *              yconst
     *   dydsigma = ------ * ( f1 * df2dsigma + f2 * df1dsigma )
     *              sigma
     *   
     *              yconst
     *   dydmu    = ------ * ( f1 * df2dmu + f2 * df1dmu )
     *              sigma
     */
    
    /**
     * calculate the derivative of the GEV w.r.t. x
     * 
     * the runtime cost is 4 transcendental functions.
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    public static double derivWRTX(float yConst, float mu, float k, float sigma, float x) {

        double z = 1. + k *( (x-mu)/sigma );

        double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        double f2 = Math.pow(z, (-1. - (1./k)) );
      
        double dzdx = k/sigma;

        //(-1-(1/k)) * z^(-2-(1/k)) * dzdx
        double df2dx = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdx;
        
        double df1dx = f1 * (1./k) * Math.pow(z, (-1. - (1./k))) * dzdx;
        
        double dydx = (yConst/sigma) * ( f1 * df2dx + f2 * df1dx );
        
        return dydx;
    }
    
    /**
     * calculate the derivative of the GEV w.r.t. k
     * 
     * the runtime cost is 5 transcendental functions.
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    /*
    public static double derivWRTK(float yConst, float mu, float k, float sigma, float x) {
        
        double z = 1. + k *( (x-mu)/sigma );
      
        double f1, f2;
        
        if (z < 0) {
            f1 = Math.exp( Math.pow(-1.*z, -1./k) );
            f2 = -1.* Math.pow(-1.*z, (-1. - (1./k)) );
        } else {
            f1 = Math.exp( -1. * Math.pow(z, -1./k) );
            f2 = Math.pow(z, (-1. - (1./k)) );
        }
        
        double dzdk = (x-mu)/sigma;
        
        double df2dk, df1dk;

        if (z == 0) {
            use alt method for df2dk.  a brute force delta using neighboring points?
            use alt method for df1dk.  a brute force delta using neighboring points?
        } else if (z > 0) {
            use alt method for df1dk.  a brute force delta using neighboring points?
            //      f2 * z^(-1-(1/k)) * ( (-1-(1/k)) * dzdk/z  +  (1/k^2) * ln(z) )
            df2dk = f2 * Math.pow(z, (-1. - (1./k))) * ( (-1-(1/k)) * (dzdk/z) + (1./k*k)*Math.log(z));
        } else {
            use alt method for df2dk.  a brute force delta using neighboring points?
            //      f1 *    -z^(-1/k)              * ( (-1/k) * (dzdk/z)   +  (1/k^2) * ln( -z ) )
            df1dk = f1 *  Math.pow(-1.*z, (-1./k)) * ( (-1./k) * (dzdk/z)  +  (1./k*k) * Math.log(-1.*z));
        }
                
        double dydk = (yConst/sigma) * ( f1 * df2dk + f2 * df1dk );
        
        return dydk;
    }*/
    
    public static Double derivWRTK(float yConst, float mu, float k, float sigma, float x) {

        return estimateDerivUsingDeltaK(mu, k, sigma, x);
    }
    
    /**
     * calculate d/dk of GEV using the difference between GEVs given minor changes in k
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    static Double estimateDerivUsingDeltaK(float mu, float k, float sigma, float x) {
        
        float deltaK = 0.001f*k;
        
        Double d0 = GeneralizedExtremeValue.generateYGEV(x, (k - deltaK), sigma, mu);
        
        Double d1 = GeneralizedExtremeValue.generateYGEV(x, (k), sigma, mu);
        
        Double d2 = GeneralizedExtremeValue.generateYGEV(x, (k + deltaK), sigma, mu);
        
        Double d = estimateDerivUsingDelta(d0, d1, d2, deltaK);
        
        return d;
    }
    
    protected static Double estimateDerivUsingDelta(Double d0, Double d1, Double d2, double delta) {
        
        if (d0 != null && d1 != null && d2 != null) {
            
            double delta0 = d1.doubleValue() - d0.doubleValue();
            
            double delta1 = d2.doubleValue() - d1.doubleValue();
            
            double d = (delta0 + delta1)/2.;
            
            return (d/delta);
        
        } else if (d1 != null && d2 != null) {
                        
            double d = d2.doubleValue() - d1.doubleValue();
                        
            return (d/delta);
            
        } else if (d0 != null && d1 != null) {
            
            double d = d1.doubleValue() - d0.doubleValue();
            
            return (d/delta);
            
        } else if (d0 != null && d2 != null) {
            
            double d = d2.doubleValue() - d0.doubleValue();
            
            return (d/delta);
            
        } else {
            
            return null;
        }
    }
    
    /**
     * calculate the derivative of the GEV w.r.t. sigma
     * 
     * the runtime cost is 5 transcendental functions.
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    public static double derivWRTSigma(float yConst, float mu, float k, float sigma, float x) {
        
        double z = 1. + k *( (x-mu)/sigma );
        
        double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        double dzdsigma = -1. * k * (x-mu) * Math.pow(sigma, -2.);
        
        //(-1-(1/k)) * z^(-2-(1/k)) * dzdsigma
        double df2dsigma = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdsigma;
        
        double df1dsigma =  f1 * (1./k) * dzdsigma;
        if (z < 0) {
            df1dsigma *= -1. * Math.pow(-1.*z, -1. - (1./k));
        } else {
            df1dsigma *= Math.pow(z, -1. - (1./k));
        }
        
        double dydSigma = (yConst/sigma) * ( f1 * df2dsigma + f2 * df1dsigma );
        
        if (Double.isNaN(dydSigma)) {
            return estimateDerivUsingDeltaSigma(mu, k, sigma, x);
        }
        
        return dydSigma;
    }
    
    /**
     * calculate d/dk of GEV using the difference between GEVs given minor changes in k
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    static Double estimateDerivUsingDeltaSigma(float mu, float k, float sigma, float x) {
        
        float deltaSigma = 0.0001f*sigma;
        
        Double d0 = GeneralizedExtremeValue.generateYGEV(x, k, (sigma - deltaSigma), mu);
        
        Double d1 = GeneralizedExtremeValue.generateYGEV(x, k, sigma, mu);
        
        Double d2 = GeneralizedExtremeValue.generateYGEV(x, k, (sigma + deltaSigma), mu);
        
        return estimateDerivUsingDelta(d0, d1, d2, deltaSigma);
    }

    /**
     * calculate the derivative of the GEV w.r.t. mu
     * 
     * the runtime cost is 5 transcendental functions.
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    public static double derivWRTMu(float yConst, float mu, float k, float sigma, float x) {
        
        double z = 1. + k *( (x-mu)/sigma );
        
        double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        double f2 = Math.pow(z, (-1. - (1./k)) );
      
        double dzdmu = -1. * k/sigma;
        
        // df2dmu = (-1-(1/k)) * z^(-2-(1/k)) * dzdmu
        double df2dmu = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdmu;
        
        double df1dmu = f1 * (1/k) * dzdmu;
        if (z < 0) {
            df1dmu *= -1. * Math.pow(-1.*z, (-1 - (1/k)));
        } else {
            df1dmu *= Math.pow(z, (-1 - (1/k)));
        }
        
        double dydmu = (yConst/sigma) * ( f1 * df2dmu + f2 * df1dmu );
        
        if (Double.isNaN(dydmu)) {
            return estimateDerivUsingDeltaMu(mu, k, sigma, x);
        }
        
        return dydmu;
    }
    
    /**
     * calculate d/dk of GEV using the difference between GEVs given minor changes in k
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    static Double estimateDerivUsingDeltaMu(float mu, float k, float sigma, float x) {
        
        float deltaMu = 0.0001f*mu;
        
        Double d0 = GeneralizedExtremeValue.generateYGEV(x, k, sigma, (mu - deltaMu));
        
        Double d1 = GeneralizedExtremeValue.generateYGEV(x, k, sigma, mu);
        
        Double d2 = GeneralizedExtremeValue.generateYGEV(x, k, sigma, (mu + deltaMu));
        
        return estimateDerivUsingDelta(d0, d1, d2, deltaMu);
    }
    
    /**
     * for given mu, k, sigma and the data x, normalizedY and normalizedYErr,
     * calculate the partial derivative of the GEV(k, sigma, mu, x[i]) with respect to 
     * k, then sigma, then mu
     * over each i in x and return the derivatives of each which minimize the chi square sum.
     * 
     * runtime cost is (x.length * 3) iterations of a function that has 5 transcendental operations
     *    + x.length iterations of a function that has 2 transcendental operations.
     *    
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param normalizedY
     * @param normalizedYErr
     * @param idx0 start of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @param idx1 stop of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @param derivs array to populate with answers.  it's given as an argument to reuse the array
     * @return
     */
    //TODO:  don't need to pass back the yconst used do I?
    public static void derivsThatMinimizeChiSqSum(float mu, float k, float sigma, float[] x, 
        float[] normalizedY, float[] normalizedYErr, float[] derivs, int idx0, int idx1) {
                
        Arrays.fill(derivs, 0);
        
        float[] chiSqMin = new float[derivs.length];
        Arrays.fill(chiSqMin, Float.MAX_VALUE);

        float[] yGEV = GeneralizedExtremeValue.genCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return;
        }
        float yMax = MiscMath.findMax(yGEV);
        for (int ii = 0; ii < yGEV.length; ii++){
            yGEV[ii] /= yMax;
        }
        float yConst = 1.f/yMax;
        
        // fGEV(x, kVar, sigmaVar, mu) = (yNorm/sigma)*exp(-1*(1+k((x-mu)/sigma))^(-1/k)) * (1+k((x-mu)/sigma))^(-1-(1/k))
        // We want to minimize sum over all data points in histogram: sum over i of (fGEV[i] - y[i]) * yerr[i]
        //    -- y is normalized to a max of 1.  make sure that fGEV is normalized too (divide it by its peak).
        //    -- deriv of sum w.r.t. sigmaVar is yerr*DerivGEV.wrtSigma... sum over i
        //    -- deriv of sum w.r.t. kVar     is yerr*DerivGEV.wrtK... sum over i
        //    -- deriv of sum w.r.t. muVar    is yerr*DerivGEV.wrtMu... sum over i
        //    ----> because yerr[i] is the same for all derivatives and we're comparing, we can drop it from the calcs
        
        float[] tmpChiSqMin = new float[derivs.length];
        float[] tmpDerivs = new float[tmpChiSqMin.length];
        
        float min = Float.MAX_VALUE;
        
        //TODO: can change this to use avg of values in spectrum, but any single value for a single x appears to be same within 
        
        for (int i = 0; i < x.length; i++) {
                        
            Arrays.fill(tmpChiSqMin, Float.MAX_VALUE);
            Arrays.fill(tmpDerivs, 0);
            
            for (int j = idx0; j <= idx1; j++) {
                Double deriv = null;
                switch (j) {
                    case 0:
                        //deriv = DerivGEV.derivWRTK(yConst, mu, k, sigma, x[i]);
                        deriv = DerivGEV.estimateDerivUsingDeltaK(mu, k, sigma, x[i]);
                        break;
                    case 1:
                        deriv = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, x[i]);
                        //deriv = DerivGEV.estimateDerivUsingDeltaSigma(mu, k, sigma, x[i]);
                        break;
                    case 2:
                        deriv = DerivGEV.derivWRTMu(yConst, mu, k, sigma, x[i]);
                        break;
                }
                if (deriv != null && !deriv.isNaN()) {
                    //TODO:  if wanted to reduce runtime, could remove the method frame load/unload here by importing the method statements
                    float chiSqSum = chiSqSum(yGEV, normalizedY, normalizedYErr);
                    if (chiSqSum < tmpChiSqMin[j]) {
                        tmpChiSqMin[j] = chiSqSum;
                        tmpDerivs[j] = deriv.floatValue();
                    }
                } else {
                    break;
                }
            }
            
            // compare chSqMins:   chose the one w/ a smallest chiSqMin OR smallest sum of all j chiSqMin?
            boolean allAreSet = true;
            boolean minIsCurrent = true;
            for (int j = idx0; j <= idx1; j++) {
                if (tmpChiSqMin[j] == Float.MAX_VALUE || Float.isNaN(tmpChiSqMin[j])) {
                    allAreSet = false;
                    break;
                }
                if (tmpChiSqMin[j] < min) {
                    min = tmpChiSqMin[j];
                    minIsCurrent = false;
                }
            }
            if (allAreSet && !minIsCurrent) {                
                 System.arraycopy(tmpChiSqMin, idx0, chiSqMin, idx0, (idx1 - idx0 + 1));
                 System.arraycopy(tmpDerivs,   idx0, derivs,   idx0, (idx1 - idx0 + 1));
            }
        }        
    }
    
    public static Float chiSqSum(float mu, float k, float sigma, float[] x, float[] normalizedY, float[] normalizedYErr) {
        float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return null;
        }
        return chiSqSum(yGEV, normalizedY, normalizedYErr);
    }
    protected static float chiSqSum(float[] normalizedYGEV, float[] normalizedY, float[] normalizedYErr) {
        float chiSqSum = 0.0f;
        for (int ii = 0; ii < normalizedYGEV.length; ii++) {
            float z = (normalizedYGEV[ii] - normalizedY[ii]);
            z *= normalizedYErr[ii];
            chiSqSum += z*z;
        }
        return chiSqSum;
    }
}
