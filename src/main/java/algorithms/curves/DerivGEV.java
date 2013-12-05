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
                
                 u(k) = -z
                 v(k) = (-1/k)
                 
     *           df(k)dk = -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
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
     *   df1dk     = f1 * -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
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
                 df2dk = f2 * ( v(k) * -------   +   dv(k)dk * ln(u(k)) )
                                         u(k)
                       = f2 * ( (-1-(1/k)) * dzdk/z  +  (1/k^2) * ln(z) )
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
     *   dydsigma:  
     *        needs to use chain rule once more
     *        
     *        f0 = (yconst/sigma)
     *        df0dsigma = -(yconst/sigma^2)
     *        
     *        f = f0 * f1 * f2
     *        
     *        dydsigma = (  df0dsigma * f1 * f2 ) + ( df1dsigma * f0 * f2 ) + (df2dsigma * f0 * f1 )
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
        
        /*if (z < 0) {
            return estimateDerivUsingDeltaK(mu, k, sigma, x);
        }*/
        
        boolean zIsNegative = (z < 0);
        
        if (zIsNegative) {
            z *= -1;
        }
        
        float a = -1.f*(float) Math.pow(z, (-1.f/k));
        
        if (zIsNegative) {
            a *= -1.f;
        }
        
        double f1 = Math.exp( a );
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        if (zIsNegative) {
            f2 *= -1.f;
        }
      
        double dzdx = k/sigma;

        double df2dx = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdx;
        
        double df1dx = f1 * (1./k) * Math.pow(z, (-1. - (1./k))) * dzdx;
        
        if (zIsNegative) {
            df2dx *= -1.f;
            df1dx *= -1.f;
        }
        
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
    public static double derivWRTK(float yConst, float mu, float k, float sigma, float x) {
        
        double z = 1. + k *( (x-mu)/sigma );
      
        /*if (z < 0) {
            return estimateDerivUsingDeltaK(mu, k, sigma, x);
        }*/
        
        boolean zIsNegative = (z < 0);
        
        if (zIsNegative) {
            z *= -1;
        }
        
        float a = -1.f*(float) Math.pow(z, (-1.f/k));
        
        if (zIsNegative) {
            a *= -1.f;
        }
        
        double f1 = Math.exp( a );
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        if (zIsNegative) {
            f2 *= -1.f;
        }
        
        // df1dk     = f1 * -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
        
        // df2dk = f2 * ( (-1-(1/k)) * dzdk/z + (1/k^2) * ln(z) )
        
        // any value of z will have trouble with ln(-z) or ln(z) because the built in logarithm doesn't
        //    use a taylor series approximation for negative values plus handling for the number of cycles
        // so df1dk and df2dk are approximated with very small deltas
        
        float deltaK = 0.0001f;
        
        double k_1 = k + deltaK;
        double z_1 = 1. + k_1 *( (x-mu)/sigma );
        if (zIsNegative) {
            z_1 *= -1;
        }
        
        float a_1 = -1.f*(float) Math.pow(z_1, (-1.f/k_1));
        
        double f2_1 = Math.pow(z_1, (-1. - (1./k_1)) );
        
        if (zIsNegative) {
            a_1 *= -1.f;
            f2_1 *= -1.f;
        }
        
        double df1dk = (Math.exp( a_1 ) - f1)/deltaK;
        
        double df2dk = (f2_1 - f2)/deltaK;
        
        // to compare to the derivative:
        /*double dzdk = (x-mu)/sigma;
        if (zIsNegative) {
            double compare_df1dk = f1 * a * ( (1.f/k)*(dzdk/z) + (1.f/(k*k))*Math.log( z ) );
            System.out.println( String.format("  df1dk   estimate=%4.6f  deriv=%4.6f ", df1dk, compare_df1dk));
        } else {
            double compare_df2dk = f2 * (  (-1.f - (1.f/k))*(dzdk/z) + (1.f/(k*k))*Math.log(z) );
            System.out.println( String.format("  df2dk   estimate=%4.6f  deriv=%4.6f   (z=%4.6f, k=%4.6f)", df2dk, compare_df2dk, z, k));
        }*/
                        
        double dydk = (yConst/sigma) * ( f1 * df2dk + f2 * df1dk );
        
        return dydk;
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
        
        //double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        //double f2 = Math.pow(z, (-1. - (1./k)) );
        
        boolean zIsNegative = (z < 0);
        
        if (zIsNegative) {
            z *= -1;
        }
        
        float a = -1.f*(float) Math.pow(z, (-1.f/k));
        
        if (zIsNegative) {
            a *= -1.f;
        }
        
        double f0 = (yConst/sigma);
        double f1 = Math.exp( a );
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        if (zIsNegative) {
            f2 *= -1.f;
        }
        
        double dzdsigma = -1. * k * (x-mu) * Math.pow(sigma, -2.);
        
        //(-1-(1/k)) * z^(-2-(1/k)) * dzdsigma
        double df2dsigma = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdsigma;
        
        //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
        double df1dsigma =  f1 * (1./k) * Math.pow(z, -1. - (1./k)) * dzdsigma;
        
        if (zIsNegative) {
            df2dsigma *= -1.f;
            df1dsigma *= -1.f;
        }
        
        double df0dsigma = -1.f/(sigma*sigma);
                
        double dydSigma = (  df0dsigma * f1 * f2 ) + ( df1dsigma * f0 * f2 ) + (df2dsigma * f0 * f1 );
        
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
        
        //double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        //double f2 = Math.pow(z, (-1. - (1./k)) );
        
        boolean zIsNegative = (z < 0);
        
        if (zIsNegative) {
            z *= -1;
        }
        
        float a = -1.f*(float) Math.pow(z, (-1.f/k));
        
        if (zIsNegative) {
            a *= -1.f;
        }
        
        double f1 = Math.exp( a );
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        if (zIsNegative) {
            f2 *= -1.f;
        }
       
        double dzdmu = -1. * k/sigma;
        
        //(-1-(1/k)) * z^(-2-(1/k)) * dzdmu
        double df2dmu = (-1. - (1./k)) * Math.pow(z, (-2. - (1./k)) ) * dzdmu;
        
        //f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
        double df1dmu =  f1 * (1./k) * Math.pow(z, -1. - (1./k)) * dzdmu;
        
        if (zIsNegative) {
            df2dmu *= -1.f;
            df1dmu *= -1.f;
        }
        
        double dydmu = (yConst/sigma) * ( f1 * df2dmu + f2 * df1dmu );
                
        if (Double.isNaN(dydmu)) {
            return estimateDerivUsingDeltaMu(mu, k, sigma, x);
        }
        
        return dydmu;
    }
    
    /**
     * estimate d/dmu of GEV using the difference between GEVs given minor changes in k
     * 
     * Note that this method does not match the results of method derivWRTMu.  Prefer method derivWRTMu
     * when possible.
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
                
        for (int i = 0; i < x.length; i++) {
                        
            Arrays.fill(tmpChiSqMin, Float.MAX_VALUE);
            Arrays.fill(tmpDerivs, 0);
            
            for (int j = idx0; j <= idx1; j++) {
                Double deriv = null;
                switch (j) {
                    case 0:
                        deriv = DerivGEV.derivWRTK(yConst, mu, k, sigma, x[i]);
                        //deriv = DerivGEV.estimateDerivUsingDeltaK(mu, k, sigma, x[i]);
                        break;
                    case 1:
                        deriv = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, x[i]);
                        //deriv = DerivGEV.estimateDerivUsingDeltaSigma(mu, k, sigma, x[i]);
                        break;
                    case 2:
                        deriv = DerivGEV.derivWRTMu(yConst, mu, k, sigma, x[i]);
                        //deriv = DerivGEV.estimateDerivUsingDeltaMu(mu, k, sigma, x[i]);
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
            chiSqSum += z*z*normalizedYErr[ii];
        }
        return chiSqSum;
    }
    
    /**
     * see notes in NonQuadraticConjugateGradientSolver.java
     * 
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
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param normalizedY
     * @param normalizedYErr
     * @param idx0 start of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @param idx1 stop of index within derivs, inclusive, to use in solution.  index 0 = k, index 1 = sigma, index 2 = mu
     * @param r array of first derivatives of k, sigma, and mu that will be transformed using the ICU0 preconditioner
     */
    public static void transformRUsingICU0Preconditioner(float mu, float k, float sigma, float[] x, 
        float[] normalizedY, float[] normalizedYErr, float[] r, int idx0, int idx1) {
        
        // Need 2nd deriv formulas.
        
        // The comparison of manually calculating the first derivative formulas versus 
        // estimating the first deriv at the points plus and minus a small delta
        // shows that an estimate is fine for all partial derivatives except for mu.
        //   either there's a bug that I haven't found yet for the mu formula or the derivatives w.r.t. mu aren't
        //   sensitive to changes (presumably because mu is related to x which is an independent variable in the eqn).
        // For that reason, can use the estimate method for 2nd derivatives which don't use mu
        //   but for mu, need to derive the 2nd derivatives
        //
        // Need to manually calculate derivatives for the following:
        // ∂^2f/∂k∂mu
        // ∂^2f/∂mu∂mu
        // ∂^2f/∂mu∂k
        // ∂^2f/∂mu∂sigma
        
        // And estimates for these:  ∂^2f/∂k∂k, ∂^2f/∂k∂sigma, ∂^2f/∂sigma∂sigma, ∂^2f/∂sigma∂k
        
    }
    
    /* 
     df1dk     = f1 * -z^(-1/k) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
     df2dk     = f2 * ( (-1-(1/k)) * dzdk/z  +  (1/k^2) * ln(z) )
     dzdk = (x-mu)/sigma
     
     ∂^2f/∂k∂mu
         = ∂/∂mu ( (yconst/sigma) * ( f1 * df2dk + f2 * df1dk )  )
         = (yconst/sigma) * ∂/∂mu( f1 * df2dk )  + (yconst/sigma) * ∂/∂mu(f2 * df1dk )
         = (yconst/sigma) *  pt1                 + (yconst/sigma) * pt2   
       
         pt1 = ∂/∂mu( f1 * df2dk )
             = f1 * ∂/∂mu( df2dk ) + df2dk * ∂/∂mu(f1)
             = f1 * pt1_0          + df2dk * df1dmu
             
             pt1_0 = ∂/∂mu( df2dk )
                   = ∂/∂mu( f2 * ( (-1-(1/k)) * dzdk * (1/z)  +  (1/k^2) * ln(z) ) )
                      
                   for a product of more than 2 factors, multiply deriv of one times all 
                   other factors then add next deriv times all other factors...
       
                   = ∂/∂mu( f2 * (-1-(1/k)) * dzdk * (1/z) )  +  ∂/∂mu( f2 * (1/k^2) * ln(z) )
                   
                   = pt1_0_0                                  + pt1_0_1
                   
                 pt1_0_0 = ∂/∂mu( f2 ) * (-1-(1/k)) * dzdk * (1/z)  +  0  
                             + ∂/∂mu(dzdk) * f2 * (-1-(1/k)) * (1/z)  
                             + ∂/∂mu((1/z)) * f2 * (-1-(1/k)) * dzdk
                         = df2dmu * (-1-(1/k)) * dzdk * (1/z)   
                             + pt1_0_0_0 * f2 * (-1-(1/k)) * (1/z)
                             + pt1_0_0_1 * f2 * (-1-(1/k)) * dzdk
                                 
                     pt1_0_0_0 = ∂/∂mu(dzdk)
                               = ∂/∂mu( (x-mu)/sigma )
                               = -mu/sigma
                     
                     pt1_0_0_1 = ∂/∂mu(1/z)  
                     
                                    z = (1 + k*( (x-mu)/sigma ) = (sigma + k*(x-mu))/sigma
                                    1/z = sigma/(sigma + k*(x-mu))
                                    
                               = ∂/∂mu(sigma/(sigma + k*(x-mu)))
                               = ∂/∂mu( sigma * (sigma + k*(x-mu))^-1 )
                               = -1 * sigma * (sigma + k*(x-mu))^-2  * (-k)
                               = k * sigma * (sigma + k*(x-mu))^-2
                               
                 pt1_0_1 = ∂/∂mu( f2 * (1/k^2) * ln(z) ) 
                         = df2dmu * (1/k^2) * ln(z)
                            + 0
                            + ∂/∂mu( ln(z) ) * f2 * (1/k^2)
                         = df2dmu * (1/k^2) * ln(z)
                            + pt1_0_1_0 * f2 * (1/k^2)
                     
                     pt1_0_1_0 = ∂/∂mu( ln(z) )   where logarithm( u(mu) ) = (1/u(mu)) * d*u(mu)/dmu
                               = ∂/∂mu( ln ( (sigma + k*(x-mu))/sigma ) 
                               = ( sigma / (sigma + k*(x-mu)) ) * (-k)
                               = (-k * sigma) / ( k * (sigma + k*(x-mu)) )
         
         pt2 = ∂/∂mu(f2 * df1dk )
             = f2 * ∂/∂mu( df1dk ) + df1dk * df2dmu
             = f2 * pt2_0 + df1dk * df2dmu
         
             pt2_0 = ∂/∂mu( df1dk )
                   = ∂/∂mu( f1 * -z^(-1/k) * ( -1*(-1/k) * dzdk * (1/z)  +  (1/k^2) * ln( -z ) ) )
                   = ∂/∂mu( f1 * (-z^(-1/k)) * (1/k) * dzdk * (1/z) )  +  ∂/∂mu( f1 * (-z^(-1/k)) * (1/k^2) * ln( -z ) )
                   = pt2_0_0   +   pt2_0_1
                 
                 pt2_0_0 = ∂/∂mu( f1 * (-z^(-1/k)) * (1/k) * dzdk * (1/z) )
                         = df1dmu * (-z^(-1/k)) * (1/k) * dzdk * (1/z)   
                             + ∂/∂mu(-z^(-1/k)) * f1 * (1/k) * dzdk * (1/z)
                             + 0
                             + ∂/∂mu(dzdk) * f1 * (-z^(-1/k)) * (1/k) * (1/z)
                             + ∂/∂mu(1/z) * f1 * (-z^(-1/k)) * (1/k) * dzdk
                         = df1dmu * (-z^(-1/k)) * (1/k) * dzdk * (1/z) 
                             + pt2_0_0_0 * f1 * (1/k) * dzdk * (1/z)
                             + pt2_0_0_1 * f1 * (-z^(-1/k)) * (1/k) * (1/z)
                             + pt2_0_0_2 * f1 * (-z^(-1/k)) * (1/k) * dzdk
                             
                     pt2_0_0_0 = ∂/∂mu(-z^(-1/k))
                               = (1/k) * z^(-1 - (1/k)) * dzdmu  <== from class level comments
                     
                     pt2_0_0_1 = ∂/∂mu(dzdk)
                               = ∂/∂mu( (x-mu)/sigma )
                               = -1/sigm
                             
                     pt2_0_0_2 = ∂/∂mu(1/z)
                               = pt1_0_0_1
                               
                 pt2_0_1 = ∂/∂mu( f1 * (-z^(-1/k)) * (1/k^2) * ln( -z ) )
                         = df1dmu * (-z^(-1/k)) * (1/k^2) * ln( -z )
                            + ∂/∂mu(-z^(-1/k)) * f1 * (1/k^2) * ln( -z )
                            + 0
                            + ∂/∂mu( ln( -z ) ) * f1 * (-z^(-1/k)) * (1/k^2)
                         = df1dmu * (-z^(-1/k)) * (1/k^2) * ln( -z )
                            + pt2_0_0_0 * f1 * (1/k^2) * ln( -z )
                            + 0
                            + pt2_0_1_0 * f1 * (-z^(-1/k)) * (1/k^2)
                        
                     pt2_0_1_0 = ∂/∂mu( ln ( -1*(sigma + k*(x-mu))/sigma ) 
                               = (k*sigma)/(sigma + k*(x-mu))
                               

*/
}
