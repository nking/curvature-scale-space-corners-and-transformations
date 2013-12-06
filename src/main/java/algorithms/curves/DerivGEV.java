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
        
        float deltaK = 0.0001f*k;
        
        Double d0 = GeneralizedExtremeValue.generateYGEV(x, (k - deltaK), sigma, mu);
        
        Double d1 = GeneralizedExtremeValue.generateYGEV(x, (k), sigma, mu);
        
        Double d2 = GeneralizedExtremeValue.generateYGEV(x, (k + deltaK), sigma, mu);
        
        Double d = estimateDerivUsingDelta(d0, d1, d2, deltaK);
        
        return d;
    }
    
    /**
     * calculate the derivative giving the 3 values which are separated by delta.
     * The first, d0, was computed with param - delta.
     * The second, d1, was computed with param.
     * The third, d2, was computed with param + delta.
     * @param d0
     * @param d1
     * @param d2
     * @param delta
     * @return
     */
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
        
        double df0dsigma = -1.f*yConst/(sigma*sigma);
                
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
        
        float delta = 0.0001f*sigma;
        
        Double d0 = GeneralizedExtremeValue.generateYGEV(x, k, (sigma - delta), mu);
        
        Double d1 = GeneralizedExtremeValue.generateYGEV(x, k, sigma, mu);
        
        Double d2 = GeneralizedExtremeValue.generateYGEV(x, k, (sigma + delta), mu);
        
        return estimateDerivUsingDelta(d0, d1, d2, delta);
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
     * calculate the changes in k, sigma, and mu which would reduce the chi sq sum.
     * 
     * Note that internally, the first and second derivatives of the curve GEV(k, sigma, mu)
     * are used to calculate what the smallest step in k, sigma, or mu would be in
     * order to create a significant change in the GEV curve.
     * 
     * The suggested change for k is calculated from d/dk modified by the preconditioner
     * at the point right after the model peak.
     * 
     * Note that the suggested changes might be applied by the NonQuadraticConguteSolver as
     * a fraction of the suggested change.
     * 
     * runtime cost is 
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
    public static void derivsThatMinimizeChiSqSum(float mu, float k, float sigma, float[] x, 
        float[] normalizedY, float[] normalizedYErr, float[] derivs, int idx0, int idx1) {
                
        Arrays.fill(derivs, 0);
        
        // to determine the suggested step for k, sigma, or mu,
        //   we look at the first derivatives of the point in the model GEV right after the maximum
        //   and then modify those derivatives using a preconditioner from ICU0 matrix.
        // that suggested step is tried as plus and minus of current variable and the one
        //   with the smallest resulting chisqsum from the curve produced by the change is the
        //   result returned in derivs for that variable.
        
        float[] yGEV = GeneralizedExtremeValue.genCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return;
        }
        int yMaxIdx = MiscMath.findYMaxIndex(yGEV);
        float yMax = yGEV[yMaxIdx];
        for (int ii = 0; ii < yGEV.length; ii++){
            yGEV[ii] /= yMax;
        }
        float yConst = 1.f/yMax;
        
        float xPoint = x[yMaxIdx + 0];
                
        for (int i = idx0; i < (idx1 + 1); i++) {

            float[] yGEVPlus = null;
            float[] yGEVMinus = null;
            double rModified = 0;
            switch(i) {
                case 0: {
                    // k
                    // calculate a step size that would affect a change in GEV by using the 1st and 2nd partial derivatives
                    rModified = calculatePreconditionerModifiedResidualK(yConst, mu, k, sigma, xPoint);
                    
                    // test whether adding or subtracting the residual results in a reduced chisqsum
                    if (rModified != 0) {
                        yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, (float)(k + rModified), sigma, mu);
                        yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, (float)(k - rModified), sigma, mu);
                    }
                    
                    break;
                }
                case 1: {
                    // sigma
                    // calculate a step size that would affect a change in GEV by using the 1st and 2nd partial derivatives
                    rModified = calculatePreconditionerModifiedResidualSigma(yConst, mu, k, sigma, xPoint);
                    // test whether adding or subtracting the residual results in a reduced chisqsum
                    if (rModified != 0) {
                        yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, (float)(sigma + rModified), mu);
                        yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, (float)(sigma - rModified), mu);
                    }
                    
                    break;
                }
                case 2: {
                    // mu
                    // calculate a step size that would affect a change in GEV by using the 1st and 2nd partial derivatives
                    rModified = calculatePreconditionerModifiedResidualMu(yConst, mu, k, sigma, xPoint);
                    
                    // test whether adding or subtracting the residual results in a reduced chisqsum
                    if (rModified != 0) {
                        yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, (float)(mu + rModified));
                        yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, (float)(mu - rModified));
                    }
                    
                    break;
                }
            }
            if (rModified != 0 && yGEVPlus != null && yGEVMinus != null) {
                float chiSqSumPlus = chiSqSum(yGEVPlus, normalizedY, normalizedYErr);
                float chiSqSumMinus = chiSqSum(yGEVMinus, normalizedY, normalizedYErr);
                float best = (chiSqSumPlus <= chiSqSumMinus) ? (float)rModified : (float)(-1.f * rModified);
                derivs[i] = best;
            }
            System.out.println(" derivs[" + i + "]=" + derivs[i]);
        }
    }

    /**
     * when derivsThatMinimizeChiSqSum is not successfully finding acceptable changes for the
     * curve to improve the minimization of the chi - square sum, use this method
     * to test a range of changes and return them in the array r if successfully found 
     * a change.
     * 
     * @param vars
     * @param x
     * @param y
     * @param ye
     * @param r
     * @param chiSqSumReturn item [0] is the best current chiSqSum for given vars.
     *    item [1] is to return the value of the best chiSqSum here to the invoker
     *    if r was set to non-zero values.
     * @param idx the index in vars to start with to help rotate which changes are
     *   set in r first and accepted.
     */
    public static void exploreChangeInVars(float[] vars, float[] x, float[] normalizedY,
        float[] normalizedYErr, float[] r, float[] chiSqSumReturn, int idx) {
        
        float defaultChange = 0.6f;
        
        Arrays.fill(r, 0);
              
        float bestChiSqSum = chiSqSumReturn[0];
                       
        // apply mu changes first by starting with i=2
        for (int c = 0; c < 3; c++) {
            
            int i = idx + c;
            if (i > 2) {
                i = i - 3;
            }
            
            float k = vars[0] + r[0];
            float sigma = vars[1] + r[1];
            float mu = vars[2] + r[2];
            
            float rModified = defaultChange;
            
            if (i == 2) {
                rModified = 0.5f - mu;
            }
            
            float[] yGEVPlus = null;
            float[] yGEVMinus = null;
            
            int nMaxIter = 10;
            int nIter = 0;
            
            boolean hasRetried = false;
            boolean useDoubleChange = false;
            
            while (nIter < nMaxIter) {
                switch(i) {
                    case 0: {
                        // k                    
                        yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, (float)(k + rModified), sigma, mu);
                        yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, (float)(k - rModified), sigma, mu);
                    
                        break;
                    }
                    case 1: {
                        // sigma
                        yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, (float)(sigma + rModified), mu);
                        yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, (float)(sigma - rModified), mu);
                        
                        break;
                    }
                    case 2: {
                        // mu
                        
                        yGEVPlus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, (float)(mu + rModified));
                        yGEVMinus = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, (float)(mu - rModified));
                                                
                        break;
                    }
                }
                        
                float chiSqSumPlus = chiSqSum(yGEVPlus, normalizedY, normalizedYErr);
                float chiSqSumMinus = chiSqSum(yGEVMinus, normalizedY, normalizedYErr);
                boolean bestIsPlus = (chiSqSumPlus < chiSqSumMinus);
                    
                if (bestIsPlus && (chiSqSumPlus < bestChiSqSum)) {
                    r[i] = rModified;
                    bestChiSqSum = chiSqSumPlus;
                    // continue to next variable
                    break;
                } else if (chiSqSumMinus < bestChiSqSum) {
                    r[i] = -1.f * rModified;
                    bestChiSqSum = chiSqSumMinus;
                    // continue to next variable
                    break;
                } else {
                    
                    if (hasRetried) {
                        if (useDoubleChange) {
                            rModified *= 2.f;
                        } else {
                            rModified /= 2.f;
                        }
                    } else {
                        rModified /= 2.f;
                        hasRetried = true;
                    }
                    
                    nIter++;
                }                
            }
        }
        chiSqSumReturn[1] = bestChiSqSum;
    }
    
    public static double calculatePreconditionerModifiedResidualK(
        float yConst, float mu, float k, float sigma, float x) {

        // using Incomplete Cholesky factorization with fill 0 (ICU0) to apply preconditioning
        // to the first derivative
        // 
        // k component to residuals = d(1,1) * (∂f/∂k)
        //       where d(1,1) is 1./(∂^2f/∂k∂k)        
        
        Double dydk = DerivGEV.derivWRTK(yConst, mu, k, sigma, x);
        
        if (dydk == null) {
            return 0;
        }
        
        Double d2ydkdk = estimateDY2DKDK(yConst, mu, k, sigma, x, dydk);
        
        double resid = dydk/d2ydkdk;
        
        return resid;
    }
    
    /**
     * estimate ∂^2f/∂k∂k empirically.
     * 
     * The method uses the tested DerivGEV.derivWRTK() for ∂f/∂k and plugs in different k's 
     * to estimate ∂^2f/∂k∂k.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    static double estimateDY2DKDK(float yConst, float mu, float k, float sigma, float x) {
     
        // ∂^2f/∂k∂k 
        
        Double dydk = DerivGEV.derivWRTK(yConst, mu, k, sigma, x);
        
        if (dydk == null) {
            return 0;
        }
        
        return estimateDY2DKDK(yConst, mu, k, sigma, x, dydk.doubleValue());
    }
    
    /**
     * estimate ∂^2f/∂k∂k empirically.  the method accepts dydk as a given to allow easier
     * resuse in other equations, but has to trust that dydk was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydk
     * @return
     */
    static double estimateDY2DKDK(float yConst, float mu, float k, float sigma, float x, double dydk) {

        // ∂^2f/∂k∂k = estimate as (dydk_2 - dydk)/dk
        
        float factor = 0.0001f;
        
        double delta = (k*factor);
        
        Double dydk_0 = DerivGEV.derivWRTK(yConst, mu, (float)(k - delta), sigma, x);
        
        Double dydk_2 = DerivGEV.derivWRTK(yConst, mu, (float)(k + delta), sigma, x);
        
        return estimateDerivUsingDelta(dydk_0, dydk, dydk_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂k∂sigma empirically.  the method accepts dydk as a given to allow easier
     * resuse in other equations, but has to trust that dydk was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydk
     * @return
     */
    static double estimateDY2DKDSigma(float yConst, float mu, float k, float sigma, float x, double dydk) {

        // ∂^2f/∂k∂sigma 
        
        float factor = 0.0001f;
        
        double delta = (sigma*factor);
        
        Double dydk_0 = DerivGEV.derivWRTK(yConst, mu, k, (float)(sigma - delta), x);
        
        Double dydk_2 = DerivGEV.derivWRTK(yConst, mu, k, (float)(sigma + delta), x);
        
        return estimateDerivUsingDelta(dydk_0, dydk, dydk_2, delta);
    }

    /**
     * estimate ∂^2f/∂k∂mu empirically.  the method accepts dydk as a given to allow easier
     * resuse in other equations, but has to trust that dydk was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydk
     * @return
     */
    static double estimateDY2DKDMu(float yConst, float mu, float k, float sigma, float x, double dydk) {

        // ∂^2f/∂k∂mu 
        
        float factor = 0.0001f;
        
        double delta = (mu*factor);
        
        Double dydk_0 = DerivGEV.derivWRTK(yConst, (float)(mu - delta), k, sigma, x);
        
        Double dydk_2 = DerivGEV.derivWRTK(yConst, (float)(mu + delta), k, sigma, x);
       
        return estimateDerivUsingDelta(dydk_0, dydk, dydk_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂sigma∂sigma empirically.  the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    static double estimateDY2DSigmaDSigma(float yConst, float mu, float k, float sigma, float x, double dydsigma) {

        // ∂^2f/∂sigma∂sigma = estimate as (dyds_2 - dyds)/ds
        
        float factor = 0.0001f;
        
        double delta = (sigma*factor);
        
        Double dyds_0 = DerivGEV.derivWRTSigma(yConst, mu, k, (float)(sigma - delta), x);
        
        Double dyds_2 = DerivGEV.derivWRTSigma(yConst, mu, k, (float)(sigma + delta), x);
                
        return estimateDerivUsingDelta(dyds_0, dydsigma, dyds_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂mu∂mu empirically.  the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    static double estimateDY2DMuDMu(float yConst, float mu, float k, float sigma, float x, double dydmu) {

        // ∂^2f/∂mu∂mu
        
        float factor = 0.0001f;
        
        double delta = (mu*factor);
        
        Double dydm_0 = DerivGEV.derivWRTMu(yConst, (float)(mu - delta), k, sigma, x);
        
        Double dydm_2 = DerivGEV.derivWRTMu(yConst, (float)(mu + delta), k, sigma, x);
                
        return estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂mu∂k empirically.  the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    static double estimateDY2DMuDK(float yConst, float mu, float k, float sigma, float x, double dydmu) {

        // ∂^2f/∂mu∂k
        
        float factor = 0.0001f;
        
        double delta = (k*factor);
        
        Double dydm_0 = DerivGEV.derivWRTMu(yConst, mu, (float)(k - delta), sigma, x);
        
        Double dydm_2 = DerivGEV.derivWRTMu(yConst, mu, (float)(k + delta), sigma, x);
                
        return estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂mu∂sigma empirically.  the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydmu was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    static double estimateDY2DMuDSigma(float yConst, float mu, float k, float sigma, float x, double dydmu) {

        // ∂^2f/∂mu∂sigma
        
        float factor = 0.0001f;
        
        double delta = (sigma*factor);
        
        Double dydm_0 = DerivGEV.derivWRTMu(yConst, mu, k, (float)(sigma - delta), x);
        
        Double dydm_2 = DerivGEV.derivWRTMu(yConst, mu, k, (float)(sigma + delta), x);
        
        return estimateDerivUsingDelta(dydm_0, dydmu, dydm_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂sigma∂k empirically.  the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    static double estimateDY2DSigmaDK(float yConst, float mu, float k, float sigma, float x, double dydsigma) {

        // ∂^2f/∂sigma∂k = estimate as (dy2_ds_dk - dyds)/dk
        
        float factor = 0.0001f;
        
        double delta = (k*factor);
        
        Double dydk_0 = DerivGEV.derivWRTSigma(yConst, mu, (float)(k - delta), sigma, x);
        
        Double dydk_2 = DerivGEV.derivWRTSigma(yConst, mu, (float)(k + delta), sigma, x);
        
        return estimateDerivUsingDelta(dydk_0, dydsigma, dydk_2, delta);
    }
    
    /**
     * estimate ∂^2f/∂sigma∂mu empirically.  the method accepts dydsigma as a given to allow easier
     * resuse in other equations, but has to trust that dydsigma was derived with the
     * same k, sigma, mu, and x.
     * 
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @param dydsigma
     * @return
     */
    static double estimateDY2DSigmaDMu(float yConst, float mu, float k, float sigma, float x, double dydsigma) {

        // ∂^2f/∂sigma∂mu = estimate as (dy2_ds_dk - dyds)/dk
        
        float factor = 0.0001f;
        
        double delta = (sigma*factor);
        
        Double dyds_0 = DerivGEV.derivWRTSigma(yConst, (float)(mu - delta), k, sigma, x);
        
        Double dyds_2 = DerivGEV.derivWRTSigma(yConst, (float)(mu + delta), k, sigma, x);
        
        return estimateDerivUsingDelta(dyds_0, dydsigma, dyds_2, delta);
    }
    
    public static double calculatePreconditionerModifiedResidualSigma(
        float yConst, float mu, float k, float sigma, float x) {

        // using Incomplete Cholesky factorization with fill 0 (ICU0) to apply preconditioning
        // to the first derivative
        // 
        // sigma component to residuals = d(2,2) * (∂f/∂sigma)
        //       where d(2,2) is ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
        
        // ∂f/∂sigma
        Double dyds = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, x);
        
        if (dyds == null) {
            return 0;
        }
        
        
        // ∂^2f/∂sigma∂sigma
        double d2ydsds = estimateDY2DSigmaDSigma(yConst, mu, k, sigma, x, dyds);
            
        
        // (∂^2f/∂sigma∂k) = estimate as (dydsigma_1 - dydsigma)/dsigma        
        double d2ydsdk = estimateDY2DSigmaDK(yConst, mu, k, sigma, x, dyds);
        
        // ∂f/∂k        
        Double dydk = DerivGEV.derivWRTK(yConst, mu, k, sigma, x);
        
        if (dydk != null) {
            
            // ∂^2f/∂k∂k
            double d2ydkdk = estimateDY2DKDK(yConst, mu, k, sigma, x, dydk.doubleValue());

            // ∂^2f/∂k∂sigma
            double d2ydkds = estimateDY2DKDSigma(yConst, mu, k, sigma, x, dydk.doubleValue());

            double modification = d2ydsds - (d2ydsdk / d2ydkdk) * d2ydkds;

            double resid = dyds / modification;

            return resid;
            
        } else {
            
            double modification = d2ydsds;
            
            double resid = dyds / modification;

            return resid;
        }
    }

    public static double calculatePreconditionerModifiedResidualMu(
        float yConst, float mu, float k, float sigma, float x) {
        
        // using Incomplete Cholesky factorization with fill 0 (ICU0) to apply preconditioning
        // to the first derivative
        // 
        // sigma component to residuals = d(3,3) * (∂f/∂mu)
        /*
           where  d(3,3) is 1./(
                     ( (∂^2f/∂mu∂mu) - (∂^2f/∂mu∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂mu) )
                     -
                     (
                        ∂^2f/∂mu∂sigma
                        *
                        ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
                        *
                        ∂^2f/∂sigma∂mu
                     )
                 )
                 
                Let pt1 = (∂^2f/∂mu∂mu) - (∂^2f/∂mu∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂mu)
            
                Let pt2_1 = ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
            
                Let pt2 = ( ∂^2f/∂mu∂sigma * pt2_1 * ∂^2f/∂sigma∂mu )
            
            d(3,3) = 1./( pt1 - pt2) 
        */
        
        // ∂f/∂mu
        Double dydmu = DerivGEV.derivWRTMu(yConst, mu, k, sigma, x);
        
        if (dydmu == null) {
            return 0;
        }
        
        // ∂f/∂sigma
        Double dyds = DerivGEV.derivWRTSigma(yConst, mu, k, sigma, x);
        if (dyds == null) {
            return 0;
        }
               
        // ∂^2f/∂sigma∂sigma
        double d2ydsds = estimateDY2DSigmaDSigma(yConst, mu, k, sigma, x, dyds.doubleValue());
        
        if (d2ydsds == 0) {
            return 0;
        }
        
                
        // ∂^2f/∂mu∂mu
        double d2ydmdm = estimateDY2DMuDMu(yConst, mu, k, sigma, x, dydmu.doubleValue());
        
        // ∂f/∂k        
        Double dydk = DerivGEV.derivWRTK(yConst, mu, k, sigma, x);
        
        double pt1;
        
        if (dydk != null) {
            
            // ∂^2f/∂mu∂k
            double d2ydmdk = estimateDY2DMuDK(yConst, mu, k, sigma, x, dydmu.doubleValue());

            // ∂^2f/∂k∂k 
            double d2ydkdk = estimateDY2DKDK(yConst, mu, k, sigma, x, dydk.doubleValue());
            
            // ∂^2f/∂k∂mu
            double d2ydkdm = estimateDY2DKDMu(yConst, mu, k, sigma, x, dydk.doubleValue());
            
            pt1 = d2ydmdm - (d2ydmdk*d2ydkdm / d2ydkdk);
            
        } else {
            
            pt1 = d2ydmdm;
        }        
            
        // (∂^2f/∂sigma∂k) = estimate as (dydsigma_1 - dydsigma)/dsigma        
        double d2ydsdk = estimateDY2DSigmaDK(yConst, mu, k, sigma, x, dyds.doubleValue());
        
        double d2ydmds = estimateDY2DMuDSigma(yConst, mu, k, sigma, x, dydmu.doubleValue());
        
        // ∂^2f/∂sigma∂mu
        double d2ydsdm = estimateDY2DSigmaDMu(yConst, mu, k, sigma, x, dyds.doubleValue());
        
        // d(3,3) = 1./( (pt1) - ( ∂^2f/∂mu∂sigma * pt2_1 * ∂^2f/∂sigma∂mu )) 
        // 
        // pt2_1 = ( 1./ ( (∂^2f/∂sigma∂sigma) - (∂^2f/∂sigma∂k)*( 1/(∂^2f/∂k∂k) ) * (∂^2f/∂k∂sigma) ) )
        // Let pt2 = ( ∂^2f/∂mu∂sigma * pt2_1 * ∂^2f/∂sigma∂mu )
        
        double pt2, pt2_1;
        
        if (d2ydmds == 0) {
            pt2 = 0;
            pt2_1 = 0;
            
        } else if (d2ydsdk != 0 && dydk != null) {
            
            // ∂^2f/∂sigma∂k
            double d2ydkdk = estimateDY2DKDK(yConst, mu, k, sigma, x, dydk.longValue());
                        
            // ∂^2f/∂k∂sigma
            double d2ydkds = estimateDY2DKDSigma(yConst, mu, k, sigma, x, dydk.doubleValue());
            
            pt2_1 = d2ydsds - (d2ydsdk * d2ydkds / d2ydkdk);
            
            pt2_1 = 1./pt2_1;
                
            pt2 = d2ydmds * pt2_1 * d2ydsdm;
            
        } else {
            
            pt2_1 = 1./d2ydsds;
            
            pt2 = d2ydmds * pt2_1 * d2ydsdm;
        }
        
        double modification = pt1 - pt2;
                                
        double resid = dydmu/modification;
        
        return resid;
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
            //TODO:  setup tests with realistic errors
            chiSqSum += z*z*normalizedYErr[ii];
        }
        return chiSqSum;
    }
    
    /**
     * retrieve the second deriv of GEV ∂^2f/∂k∂mu.
     * 
     * Note:  the method hasn't been tested yet and will be compared to the results from method estimateDY2DKDMu.
     * 
     * @param yConst
     * @param mu
     * @param k
     * @param sigma
     * @param x
     * @return
     */
    public static double secondDerivKDerivMu(float yConst, float mu, float k, float sigma, float x) {
        
        double z = 1. + k *( (x-mu)/sigma );
        
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
                
        float deltaK = 0.0001f;
        float deltaMu = 0.0001f;
        double k_1 = k + deltaK;
        double mu_1 = mu + deltaMu;
        double z_1 = 1. + k_1 *( (x-mu)/sigma );
        double z_1_1 = 1. + k_1 *( (x-mu_1)/sigma );
        
        if (zIsNegative) {
            df2dmu *= -1.f;
            df1dmu *= -1.f;
            z_1 *= -1.f;
            z_1_1 *= -1.f;
        }
                
        // approximating  df1dk and df2dk  to avoid the correction needed for log(negative number)
        float a_1 = -1.f*(float) Math.pow(z_1, (-1.f/k_1)); 
        float a_1_1 = -1.f*(float) Math.pow(z_1_1, (-1.f/k_1)); 
        double f2_1 = Math.pow(z_1, (-1. - (1./k_1)) );
        double f2_1_1 = Math.pow(z_1_1, (-1. - (1./k_1)) );
        if (zIsNegative) {
            a_1 *= -1.f;
            a_1_1 *= -1.f;
            f2_1 *= -1.f;
            f2_1_1 *= -1.f;
        }
        double f1_1 = Math.exp( a_1 );
        double df1dk = (Math.exp( a_1 ) - f1)/deltaK;
        double df1dk_1 = (Math.exp( a_1_1 ) - f1_1)/deltaK;
        double df2dk = (f2_1 - f2)/deltaK;
        double df2dk_1 = (f2_1_1 - f2_1)/deltaK;
        
        
        //pt1_0 = ∂/∂mu( df2dk )  the derivative involves potential log(negative number) so work around: 
        double pt1_0 = (df2dk_1 - df2dk)/deltaMu;
        
        // pt2_0 = ∂/∂mu( df1dk ) the derivative involves potential log(negative number) so work around:
        double pt2_0 = (df1dk_1 - df1dk)/deltaMu;
        
        
        double pt2 = f2 * pt2_0 + df1dk * df2dmu;       
        
        double pt1 = (f1 * pt1_0) + (df2dk * df1dmu);
        
        double d2fdkdmu = (yConst/sigma) * (pt1 + pt2);
        
        return d2fdkdmu;
    }

    /* 
     Can see that calculating the formula for the partial derivative ∂^2f/∂k∂mu
     is time consuming and error prone so will only implement the "estimate" methods hereafter.
     
     
     df1dk     = f1 * (-1*z^(-1/k)) * ( -1*(-1/k) * (dzdk/z)  +  (1/k^2) * ln( -z ) )
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
                   = ∂/∂mu( f1 * (-1*z^(-1/k)) * ( -1*(-1/k) * dzdk * (1/z)  +  (1/k^2) * ln( -z ) ) )
                   = ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z) )  +  ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k^2) * ln( -z ) )
                   = pt2_0_0   +   pt2_0_1
                 
                 pt2_0_0 = ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z) )
                         = df1dmu * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z)   
                             + ∂/∂mu(-1*z^(-1/k)) * f1 * (1/k) * dzdk * (1/z)
                             + 0
                             + ∂/∂mu(dzdk) * f1 * (-1*z^(-1/k)) * (1/k) * (1/z)
                             + ∂/∂mu(1/z) * f1 * (-1*z^(-1/k)) * (1/k) * dzdk
                         = df1dmu * (-1*z^(-1/k)) * (1/k) * dzdk * (1/z) 
                             + pt2_0_0_0 * f1 * (1/k) * dzdk * (1/z)
                             + pt2_0_0_1 * f1 * (-1*z^(-1/k)) * (1/k) * (1/z)
                             + pt2_0_0_2 * f1 * (-1*z^(-1/k)) * (1/k) * dzdk
                             
                     pt2_0_0_0 = ∂/∂mu(-1*z^(-1/k))
                               = (1/k) * z^(-1 - (1/k)) * dzdmu  <== from class level comments
                     
                     pt2_0_0_1 = ∂/∂mu(dzdk)
                               = ∂/∂mu( (x-mu)/sigma )
                               = -1/sigm
                             
                     pt2_0_0_2 = ∂/∂mu(1/z)
                               = pt1_0_0_1
                               
                 pt2_0_1 = ∂/∂mu( f1 * (-1*z^(-1/k)) * (1/k^2) * ln( -z ) )
                         = df1dmu * (-1*z^(-1/k)) * (1/k^2) * ln( -z )
                            + ∂/∂mu(-1*z^(-1/k)) * f1 * (1/k^2) * ln( -z )
                            + 0
                            + ∂/∂mu( ln( -z ) ) * f1 * (-1*z^(-1/k)) * (1/k^2)
                         = df1dmu * (-1*z^(-1/k)) * (1/k^2) * ln( -z )
                            + pt2_0_0_0 * f1 * (1/k^2) * ln( -z )
                            + 0
                            + pt2_0_1_0 * f1 * (-1*z^(-1/k)) * (1/k^2)
                        
                     pt2_0_1_0 = ∂/∂mu( ln ( -1*(sigma + k*(x-mu))/sigma ) 
                               = (k*sigma)/(sigma + k*(x-mu))

*/
}
