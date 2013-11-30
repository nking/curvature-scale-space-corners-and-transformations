package algorithms.curves;

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
     *                                dzdx
     *   df2dx = f2 * ( (-1-(1/k)) * ------  +  0 ) 
     *                                 z 
     *                                  
     *   df2dk = f2 * z^(-1-(1/k)) * ( (-1-(1/k)) * dzdk/z  +  (1/k^2) * ln(z) )
     *       
     *                                   dzdsigma
     *   df2dsigma = f2 * ( (-1-(1/k)) * --------  +  0 )
     *                                      z
     *   
     *                                 dzdmu
     *   df2dmu = f2 * ( (-1-(1/k)) * --------  +  0 )
     *                                   z
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

        // f1 = exp(-1*(z^(-1/k)))
        // f2 = z^(-1-(1/k))
        // df1dx     = f1 * (1/k) * z^(-1 - (1/k)) * dzdx
        // df2dx     = f2 * ( (-1-(1/k)) * (dzdx/z) )

        double f1 = Math.exp( -1.*Math.pow(z, -1./k) );
        double f2 = Math.pow(z, (-1. - (1./k)) );
      
        double dzdx = k/sigma;

        double df2dx = f2 * (  (-1. - (1./k)) * ( dzdx/ z) );
        
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
        
        // TODO:  recheck the math in this one
        //        and decide on tests.
        
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
    
    public static double derivWRTK(float yConst, float mu, float k, float sigma, float x) {

        return calculateDerivUsingDeltaK(yConst, mu, k, sigma, x);
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
    protected static Double calculateDerivUsingDeltaK(float yConst, float mu, float k, float sigma, float x) {
        
        float deltaK = 0.01f*k;
        
        Double d0 = GeneralizedExtremeValue.generateYGEV(x, (k - deltaK), sigma, mu);
        
        Double d1 = GeneralizedExtremeValue.generateYGEV(x, (k), sigma, mu);
        
        Double d2 = GeneralizedExtremeValue.generateYGEV(x, (k + deltaK), sigma, mu);
        
        if (d0 != null && d1 != null && d2 != null) {
            
            double delta0 = d1.doubleValue() - d0.doubleValue();
            
            double delta1 = d2.doubleValue() - d1.doubleValue();
            
            double d = (delta0 + delta1)/2.;
            
            return (d/deltaK);
        
        } else if (d1 != null && d2 != null) {
                        
            double delta = d2.doubleValue() - d1.doubleValue();
                        
            return (delta/deltaK);
            
        } else if (d0 != null && d1 != null) {
            
            double delta = d1.doubleValue() - d0.doubleValue();
            
            return (delta/deltaK);
            
        } else if (d0 != null && d2 != null) {
            
            double delta = d2.doubleValue() - d0.doubleValue();
            
            return (delta/deltaK);
            
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
        
        double df2dsigma = f2 * (  (-1. - (1./k)) * ( dzdsigma/ z)  );
        
        double df1dsigma =  f1 * (1./k) * dzdsigma;
        if (z < 0) {
            df1dsigma *= -1. * Math.pow(-1.*z, -1. - (1./k));
        } else {
            df1dsigma *= Math.pow(z, -1. - (1./k));
        }
        
        double dydk = (yConst/sigma) * ( f1 * df2dsigma + f2 * df1dsigma );
        
        return dydk;
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
        
        double df2dmu = f2 * (  (-1. - (1./k)) * ( dzdmu / z)  );
        
        double df1dmu = f1 * (1/k) * dzdmu;
        if (z < 0) {
            df1dmu *= -1. * Math.pow(-1.*z, (-1 - (1/k)));
        } else {
            df1dmu *= Math.pow(z, (-1 - (1/k)));
        }
        
        double dydmu = (yConst/sigma) * ( f1 * df2dmu + f2 * df1dmu );
        
        return dydmu;
    }
    
}
