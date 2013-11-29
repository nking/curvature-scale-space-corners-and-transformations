package algorithms.curves;

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
                
                                                  -1*dzdk       
                 df(k)dk = -z^(-1/k) * ( (-1/k) * -------  +  (1/k^2) * ln(-z) )
     *                                               z
     *   
     *   deriv of -z^(-1/k) w.r.t. sigma is
     *                                                  -1*dzdsigma
     *           df(sigma)dsigma = -z^(-1/k) * (-1/k) * ----------- + 0
     *                                                       z
     *   
     *   deriv of -z^(-1/k) w.r.t. mu is            -1*dzdmu
     *           df(mu)mu =  -z^(-1/k) * ( (-1/k) * ----------- + 0 )
     *                                                  z
     *   
     *   
     *   Let f1 = the first exponential in y
     *          = exp(-z^(1/k))
     *   
     *   df1dx = f1 * d(-z^(1/k))dx = f1 * (1/k) * z^(-1 - (1/k)) * dzdx
     *   
     *                                      -1*dzdk
     *   df1dk = f1 * -z^(-1/k) * ( (-1/k) * -------  +  (1/k^2) * ln(-z) )
     *                                          z
     *                                  
     *   df1dsigma = f1 * (1/k) * z^(-1 - (1/k)) * dzdsigma
     *                                             
     *   
     *   df1dmu = f1 * (1/k) * z^(-1 - (1/k)) * dzdmu
     *   
     *   
     *   
     *   Let f2 = z^(-1-(1/k))
     *                                dzdx
     *   df2dx = f2 * ( (-1-(1/k)) * ------  +  0 ) 
     *                                 z 
     *                                     
     *                                dzdk
     *   df2dk = f2 * ( (-1-(1/k)) * -------  +  (-1-(1/k)) * ln(z) )
     *                                  z
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
     *   
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
        
        double f1 = Math.exp( Math.pow(-1.*z, 1./k) );
        
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
    public static double derivWRTK(float yConst, float mu, float k, float sigma, float x) {
        
        double z = 1. + k *( (x-mu)/sigma );
        
        double f1 = Math.exp( Math.pow(-1.*z, 1./k) );
        
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        double dzdk = (x-mu)/sigma;
        
        double df2dk = f2 * (  (-1. - (1./k)) * ( dzdk/ z)  +  (-1. - (1./k)) * Math.log(z) );
        
        double df1dk = f1 * Math.pow(-1*z, (1./k)) * ( (-1./k) * (-1.*dzdk/z)  +  (1./k*k) * Math.log(-z) );
        
        double dydk = (yConst/sigma) * ( f1 * df2dk + f2 * df1dk );
        
        return dydk;
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
        
        double f1 = Math.exp( Math.pow(-1.*z, 1./k) );
        
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        double dzdsigma = -1. * k * (x-mu) * Math.pow(sigma, -2.);
        
        double df2dsigma = f2 * (  (-1. - (1./k)) * ( dzdsigma/ z)  );
        
        double df1dsigma = f1 * (-1./k) * Math.pow(-1*z, -1. - (1./k)) * dzdsigma;
        
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
        
        double f1 = Math.exp( Math.pow(-1.*z, 1./k) );
        
        double f2 = Math.pow(z, (-1. - (1./k)) );
        
        double dzdmu = -1. * k/sigma;
        
        double df2dmu = f2 * (  (-1. - (1./k)) * ( dzdmu/ z)  );
        
        double df1dmu = f1 * Math.pow(-1*z, (1./k)) * ( (-1./k) * (-1.*dzdmu/z) );
        
        double dydmu = (yConst/sigma) * ( f1 * df2dmu + f2 * df1dmu );
        
        return dydmu;
    }
    
}
