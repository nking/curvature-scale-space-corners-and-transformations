package algorithms.curves;

public class DerivGEV {

    /*
     *
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
     *   dzdsigma =   d( k * (x-mu) * sigma^-1)  =  -1 * k * (x-mu) * (sigma^-2) 
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
     *                                         -1*dzdsigma
     *   df1dsigma = f1 * -z^(-1/k) * (-1/k) * -----------
     *                                              z
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
     *   
               
     */
}
