package thirdparty.brendano.LBFGS;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import java.util.Arrays;
import thirdparty.brendano.LBFGS.LBFGS.Function;

/**
 *
 * @author nichole
 */
public class Helper {
   
    /**
     * NOTE: NOT READY FOR USE.  THE Gradient calculation needs to be
     * improved.
     * a function that calculates the negative of the log likelihood
     * useful for 2D curve fitting.
     */
    public static class FunctionPolyML implements Function {

        final double[] xp;
        final double[] yp;
        final double[] initVars;

        //TODO: consider a constructor that accepts errors for the points
        
        public FunctionPolyML(final double[] xPoints, double[] yPoints,
            double[] init) {
            if (xPoints.length != yPoints.length) {
                throw new IllegalArgumentException(
                    "xPoints and yPoints must be same length");
            }
            this.xp = xPoints;
            this.yp = yPoints;
            this.initVars = init;
        }

        @Override
        public double evaluate(double[] vars,
            double[] outputGradient, int nVars, double step) {
            
            double[] y = new double[xp.length];

            Arrays.fill(outputGradient, 0.);
            
  //TODO: there is still an error here          
            
            /*
            calculating the gradient of the polynomial coefficients
           
            example function generated by 3 coefficients:
               gen = c0 * x^2 + c1 * x^1 + c2

               d(gen)/d(c0) = x^2
               d(gen)/d(c1) = x^1
               d(gen)/d(c2) = 1

               d(c0)= d(gen)/x^2
               d(c1)= d(gen)/x
               d(c2)= d(gen)

               gen = c0 * x^2 + c1 * x^1 + c2
               d(gen)/d(x) = 2 * c0 * x + c1
               d(gen) = d(x) * (2 * c0 * x + c1)

               d(c0)= d(gen)/x^2
                    = d(x) * (2 * c0 * x + c1) / x^2
               d(c1)= d(gen)/x
                    = d(x) * (2 * c0 * x + c1) / x
               d(c2)= d(gen)/x
                    = d(x) * (2 * c0 * x + c1)
            */

            double[] gen = Misc.generate(vars, xp);

            double[] diff = MatrixUtil.subtract(gen, yp);

            double sumDiff = 0;
            
            for (int i = 0; i < xp.length; ++i) {
                double x2 = 1;
                double dx = diff[i];
                double d = dPolydX(vars, xp[i]);
                for (int j = vars.length - 1; j > -1; j--) {
                    int varIdx = vars.length - j - 1;
                    outputGradient[varIdx] += (dx * d / x2);                      
                    x2 *= xp[i];
                    if (x2 == 0.0) {
                        break;
                    }
                }
                sumDiff += (dx * dx);
            }
            sumDiff = Math.sqrt(sumDiff);
            
            for (int j = 0; j < vars.length; ++j) {
                outputGradient[j] /= xp.length;
            }
                
            // TODO: revisit this to consider including errors given
            //   to the code for each point.
            
            double[] mnAndStDv = MiscMath.getAvgAndStDev(diff);
            double sigma = mnAndStDv[1] * mnAndStDv[1];
            
            double f = Math.pow((1.0/(2.0*Math.PI*sigma)), (vars.length/2))
                * Math.exp(-1. * sumDiff/(2. * sigma));
    
            double lnf = -2. * Math.log(f);
            
            //System.out.println(" lnf=" + lnf);
         
            return lnf;
        }
    }
    
    /**
     * assuming that the polynomial coefficients coeff are given from
     * highest order to lowest, return the derivative of y with
     * respect to x.
     * e.g. for y = coeff[0] * x^2 + coeff[1] * x^1 + coeff[2]
     * it returns dydx = 2 * coeff[0] * x + coeff[1]
     * @param coeff
     * @param x
     * @return 
     */
    public static double dPolydX(double[] coeff, double x) {
        /*
        y = c0 * x^2 + c1 * x^1 + c2
        */
        double sum = 0;
        for (int order = (coeff.length - 1); order > 0; --order) {
            sum += dPolydX(order, coeff[order], x);
        }
        
        return sum;
    }
    
    private static double dPolydX(int order, double coeff, double x) {
        
        if (order == 0) {
            return 0;
        }
        
        double dydx = coeff * (double)order * Math.pow(x, order - 1);
        
        return dydx;
    }
}
