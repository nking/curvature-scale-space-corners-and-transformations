package thirdparty.brendano.LBFGS;

import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import java.util.Arrays;
import thirdparty.brendano.LBFGS.LBFGS.Function;
import thirdparty.brendano.LBFGS.LBFGS.Result;

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

            //TODO: reread section on vars and nVars
            
            double[] y = new double[xp.length];

            Arrays.fill(outputGradient, 9.);
            
            int n0 = vars.length;
            int n = n0 - 1;
            
            double[] xsub = Arrays.copyOf(xp, xp.length);
            for (int i = 0; i < n; ++i) {
                int varIdx = n - i - 1;
                double c = vars[varIdx];
                for (int j = 0; j < xp.length; ++j) {
                    double t = c * xsub[j];
                    y[j] += t;
                    
             //TODO: correct this
                    outputGradient[varIdx] += ((t - yp[j]) * xsub[j]);
                    
                    xsub[j] *= xp[j];
                }
                outputGradient[varIdx] /= xp.length;
                outputGradient[varIdx] -= c;
            }
            for (int j = 0; j < xp.length; ++j) {
                y[j] += vars[n];
                outputGradient[n] += (vars[n] - yp[j]);
            }
            outputGradient[n] /= xp.length;
            outputGradient[n] -= vars[n];
            
            double[] diffByData = new double[xp.length];
            double sumDiff = 0; 
            for (int j = 0; j < xp.length; ++j) {
                diffByData[j] = yp[j] - y[j];
                sumDiff += (diffByData[j] * diffByData[j]);
            }
            
            // TODO: revisit this to consider including errors given
            //   to the code for each point.
            
            double[] mnAndStDv = MiscMath.getAvgAndStDev(diffByData);
            double sigma = mnAndStDv[1] * mnAndStDv[1];
            
            double f = Math.pow((1.0/(2.0*Math.PI*sigma)), (vars.length/2))
                * Math.exp(-1. * sumDiff/(2. * sigma));
    
            double lnf = -2. * Math.log(f);
  System.out.println("vars=" + Arrays.toString(vars) + 
      " lnf=" + lnf + " sd=" + sumDiff
  + " \ngradient=" + Arrays.toString(outputGradient));    
            return lnf;
        }
    }
    
    /*
    public static class Result2 extends Result {

        public double[] coeff;

        public Result2(Result r, double[] finalCoeff) {
            super(r.status);
            this.coeff = finalCoeff;
        }

        public double[] getCoeff() {
            return coeff;
        }
        
    }*/
}
