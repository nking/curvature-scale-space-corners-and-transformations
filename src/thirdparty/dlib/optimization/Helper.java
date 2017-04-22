package thirdparty.dlib.optimization;

import java.util.Arrays;
import thirdparty.dlib.optimization.LBFGSOptimization.IFunction;

/**
 *
 * @author nichole
 */
public class Helper {
   
    // NOT READY FOR USE
    public static class FunctionPoly implements IFunction {

        final double[] xp;
        final double[] yp;
        
        public FunctionPoly(double[] xData, double[] yData) {
            if (xData.length != yData.length) {
                throw new IllegalArgumentException(
                    "xData and yData must be same length");
            }
            this.xp = Arrays.copyOf(xData, xData.length);
            this.yp = Arrays.copyOf(yData, yData.length);
        }
        
        @Override
        public double f(double[] coeffs) {
    
            System.out.println("aa0");
            
            /*
            double[] gen = new double[xp.length];
            generatePolynomial(coeffs, xp, gen);
            double sumDiff = 0;
            for (int i = 0; i < 11; ++i) {
                double diff = gen[i] - yp[i];
                sumDiff += (diff * diff);
            }
            sumDiff = Math.sqrt(sumDiff);
            */
            
            double[] gen = new double[xp.length];
            double[] gradient = new double[coeffs.length];
            double[] diffY = new double[xp.length];
       
            double sumDiff = calcGradient(coeffs, gen, gradient, diffY);
            
            System.out.println("poly coeffs=" + Arrays.toString(coeffs));
            System.out.println("  diff=" + sumDiff);

            return sumDiff;
        }

        @Override
        public double f(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public double[] der(double[] coeffs) {
            
            System.out.println("a0 der");
            
            double[] gen = new double[xp.length];
            double[] gradient = new double[coeffs.length];
            double[] diffY = new double[xp.length];
       
            double diffSum = calcGradient(coeffs, gen, gradient, diffY);
        
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            //System.out.print("==>diff="); printFormattedArray(outputDiffY);
            System.out.println("==>sumDiff=" + diffSum);
            System.out.print("==>gradient="); printFormattedArray(gradient);
            
            return gradient;
        }

        @Override
        public double der(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        // ----------------------
        // NOT READY for use.  need to revisit to make sure
        //   that result is convex 
        double calcGradient(double[] coeffs, 
            double[] gen, double[] outputCoeffGrad,
            double[] outputDiffY) {

            generatePolynomial(coeffs, xp, gen);

            double sumDiff = 0;
            for (int i = 0; i < 11; ++i) {
                outputDiffY[i] = gen[i] - yp[i];
                sumDiff += (outputDiffY[i] * outputDiffY[i]);
            }
            sumDiff = Math.sqrt(sumDiff);

            for (int j = 0; j < 3; ++j) {
                outputCoeffGrad[j] = 0.;
            }

            for (int i = 0; i < 11; ++i) {
                double x2 = 1;
                double dyAtX = outputDiffY[i];
                double dydx = dPolydXHL(coeffs, xp[i]);
                for (int j = 2; j > -1; j--) {
                    int varIdx = 3 - j - 1;

                    //dc/dx = dy/dx * dc/dy
                    //dc = dx * dy/dx * dc/dy
                    //outputCoeffGrad(varIdx) += (dyAtX * dydx / x2);

                    double dx = dyAtX * (xp[i]/gen[i]);
                    outputCoeffGrad[varIdx] +=
                        ( dx * dydx / x2);

                    x2 *= xp[i];

                    if (x2 == 0.0) {
                        break;
                    }
                }
            }

            for (int j = 0; j < 3; ++j) {
                outputCoeffGrad[j] /= 11.;
            }
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            //System.out.print("==>diff="); printFormattedArray(outputDiffY);
            System.out.println("==>sumDiff=" + sumDiff);
            System.out.print("==>gradient="); printFormattedArray(outputCoeffGrad);
            */
            
            return sumDiff;
        }

        void generatePolynomial(double[] coeffs, 
            double[] xPoly, double[] gen) {

            for (int i = 0; i < 11; ++i) {
                gen[i] = 0.;
            }        

            for (int i = 0; i < 11; ++i) {
                double x2 = 1.;
                for (int j = 2; j > -1; j--) {
                    double c = coeffs[j];
                    gen[i] += (c * x2);
                    x2 *= xPoly[i];
                }
            }
        }

        double dPolydX(int order, double coeff, double x) {
            if (order == 0) {
                return 0;
            }
            double dydx = coeff * (double)order * Math.pow(x, order - 1);
            return dydx;
        }
    
        double dPolydXHL(double[] coeffs,  double x) {
            double sum = 0;
            for (int i = 0; i < 3; ++i) {
                int order = 3 - i - 1;
                if (order == 0) continue;
                sum += dPolydX(order, coeffs[i], x);
            }
            return sum;
        }

        double calcStDev(double[] x) {
            double sumX = 0;
            for (int i = 0; i < 11; i++) {
                sumX += x[i];
            }
            double avgX = sumX/11.;
            sumX = 0;
            for (int i = 0; i < 11; i++) {
                double diffX = x[i] - avgX;
                sumX += (diffX * diffX);
            }
            double stdDevX = Math.sqrt(sumX/10.);
            return stdDevX;
        }
    }

    private static void printFormattedArray(double[] a) {
        System.out.print("[");
        for (double m : a) {
            System.out.format("%.3f, ", (float)m);
        }
        System.out.println("]");
    }
    
    // add the dlib reference here.  this is adapted
    //    from their optimization tests
    public static class CentralDifferences implements IFunction {
        
        private final IFunction f;
        private final double eps;
   
        public CentralDifferences(IFunction f) {
            this.f = f;
            this.eps = 1.e-7;
        }
        public CentralDifferences(IFunction f, double eps) {
            this.f = f;
            this.eps = eps;
        }
    
        @Override
        public double f(double[] coeffs) {
            return f.f(coeffs);
        }

        @Override
        public double f(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public double[] der(double[] coeffs) {

            return derivative(coeffs);
        }

        @Override
        public double der(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        
        /**
        adapted from dlib optimization.h
        Copyright (C) 2008  Davis E. King (davis@dlib.net)
        License: Boost Software License   See LICENSE.txt for the full license.
        */
        private double[] derivative(double[] coeffs) {
        
            System.out.println("a1  x.size=" + coeffs.length);
            
            int n = coeffs.length;
            
            double[] der = new double[n];
            double[] e = Arrays.copyOf(coeffs, n);
            
            for (int i = 0; i < n; ++i) {
                final double old_val = e[i];
                e[i] += eps;
                final double delta_plus = f(e);
                e[i] = old_val - eps;
                final double delta_minus = f(e);

                // finite difference:  this is the approx jacobian
                der[i] = (delta_plus - delta_minus)/(2.*eps); 

                //NOTE: newtons method would continue with:
                // x_(i+1) = x_i - (delta_plus/der(i))

                // and finally restore the old value of this element
                e[i] = old_val;
            }
            
            //NLK
            System.out.println("derivitave=" + Arrays.toString(der));
            

            return der;
        }
    };
}
