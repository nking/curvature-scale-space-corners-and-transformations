package algorithms.search.global;

import algorithms.random.AFunction;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class QMCLBFGSTest extends TestCase {
    
    public QMCLBFGSTest() {
    }

    public void test0() {
        
        int nPoints = 10;
        double[] bounds0 = new double[]{0, 0, 0};
        double[] bounds1 = new double[]{5, 5, 5};
        
        AFunction function = new Rastrigin(3, 3);
        
        QMCLBFGS srch = new QMCLBFGS();
        
        double[] r = srch.search(function, nPoints, bounds0, bounds1);
        
        System.out.println("final results=" + Arrays.toString(r));
    }
    
    public class Rastrigin extends AFunction {
        public final double axisratio = 1;
        public final double amplitude = 10;
        public final double eps = 1.e-7;
        
        public Rastrigin(int nParams1, int nParams2) {
            super(nParams1, nParams2);
        }

        /**
         function is either from COCO or CMA-es
         need to look up the reference.
          
         COCO Numerical Black-Box Optimization
         Benchmarking Framework.
         https://github.com/numbbo/coco/blob/master/code-experiments/src
         which has
         Copyright (c) 2013-2016 by the NumBBO/CoCO team. See AUTHORS file and
        remarks below for exceptions and more details.

         and CMAES is
             https://www.lri.fr/~hansen/CMAExample1.java
        
         * @param coeffs
         * @return 
         */
        @Override
        public double f(double[] coeffs) {
            
            if (coco_vector_contains_nan(coeffs)) {
                return Double.NaN;
            }

            int d = coeffs.length;
            double sum = 0;
            for (int ii = 1; ii < d; ++ii) {
                double xi = coeffs[ii];
                sum += (xi*xi - amplitude * Math.cos(2 * Math.PI * xi));
            }
            sum += amplitude * d;
            return sum;
        }

        /**
        finite difference method to calculate the local gradient.
        adapted from dlib optimization.h
        Copyright (C) 2008  Davis E. King (davis@dlib.net)
        License: Boost Software License   See LICENSE.txt for the full license.
        */
        @Override
        public double[] der(double[] coeffs) {
            
            //System.out.println("a1  x.size=" + coeffs.length);
            
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

                //NOTE: newton's method would continue with:
                // x_(i+1) = x_i - (delta_plus/der(i))

                // and finally restore the old value of this element
                e[i] = old_val;
            }
            
            return der;
        }
    }
    
    static boolean coco_vector_contains_nan(double[] aa) {

        for (double a : aa) {
            if (Double.isNaN(a)) {
                return true;
            }
        }
        
        return false;
    }
}
