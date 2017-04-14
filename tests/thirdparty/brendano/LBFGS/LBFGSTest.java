package thirdparty.brendano.LBFGS;

import java.util.Arrays;
import junit.framework.TestCase;
import thirdparty.brendano.LBFGS.LBFGS.Function;
import thirdparty.brendano.LBFGS.LBFGS.Result;

/**
 *
 * @author nichole
 */
public class LBFGSTest extends TestCase {

    /**
     * unit test ported from gsl project file
     * multimin/test.c
     * donloaded from 
     * http://mirror.team-cymru.org/gnu/gsl/
     * the code has license 
     * GNU General Public License.
     * https://www.gnu.org/software/gsl/
     * http://www.gnu.org/copyleft/gpl.html
     * 
     */
    public void test0() {
        
        //TODO: port the minimizer test of the GSL project method
        // gsl_multimin_fdfminimizer_vector_bfgs
        // in multimin/test.c
        
        // check BFGS
        rosenbrock();
        rosenbrock1();
        
        //roth();
        //wood();
        //simpleAbs();
       
    }
   
    private void roth() {
        
        double[] coeffs = new double[] {4.5, 3.5};
        
        Roth f = new Roth();
        
        Result r = LBFGS.lbfgs(coeffs, f);
        
        System.out.println("roth coeffs=" + Arrays.toString(coeffs));
        
        //assertTrue(Math.abs(coeffs[0] - 2.) < 0.1);
        //assertTrue(Math.abs(coeffs[1] - 0.) < 0.1);
        
    }
    
    private void wood() {
        
        double[] coeffs = new double[] {-3.0, -1.0, -3.0, -1.0};
        
        //4, 0
        Wood f = new Wood();
        
        Result r = LBFGS.lbfgs(coeffs, f);
        
        System.out.println("wood coeffs=" + Arrays.toString(coeffs));
        
        //assertTrue(Math.abs(coeffs[0] - 4.) < 0.1);
        //assertTrue(Math.abs(coeffs[1] - 0.) < 0.1);
    }
    
    private void rosenbrock() {
        
        double[] coeffs = new double[] {-1.2, 1.0};
        
        Rosenbrock f = new Rosenbrock();
        
        Result r = LBFGS.lbfgs(coeffs, f);
        
        System.out.println("rb coeffs=" + Arrays.toString(coeffs));
        
        assertTrue(Math.abs(coeffs[0] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 1.) < 0.1);
    }
    
    private void rosenbrock1() {
        double[] coeffs = new double[] {2., 2.0};
        
        Rosenbrock f = new Rosenbrock();
        
        Result r = LBFGS.lbfgs(coeffs, f);
        
        System.out.println("rb1 coeffs=" + Arrays.toString(coeffs));
        
        assertTrue(Math.abs(coeffs[0] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 1.) < 0.1);
    }
    
    private void simpleAbs() {
        
        double[] coeffs = new double[] {1., 2.0};
        
        //2, 0
        SimpleAbs f = new SimpleAbs();
        
        Result r = LBFGS.lbfgs(coeffs, f);
        
        System.out.println("sa coeffs=" + Arrays.toString(coeffs));
        
        //assertTrue(Math.abs(coeffs[0] - 2.) < 0.1);
        //assertTrue(Math.abs(coeffs[1] - 0.) < 0.1);
    }

    private static class Roth implements Function {
        
        public Roth() {}

        @Override
        public double evaluate(final double[] coeffs,
            final double[] outputGradient, int nCoeffs, double step) {
        
            // populate gradient and return the evaluation of min
            
            double u = coeffs[0];
            double v = coeffs[1];
            double a = -13.0 + u + ((5.0 - v) * v - 2.0) * v;
            double b = -29.0 + u + ((v + 1.0) * v - 14.0) * v;
            
            double c = -2 + v * (10 - 3 * v);
            double d = -14 + v * (2 + 3 * v);
            
            outputGradient[0] = 2 * a + 2 * b;
            outputGradient[1] = 2 * a * c + 2 * b * d;
            
            return a * a + b * b;
        }
    }
    
    private static class Wood implements Function {
        
        public Wood() {}

        @Override
        public double evaluate(final double[] coeffs,
            final double[] outputGradient, int nCoeffs, double step) {
        
            // populate gradient and return the evaluation of min
            
            double u1 = coeffs[0];
            double u2 = coeffs[1];
            double u3 = coeffs[2];
            double u4 = coeffs[3];

            double t1 = u1 * u1 - u2;
            double t2 = u3 * u3 - u4;
            outputGradient[0] = 400 * u1 * t1 - 2 * (1 - u1);
            outputGradient[1] = -200 * t1 - 20.2 * (1 - u2) - 19.8 * (1 - u4);
            outputGradient[2] = 360 * u3 * t2 - 2 * (1 - u3);
            outputGradient[3] = -180 * t2 - 20.2 * (1 - u4) - 19.8 * (1 - u2);
        
            return 100 * t1 * t1 + (1 - u1) * (1 - u1)
                + 90 * t2 * t2 + (1 - u3) * (1 - u3)
                + 10.1 * ((1 - u2) * (1 - u2) + (1 - u4) * (1 - u4))
                + 19.8 * (1 - u2) * (1 - u4);
        }
    }
    
    private static class Rosenbrock implements Function {
        
        public Rosenbrock() {}

        @Override
        public double evaluate(final double[] coeffs,
            final double[] outputGradient, int nCoeffs, double step) {
        
            // populate gradient and return the evaluation of min
            
            double u = coeffs[0];
            double v = coeffs[1];
            double b = u * u - v;
            outputGradient[0] = 2 * (u - 1) + 40 * u * b;
            outputGradient[1] = -20 * b;
            
            double a = u - 1;
            return a * a + 10 * b * b;
        }
    }
    
    private static class SimpleAbs implements Function {
        
        public SimpleAbs() {}

        @Override
        public double evaluate(final double[] coeffs,
            final double[] outputGradient, int nCoeffs, double step) {
        
            // populate gradient and return the evaluation of min
            
            double u = coeffs[0];
            double v = coeffs[1];
            double a = u - 1;
            double b = v - 2;
  
            double sign0 = u - 1;
            if (sign0 >= 0) {
                sign0 = 1;
            } else {
                sign0 = -1;
            }
            
            double sign1 = v - 2;
            if (sign1 >= 0) {
                sign1 = 1;
            } else {
                sign1 = -1;
            }
            
            outputGradient[0] = sign0;
            outputGradient[1] = sign1;
            
            return Math.abs(a) + Math.abs(b);
        }
    }
}
