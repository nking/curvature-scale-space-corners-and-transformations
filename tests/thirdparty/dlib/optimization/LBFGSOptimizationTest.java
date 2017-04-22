package thirdparty.dlib.optimization;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LBFGSOptimizationTest extends TestCase {
    
    boolean debug = true;
    
    public LBFGSOptimizationTest() {
    }
    
    public void testLBFGS_poly() {
        
        double[] x2 = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double[] y2 = new double[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
       
        double[] init = new double[]{1, 1, 1};
        double[] expected = new double[]{3., 2., 1.};
        
        Helper.FunctionPoly f = new Helper.FunctionPoly(x2, y2);
        
        // use an empirical derivative derived from finite differences method
        Helper.CentralDifferences f2 = new Helper.CentralDifferences(f, 1.e-7);
                
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(1.e-5, 100);
        
        double fLower = -10;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, 
            stopStrategy, f2, init, fLower);
       
        System.out.println("min=" + min + " \n   coeffs=" +
            Arrays.toString(init));
        
    }
    
    /**
     * unit test ported from gsl project file
     * multimin/test.c
     * downloaded from 
     * http://mirror.team-cymru.org/gnu/gsl/
     * the code has license 
     * GNU General Public License.
     * https://www.gnu.org/software/gsl/
     * http://www.gnu.org/copyleft/gpl.html
     * 
     */
    public void test2() {
        
        rosenbrock();
        rosenbrock1();
        
        roth();
        wood();
        simpleAbs();
       
    }
   
    private void roth() {
        
        double[] coeffs = new double[] {4.5, 3.5};
        
        Roth f = new Roth();
                
        // use an empirical derivative derived from finite differences method
        //Helper.CentralDifferences f2 = new Helper.CentralDifferences(f, 1.e-7);
                
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(1.e-5, 100);
        
        double fLower = -10;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, 
            stopStrategy, f, coeffs, fLower);
        
        if (debug)
        System.out.println("roth coeffs=" + Arrays.toString(coeffs));
        
        assertTrue(Math.abs(coeffs[0] - 5.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 4.) < 0.1);
        
    }
    
    private void wood() {
        
        double[] coeffs = new double[] {-3.0, -1.0, -3.0, -1.0};
        
        Wood f = new Wood();
        
        // use an empirical derivative derived from finite differences method
        //Helper.CentralDifferences f2 = new Helper.CentralDifferences(f, 1.e-7);
                
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(1.e-5, 100);
        
        double fLower = -10;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, 
            stopStrategy, f, coeffs, fLower);
        
        
        if (debug)
        System.out.println("wood coeffs=" + Arrays.toString(coeffs));
        
        // CHECK THESE
        //assertTrue(Math.abs(coeffs[0] - -1.) < 0.1);
        //assertTrue(Math.abs(coeffs[1] - 1.) < 0.1);
        //assertTrue(Math.abs(coeffs[2] - -1.) < 0.1);
        //assertTrue(Math.abs(coeffs[3] - 1.) < 0.1);
    
        assertTrue(Math.abs(coeffs[0] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[2] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[3] - 1.) < 0.1);
    
    }
    
    private void rosenbrock() {
        
        double[] coeffs = new double[] {-1.2, 1.0};
        
        Rosenbrock f = new Rosenbrock();
        
        // use an empirical derivative derived from finite differences method
        //Helper.CentralDifferences f2 = new Helper.CentralDifferences(f, 1.e-7);
                
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(1.e-5, 100);
        
        double fLower = -10;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, 
            stopStrategy, f, coeffs, fLower);
        
        
        if (debug)
        System.out.println("rb coeffs=" + Arrays.toString(coeffs));
        
        assertTrue(Math.abs(coeffs[0] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 1.) < 0.1);
    }
    
    private void rosenbrock1() {
        double[] coeffs = new double[] {2., 2.0};
        
        Rosenbrock f = new Rosenbrock();
        
        // use an empirical derivative derived from finite differences method
        //Helper.CentralDifferences f2 = new Helper.CentralDifferences(f, 1.e-7);
                
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(1.e-5, 100);
        
        double fLower = -10;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, 
            stopStrategy, f, coeffs, fLower);
        
        
        if (debug)
        System.out.println("rb1 coeffs=" + Arrays.toString(coeffs));
        
        assertTrue(Math.abs(coeffs[0] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 1.) < 0.1);
    }
    
    private void simpleAbs() {
        
        double[] coeffs = new double[] {1., 2.0};
        
        SimpleAbs f = new SimpleAbs();
        
        // use an empirical derivative derived from finite differences method
        //Helper.CentralDifferences f2 = new Helper.CentralDifferences(f, 1.e-7);
                
        LBFGSSearchStrategy searchStrategy = new LBFGSSearchStrategy(5);
        ObjectiveDeltaStopStrategy stopStrategy 
            = new ObjectiveDeltaStopStrategy(1.e-5, 100);
        
        double fLower = -10;
        
        LBFGSOptimization opt = new LBFGSOptimization();
        double min = opt.findMin(searchStrategy, 
            stopStrategy, f, coeffs, fLower);
        
        
        if (debug)
        System.out.println("sa coeffs=" + Arrays.toString(coeffs));
        
        assertTrue(Math.abs(coeffs[0] - 1.) < 0.1);
        assertTrue(Math.abs(coeffs[1] - 2.) < 0.1);
    }

    private static class Roth implements LBFGSOptimization.IFunction {
        
        public Roth() {}

        @Override
        public double f(double[] coeffs) {
            
            double u = coeffs[0];
            double v = coeffs[1];
            double a = -13.0 + u + ((5.0 - v) * v - 2.0) * v;
            double b = -29.0 + u + ((v + 1.0) * v - 14.0) * v;
            
            double c = -2 + v * (10 - 3 * v);
            double d = -14 + v * (2 + 3 * v);
            
            double m = a * a + b * b;
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return m;
        }

        @Override
        public double f(double d) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public double[] der(double[] coeffs) {
            
            double[] outputGradient = new double[coeffs.length];
            
            double u = coeffs[0];
            double v = coeffs[1];
            double a = -13.0 + u + ((5.0 - v) * v - 2.0) * v;
            double b = -29.0 + u + ((v + 1.0) * v - 14.0) * v;
            
            double c = -2 + v * (10 - 3 * v);
            double d = -14 + v * (2 + 3 * v);
            
            outputGradient[0] = 2 * a + 2 * b;
            outputGradient[1] = 2 * a * c + 2 * b * d;
       
            double m = a * a + b * b;
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return outputGradient;
        }

        @Override
        public double der(double d) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }
    
    private static class Wood implements LBFGSOptimization.IFunction {
        
        public Wood() {}

        @Override
        public double f(double[] coeffs) {
            
            double u1 = coeffs[0];
            double u2 = coeffs[1];
            double u3 = coeffs[2];
            double u4 = coeffs[3];

            double t1 = u1 * u1 - u2;
            double t2 = u3 * u3 - u4;
        
            double m = 100 * t1 * t1 + (1 - u1) * (1 - u1)
                + 90 * t2 * t2 + (1 - u3) * (1 - u3)
                + 10.1 * ((1 - u2) * (1 - u2) + (1 - u4) * (1 - u4))
                + 19.8 * (1 - u2) * (1 - u4);
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return m;
        }

        @Override
        public double f(double d) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public double[] der(double[] coeffs) {
            
            double[] outputGradient = new double[coeffs.length];
         
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
        
            double m = 100 * t1 * t1 + (1 - u1) * (1 - u1)
                + 90 * t2 * t2 + (1 - u3) * (1 - u3)
                + 10.1 * ((1 - u2) * (1 - u2) + (1 - u4) * (1 - u4))
                + 19.8 * (1 - u2) * (1 - u4);
            
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return outputGradient;
        }

        @Override
        public double der(double d) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    }
    
    private static class Rosenbrock implements LBFGSOptimization.IFunction {
        
        public Rosenbrock() {}

        @Override
        public double f(double[] coeffs) {
            double u = coeffs[0];
            double v = coeffs[1];
            double b = u * u - v;            
            double a = u - 1;
            return a * a + 10 * b * b;
        }
        
        @Override
        public double f(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public double[] der(double[] coeffs) {
            
            double[] der = new double[coeffs.length];
            
            double u = coeffs[0];
            double v = coeffs[1];
            double b = u * u - v;
            der[0] = 2 * (u - 1) + 40 * u * b;
            der[1] = -20 * b;
            
            return der;
        }
        
        @Override
        public double der(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
       
    }
    
    private static class SimpleAbs implements LBFGSOptimization.IFunction {
        
        public SimpleAbs() {}

        @Override
        public double f(double[] coeffs) {
            
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
            
            double m = Math.abs(a) + Math.abs(b);
        
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return m;
        }
        
        @Override
        public double f(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public double[] der(double[] coeffs) {
            
            double[] der = new double[coeffs.length];
            
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
            
            der[0] = sign0;
            der[1] = sign1;
            
            double m = Math.abs(a) + Math.abs(b);
        
            /*
            System.out.println("==>vars=" + Arrays.toString(coeffs));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>gradient=" + Arrays.toString(outputGradient));
            System.out.println("==>m=" + m);
            */
            
            return der;
        }
        
        @Override
        public double der(double a) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        
    }
}
