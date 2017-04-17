package algorithms.util;

import algorithms.misc.Misc;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import thirdparty.brendano.LBFGS.Helper.FunctionPolyML;
import thirdparty.brendano.LBFGS.LBFGS;
import thirdparty.brendano.LBFGS.LBFGS.Function;
import thirdparty.brendano.LBFGS.LBFGS.Params;
import thirdparty.brendano.LBFGS.LBFGS.ProgressCallback;
import thirdparty.brendano.LBFGS.LBFGS.Result;
import thirdparty.brendano.LBFGS.LBFGS.Status;

/**
 *
 * @author nichole
 */
public class PolynomialFitterTest extends TestCase {
    
    public PolynomialFitterTest() {
    }
    
    public void testMisc_generate() {

        double[] x2 = new double[]{0, 1, 2,  3,  4,  5, 
            6,    7,   8,  9,   10};
        double[] y2 = new double[]{1, 6, 17, 34, 57, 86, 
            121, 162, 209, 262, 321};

        double[] yGen2 = Misc.generate(new double[]{3., 2., 1.}, x2);

        //System.out.println("yGen=" + Arrays.toString(yGen2));
        
        for (int i = 0; i < y2.length; ++i) {
            double a = y2[i];
            double b = yGen2[i];
            double d = a - b;
            //System.out.format("%.3f - %.3f = %.3f\n", (float) a, (float) b, (float) d);
            assertTrue(Math.abs(d) < 0.0001);
        }
        
        // ------
        
        float[] x = new float[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        float[] y = new float[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
        float[] yGen = Misc.generate(new float[]{3.f, 2.f, 1.f}, x);

        //System.out.println("yGen=" + Arrays.toString(yGen));
        
        for (int i = 0; i < y2.length; ++i) {
            float a = y[i];
            float b = yGen[i];
            float d = a - b;
            //System.out.format("%.3f - %.3f = %.3f\n", (float) a, (float) b, (float) d);
            assertTrue(Math.abs(d) < 0.0001);
        }
    }

    public void test0() throws Exception {
             
        //NOTE: running this test w/ the dlib versions of lbfgs shows
        //   that it passes for search strategy 5,
        //   but not search strategy 20, and the results
        //   below are the same as the later.
        
        // test from: http://rosettacode.org/wiki/Polynomial_Fitting
        //3 x2 + 2 x + 1
        float[] x = new float[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        float[] y = new float[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
        
        PolynomialFitter polyFitter = new PolynomialFitter();
        float[] coef = polyFitter.solveOLS(x, y, 2);
                
        assertNotNull(coef);
        
        //System.out.println("coef=" + Arrays.toString(coef));
        
        assertTrue(Math.abs(coef[0] - 1) < 0.01);
        assertTrue(Math.abs(coef[1] - 2) < 0.01);
        assertTrue(Math.abs(coef[2] - 3) < 0.01);
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < x.length; i++) {
            PairInt p = new PairInt((int)x[i], (int)y[i]);
            points.add(p);
        }
        
        PolynomialFitter.plotFit(coef, points, 20, 400, 12, "test");
        
        double resid = PolynomialFitter.calcResiduals(coef, points);
        
        assertTrue(Math.abs(resid) < 0.01);
        
        // -----------
       
        double[] x2 = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double[] y2 = new double[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
       
        //NOTE: this is a local search method. 
        // -- need to add standard tests for LBFGS
        
        // would be good to compare to a global search method
        //   random search, a genetic algorithm, simulated annealing, 
        //   or particle swarm)
        
        //double[] init = new double[]{3, 2, 1};// ans is 3, 2, 1
        //double[] init = new double[]{3.1, 1.3, 1.06};
        double[] init = new double[]{1, 1, 1};
        Function f = new FunctionPolyML(x2, y2, init);
        
		Params p = new Params();
		ProgressCallback cb = new ProgressCallback() {
			@Override
			public int apply(double[] x, double[] g, double fx, double xnorm,
					double gnorm, double step, int n, int k, Status ls) {
				System.out.printf("ITER %d obj=%g sol=%.6g\n", k, fx, x[0]);
				return 0;
			}
		};
       
        //Result r = LBFGS.lbfgsNice(init, 1000, f, cb);
        
        //Result r = LBFGS.lbfgs(init, f);
        Result r = LBFGS.lbfgs(init, f, cb);
        
        System.out.println("lbfgs coef = " + 
            Arrays.toString(init));
        
        // 3 x2 + 2 x + 1
        //assertTrue(Math.abs(init[0] - 3) < 0.1);
        //assertTrue(Math.abs(init[1] - 2) < 0.2);
        //assertTrue(Math.abs(init[2] - 1) < 0.3);
        
    }
}
