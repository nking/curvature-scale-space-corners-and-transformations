package algorithms.util;

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

    public void test0() throws Exception {
                
        // test from: http://rosettacode.org/wiki/Polynomial_Fitting
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
        
        polyFitter.plotFit(coef, points, 20, 400, 12, "test");
        
        double resid = polyFitter.calcResiduals(coef, points);
        
        assertTrue(Math.abs(resid) < 0.01);
        
        // -----------
       
        double[] x2 = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double[] y2 = new double[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
       
        double[] init = new double[]{1, 1, 1};// ans is 1, 2, 3
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
        
        Result r = LBFGS.lbfgs(init, f);
        
        System.out.println("lbfgs coef = " + 
            Arrays.toString(init));
        
        // gradient should be imprived...
        assertTrue(Math.abs(init[2] - 1) < 1.5*0.1);
        assertTrue(Math.abs(init[1] - 2) < 1.5*0.2);
        assertTrue(Math.abs(init[0] - 3) < 1.5*0.3);
        
    }
}
