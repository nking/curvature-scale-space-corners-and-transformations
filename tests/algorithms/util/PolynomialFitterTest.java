package algorithms.util;

import algorithms.misc.Misc;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PolynomialFitterTest extends TestCase {
    
    public PolynomialFitterTest() {
    }
    
    public void estMisc_generate() {

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
       
    }
}
