package algorithms.util;

import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

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
        float[] coef = polyFitter.solve(x, y);
                
        assertNotNull(coef);
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
                
    }
}
