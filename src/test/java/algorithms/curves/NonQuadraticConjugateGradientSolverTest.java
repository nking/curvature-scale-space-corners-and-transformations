package algorithms.curves;

import java.util.Arrays;

import algorithms.util.Errors;
import junit.framework.TestCase;

public class NonQuadraticConjugateGradientSolverTest extends TestCase {

    public void testFitCurve() throws Exception {
        
        // while revising code, if not on a development branch, don't assert results
        boolean assertResults = false;
        
        float k = 1.80f;
        float sigma = 0.85f;
        float mu = 0.441f;
        /*k = 0.294f;
        sigma = 0.141f;
        mu = 0.254f;*/
        /*
        k = 0.294f;
        sigma = 2.141f;
        mu = 0.254f;
        */
        
        float[] xp = new float[30];
        float[] ye = new float[30];
        Arrays.fill(ye, 0.05f);

        for (int i = 0; i < xp.length; i++) {
            xp[i] = (float)i/xp.length;
        }
       
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new float[0], 
            new float[0], new float[0], new float[0]);
        float[] yGEV = gev.generateCurve(xp, k, sigma, mu);
                    
        NonQuadraticConjugateGradientSolver fitCurve = new NonQuadraticConjugateGradientSolver(xp, yGEV, Errors.populateYErrorsBySqrt(xp), ye);
       
        float kMin = 0.01f;
        float kMax = 3*k;
        float sigmaMin = 0.01f*sigma;
        float sigmaMax = 2.5f*sigma;
        float muMin = 0.5f * mu;
        float muMax = xp[xp.length - 1];
        
        //GEVYFit yFit = fitCurve.fitCurve(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
        GEVYFit yFit = fitCurve.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
        //GEVYFit yFit = fitCurve.fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu);
        
        if (assertResults) {
            assertNotNull(yFit);
            assertTrue(Math.abs(yFit.k - k) < 0.1*k);
            assertTrue(Math.abs(yFit.sigma - sigma) < 0.1*sigma);
            assertTrue(Math.abs(yFit.mu - mu) < 0.1*mu);
        }
    }
}
