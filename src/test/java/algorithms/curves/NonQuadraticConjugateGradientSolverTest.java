package algorithms.curves;

import java.security.SecureRandom;
import java.util.Arrays;

import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

public class NonQuadraticConjugateGradientSolverTest extends TestCase {

    public void estFitCurve() throws Exception {
        
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
       
        // problem fiting curve: k=0.91, s=0.22, m=.218
        
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
    
    public void testFitRandomCurves() throws Exception {
        
        System.out.println("testFitRandomCurves");
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
        
        int n = 40;
        int nX = 40;
        
        // do for number of x points being 40, 30, 20, 10
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed );
        
        System.out.println("SEED=" + seed);
        
        float[] kRange = new float[]{0.00001f, 2.0f};
        float[] sRange = new float[]{0.025f, 0.5f};
        
        String path = null;
        
        while (nX > 9) {
            
            float[] xp = new float[nX];
            float[] ye = new float[nX];
            Arrays.fill(ye, 0.05f);

            for (int i = 0; i < xp.length; i++) {
                xp[i] = (float)i/xp.length;
            }
            float[] xe = Errors.populateYErrorsBySqrt(xp);
            
            float[] mRange = new float[]{xp[0], 0.3f };
            
            for (int i = 0; i < n; i++) {
                
                float k = kRange[0] + sr.nextFloat()*(kRange[1] - kRange[0]);
                float s = sRange[0] + sr.nextFloat()*(sRange[1] - sRange[0]);
                float m = mRange[0] + sr.nextFloat()*(mRange[1] - mRange[0]);
                
                float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, s, m);
                
                NonQuadraticConjugateGradientSolver solver = 
                    new NonQuadraticConjugateGradientSolver(xp, yGEV, xe, ye);

                GEVYFit fit = solver.fitCurveParametersSeparately(kRange[0], kRange[1], sRange[0], sRange[1], mRange[0], mRange[1]);
                
                assertNotNull(fit);
                
                String label = String.format(
                    "k=%4.4f <%2.2f>  s=%4.4f <%2.2f>  m=%4.4f <%2.3f>  (nx=%d,i=%d)  chi=%4.6f", 
                    fit.getK(), k, 
                    fit.getSigma(), s, 
                    fit.getMu(), m, 
                    nX, i, fit.getChiSqSum());
                
                plotter.addPlot(xp, yGEV, xe, ye, fit.getX(), fit.getYFit(), label);
                    
                path = plotter.writeFile2();
            }
            
            nX >>= 1;
        }
        System.out.println(" plot is at path=" + path);
    }
}
