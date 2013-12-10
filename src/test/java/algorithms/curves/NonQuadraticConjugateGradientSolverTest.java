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
       
        // TODO:  problem fitting parameters in <>  k=1.4680 <0.09> s=0.1307 <0.05> m=0.1543 <0.106> (nx=40,i=21) chi=0.013889
        
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
        
        int n = 50;
        int nX = 40;
        
        // do for number of x points being 40, 30, 20, 10
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1386620575944l;
        sr.setSeed( seed );
        
        System.out.println("SEED=" + seed);
        
        float[] kRange = new float[]{0.001f, 2.0f};
        float[] sRange = new float[]{0.025f, 1.0f};
        
        String path = null;
        
        while (nX > 9) {
            
            float[] xp = new float[nX];
            float[] ye = new float[nX];
            Arrays.fill(ye, 0.04f);

            for (int i = 0; i < xp.length; i++) {
                xp[i] = (float)i/xp.length;
            }
            float[] xe = Errors.populateYErrorsBySqrt(xp);
            
            float[] mRange = new float[]{xp[0], 0.35f };
            
            float kDiff = kRange[1] - kRange[0];
            float sDiff = sRange[1] - sRange[0];
            float mDiff = mRange[1] - mRange[0];
            
            for (int i = 0; i < n; i++) {
                
                float k = kRange[0] + sr.nextFloat()*kDiff;
                float s = sRange[0] + sr.nextFloat()*sDiff;
                float m = mRange[0] + sr.nextFloat()*mDiff;
                
                float[] yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, k, s, m);
                
                //if (nX == 40 && i == 45) {//40,17 40,20  40,24  40,25   40,44 40,45   20,3
                NonQuadraticConjugateGradientSolver solver = 
                    new NonQuadraticConjugateGradientSolver(xp, yGEV, xe, ye);

                    solver.setDebug(true);
                
                //GEVYFit fit = solver.fitCurveParametersSeparately(kRange[0], kRange[1], sRange[0], sRange[1], mRange[0], mRange[1]);
                
                GEVYFit fit = solver.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
                assertNotNull(fit);
                //}
                
                
                String label = String.format(
                    "k=%4.4f <%2.2f>  s=%4.4f <%2.2f>  m=%4.4f <%2.3f>  (nx=%d,i=%d)  chi=%4.6f", 
                    fit.getK(), k, 
                    fit.getSigma(), s, 
                    fit.getMu(), m, 
                    nX, i, fit.getChiSqSum());
                
                plotter.addPlot(xp, yGEV, xe, ye, fit.getX(), fit.getYFit(), label);
                    
                path = plotter.writeFile2();
                
                
                //assertTrue(fit.getChiSqSum() < 0.05f);
                //assertTrue(Math.abs(fit.getK() - k)     < k*0.4);
                //assertTrue(Math.abs(fit.getSigma() - s) < s*0.25);
                //assertTrue(Math.abs(fit.getMu() - m)    < 0.3);
  
                if (false) { // for print out to improve fit using NonQuadraticConjugateGradientSolverTest
                    if (nX == 40 && i == 36) {
                        StringBuilder xsb = new StringBuilder();
                        StringBuilder ysb = new StringBuilder();
                        StringBuilder xesb = new StringBuilder();
                        StringBuilder yesb = new StringBuilder();
    
                        for (int z = 0; z < xp.length; z++) {
                            if (z > 0) {
                                xsb.append("f, ");
                                ysb.append("f, ");
                                xesb.append("f, ");
                                yesb.append("f, ");
                            }
                            xsb.append(xp[z]);
                            ysb.append(yGEV[z]);
                            xesb.append(xe[z]);
                            yesb.append(ye[z]);
                        }
                        System.out.println("float[] x = new float[]{"  + xsb.append("f").toString() + "};");
                        System.out.println("float[] y = new float[]{"  + ysb.append("f").toString() + "};");
                        System.out.println("float[] xe = new float[]{" + xesb.append("f").toString() + "};");
                        System.out.println("float[] ye = new float[]{" + yesb.append("f").toString() + "};");
                    }
                }
            }
            
            nX >>= 1;
        }
        System.out.println(" plot is at path=" + path);
    }
    
    public void estAFit() throws Exception {
        
        float[] x = new float[]{0.010065701f, 0.030197103f, 0.050328504f, 0.0704599f, 0.09059131f, 0.110722706f, 0.13085411f, 0.15098551f, 0.17111692f, 0.19124833f, 0.21137972f};
        float[] y = new float[]{280.f, 328.f, 225.f, 177.f, 149.f, 123.f, 78.f, 71.f, 66.f, 40.f, 26.f};
        float[] xe = new float[]{5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f, 5.768689f}; 
        float[] ye = new float[]{0.026789403f, 0.039247185f, 0.06930564f, 0.10465008f, 0.14083745f, 0.18868062f, 0.2735232f, 0.31728983f, 0.37381834f, 0.5396733f, 0.76859f}; 
        
        NonQuadraticConjugateGradientSolver solver = 
            new NonQuadraticConjugateGradientSolver(x, y, xe, ye);
        
        solver.setDebug(true);

        GEVYFit fit = solver.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
        //GEVYFit fit = solver.fitCurveKGreaterThanZeroAllAtOnce(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
        
    }
}
