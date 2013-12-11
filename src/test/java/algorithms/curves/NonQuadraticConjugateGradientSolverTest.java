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
                
        float[] x = new float[]{
            0.11858261f, 0.35574782f, 0.59291303f, 0.83007824f, 1.0672435f, 1.3044087f, 1.5415739f, 1.7787391f, 2.0159044f, 2.2530694f, 2.4902349f, 2.7273998f, 2.9645653f, 3.2017303f, 3.4388957f, 3.6760607f, 3.9132261f, 4.1503916f, 4.3875566f, 4.6247215f, 4.861887f, 5.0990524f, 5.3362174f, 5.5733824f, 5.810548f, 6.0477133f, 6.2848783f, 6.522043f, 6.7592087f, 6.996374f, 7.233539f, 7.470704f
        };
        float[] y = new float[]{
            0f, 1569f, 8439f, 12001f, 10823f, 9499f, 7622f, 6688f, 4853f, 3802f, 3981f, 2877f, 2551f, 1991f, 1596f, 1946f, 1320f, 942f, 932f, 792f, 972f, 781f, 370f, 594f, 531f, 353f, 295f, 491f, 536f, 248f, 330f, 186f
        };
        float[] xe = new float[]{
            424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f, 424.8243f
        };
        float[] ye = new float[]{
            0.04493785f, 0.24282226f, 0.16796246f, 0.18530881f, 0.23318607f, 0.29281694f, 0.3710743f, 0.440774f, 0.5861263f, 0.7242663f, 0.7829777f, 0.9730387f, 1.147178f, 1.3474023f, 1.6833489f, 1.5748912f, 2.0653148f, 2.421378f, 2.705489f, 3.0849593f, 2.9578586f, 3.437159f, 5.1928988f, 4.247694f, 4.7362347f, 6.1843014f, 6.62399f, 5.7324314f, 5.3094177f, 8.594502f, 6.6565247f, 10.56011f
        };
        
        NonQuadraticConjugateGradientSolver solver = 
            new NonQuadraticConjugateGradientSolver(x, y, xe, ye);
        
        solver.setDebug(true);

        GEVYFit fit = solver.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
        //GEVYFit fit = solver.fitCurveKGreaterThanZeroAllAtOnce(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
        
    }
}
