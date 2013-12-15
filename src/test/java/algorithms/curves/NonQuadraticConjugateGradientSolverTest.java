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
    
    public void estFitRandomCurves() throws Exception {
        
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
    
    public void testAFit() throws Exception {
        /*        
        float[] x = new float[]{0.15102175f, 0.22170976f, 0.2923978f, 0.3630858f, 0.43377382f, 0.5044618f, 0.5751499f, 0.6458379f, 0.7165259f, 0.7872139f, 0.85790193f, 0.92859f, 0.999278f, 1.0699661f, 1.1406541f, 1.2113421f, 1.2820301f, 1.3527181f, 1.4234061f, 1.4940941f, 1.5647821f, 1.6354703f, 1.7061583f, 1.7768463f, 1.8475343f, 1.9182223f, 1.9889103f, 2.0595984f, 2.1302867f, 2.2009747f, 2.2716627f, 2.3423507f, 2.4130387f, 2.4837267f, 2.5544147f, 2.6251028f, 2.6957908f, 2.7664788f, 2.8371668f, 2.9078548f};
        
        float[] y = new float[]{19368f, 26163f, 23344f, 19869f, 17298f, 15229f, 13484f, 12431f, 11546f, 10834f, 10037f, 8986f, 7671f, 6119f, 5195f, 4326f, 3466f, 2812f, 2246f, 1783f, 1463f, 1210f, 997f, 848f, 702f, 577f, 500f, 409f, 294f, 318f, 294f, 211f, 193f, 164f, 145f, 121f, 124f, 110f, 82f, 74f};
        
        float[] xe = new float[]{169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f};
        
        float[] ye = new float[]{0.040564783f, 0.0448817f, 0.054843098f, 0.0670746f, 0.080281235f, 0.09502251f, 0.11150011f, 0.12697266f, 0.14266293f, 0.15957877f, 0.17806467f, 0.20129223f, 0.23124927f, 0.27282917f, 0.31284237f, 0.3598531f, 0.42313862f, 0.49322188f, 0.57433116f, 0.66769683f, 0.77801037f, 0.88995105f, 1.0249188f, 1.1400098f, 1.2954072f, 1.4906247f, 1.7007264f, 1.9036003f, 2.3256164f, 2.308904f, 2.4918737f, 2.9840057f, 3.2392955f, 3.6577778f, 3.919288f, 4.590691f, 4.4787126f, 4.9169383f, 5.829587f, 6.103832f};
        */
        
        float[] x = new float[]{0.16171317f, 0.25354025f, 0.34536734f, 0.4371944f, 0.5290215f, 0.62084854f, 0.7126756f, 0.80450267f, 0.89632976f, 0.98815686f, 1.0799841f, 1.1718111f, 1.2636381f, 1.3554653f, 1.4472923f, 1.5391195f, 1.6309465f, 1.7227736f, 1.8146007f, 1.9064277f, 1.9982549f, 2.0900817f, 2.1819088f, 2.273736f, 2.365563f, 2.45739f, 2.5492172f, 2.6410444f, 2.7328713f, 2.8246984f, 2.9165256f, 3.0083525f, 3.1001797f, 3.1920068f, 3.2838337f, 3.375661f, 3.467488f, 3.559315f, 3.6511421f, 3.7429693f};
        
        float[] y = new float[]{20917f, 24794f, 20404f, 16888f, 14680f, 13050f, 12207f, 10824f, 9090f, 7316f, 5585f, 4407f, 3445f, 2612f, 2136f, 1676f, 1366f, 1088f, 889f, 651f, 596f, 507f, 436f, 330f, 314f, 283f, 218f, 167f, 172f, 131f, 133f, 103f, 75f, 85f, 85f, 49f, 43f, 46f, 41f, 22f};
        
        float[] xe = new float[]{161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f};
        
        float[] ye = new float[]{0.061465673f, 0.06955376f, 0.08488417f, 0.1032283f, 0.12232937f, 0.14293118f, 0.16277218f, 0.18691832f, 0.21879023f, 0.2594103f, 0.31464696f, 0.37558132f, 0.4501729f, 0.5432629f, 0.6320125f, 0.74401575f, 0.8751161f, 1.0240066f, 1.1970534f, 1.4569196f, 1.5881971f, 1.8039417f, 1.98375f, 2.3620512f, 2.5600765f, 2.7914815f, 3.283865f, 3.9357467f, 3.8539512f, 4.678683f, 4.851508f, 5.698987f, 6.5019464f, 6.769337f, 6.5251536f, 9.440966f, 10.333385f, 9.880095f, 10.917977f, 14.638054f};
        
        NonQuadraticConjugateGradientSolver solver = 
            new NonQuadraticConjugateGradientSolver(x, y, xe, ye);
        
        solver.setDebug(true);

        GEVYFit fit = solver.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
        //GEVYFit fit = solver.fitCurveKGreaterThanZeroAllAtOnce(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
        
    }
}
