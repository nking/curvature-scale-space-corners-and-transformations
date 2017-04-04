package algorithms.curves;

import java.io.IOException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NonQuadraticConjugateGradientSolverTest extends TestCase {

    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     *
     * @throws Exception
     */
    public void testFitCurve() throws Exception {

        // while revising code, if not on a development branch, don't assert results
        boolean assertResults = false;

        float k = 1.80f;
        float sigma = 0.85f;
        float mu = 0.441f;

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
        fitCurve.setMaximumNumberOfIterations(100);

        // TODO:  problem fitting parameters in <>  k=1.4680 <0.09> s=0.1307 <0.05> m=0.1543 <0.106> (nx=40,i=21) chi=0.013889

        float kMin = 0.01f;
        float kMax = 3*k;
        float sigmaMin = 0.01f*sigma;
        float sigmaMax = 2.5f*sigma;
        float muMin = 0.5f * mu;
        float muMax = xp[xp.length - 1];

        //GEVYFit yFit = fitCurve.fitCurve(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
        //GEVYFit yFit = fitCurve.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
        GEVYFit yFit = fitCurve.fitCurveParametersAllAtOnce(kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
        //GEVYFit yFit = fitCurve.fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu);

        if (assertResults) {
            assertNotNull(yFit);
            assertTrue(Math.abs(yFit.k - k) < 0.1*k);
            assertTrue(Math.abs(yFit.sigma - sigma) < 0.1*sigma);
            assertTrue(Math.abs(yFit.mu - mu) < 0.1*mu);
        }
    }

    /**
     *
     * @throws Exception
     */
    public void estFitRandomCurves() throws Exception {

        log.info("testFitRandomCurves");

        boolean createPlots = true;
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);

        int n = 50;
        int nX = 40;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1386620575944l;
        sr.setSeed( seed );

        log.info("SEED=" + seed);

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

                NonQuadraticConjugateGradientSolver solver =
                    new NonQuadraticConjugateGradientSolver(xp, yGEV, xe, ye);

                GEVYFit fit = solver.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
                assertNotNull(fit);

                if (createPlots) {
                    
                    String label = String.format(
                        "k=%4.4f <%2.2f>  s=%4.4f <%2.2f>  m=%4.4f <%2.3f>  (nx=%d,i=%d)  chi=%4.6f",
                        fit.getK(), k,
                        fit.getSigma(), s,
                        fit.getMu(), m,
                        nX, i, fit.getChiSqSum());

                    plotter.addPlot(xp, yGEV, xe, ye, fit.getX(), fit.getYFit(), label);

                    path = plotter.writeFile2();
                }

                System.out.println(fit.getChiSqStatistic() + "(" + nX + ", " + i + ")");

                assertTrue(fit.getChiSqStatistic() < 1e-4);

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
                        log.info("float[] x = new float[]{"  + xsb.append("f").toString() + "};");
                        log.info("float[] y = new float[]{"  + ysb.append("f").toString() + "};");
                        log.info("float[] xe = new float[]{" + xesb.append("f").toString() + "};");
                        log.info("float[] ye = new float[]{" + yesb.append("f").toString() + "};");
                    }
                }
            }

            nX >>= 1;
        }
        
        if (createPlots) {
            log.fine(" plot is at path=" + path);
        }
    }

    /**
     *
     * @throws Exception
     */
    public void testAFit() throws Exception {
        /*
        float[] x = new float[]{0.15102175f, 0.22170976f, 0.2923978f, 0.3630858f, 0.43377382f, 0.5044618f, 0.5751499f, 0.6458379f, 0.7165259f, 0.7872139f, 0.85790193f, 0.92859f, 0.999278f, 1.0699661f, 1.1406541f, 1.2113421f, 1.2820301f, 1.3527181f, 1.4234061f, 1.4940941f, 1.5647821f, 1.6354703f, 1.7061583f, 1.7768463f, 1.8475343f, 1.9182223f, 1.9889103f, 2.0595984f, 2.1302867f, 2.2009747f, 2.2716627f, 2.3423507f, 2.4130387f, 2.4837267f, 2.5544147f, 2.6251028f, 2.6957908f, 2.7664788f, 2.8371668f, 2.9078548f};

        float[] y = new float[]{19368f, 26163f, 23344f, 19869f, 17298f, 15229f, 13484f, 12431f, 11546f, 10834f, 10037f, 8986f, 7671f, 6119f, 5195f, 4326f, 3466f, 2812f, 2246f, 1783f, 1463f, 1210f, 997f, 848f, 702f, 577f, 500f, 409f, 294f, 318f, 294f, 211f, 193f, 164f, 145f, 121f, 124f, 110f, 82f, 74f};

        float[] xe = new float[]{169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f, 169.9074f};

        float[] ye = new float[]{0.040564783f, 0.0448817f, 0.054843098f, 0.0670746f, 0.080281235f, 0.09502251f, 0.11150011f, 0.12697266f, 0.14266293f, 0.15957877f, 0.17806467f, 0.20129223f, 0.23124927f, 0.27282917f, 0.31284237f, 0.3598531f, 0.42313862f, 0.49322188f, 0.57433116f, 0.66769683f, 0.77801037f, 0.88995105f, 1.0249188f, 1.1400098f, 1.2954072f, 1.4906247f, 1.7007264f, 1.9036003f, 2.3256164f, 2.308904f, 2.4918737f, 2.9840057f, 3.2392955f, 3.6577778f, 3.919288f, 4.590691f, 4.4787126f, 4.9169383f, 5.829587f, 6.103832f};
        */

        /*
        float[] x = new float[]{0.16171317f, 0.25354025f, 0.34536734f, 0.4371944f, 0.5290215f, 0.62084854f, 0.7126756f, 0.80450267f, 0.89632976f, 0.98815686f, 1.0799841f, 1.1718111f, 1.2636381f, 1.3554653f, 1.4472923f, 1.5391195f, 1.6309465f, 1.7227736f, 1.8146007f, 1.9064277f, 1.9982549f, 2.0900817f, 2.1819088f, 2.273736f, 2.365563f, 2.45739f, 2.5492172f, 2.6410444f, 2.7328713f, 2.8246984f, 2.9165256f, 3.0083525f, 3.1001797f, 3.1920068f, 3.2838337f, 3.375661f, 3.467488f, 3.559315f, 3.6511421f, 3.7429693f};

        float[] y = new float[]{20917f, 24794f, 20404f, 16888f, 14680f, 13050f, 12207f, 10824f, 9090f, 7316f, 5585f, 4407f, 3445f, 2612f, 2136f, 1676f, 1366f, 1088f, 889f, 651f, 596f, 507f, 436f, 330f, 314f, 283f, 218f, 167f, 172f, 131f, 133f, 103f, 75f, 85f, 85f, 49f, 43f, 46f, 41f, 22f};

        float[] ye = new float[]{161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f, 161.12688f};

        float[] xe = new float[]{0.061465673f, 0.06955376f, 0.08488417f, 0.1032283f, 0.12232937f, 0.14293118f, 0.16277218f, 0.18691832f, 0.21879023f, 0.2594103f, 0.31464696f, 0.37558132f, 0.4501729f, 0.5432629f, 0.6320125f, 0.74401575f, 0.8751161f, 1.0240066f, 1.1970534f, 1.4569196f, 1.5881971f, 1.8039417f, 1.98375f, 2.3620512f, 2.5600765f, 2.7914815f, 3.283865f, 3.9357467f, 3.8539512f, 4.678683f, 4.851508f, 5.698987f, 6.5019464f, 6.769337f, 6.5251536f, 9.440966f, 10.333385f, 9.880095f, 10.917977f, 14.638054f};
        */
        
        float[] x = new float[]{0.07006695f, 0.08293f, 0.09579305f, 0.1086561f, 0.121519156f, 0.13438219f, 0.14724524f, 0.1601083f, 0.17297134f, 0.1858344f, 0.19869745f, 0.2115605f, 0.22442354f, 0.2372866f, 0.25014967f, 0.2630127f, 0.27587575f, 0.2887388f, 0.30160183f, 0.3144649f, 0.32732797f, 0.340191f, 0.35305405f, 0.3659171f, 0.37878013f, 0.39164323f, 0.40450627f, 0.4173693f, 0.43023235f, 0.4430954f, 0.45595843f, 0.46882153f, 0.48168457f, 0.4945476f, 0.50741065f, 0.5202737f, 0.5331367f, 0.5459998f, 0.55886286f, 0.5717259f};
        float[] y = new float[]{2249f, 4722f, 5745f, 6297f, 6070f, 5660f, 5349f, 4968f, 4765f, 4547f, 4358f, 4091f, 3894f, 3695f, 3639f, 3409f, 3417f, 3302f, 3151f, 3052f, 3040f, 2953f, 2893f, 2764f, 2756f, 2759f, 2656f, 2593f, 2497f, 2281f, 2236f, 2103f, 2011f, 1942f, 1797f, 1646f, 1573f, 1450f, 1342f, 1330f};
        float[] xe = new float[]{0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f, 0.0064315237f};
        float[] ye = new float[]{0.06368805f, 0.05134599f, 0.053022746f, 0.056779608f, 0.06454764f, 0.07324291f, 0.08237399f, 0.09312447f, 0.102609545f, 0.112454824f, 0.122811265f, 0.1344264f, 0.14482236f, 0.15727329f, 0.16736402f, 0.18021643f, 0.18822077f, 0.20053904f, 0.21431826f, 0.22658357f, 0.23538195f, 0.2489218f, 0.25953734f, 0.27523327f, 0.28748488f, 0.2942547f, 0.30991468f, 0.32356742f, 0.34100977f, 0.36846435f, 0.3814147f, 0.40344337f, 0.42271626f, 0.44184867f, 0.46721852f, 0.5086379f, 0.5280288f, 0.56130296f, 0.59482855f, 0.6136231f};
        
        for (int i = 0; i < 2; i++) {
            
            GEVYFit fit;
            
            if (i == 0) {
                
                NonQuadraticConjugateGradientSolver solver = new NonQuadraticConjugateGradientSolver(
                    x, y, xe, ye);
                
                fit = solver.fitCurveKGreaterThanZeroAllAtOnce(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
                
                if (fit != null) {

                    String label = String.format("chisq=%.8f k=%.4e s=%.4e m=%.4e",
                        fit.getChiSqStatistic(), fit.getK(), fit.getSigma(), fit.getMu());

                    plotFit2(fit, label, solver);
                }
                
            } else {

                GEVChiSquareMinimization solver = new GEVChiSquareMinimization(
                    x, y, xe, ye);
            
                fit = solver.fitCurveKGreaterThanZero(
                    GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
                
                if (fit != null) {

                    String label = String.format("chisq=%.8f k=%.4e s=%.4e m=%.4e",
                        fit.getChiSqStatistic(), fit.getK(), fit.getSigma(), fit.getMu());

                    plotFit(fit, label, solver);
                }
            }
        }
    }
    
    /**
     *
     * @param yfit
     * @param label
     * @param solver
     * @throws IOException
     */
    protected void plotFit(GEVYFit yfit, String label, GEVChiSquareMinimization solver) throws IOException {

        float xIntervalHalf = (yfit.getX()[1] - yfit.getX()[0]) / 2.0f;
        float xmin = yfit.getX()[0] - xIntervalHalf;
        float xmax = yfit.getX()[ yfit.getX().length - 1];
        float ymin = 0.0f;
        float ymax = MiscMath.findMax(solver.y);

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);

        try {
            plotter.addPlot(solver.x, solver.y, solver.xe, solver.ye, yfit.getX(), yfit.getYFit(), label);
            String filePath = plotter.writeFile3();
            System.out.println("filePath=" + filePath);
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }
    }
    
    /**
     *
     * @param yfit
     * @param label
     * @param solver
     * @throws IOException
     */
    protected void plotFit2(GEVYFit yfit, String label, NonQuadraticConjugateGradientSolver solver) throws IOException {

        float xIntervalHalf = (yfit.getX()[1] - yfit.getX()[0]) / 2.0f;
        float xmin = yfit.getX()[0] - xIntervalHalf;
        float xmax = yfit.getX()[ yfit.getX().length - 1];
        float ymin = 0.0f;
        float ymax = MiscMath.findMax(solver.y);

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);

        try {
            plotter.addPlot(solver.x, solver.y, solver.xe, solver.ye, yfit.getX(), yfit.getYFit(), label);
            String filePath = plotter.writeFile3();
            System.out.println("filePath=" + filePath);
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }
    }
}