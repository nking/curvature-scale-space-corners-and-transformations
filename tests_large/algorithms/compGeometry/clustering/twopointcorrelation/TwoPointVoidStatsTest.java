package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class TwoPointVoidStatsTest extends BaseTwoPointTest {

    protected boolean debug = true;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testCalc_Compare_Complete_toPartial() throws Exception {

        log.info("testCalc_Compare_Complete_toPartial()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        try {
            TwoPointVoidStatsPlotter plotter = new TwoPointVoidStatsPlotter();

            PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter();

            SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
            srr.setSeed(System.currentTimeMillis());
            long seed = srr.nextLong();

            SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(seed);

            int nIter = 5;

            for (int i = 0; i < nIter; i++) {

                // compare the distributions for 'SPARSELY' populated background
                DoubleAxisIndexer indexer = indexer = createIndexerWithRandomPoints(
                    sr, xmin, xmax, ymin, ymax, 3, 35, 60, 0.1f);

                TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
                stats.setDebug(debug);
                stats.setUseCompleteSampling(true);
                if (debug) {
                    log.info("run for complete sampling");
                }
                stats.calc();

                stats.plot(plotter, xmin, xmax, ymin, ymax);
                stats.plotFit(plotter2);

                float bSurfDensC = stats.getBackgroundSurfaceDensity();
                float bSurfDensCErr = stats.getBackgroundSurfaceDensityError();

                assertTrue(bSurfDensC > 0.0f);
                assertTrue(bSurfDensCErr > 0.0f);


                stats = new TwoPointVoidStats(indexer);
                stats.setDebug(debug);
                stats.setUseCompleteSampling(false);
                if (debug) {
                    log.info("run for default incomplete sampling");
                }
                stats.calc();

                stats.plot(plotter, xmin, xmax, ymin, ymax);
                stats.plotFit(plotter2);

                float bSurfDensI = stats.getBackgroundSurfaceDensity();
                float bSurfDensIErr = stats.getBackgroundSurfaceDensityError();

                assertTrue(bSurfDensI > 0.0f);
                assertTrue(bSurfDensIErr > 0.0f);

                float diff = bSurfDensC - bSurfDensI;

                boolean sameWithinFivePercent = (Math.abs(diff) < 0.05*bSurfDensC);
                boolean sameWithinOneSigmaError = (Math.abs(diff) < bSurfDensCErr);
                assertTrue(sameWithinOneSigmaError);
            }

        } catch (IOException e) {}
    }

    public void checkRuntimes() throws Exception {

        int[] nn = new int[]{10, 1000, 1000000};

        for (int i = 0; i < nn.length; i++) {
            int n = nn[i];
            // complete sampling:
            double permutations = computeDiv(n, 2)/2.f;
            log.info( "[0] N=" + n + " niterations=" + (permutations*permutations));
        }

        for (int i = 0; i < nn.length; i++) {
            int n = nn[i];
            // incomplete, range step sampling:
            double permutations = 0.375*n*n - n;
            log.info( "[1] N=" + n + " niterations=" + permutations);
        }
    }

    protected static long computeDiv(int n, int k) {

        long result = 1;
        for (int i = n; i > (n-k); i--) {
            result *= i;
        }
        return result;
    }
}
