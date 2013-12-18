package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import java.io.File;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

/**
 *
 * @author nichole
 */
public class TwoPointVoidStatsTest extends BaseTwoPointTest {

    boolean debug = true;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testCalc0() throws Exception {

        log.info("testCalc0()");

        /* Create a simple rectangle of values with 4 minima areas to find
         *
         * 100 100 100 100 100 100 100 100 100  10  <== [90:98]
         * 100 100  10  10 100 100 100 100 100  9   <== [81:89]
         * 100 100  10  10 100 100 100 100 100  8   <== [72:80]
         * 100 100 100 100 100  10  10  10 100  7   <== [63:71]
         * 100 100 100 100 100  10  10  10 100  6   <== [54:62]
         * 100 100 100 100 100  10  10  10 100  5   <== [45:53]
         * 100 100 100 100 100 100 100 100 100  4   <== [36:44]
         * 100  10  10  10 100 100 100 100 100  3   <== [27:35]
         * 100  10  10  10 100 100 100 100 100  2   <== [18:26]
         * 100  10  10  10 100 100 100 100 100  1   <== [ 9:17]
         * 100 100 100 100 100 100 100 100 100  0   <== [ 0: 8]
         *  0   1   2   3   4   5   6   7   8
         *
         * where the 100's are where points are and the 10's are where points are not (==voids)
         */

        generator.x = new float[] {
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,              8,
            0,  1,  2,  3,  4,              8,
            0,  1,  2,  3,  4,              8,
            0,  1,          4,  5,  6,  7,  8,
            0,  1,          4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8
        };
        generator.y = new float[] {
            0,  0,  0,  0,  0,  0,  0,  0,  0,
            1,              1,  1,  1,  1,  1,
            2,              2,  2,  2,  2,  2,
            3,              3,  3,  3,  3,  3,
            4,  4,  4,  4,  4,  4,  4,  4,  4,
            5,  5,  5,  5,  5,              5,
            6,  6,  6,  6,  6,              6,
            7,  7,  7,  7,  7,              7,
            8,  8,          8,  8,  8,  8,  8,
            9,  9,          9,  9,  9,  9,  9,
           10, 10, 10, 10, 10, 10, 10, 10, 10
        };

        // make uniform errors for x and y
        generator.xErrors = new float[generator.x.length];
        generator.yErrors = new float[generator.x.length];
        Arrays.fill(generator.xErrors, 0.1f);
        Arrays.fill(generator.yErrors, 0.1f);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(debug);
        stats.setGEVRangeParameters(0.0001f, 1.0f, 0.00001f, 0.002f);
        //stats.setStandardDeviationFactor(2.5f);        
        stats.calc();

        float bSurfDens    = stats.getBackgroundDensity();
        float bSurfDensErr = stats.getBackgroundDensityError();

        String filePath = stats.persistTwoPointBackground();
        assertNotNull(filePath);

        File fl = new File(filePath);
        assertTrue(fl.exists());

        boolean didRead = stats.readTwoPointBackground(filePath);
        assertTrue(didRead);
    }

    public void testCalulateAndHistogramBinErrors() throws Exception {

        log.info("testCalulateAndHistogramBinErrors()");

        // note, this dataset needs a GEV with k < 0 which is not allowed by the fitting code

        // make evenly distributed data with slight offset
        /*
         *
         *
         *   1|  1   3
         *    |
         *   0|  0   2
         *     --------
         *       0   1
         */
        int nX = 2;
        generator.x = new float[nX*nX];
        generator.y = new float[nX*nX];
        int count = 0;
        for (int i = 0 ; i < nX; i++) {
            for (int j = 0; j < nX; j++) {
                generator.x[count] = i*100;
                generator.y[count] = j*100;
                count++;
            }
        }

        // make uniform errors for x and y
        generator.xErrors = Errors.populateXErrorsByPointSeparation(generator.x);
        generator.yErrors = Errors.populateYErrorsBySqrt(generator.y);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(debug);
        //stats.setStandardDeviationFactor(2.5f);
        stats.calc();

        assertNotNull(stats.getBestFit());
        assertNotNull(stats.getStatsHistogram());

        HistogramHolder hist = stats.getStatsHistogram();
        float[] xh = hist.getXHist();
        float[] yh = hist.getYHistFloat();
        float[] xhe = hist.getXErrors();
        float[] yhe = hist.getYErrors();

    }

    public void testProcessIndexedRegion() throws Exception {
    }
    public void testCalculateStatsForBackground() throws Exception {
    }

    public void testCalulateHistogramBinErrors_2() throws Exception {

        // note, this dataset needs a GEV with k < 0 which is not allowed by the fitting code

        // make evenly distributed data with slight offset
        /*
         *
         *
         *   1|  1   3
         *    |
         *   0|  0   2
         *     --------
         *       0   1
         */
        int nX = 2;
        generator.x = new float[nX*nX];
        generator.y = new float[nX*nX];
        int count = 0;
        for (int i = 0 ; i < nX; i++) {
            for (int j = 0; j < nX; j++) {
                generator.x[count] = i ;//+ count*0.01f;
                generator.y[count] = j ;//+ count*0.01f;
                count++;
            }
        }

        // make uniform errors for x and y
        generator.xErrors = new float[generator.x.length];
        generator.yErrors = new float[generator.y.length];
        Arrays.fill(generator.xErrors, 0.1f);
        Arrays.fill(generator.yErrors, 0.1f);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(debug);
        stats.setUseCompleteSampling();
        //stats.setStandardDeviationFactor(2.5f);

        stats.calc();

        assertNotNull(stats.getBestFit());

    }

    public void testCalc3() throws Exception {

        log.info("testCalc3");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        //sr.setSeed( seed );
        sr.setSeed(-2384802679227907254l);

        log.info("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        createPoints((30+40+60), new int[]{30, 40, 60},
            RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION.LARGE,
            xmin, xmax, ymin, ymax, sr, false);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(debug);
        //stats.setStandardDeviationFactor(2.5f);
        
        stats.calc();

        assertNotNull(stats.getBestFit());
    }

    public void testCalc4() throws Exception {

        log.info("testCalc4");

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        //sr.setSeed( seed );
        sr.setSeed(-2384802679227907254l);

        log.info("using SEED=" + seed);

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        createPoints((3), new int[]{3},
            RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION.LARGE,
            xmin, xmax, ymin, ymax, sr, false);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(generator.x, generator.y,
            generator.xErrors, generator.yErrors, generator.x.length);

        TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
        stats.setDebug(debug);
        //stats.setStandardDeviationFactor(2.5f);
        
        stats.calc();

        assertTrue(stats.voidFinder.getTwoPointDensities().length == 10);
        
        assertNotNull(stats.getBestFit());
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
