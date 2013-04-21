package algorithms.curves;

import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

public class GEVChiSquareMinimizationTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    protected float[] x = null;
    protected float[] y = null;
    protected float[] dx = null;
    protected float[] dy = null;

    protected boolean debug = true;

    protected GEVChiSquareMinimization chiSqMin = null;

    protected boolean enable = true;

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testSortFromMinToMax() throws Exception {

        if (!enable) {
            return;
        }

        // placeholders, they don't resemble real data
        x = new float[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        y = new float[]{9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
        dx = new float[x.length];
        dy = new float[x.length];

        GEVYFit[] yfits = new GEVYFit[3];
        yfits[0] = new GEVYFit();
        yfits[0].setChiSqSum(10.0f);
        yfits[0].setX(x);
        yfits[0].setYFit(y);

        yfits[1] = new GEVYFit();
        yfits[1].setChiSqSum(5.0f);
        yfits[1].setX(x);
        yfits[1].setYFit(y);

        yfits[2] = new GEVYFit();
        yfits[2].setChiSqSum(1.0f);
        yfits[2].setX(x);
        yfits[2].setYFit(y);

        chiSqMin = new GEVChiSquareMinimization(x, y, dx, dx);
        chiSqMin.sortFromMinToMax(yfits, 0, 2);

        assertTrue(yfits[0].getChiSqSum() == 1.0f);
        assertTrue(yfits[1].getChiSqSum() == 5.0f);
        assertTrue(yfits[2].getChiSqSum() == 10.0f);
    }

    public void testGetSectionBoundariesFromGrid() throws Exception {

        if (!enable) {
            return;
        }

        float[] minMaxOut = new float[4];
        float kMin = 0;
        float kMax = 3.0f;
        float sigmaMin = 100.0f;
        float sigmaMax = 400.0f;

        float nDimension = 10;

        float kInterval = (kMax - kMin)/nDimension;
        float sInterval = (sigmaMax - sigmaMin)/nDimension;

        Set<String> uniqueSets = new HashSet<String>();

        // placeholders, they don't resemble real data
        x = new float[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        y = new float[]{9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
        dx = new float[x.length];
        dy = new float[x.length];

        GEVChiSquareMinimization chiSqMin = new GEVChiSquareMinimization(x, y, dx, dx);

        for (int position = 0; position < nDimension*nDimension; position++) {

            chiSqMin.getSectionBoundariesFromGrid(kMin, kMax, sigmaMin, sigmaMax, position, minMaxOut, nDimension);

            int row = (int)(position/nDimension);// rounds down

            int col = (position % (int)nDimension);

            float k0 = kMax - ((row + 1)*kInterval);
            float k1 = kMax - (row*kInterval);

            float s0 = (sigmaMin + (col)*sInterval);
            float s1 = (sigmaMin + (col + 1)*sInterval);

            assertTrue(minMaxOut[0] == k0);
            assertTrue(minMaxOut[1] == k1);

            assertTrue(minMaxOut[2] == s0);
            assertTrue(minMaxOut[3] == s1);

            StringBuffer sb = new StringBuffer();
            sb.append(Float.toString(k0)).append(":").append(Float.toString(k1)).append(" ")
                .append(Float.toString(s0)).append(":").append(Float.toString(s1));

            uniqueSets.add(sb.toString());

            //log.info(sb);
        }
        assertTrue(uniqueSets.size() == nDimension*nDimension);
    }

    public void testGetSectionParametersFromGrid() throws Exception {

        if (!enable) {
            return;
        }

        float[] parametersOut = new float[2];
        float[] minMaxOut = new float[4];

        float kMin = 0;
        float kMax = 3.0f;
        float sigmaMin = 100.0f;
        float sigmaMax = 400.0f;

        Set<String> uniqueSets = new HashSet<String>();

        // placeholders, they don't resemble real data
        x = new float[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        y = new float[]{9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
        dx = new float[x.length];
        dy = new float[x.length];

        GEVChiSquareMinimization chiSqMin = new GEVChiSquareMinimization(x, y, dx, dx);

        float nDimension = 10;

        float kInterval = (kMax - kMin)/nDimension;
        float sInterval = (sigmaMax - sigmaMin)/nDimension;

        for (int position = 0; position < nDimension*nDimension; position++) {

            chiSqMin.getSectionParametersFromGrid(kMin, kMax, sigmaMin, sigmaMax, position, parametersOut, nDimension);

            int row = (int)(position/nDimension);// rounds down

            int col = (position % (int)nDimension);

            float k = kMax - (row * kInterval) - (kInterval/2.0f);
            float s = sigmaMin + (col * sInterval) + (sInterval/2.0f);

            StringBuffer sb = new StringBuffer();
            sb.append(Float.toString(k)).append("=").append(Float.toString(parametersOut[0]))
                .append(" ").append(Float.toString(s)).append("=").append(Float.toString(parametersOut[1]));
            uniqueSets.add(sb.toString());

            //log.info(sb);

            assertTrue(parametersOut[0] == k);
            assertTrue(parametersOut[1] == s);

            // compare consistency of ranges

            chiSqMin.getSectionBoundariesFromGrid(kMin, kMax, sigmaMin, sigmaMax, position, minMaxOut, nDimension);
            float k0 = minMaxOut[0];
            float k1 = minMaxOut[1];

            float s0 = minMaxOut[2];
            float s1 = minMaxOut[3];

            sb = new StringBuffer();
            sb.append(Float.toString(k0)).append("=>").append(k).append("<=").append(Float.toString(k1)).append(" ")
                .append(Float.toString(s0)).append("=>").append(s).append("<=").append(Float.toString(s1));

            //log.info(sb);

            assertTrue(k0 < k);
            assertTrue(k < k1);
            assertTrue(s0 < s);
            assertTrue(s < s1);

        }
        assertTrue(uniqueSets.size() == nDimension*nDimension);
    }

    public void testCalculateChiSqSumForCurve_0() throws Exception {

        if (!enable) {
            return;
        }

        useTestData1();

        float k = 0.000198f;    // 1E-5 to 1E-3
        float sigma = 0.115f;   // 0.025 to 0.5
        float mu = 0.17647058f; // normalized to 0:1 scale

        GEVYFit yfit = chiSqMin.calculateChiSqSumAndCurve(k, sigma, mu,
            GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        if (debug) {

            float xmin = 0;
            float xmax = x[x.length - 1];
            float ymin = 0;
            float ymax = MiscMath.findMax(y);

            float[] xf = yfit.getOriginalScaleX();
            float[] yf = yfit.getOriginalScaleYFit();

            float ymax2 = MiscMath.findMax(yf);
            if (ymax2 > ymax) {
                ymax = ymax2;
            }

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
            plotter.addPlot(x, y, xf, yf, "");
            plotter.writeFile();

            log.info(yfit.toString());
        }
    }

    public void testFitCurve_WEIGHTED_BY_ERRORS_SIM_DATA_00() throws Exception {

        if (!enable) {
            return;
        }

        useTestData1();

        float k = 0.000198f;    // 1E-5 to 1E-3
        float sigma = 0.115f;   // 0.025 to 0.5
        float mu = 0.17647058f; // normalized to 0:1 scale

        GEVYFit yfit = chiSqMin.fitCurve(k/2.f, k*2.f, sigma/2.f, sigma*2.f, mu, 808.0f,
            GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        assertNotNull(yfit);

        if (debug) {

            float xmin = 0;
            float xmax = x[x.length - 1];
            float ymin = 0;
            float ymax = MiscMath.findMax(y);

            float[] xf = yfit.getOriginalScaleX();
            float[] yf = yfit.getOriginalScaleYFit();

            float ymax2 = MiscMath.findMax(yf);
            if (ymax2 > ymax) {
                ymax = ymax2;
            }

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
            plotter.addPlot(x, y, xf, yf, "");
            plotter.writeFile();

            log.info(yfit.toString());
        }

        float kdiff = Math.abs(yfit.getK() - k);
        float sdiff = Math.abs(yfit.getSigma() - sigma);
        assertTrue(kdiff < k/10);
        assertTrue(sdiff < sigma/10);
    }

    public void testFitCurve_WEIGHTED_BY_ERRORS_RANDOM_DATA_00() throws Exception {

        if (!enable) {
            return;
        }

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        //sr.setSeed( seed );
        sr.setSeed(9058936715185094572l);


        List<Long> seedsOfFailedTests = new ArrayList<Long>();

        log.info("using SEED=" + seed);
        int nIter = 4;

        for (int i = 0; i < nIter; i++) {

            float[] parameters = createTestDataUsingRandomParameters(sr);

            if (debug) {
                StringBuffer s0 = new StringBuffer("\n\nx = new float[]{");
                StringBuffer s1 = new StringBuffer("y = new float[]{");
                for (int ii = 0; ii < x.length; ii++) {
                    if (ii > 0) {
                        s0.append(", ");
                        s1.append(", ");
                    }
                    s0.append(x[ii]).append("f");
                    s1.append(y[ii]).append("f");
                }
                s0.append("};");
                s1.append("};");
                log.info(s0.toString());
                log.info(s1.toString());
            }

            float k = parameters[0];
            float sigma = parameters[1];
            float mu = parameters[2];
            float yConst = parameters[3];

            log.info("           k=" + k + "                 sigma=" + sigma + " mu=" + mu);

            GeneralizedExtremeValue gev = new GeneralizedExtremeValue(x, y, dx, dy);
            float[] yGEV = gev.generateNormalizedCurve(new float[]{k, sigma, mu});

            chiSqMin = new GEVChiSquareMinimization(x, yGEV, dx, dy);

            GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            if (yfit == null) {

                seedsOfFailedTests.add(seed);

                log.info("ERROR:  curve was not fit");

            } else if (debug) {

                float xmin = 0;
                float xmax = x[x.length - 1];
                float ymin = 0;
                float ymax = MiscMath.findMax(y);

                PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 1, 0, 1);
                plotter.addPlot(x, yGEV, dx, dy, yfit.getX(), yfit.getYFit(), "");
                plotter.writeFile();

                log.info("expecting: k=" + k + " sigma=" + sigma + " mu=" + mu);

                log.info(yfit.toString());

                // plot 2 fits.  the 2nd should be closer by number, but chisqsum is smaller for first.
                //               errors incorect?   <====
                /*
                 float[] yGEV1 = gev.generateNormalizedCurve(new float[]{50.0000076f, 422.1334229f, 0.067541346f});
                 float[] yGEV2 = gev.generateNormalizedCurve(new float[]{16.6666775f, 84.4266815f, 0.067541346f});
                 PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
                 plotter.addPlot(x, yGEV1, x, yfit.yfit, "0");
                 plotter.addPlot(x, yGEV2, x, yfit.yfit, "1");
                 plotter.writeFile();*/

                log.info("see plot");
                try {
                    Thread.sleep(5000);
                } catch (InterruptedException e) {
                }

                boolean fitIsBetterThanErrors = yfit.getChiSqSum() < yfit.getYDataErrSq();

                if (yfit.getChiSqStatistic() < 5) {
                    assertTrue(yfit.getChiSqStatistic() < 5);
                } else {
                    seedsOfFailedTests.add(seed);
                }
            }
        }

        if (!seedsOfFailedTests.isEmpty()) {
            log.info("Failed tests: ");
            for (Long s : seedsOfFailedTests) {
                log.info("seed: " + s);
            }
        }
    }

    protected void useTestData1() {

        x = new float[]{0.0005f, 0.0015f, 0.0025f, 0.0035f, 0.0045f, 0.0055f, 0.0065f, 0.0075f, 0.0085f};

        y = new float[]{0.46f, 1.0f, 0.517f, 0.265f, 0.148f, 0.083f, 0.078f, 0.048f, 0.0043f};

        for (int i = 0; i < y.length; i++) {
            y[i] *= 230;
        }

        dy = Errors.populateYErrorsBySqrt(y);

        dx = Errors.populateXErrorsByPointSeparation(x);

        chiSqMin = new GEVChiSquareMinimization(x, y, dx, dy);
    }

    protected float[] createTestDataUsingRandomParameters(SecureRandom sr) throws NoSuchAlgorithmException {
        // generate histogram
        float xmin = (float) (sr.nextFloat() * Math.pow(10, -1*sr.nextInt(2)));
        float xmax = xmin * (float) (Math.pow(10, sr.nextInt(2)));
        while (xmax == xmin) {
            xmax = xmin * (float) (Math.pow(10, sr.nextInt(2)));
        }
        return createTestDataUsingRandomParameters(sr, xmin, xmax);
    }

    protected float[] createTestDataUsingRandomParameters(SecureRandom sr, float xmin, float xmax) throws NoSuchAlgorithmException {

        int nPoints = sr.nextInt(30);
        while (nPoints < 8) {
            nPoints = sr.nextInt(30);
        }
        x = new float[nPoints];
        y = new float[nPoints];
        float deltaX = (xmax - xmin)/(float)nPoints;
        for (int i = 0; i < nPoints; i++) {
            x[i] = xmin + i*deltaX;
        }

        xmin = x[0];
        xmax = x[x.length - 1];
        float mu = x[sr.nextInt(3)];

        /*
        float kMin = 0.00001f;
        float kMax = 0.001f;
        float mu = x[1];
        float sigmaMin = 0.025f;
        float sMax = 20.0f*sigmaMin;
         */
        float[] parameters = GEVChiSquareMinimization.generateRandomParameters(sr, mu);
        float k = parameters[0];
        float sigma = parameters[1];

        log.info("nPoints=" + nPoints + " xmin=" + xmin + " xmax=" + xmax + " mu=" + mu);

        dx = Errors.populateXErrorsByPointSeparation(x);
        dy = new float[nPoints];

        float normFactor = sr.nextInt(100);

        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(x, y, dx, dy);
        float[] yGEV = gev.generateNormalizedCurve(new float[]{k, sigma, mu}, normFactor);

        if (yGEV == null) {
            return createTestDataUsingRandomParameters(sr);
        }

        for (int i = 0; i < nPoints; i++) {
            y[i] = yGEV[i] * normFactor;
        }

        dy = Errors.populateYErrorsBySqrt(y);

        // simulate 1/10th the shot noise
        /*for (int i = 0; i < dy.length; i++) {
            dy[i] = 0.1f * dy[i];
        }*/

        chiSqMin = new GEVChiSquareMinimization(x, y, dx, dy);

        /*
        try {
            xmin = x[0];
            xmax = x[x.length - 1];
            float ymin = MiscMath.findMin(y);
            float ymax = MiscMath.findMax(y);

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
            plotter.addPlot(x, yGEV, x, y, "");
            plotter.writeFile();

            log.info("nPoints=" + nPoints + " k=" + k + " sigma=" + sigma + " mu=" + mu);
        } catch (Exception e) {
        }*/

        return new float[]{k, sigma, mu, normFactor};
    }

    public void testUsingRandomClusters() throws Exception {

        log.info("testUsingRandomClusters()");

        if (enable) {
            return;
        }

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        String[] filePaths = CreateClusterDataTest.getHistogramFilePaths();

        for (int i = 0; i < filePaths.length; i++) {
        //for (int i = 17; i < 18; i++) {

            //String filePath = filePaths[i];
            String filePath =
                "/Volumes/portable/data/projects/github/algorithms/algorithms/algorithms_in_java/submodules/two-point-correlation/bin/test-classes/../../tmpdata/"
                + "histogram_random_background_with_0_clusters_017.dat";

            HistogramHolder histogram = CreateClusterDataTest.readHistogram(filePath);

            chiSqMin = new GEVChiSquareMinimization(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors());

            //float kMin, float kMax, float sigmaMin, float sigmaMax, float mu, float yErrSquareSum,
            //     WEIGHTS_DURING_CHISQSUM weightMethod, float yNorm
            float kMin = 0.00011246394f/1.1f;
            float kMax = 0.00011246394f*1.1f;
            float sMin = 0.07416199f/2;
            float sMax = 2.0f*0.07416199f;

            float mu = 0.10344828f;
            float yNorm = 1.0f;
            float yErrSqSum = chiSqMin.calcYErrSquareSum();

            float k = 1.4f;// changes sharpness of left side slope
            float s = 0.35f;
            mu = 0.26f;

            //k = 0.0000811f;
            //mu = 0.1176f;
            //s = 0.059f;

            // k=8.110986E-5 sigma=0.059118375 mu=0.11764706 chiSqSum=11968.323 chiSqStatistic=544.0147
            // k=1.4         sigma=0.35        mu=0.26       chiSqSum=5268.0234 chiSqStatistic=239.45561

            //GEVYFit yfit = chiSqMin.fitCurve(kMin, kMax, sMin, sMax, mu, yErrSqSum, GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, k, s, mu, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            chiSqMin.setDebug(true);
            GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            float[] xf = yfit.getOriginalScaleX();
            float[] yf = yfit.getOriginalScaleYFit();

            plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors(), xf, yf, String.valueOf(i));
            plotter.writeFile();

            System.out.println(yfit.toString());

            log.info("see plot for " + i);
/*
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) {
            } finally {
                log.info("next");
            }*/
        }
    }
}
