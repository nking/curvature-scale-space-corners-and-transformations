package algorithms.curves;

import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.curves.GEVChiSquareMinimization;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.PolygonAndPointPlotter;
import java.util.logging.Logger;
import junit.framework.TestCase;

public class GEVChiSquareMinimization2Test extends TestCase {

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

    public void testUsingRandomClusters() throws Exception {

        log.info("testUsingRandomClusters()");

        if (!enable) {
            return;
        }

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        String[] filePaths = CreateClusterDataTest.getHistogramFilePaths();

        //for (int i = 0; i < filePaths.length; i++) {
        for (int i = 61; i < 62; i++) {//120

            String filePath = filePaths[i];//ndexer_random_background_with_3_clusters_130.dat"

            HistogramHolder histogram = CreateClusterDataTest.readHistogram(filePath);

            chiSqMin = new GEVChiSquareMinimization(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors());


            float k = 1.1f;// changes sharpness of left side slope
            float s = 0.27f;
            float kMin = k/2.f;
            float kMax = k*2.f;
            float sMin = s/2;
            float sMax = s*2.0f;

            float mu = 0.26f;
            float yNorm = 1.0f;
            float yErrSqSum = chiSqMin.calcYErrSquareSum();

            chiSqMin.setDebug(true);

            // k=1.798601 sigma=0.45760158 mu=0.010309278 chiSqSum=286.08374 chiSqStatistic=6.3574166
            //  trial and error:
            //    k=1.1 sigma=0.23 mu=0.25 chiSqSum=372.00134 chiSqStatistic=8.266697

            //GEVYFit yfit = chiSqMin.fitCurve(kMin, kMax, sMin, sMax, mu, yErrSqSum, GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, k, s, mu, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            float[] xf = yfit.getOriginalScaleX();
            float[] yf = yfit.getOriginalScaleYFit();

            plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors(), xf, yf, String.valueOf(i));
            plotter.writeFile();

            log.fine(yfit.toString());

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
