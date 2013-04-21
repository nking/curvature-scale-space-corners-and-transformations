package algorithms.curves;

import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
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
        for (int i = 80; i < 81; i++) {

            String filePath = filePaths[i];

            HistogramHolder histogram = CreateClusterDataTest.readHistogram(filePath);

            chiSqMin = new GEVChiSquareMinimization(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors());


            float k = 0.5f;// changes sharpness of left side slope
            float s = 0.095f;
            float kMin = k/2.f;
            float kMax = k*2.f;
            float sMin = s/2;
            float sMax = s*2.0f;

            float mu = 0.1f;
            float yNorm = 1.0f;
            float yErrSqSum = chiSqMin.calcYErrSquareSum();

            chiSqMin.setDebug(true);

            // k=0.1842273 sigma=0.046056826 mu=0.07563025 chiSqSum=1077158.6 chiSqStatistic=19234.97
            // found a better fit by adjustments:
            //     k=0.5 sigma=0.095 mu=0.1                chiSqSum=7533744.0 chiSqStatistic=134531.14
            // k=0.37730515 sigma=0.07168771 mu=0.1        chiSqSum=1130325.2 chiSqStatistic= 20184.379

            //GEVYFit yfit = chiSqMin.fitCurve(kMin, kMax, sMin, sMax, mu, yErrSqSum, GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, k, s, mu, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

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
