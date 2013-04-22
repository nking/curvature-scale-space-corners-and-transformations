package algorithms.curves;

import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.misc.HistogramHolder;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;

public class GEVChiSquareMinimization3Test extends TestCase {

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

        String fileNamePostfix = "wikipedia_dbscan.dat";
        String fileName = CreateClusterDataTest.histogramFileNamePrefix + fileNamePostfix;
        String filePath = ResourceFinder.getAFilePathInTmpData(fileName);


            HistogramHolder histogram = CreateClusterDataTest.readHistogram(filePath);

            chiSqMin = new GEVChiSquareMinimization(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors());


            float k = 1.04f;// changes sharpness of left side slope
            float s = 0.195f;
            float kMin = k/2.f;
            float kMax = k*2.f;
            float sMin = s/2;
            float sMax = s*2.0f;

            float mu = 0.06f;
            float yNorm = 1.0f;
            float yErrSqSum = chiSqMin.calcYErrSquareSum();

            chiSqMin.setDebug(true);

            // k=0.67109966 sigma=0.16777492 mu=0.016949153 chiSqSum=96.46002 chiSqStatistic=3.7100008
            //  trial and error:
            //   k=1.0402603 sigma=0.19504745 mu=0.06 chiSqSum=86.89037 chiSqStatistic=3.3419375

            GEVYFit yfit = chiSqMin.fitCurve(kMin, kMax, sMin, sMax, mu, yErrSqSum, GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, k, s, mu, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            float[] xf = yfit.getOriginalScaleX();
            float[] yf = yfit.getOriginalScaleYFit();

            plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors(), xf, yf, "");
            plotter.writeFile();

            System.out.println(yfit.toString());

/*
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) {
            } finally {
                log.info("next");
            }*/
    }
}
