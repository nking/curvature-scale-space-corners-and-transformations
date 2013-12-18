package algorithms.curves;

import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.util.logging.Logger;
import junit.framework.TestCase;

public class GEVChiSquareMinimization3Test extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());

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

        if (!(new File(filePath)).exists()) {
            return;
        }

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

            log.info(yfit.toString());

/*
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) {
            } finally {
                log.info("next");
            }*/
    }
    
    public void estAFit() throws Exception {

        log.info("testAFit()");

        float[] xx = new float[]{1.1976953f, 1.4704535f, 1.7432117f, 2.01597f, 2.2887282f, 2.5614867f, 
            2.8342447f, 3.1070032f, 3.3797612f, 3.6525197f, 3.9252777f, 4.198036f, 4.470794f, 4.7435527f, 
            5.0163107f, 5.2890687f, 5.561827f, 5.8345857f, 6.1073437f, 6.3801017f, 6.65286f, 6.9256186f, 
            7.1983767f, 7.4711347f, 7.743893f};
        float[] yy = new float[]{19f, 34f, 70f, 86f, 85f, 79f, 52f, 59f, 42f, 29f, 25f, 
            17f, 17f, 14f, 7f, 6f, 5f, 3f, 2f, 4f, 4f, 1f, 2f, 1f, 0f};
        float[] xxe = Errors.populateXErrorsByPointSeparation(xx);
        float[] yye = Errors.populateYErrorsBySqrt(yy);
        
        chiSqMin = new GEVChiSquareMinimization(xx, yy, xxe, yye);

        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        float[] xf = yfit.getOriginalScaleX();
        float[] yf = yfit.getOriginalScaleYFit();

        plotter.addPlot(xx, yy, xxe, yye, xf, yf, "");
        plotter.writeFile();

        log.info(yfit.toString());

    }
    
}
