package algorithms.curves;

import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GEVChiSquareMinimization3Test extends TestCase {

    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     *
     */
    protected float[] x = null;

    /**
     *
     */
    protected float[] y = null;

    /**
     *
     */
    protected float[] dx = null;

    /**
     *
     */
    protected float[] dy = null;

    /**
     *
     */
    protected boolean debug = true;

    /**
     *
     */
    protected GEVChiSquareMinimization chiSqMin = null;

    /**
     *
     */
    protected boolean enable = true;

    /**
     *
     * @throws Exception
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     *
     * @throws Exception
     */
    public void testAFit() throws Exception {

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

        GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(
            GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        float[] xf = yfit.getOriginalScaleX();
        float[] yf = yfit.getOriginalScaleYFit();

        plotter.addPlot(xx, yy, xxe, yye, xf, yf, "");
        plotter.writeFile();

        log.info(yfit.toString());

    }
    
}
