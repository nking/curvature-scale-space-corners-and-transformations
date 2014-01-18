package algorithms.curves;

import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GEVYFitTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of calculateArea method, of class GEVYFit.
     */
    public void testCalculateArea() {

        log.info("calculateArea");

        GEVYFit yfit = new GEVYFit();

        float[] x = new float[]{1,  2, 3, 4, 5, 6, 7, 8, 9};
        float[] y = new float[]{9,  8, 7, 6, 5, 4, 3, 2, 1};
        /*
         *   *
         *      *
         *         *
         *            *
         */
        yfit.setX(x);
        yfit.setYFit(y);

        float[] expectedArea = new float[]{
            9, 8, 7, 6, 5, 4, 3, 2, 1
        };

        for (int i = 0; i < x.length; i++) {

            float histBinArea = yfit.calculateArea(i, true);
            float curvePointArea = yfit.calculateArea(i, false);

            assertTrue(histBinArea == expectedArea[i]);
            assertTrue(curvePointArea == expectedArea[i]);
        }

        for (int i = 0; i < x.length; i++) {
            assertTrue(yfit.getX(i) == x[i]);
        }

        float xScale = 2.0f;
        float yScale = 3.0f;
        yfit.setXScale(xScale);
        yfit.setYScale(yScale);

        for (int i = 0; i < x.length; i++) {
            yfit.getX()[i] /= xScale;
            yfit.getYFit()[i] /= yScale;
        }

        for (int i = 0; i < x.length; i++) {
            float histBinArea = yfit.calculateArea(i, true);
            float curvePointArea = yfit.calculateArea(i, false);

            assertTrue(histBinArea == expectedArea[i]);
            assertTrue(curvePointArea == expectedArea[i]);
        }
        
        assertTrue(yfit.getXPeak() == 1);
        assertTrue(yfit.getXScale() == xScale);
        assertTrue(yfit.getYScale() == yScale);
    }
    
    public void test0() throws Exception {
        // for coverage tool:
        CurveMisc cm = new CurveMisc();
    }
}
