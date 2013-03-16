package algorithms.misc;

import algorithms.util.Errors;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class HistogramTest extends TestCase {

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
     * Test of createHistogram method, of class Histogram.
     */
    public void testCreateHistogram_4args() {

        log.info("testCreateHistogram_4args");

        //  1  2  3  4  5
        //    0  1  2  3  4
        float[] a = new float[]{1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5};
        int nBins = 5;

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];

        Histogram.createHistogram(a, nBins, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 1));
            assertTrue(xHist[i] == (i + 1 + 0.5));
        }

        a = new float[]{0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4};
        xHist = new float[nBins];
        yHist = new int[nBins];

        Histogram.createHistogram(a, nBins, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 1));
            assertTrue(xHist[i] == (i + 0.5));
        }


        //  .1  .2  .3  .4  .5
        //   |   |   |   |   |
        //     0   1   2   3   4
        a = new float[]{0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.4f, 0.4f, 0.4f, 0.4f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};

        xHist = new float[nBins];
        yHist = new int[nBins];

        Histogram.createHistogram(a, nBins, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            float expected = (i + 1);
            float found = yHist[i];
            assertTrue(expected == found);

            expected = (float) ((i + 1 + 0.5)*0.1f);
            found = xHist[i];
            assertTrue( Math.abs(expected - found) < 0.01);
        }

        //  1  2  3  4  5
        //    0  1  2  3  4
        a = new float[]{
            100, 100, 100, 100, 100,
            200, 200, 200, 200, 200, 200,
            300, 300, 300, 300, 300, 300, 300,
            400, 400, 400, 400, 400, 400, 400, 400,
            500, 500, 500, 500, 500, 500, 500, 500, 500,
        };

        float[] ae = new float[a.length];
        for (int i = 0; i < a.length; i++) {
            ae[i] = 0.1f*(float)Math.sqrt(a[i]);
        }

        xHist = new float[nBins];
        yHist = new int[nBins];

        Histogram.createHistogram(a, nBins, xHist, yHist);
        for (int i = 0; i < yHist.length; i++) {
            float expected = (i + 5);
            float found = yHist[i];
            assertTrue(expected == found);

            expected = (float) ((i + 1 + 0.5)*100.f);
            found = xHist[i];
            assertTrue( Math.abs(expected - found) < 0.01);
        }

        float[] xHistErrorsOutput = new float[xHist.length];
        float[] yHistErrorsOutput = new float[xHist.length];

        Histogram.calulateHistogramBinErrors(xHist, yHist, a, ae, xHistErrorsOutput, yHistErrorsOutput);
        for (int i = 0; i < yHist.length; i++) {
            float yh = yHist[i];
            float xh = xHist[i];
            float xhe = xHistErrorsOutput[i];
            float yhe = yHistErrorsOutput[i];
            //TODO:  add assert of rough expected values
            assertTrue(xhe < xh);
            assertTrue(yhe < yh);
        }

        HistogramHolder hist = Histogram.createHistogramForSkewedData(nBins, a, ae);
        assertTrue(hist.xHist.length == nBins);
        for (int i = 0; i < hist.yHist.length; i++) {
            float expected = (i + 5);
            float found = hist.yHist[i];
            assertTrue(expected == found);

            expected = (float) ((i + 1 + 0.5)*100.f);
            found = hist.xHist[i];
            assertTrue( Math.abs(expected - found) < 0.01);
        }
    }

    /**
     * Test of createHistogram method, of class Histogram.
     */
    public void testCreateHistogram_6args() {

        log.info("testCreateHistogram_6args");

        //  2   3   4   5   6
        //    0   1   2   3   4
        float[] a = new float[]{2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6};
        int nBins = 5;

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];

        Histogram.createHistogram(a, nBins, 2, 6, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 2));
            assertTrue(xHist[i] == (i + 2 + 0.5));
        }

        nBins = 4;
        xHist = new float[nBins];
        yHist = new int[nBins];
        Histogram.createHistogram(a, nBins, 2, 5, xHist, yHist);

        for (int i = 0; i < yHist.length; i++) {
            assertTrue(yHist[i] == (i + 2));
            assertTrue(xHist[i] == (i + 2 + 0.5));
        }
    }

    /**
     * Test of createHistogram method, of class Histogram.
     */
    public void testCreateHistogram_3args() {

        log.info("createHistogram");

        float[] a = new float[]{
            2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11,
            12, 13, 14, 15
        };
        float[] ae = Errors.populateYErrorsBySqrt(a);

        int nBins = 12;

        HistogramHolder result = Histogram.createHistogramForSkewedData(nBins, a, ae);

        assertTrue(result.getXHist().length >= 6);
    }

    /**
     * Test of calulateHistogramBinErrors method, of class Histogram.
     */
    public void testCalulateHistogramBinErrors() {

        log.info("calulateHistogramBinErrors");

        float[] values = new float[]{
            200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
            300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
            400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400,
            500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500,
            600, 600, 600, 600, 600, 600, 600, 600, 600, 600,
            700, 700, 700, 700, 700, 700, 700, 700, 700,
            800, 800, 800, 800, 800, 800, 800, 800,
            900, 900, 900, 900, 900, 900, 900,
            1000, 1000, 1000, 1000, 1000,
            1100, 1100, 1100, 1100, 1100,
            1200, 1200, 1200, 1200, 1200,
            1300, 1300, 1300, 1300, 1300,
            1400, 1400, 1400, 1400, 1400,
            1500, 1500, 1500, 1500, 1500
        };
        float[] valueErrors = new float[values.length];
        for (int i = 0; i < values.length; i++) {
            valueErrors[i] = values[i]/1000.0f;
        }

        HistogramHolder hist = Histogram.createHistogramForSkewedData(values.length, values, valueErrors);

        for (int i = 0; i < hist.xHist.length; i++) {
            float xh = hist.xHist[i];
            float yh = hist.yHist[i];
            float xhe = hist.xErrors[i];
            float yhe = hist.yErrors[i];
            assertTrue(xhe < xh);
            //assertTrue(yhe < yh);
        }
    }

    /**
     * Test of calculateHistogramAreaError method, of class Histogram.
     */
    public void testCalculateHistogramAreaError_4args() {

        log.info("calculateHistogramAreaError");

        //float result = Histogram.calculateHistogramAreaError(xHist, yHist, xErrors, yErrors);

    }
}
