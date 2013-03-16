package algorithms.curves;

import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

public class ATest extends TestCase {

    protected float[] x = null;
    protected float[] y = null;
    protected float[] dx = null;
    protected float[] dy = null;

    protected boolean debug = true;

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testCalculateChiSqSumForCurve_0() throws Exception {

        // not a test.  this is a way to fit a curve by eye if needed.

        x = new float[]{0.017533371f, 0.052600116f, 0.087666854f, 0.1227336f, 0.15780035f,
            0.19286709f, 0.22793384f, 0.26300058f, 0.2980673f, 0.33313406f, 0.36820078f,
            0.40326753f, 0.4383343f, 0.473401f, 0.5084678f, 0.5435345f, 0.57860124f,
            0.61366796f, 0.64873475f, 0.6838015f};
        y = new float[]{0.0f, 0.0f, 50.0f, 220.0f, 340.0f, 357.0f, 339.0f, 295.0f,
            217.0f, 199.0f, 163.0f, 135.0f, 113.0f, 87.0f, 85.0f, 48.0f, 47.0f,
            44.0f, 31.0f, 40.0f};

        dy = Errors.populateYErrorsBySqrt(y);
        dx = Errors.populateXErrorsByPointSeparation(x);

        /*
         * k=1.7483624E-4 sigma=0.10408807 mu=0.2820513 chiSqSum=859361.06 chiSqStatistic=53710.066
         *
         * k=4.2631186E-4 sigma=0.13125233 mu=0.28205 chiSqSum=364.20374 chiSqStatistic=22.762733
         */
        float k = 0.00034782593f;
        float sigma = 0.1070881f;
        int yConstIndex = 1;
        float mu = 0.28205f;
        float yNorm = 1.0f;

        GEVChiSquareMinimization chiSqMin = new GEVChiSquareMinimization(x, y, dx, dy);

        GEVYFit yfit = null;

        //yfit = chiSqMin.calculateChiSqSumAndCurve(k, sigma, mu, GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yConstIndex);

        yfit = chiSqMin.fitCurve(k/2.f, k*2.f, sigma/2.f, sigma*2.f, mu, 808.0f,
            GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

        //yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        if (debug) {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 0.7f, 0, 360);
            plotter.addPlot(x, y, yfit.getOriginalScaleX(), yfit.getOriginalScaleYFit(), "");
            plotter.writeFile();

            System.out.println(yfit.toString());
        }
    }

    public void estCalculateChiSqSumForCurve_1() throws Exception {

        float[] x = new float[20];
        float xDelta = 1.0f/x.length;
        for (int i = 0; i < x.length; i++) {
            x[i] = xDelta*i;
        }

        //* sigma range:    0.025 at       50000.0f*kMin * x[1];  where k={1.0E-5:100.0}  first several
        //*     and then    0.25  at      500000.0f*kMin * x[1]; for k={1.0E-5:100.0}
        //*     and then    2.5   at     5000000.0f*kMin * x[1]; for k={1.0E-5:100.0}
        //*     and then   25.    at    50000000.0f*kMin * x[1]; for k={1.0E-5:100.0}  <=== nearly horizontal lines
        //      and then  250.    at   500000000.0f*kMin * x[1]; for k={1.0E-5:100.0}         same
        //float kMin = 0.000001f;
        //float kMax = 1000.0f;
        //float mu = x[1];

        float kMin = 0.00001f;
        float kMax = 0.001f;
        float mu = x[1];

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.0f, 1.0f, 0.0f, 1.0f);
        //  1 + k*(-0.5)/sigma <=== if k == (1/(|x1-x0|))*sigma, total is zero

        float k = kMin;
        //float sigma = 50000.0f*kMin * x[1];// to kMax*x[x.length - 1];
        float sigma = 0.025f;
        float s = sigma;
        float sMax = 20.0f*s;
        System.out.println("k={" + kMin + ":" + kMax + "} sigma=" + sigma);
        while (k < kMax) {

            while (s < sMax) {
                float[] y = GeneralizedExtremeValue.generateNormalizedCurve(x, k, s, mu);

                if (y != null) {
                    String str = String.format("k=%.7f s=%.7f", k, s);
                    plotter.addPlot(x, y, x, y, str);
                    plotter.writeFile();
                    System.out.println(str);
                }
                s*=1.1f;
            }

            k*=1.1f;
        }
    }
}
