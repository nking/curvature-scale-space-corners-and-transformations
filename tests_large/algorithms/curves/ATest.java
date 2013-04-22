package algorithms.curves;

import algorithms.curves.GEVChiSquareMinimization;
import algorithms.curves.GEVYFit;
import algorithms.curves.GeneralizedExtremeValue;
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

    public void test() {}

    public void estCalculateChiSqSumForCurve_0() throws Exception {

        // not a test.  this is a way to fit a curve by eye if needed.
        x = new float[]{0.04042134f, 0.121264026f, 0.2021067f, 0.2829494f, 0.36379206f, 0.44463474f, 0.52547747f, 0.60632014f, 0.6871628f, 0.7680055f, 0.84884816f, 0.92969084f, 1.0105336f, 1.0913762f, 1.1722189f, 1.2530615f, 1.3339043f, 1.414747f, 1.4955896f, 1.5764323f, 1.657275f};
        y = new float[]{18.0f, 32.0f, 39.0f, 23.0f, 18.0f, 10.0f, 7.0f, 6.0f, 6.0f, 5.0f, 3.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f};

        dy = Errors.populateYErrorsBySqrt(y);
        dx = Errors.populateXErrorsByPointSeparation(x);

        /*
         * k=1.5e-04 s=9.2e-02 m=1.2e-01 chisq=2284.0 yerrsq=200.8
         * k=1.2529743E-4 sigma=0.0803195 mu=0.09756097 chiSqSum=16.608713 chiSqStatistic=0.97698313
         */

        //        0.00001f;
        float k = 0.000058f;
        float sigma = 0.14f;
        float mu = 0.1f;
        float yNorm = 1.0f;
        int yConstIndex = 2;

        GEVChiSquareMinimization chiSqMin = new GEVChiSquareMinimization(x, y, dx, dy);

        GEVYFit yfit = null;

        //yfit = chiSqMin.calculateChiSqSumAndCurve(k, sigma, mu,
        //    GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        //yfit = chiSqMin.fitCurve(k/2.f, k*2.f, sigma/2.f, sigma*2.f, mu, 808.0f,
        //    GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

        //yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

        chiSqMin.setDebug(true);
        yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS,
            GEVChiSquareMinimization.kMinDefault, GEVChiSquareMinimization.kMaxDefault,
            GEVChiSquareMinimization.sigmaMinDefault, GEVChiSquareMinimization.sigmaMaxDefault);

        yfit.setChiSqStatistic(chiSqMin.calculateChiSquareStatistic(yfit.getYFit(),
            GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS));

        if (debug) {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
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
        float kMax = 2.0f; // 0.001f;
        float mu = x[1];
        float muMax = x[x.length/2];

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.0f, 1.0f, 0.0f, 1.0f);
        //  1 + k*(-0.5)/sigma <=== if k == (1/(|x1-x0|))*sigma, total is zero

        float k = kMin;
        //float sigma = 50000.0f*kMin * x[1];// to kMax*x[x.length - 1];
        float sigma = 0.025f;
        float s = sigma;
        float sMax = 20.0f*s;
        System.out.println(String.format("k={%.7f : %.7f} sigma=%.7f mu=%.4f", kMin, kMax, sigma, mu));
        while (mu < muMax) {

            while (k < kMax) {

                while (s < sMax) {
                    float[] y = GeneralizedExtremeValue.generateNormalizedCurve(x, k, s, mu);

                    if (y != null) {
                        String str = String.format("k=%.7f s=%.7f m=%.4f", k, s, mu);
                        plotter.addPlot(x, y, x, y, str);
                        plotter.writeFile();
                        System.out.println(str);
                    }
                    s*=1.5f;
                }

                k *= 1.5f;
                s = sigma;
            }
            k = kMin;
            s = sigma;
            mu = mu*=1.1f;
        }
    }
}
