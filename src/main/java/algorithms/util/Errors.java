package algorithms.util;

/**
 * convenience methods for generating errors used with unit tests.
 *
 * @author nichole
 */
public class Errors {

    /**
     * make assumption that errors are dominated by shot noise and so the noise is sqrt(y[i]).
     *
     * @param y
     * @return
     */
    public static float[] populateYErrorsBySqrt(float[] y) {

        float[] dy = new float[y.length];

        float maxError = 0.f;

        for (int i = 0; i < dy.length; i++) {
            if (y[i] > 0) {
                dy[i] =(float)(Math.sqrt(y[i]));
                if (dy[i] > maxError) {
                    maxError = dy[i];
                }
            }
        }
        for (int i = 0; i < dy.length; i++) {
            if (y[i] == 0) {
                dy[i] = maxError;
            }
        }

        return dy;
    }

    /**
     *
     * @param x
     * @param y
     * @param sigma
     * @param mu
     * @return
     */
    public static float[] calculateYErrorsByGaussian(float[] x, float[] y, float sigma, float mu) {

        /*                  1
         * p(x) dx = -----------------  exp( -(x - mu)^2/(2sigma^2) ) dx
         *           sigma * sqrt(2pi)
         */
        float[] erry = new float[y.length];

        float a = (float) (1.0f/(sigma * Math.sqrt(2.0f*Math.PI)));

        for (int i = 0; i < erry.length; i++) {

            float b = (x[i] - mu)/sigma;
            b *= b;

            b *= -1.0f/2.0f;

            float dx;
            if (i < (erry.length - 1)) {
                dx = x[i + 1] - x[i];
            } else {
                dx = x[i] - x[i-1];
            }

            erry[i] = (float) (a * y[i] * Math.exp(b) * dx);
        }

        return erry;
    }


    /**
     * Make assumption that errors are as large as half the distance between one point and another.
     * This roughly generates an error usable in tests.
     *
     * Note that usually, the errors that go into creating a histogram, for example,
     * are the measurement errors and systematic errors for the points used to generate
     * the histogram bins.  Those errors are normalized and added in quadrature.
     * The method here is used to replace such errors with the data resolution.
     * It notices that the x resolution of the histogram, for example, can be no
     * finer than half the distance between a point and it's nearest neighbor and
     * returns errors based upon that.
     *
     * @param x
     * @return
     */
    public static float[] populateXErrorsByPointSeparation(float[] x) {

        float[] errx = new float[x.length];

        // take minimum of nearest neighboring point instead of average of both
        float xdiffSum = 0;
        for (int i = 1; i < (errx.length - 1); i++) {
            float diff0 = x[i] - x[i - 1];
            float diff1 = x[i+1] - x[i];

            float diffMin = (diff0 < diff1) ? diff0 : diff1;

            errx[i] = diffMin/2.f;

            xdiffSum += errx[i];
        }

        // handle the 2 endpoints:
        xdiffSum /= (errx.length - 2);
        errx[0] = xdiffSum;
        errx[errx.length - 1] = xdiffSum;

        return errx;
    }

}
