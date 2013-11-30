package algorithms.curves;

import algorithms.misc.MiscMath;

/**
 * Generate curves following the Generalized Extreme Value probability density
 * function.
 *
 *                          (   (      ( x-mu))-(1/k))
 *                          (-1*(1 + k*(-----))      )
 *                 1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
 * y = y_const * ----- * exp                           * (1 + k*(-----))
 *               sigma                                   (      (sigma))
 *
 * mu is  the location parameter
 * sigma is the scale parameter and is > 0
 * k is the shape parameter
 *
 *
 * if k != 0,
 *     1 + (k * (x-mu)/sigma) > 0
 *
 *
 * for below, use z = (x-mu)/sigma
 *
 * Extreme Value Type I:
 *                     1
 *     y = y_const * ----- * exp( -z -exp(-z))
 *                   sigma
 *     sigma > 0
 *
 * Extreme Value Type II:
 *                     k     (sigma)^(k+1)      (   (sigma)^k)
 *     y = y_const * ----- * (-----)       * exp( - (-----)  )
 *                   sigma   (  x  )            (   (  x  )  )
 *
 *     k  > 0
 *     sigma > 0
 *     x > 0
 *
 * Extreme Value Type III:
 *                     k     (x - mu)^(k+1)      (   (x - mu)^k)
 *     y = y_const * ----- * (------)       * exp( - (------)  )
 *                   sigma   (sigma )            (   (sigma )  )
 *
 *     k > 0
 *     sigma > 0
 *     x > 0
 *
 * @author nichole
 */
public class GeneralizedExtremeValue implements ICurveGenerator {

    protected final float[] x;
    protected final float[] y;

    protected final float[] dx;
    protected final float[] dy;

    public GeneralizedExtremeValue(float[] xPoints, float[] yPoints, float[] dXPoints, float[] dYPoints) {
        this.x = xPoints;
        this.y = yPoints;
        this.dx = dXPoints;
        this.dy = dYPoints;
    }

    public float[] generateNormalizedCurve(float[] parameters, float yConst) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold k, sigma, and mu");
        }

        return generateNormalizedCurve(parameters[0], parameters[1], parameters[2], yConst);
    }

    public float[] generateNormalizedCurve(float[] parameters) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold k, sigma, and mu");
        }

        return generateNormalizedCurve(parameters[0], parameters[1], parameters[2]);
    }

    public float[] generateNormalizedCurve(float k, float sigma, float mu) {

        float[] yGEV = generateCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return null;
        }
        float yConst = determineYConstant(yGEV, mu);

        return generateNormalizedCurve(k, sigma, mu, yConst);
    }

    /**
     * generate a GEV curve w/ parameters k, sigma, and mu and then multiply it
     * by yConst so that the highest peak has value yConst
     *
     * @param k
     * @param sigma
     * @param mu
     * @param yConst
     * @return
     */
    public float[] generateNormalizedCurve(float k, float sigma, float mu, float yConst) {

        float[] yGEV = generateCurve(x, k, sigma, mu);
        if (yGEV == null) {
            return null;
        }

        float yMax = MiscMath.findMax(yGEV);

        for (int i = 0; i < yGEV.length; i++) {
            yGEV[i] *= yConst/yMax;
        }
        return yGEV;
    }

    public float[] generateCurve(float[] parameters) {
        if (parameters == null) {
            throw new IllegalArgumentException("parameters cannot be null");
        }
        if (parameters.length != 3) {
            throw new IllegalArgumentException("parameters must hold k, sigma, and mu");
        }

        return generateCurve(x, parameters[0], parameters[1], parameters[2]);
    }
    
    /*
     *                          (   (      ( x-mu))-(1/k))
     *                          (-1*(1 + k*(-----))      )
     *                 1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
     * y = y_const * ----- * exp                           * (1 + k*(-----))
     *               sigma                                   (      (sigma))
     *
     * mu is  the location parameter
     * sigma is the scale parameter and is > 0
     * k is the shape parameter
     * 
     * 
     * Components needed in the derivatives:
     * 
     *   Let z = (1 + k*( (x-mu)/sigma )
     *   
     */
    public static Double generateYGEV(float xPoint, float k, float sigma, float mu) {
        
        if (sigma == 0) {
            //throw new IllegalArgumentException("sigma must be > 0");
            return null;
        }
        
        if (k == 0) {
            // use Gumbel
            return generateYEVTypeI(xPoint, sigma, mu);
        }
        
        // if k == 0, use Gumbel which is Type I
        // if k > 0, use Frechet which is Type II
        // if k < 0, use Weibull which is Type III.  Not using k < 0 in this project

        float z = 1.f + k * ((xPoint - mu)/sigma);

        float a = -1.f * (float) Math.pow(z, (-1.f/k));
       
        if (Float.isInfinite(a)) {
            // k is very small
            return generateYEVTypeI(xPoint, sigma, mu);
        }

        float b = (float) Math.pow(z, (-1.f - (1.f/k)));

        if (Float.isInfinite(b)) {
            // k is very small
            return generateYEVTypeI(xPoint, sigma, mu);
        }
        
        if (Float.isNaN(a) || Float.isNaN(b)) {
            // or return 0?
            return null;
        } else {
            float t = (float) ((1.f/sigma) * Math.exp(a) * b);
            return Double.valueOf(t);
        }
    }
    
    public static Double generateYEVTypeI(float xPoint, float sigma, float mu) {
        
        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        double z = (xPoint - mu)/sigma;

        double a = (float) (-1.f - Math.exp(-1.0f*z));
         
        double yGEV = (float) ((1.f/sigma) * Math.exp(a));

        return yGEV;
    }
    

    public float[] generateCurve(float[] x1, float k, float sigma, float mu) {

        if (sigma == 0) {
            //throw new IllegalArgumentException("sigma must be > 0");
            return null;
        }

        if (k == 0) {
            return generateEVTypeICurve(sigma, mu);
        }

        float[] yGEV = new float[x1.length];

        for (int i = 0; i < x1.length; i++) {

            float z = 1.f + k*((x1[i] - mu)/sigma);

            float a = -1.f*(float) Math.pow(z, (-1.f/k));

            if (Float.isInfinite(a)) {
                return null;
            }

            float b = (float) Math.pow(z, (-1.f - (1.f/k)));

            if (Float.isInfinite(b)) {
                return null;
            }

            if (Float.isNaN(a) || Float.isNaN(b)) {
                // return null;
                yGEV[i] = 0;
            } else {
                float t = (float) ((1.f/sigma) * Math.exp(a) * b);
                yGEV[i] = t;
            }
        }

        return yGEV;
    }

    public float[] generateEVTypeICurve(float sigma, float mu) {

        if (sigma == 0) {
            throw new IllegalArgumentException("sigma must be > 0");
        }

        float[] yGEV = new float[x.length];

        for (int i = 0; i < x.length; i++) {

            float z = (x[i] - mu)/sigma;

            float a = (float) (-1.f - Math.exp(-1.0f*z));

            yGEV[i] = (float) ((1.f/sigma) * Math.exp(a));// times Math.exp(z) too???
        }

        return yGEV;
    }

    public float determineYConstant(float[] yGEV, float mu) {

        int index = MiscMath.findYMaxIndex(y);
        if (index == -1) {
            return 0;
        }

        return determineYConstant(yGEV, mu, index);
    }

    public float determineYConstant(float[] yGEV, float mu, int index) {

        float yConst = y[index]/yGEV[index];

        return yConst;
    }

    /**
     * Errors in fitting the curve can be calculated via chain rule of derivatives of the fit:
     *
     *                                  (   (      ( x-mu))-(1/k))
     *                                  (-1*(1 + k*(-----))      )
     *                         1        (   (      (sigma))      )  (      ( x-mu))(-1-(1/k))
     *     f(y) = y_const * ----- * exp                           * (1 + k*(-----))
     *                      sigma                                   (      (sigma))
     *
     *                               | df |^2               | df |^2         df   df
     *     (sigma_f)^2 =  (sigma_x)^2|----|   +  (sigma_y)^2|----|    +  2 * -- * -- * cov_x,y
     *                               | dx |                 | dy |           dx   dy
     *
     *     For uncorrelated variables the covariance terms are zero.
     *         else, cov(x,y) = <(x - mu_x)(y - mu_y)>
     *
     *     For example, for just 2 variables, x and y the equation for f(y) = XY reduces to:
     *         sigma^2  =  (Y^2)*xError^2  +  (X^2)*yError^2 + xError^2*yError^2
     *
     *     For the GEV, the equation differentials should be x, y, k, and sigma.
     *     (Within use of this code, mu is set rather than fit, so it's contribution is different.)
     *
     *     At this time, the full errors are not calculated, and instead the chi square
     *     is used.
     *
     *     TODO: The full derivative could be calculated if needed.
     *
     * @param y1 the y data array
     * @param yGEV the generated GEV y model array
     * @return the mean error of the fit
     */
    public static float calculateTotalMeanFittingError(float[] y1, float[] yGEV) {

        float yNorm = MiscMath.findMax(y1);

        float chiSum = 0.f;

        for (int i = 0; i < yGEV.length; i++) {

            float z = (yGEV[i] - y1[i])/yNorm;
            z *= z;

            chiSum += z;
        }

        return chiSum;
    }

    public static float calculateTotalMeanFittingError(int[] y1, float[] yGEV) {

        float yNorm = MiscMath.findMax(y1);

        float sum = 0.f;

        for (int i = 0; i < yGEV.length; i++) {

            float z = (yGEV[i] - y1[i])/yNorm;
            z *= z;

            sum += z;
        }

        return sum;
    }

    public static float[] generateNormalizedCurve(float[] x1, float k, float sigma, float mu) {

        if (sigma == 0) {
            //throw new IllegalArgumentException("sigma must be > 0");
            return null;
        }

        float yMax = Float.MIN_VALUE;

        float[] yGEV = new float[x1.length];

        for (int i = 0; i < x1.length; i++) {

            float z = 1.f + k*((x1[i] - mu)/sigma);

            float a = -1.f*(float) Math.pow(z, (-1.f/k));

            if (Float.isInfinite(a)) {
                return null;
            }

            float b = (float) Math.pow(z, (-1.f - (1.f/k)));

            if (Float.isInfinite(b)) {
                return null;
            }

            if (Float.isNaN(a) || Float.isNaN(b)) {
                // return null;
                yGEV[i] = 0;
            } else {
                float t = (float) ((1.f/sigma) * Math.exp(a) * b);
                yGEV[i] = t;
            }

            if (yGEV[i] > yMax) {
                yMax = yGEV[i];
            }
        }

        for (int i = 0; i < yGEV.length; i++){
            yGEV[i] /= yMax;
        }

        return yGEV;
    }
}
