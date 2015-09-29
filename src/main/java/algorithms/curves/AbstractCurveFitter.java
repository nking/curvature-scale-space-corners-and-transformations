package algorithms.curves;

import java.io.IOException;
import java.security.SecureRandom;
import java.util.Arrays;

import algorithms.curves.GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM;
import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public abstract class AbstractCurveFitter implements ICurveFitter {

    /**
     *
     */
    protected final float[] x;

    /**
     *
     */
    protected final float[] y;

    /**
     *
     */
    protected final float[] xe;

    /**
     *
     */
    protected final float[] ye;

    /**
     *
     */
    protected final float xmin; 

    /**
     *
     */
    protected final float xmax;

    /**
     *
     */
    protected final float ymin;

    /**
     *
     */
    protected final float ymax;

    /**
     *
     */
    protected float xScale = -1;
    
    // the factor by which the y axis can be multiplied to return it to the original values

    /**
     *
     */
        protected float yScale = -1;
    
    /**
     *
     */
    protected final GeneralizedExtremeValue gev;
    
    /**
     *
     */
    protected boolean debug = false;
    
    /**
     *
     * @param xPoints
     * @param yPoints
     * @param xErrPoints
     * @param yErrPoints
     */
    public AbstractCurveFitter(float[] xPoints, float[] yPoints,
        float[] xErrPoints, float[] yErrPoints) {

        if (xPoints == null) {
            throw new IllegalArgumentException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalArgumentException("yPoints cannot be null");
        }
        if (xErrPoints == null) {
            throw new IllegalArgumentException("xErrPoints cannot be null");
        }
        if (yErrPoints == null) {
            throw new IllegalArgumentException("yErrPoints cannot be null");
        }

        float[] tmp = Arrays.copyOf(xPoints, xPoints.length);
        this.xScale = scaleDataTo1(tmp);
        this.x = tmp;

        tmp = Arrays.copyOf(yPoints, yPoints.length);
        this.yScale = scaleDataTo1(tmp);
        this.y = tmp;

        tmp = Arrays.copyOf(xErrPoints, xErrPoints.length);
        scaleDataTo1(tmp, xScale);
        this.xe = tmp;

        tmp = Arrays.copyOf(yErrPoints, yErrPoints.length);
        scaleDataTo1(tmp, yScale);
        this.ye = tmp;

        this.gev = new GeneralizedExtremeValue(x, y, xe, ye);
        
        xmin = x[0];
        xmax = x[x.length - 1];
        ymin = 0.0f;
        ymax = MiscMath.findMax(y);
      
    }
    
    /**
     *
     * @param a
     * @return
     */
    protected final float scaleDataTo1(float[] a) {
        float max = MiscMath.findMax(a);
        scaleDataTo1(a, max);
        return max;
    }

    /**
     *
     * @param a
     * @param scale
     */
    protected final void scaleDataTo1(float[] a, float scale) {
        if (a == null) {
            return;
        }
        for (int i = 0; i < a.length; i++) {
            a[i] = a[i]/scale;
        }
    }

    /**
     * turn debugging results in comments such as intermediate fit results being
     * printed to standard out and plots being generated.
     *
     * @param doUseDebug
     */
    public void setDebug(boolean doUseDebug) {
        debug = doUseDebug;
    }

    
    /**
     * calculate the chi square of the yModel using the method for weights.
     * The best results for the fit are usually from using errors for the
     * fit.
     *
     * @param yNormalizedModel
     * @param wdc
     * @return
     */
    public float calculateChiSquareSum(float[] yNormalizedModel, WEIGHTS_DURING_CHISQSUM wdc) {

        if (yNormalizedModel == null) {
            return Float.POSITIVE_INFINITY;
        }

        float[] w = calcWeights(wdc, yNormalizedModel);

        float chiSum = 0.f;

        for (int i = 0; i < yNormalizedModel.length; i++) {

            float z = yScale*(yNormalizedModel[i] - y[i]);
            z *= z*w[i];

            chiSum += z;
        }

        return chiSum;
    }

    /**
     * create the weight array used by the chi square sum method for the given
     * weight method.
     *
     * @param wdc
     * @param yModel
     * @return
     */
    float[] calcWeights(WEIGHTS_DURING_CHISQSUM wdc, float[] yModel) {

        float[] w = new float[yModel.length];

        if ((wdc != null) && (wdc.ordinal() == WEIGHTS_DURING_CHISQSUM.INVERSE_Y.ordinal())) {
            for (int i = 0; i < w.length; i++) {
                w[i] = 1.0f/(yScale*y[i]);
            }
         } else if ((wdc != null) && (wdc.ordinal() == WEIGHTS_DURING_CHISQSUM.MODEL_Y.ordinal())) {
            for (int i = 0; i < w.length; i++) {
                w[i] = yModel[i]*yScale;
            }
         } else {
            boolean hasValidValues = hasValidValues(ye);
            if (hasValidValues) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = ye[i]*yScale;
                }
            } else {
                throw new IllegalStateException("dy has invalid values");
            }
        }

        return w;
    }

    /**
     *
     * @param a
     * @return
     */
    protected boolean hasValidValues(float[] a) {
        if (a == null) {
            return false;
        }
        for (int i = 0; i < a.length; i++) {
            if (Float.isNaN(a[i])) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * given the GEV model points yModel, calculate the chi-square statistic.
     *
     * @param yModel
     * @param wdc
     * @return
     */
    public float calculateChiSquareStatistic(float[] yModel, WEIGHTS_DURING_CHISQSUM wdc) {

        float chiSquareSum = calculateChiSquareSum(yModel, wdc);

        // minus one due to having to calculate the mean by using the data
        float degreesOfFreedom = yModel.length - 3 - 1;

        return chiSquareSum/degreesOfFreedom;
    }

    /**
     * return the y array errors summed in quadrature
     * @return
     */
    public float calcYErrSquareSum() {
        float sum = 0.f;
        for (int i = 0; i < y.length; i++) {
            float z = yScale*ye[i];
            z *= z;
            sum += z;
        }
        return sum;
    }
    
    /**
    generate randomly values for k and sigma between a default of

    float kMin = 0.00001f;
    float kMax = 0.001f;
    float sigmaMin = 0.025f;
    float sMax = 20.0f*sigmaMin;

    * @param sr
    * @param mu
    * @return an array of k and sigma, respectively
    */
   public static float[] generateRandomParameters(SecureRandom sr, float mu) {
       float kMin = 0.00001f;
       float kMax = 0.001f;
       float sigmaMin = 0.025f;
       float sigmaMax = 20.0f*sigmaMin;
       return generateRandomParameters(sr, kMin, kMax, sigmaMin, sigmaMax, mu);
   }

   /**
    * generate randomly values for k and sigma
    *
    * @param sr
    * @param kMin
    * @param kMax
    * @param sigmaMin
    * @param sigmaMax
    * @param mu
    * @return an array of k and sigma, respectively
    */
   public static float[] generateRandomParameters(SecureRandom sr, float kMin, float kMax, float sigmaMin, float sigmaMax, float mu) {

       float sDiff = sigmaMax - sigmaMin;
       float kDiff = kMax - kMin;

       float sigma = sigmaMin + sr.nextFloat() * sDiff;

       float k = kMin + sr.nextFloat() * kDiff;

       boolean go = true;
       while (go) {
           if ( (k * mu) > -1*sigma) {

               go = false;

               k = kMin + sr.nextFloat() * kDiff;
           }
       }

       return new float[]{k, sigma};
   }

    /**
     *
     * @param weightMethod
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public abstract GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) 
        throws FailedToConvergeException, IOException;
       
}
