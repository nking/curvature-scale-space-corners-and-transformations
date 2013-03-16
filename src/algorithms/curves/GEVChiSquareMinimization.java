package algorithms.curves;

import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * A chi square minimization routine for the Generalized Extreme Value function.
 * It has been tailored for use with the two-point correlation algorithm,
 * that is for k > 0.
 *
 * If one knows the range of parameter space for k and sigma, it is better to
 * use methods which calculate the fit with that better knowledge.
 *
 * mu is  the location parameter
 * sigma is the scale parameter and is > 0
 * k is the shape parameter
 *
 *                          (   (      ( x-mu))-(1/k))
 *                          (-1*(1 + k*(-----))      )
 *                 1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
 * y = y_const * ----- * exp                           * (1 + k*(-----))
 *               sigma                                   (      (sigma))
 *
 * @author nichole
 */
public class GEVChiSquareMinimization {

    protected final float[] x;
    protected final float[] y;

    protected final float[] dx;
    protected final float[] dy;

    protected float xScale = -1;
    protected float yScale = -1;

    public static final float kMinDefault = 0.00001f;
    public static final float kMaxDefault = 0.001f;
    public static final float sigmaMinDefault = 0.025f;
    public static final float sigmaMaxDefault = 0.5f;

    protected final GeneralizedExtremeValue gev;

    protected boolean debug = false;

    public enum WEIGHTS_DURING_CHISQSUM {
        ERRORS, INVERSE_Y, MODEL_Y
    }

    public GEVChiSquareMinimization(float[] xPoints, float[] yPoints,
        float[] dXPoints, float[] dYPoints) {

        if (xPoints == null) {
            throw new IllegalArgumentException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalArgumentException("yPoints cannot be null");
        }
        if (dXPoints == null) {
            throw new IllegalArgumentException("dXPoints cannot be null");
        }
        if (dYPoints == null) {
            throw new IllegalArgumentException("dYPoints cannot be null");
        }

        float[] tmp = Arrays.copyOf(xPoints, xPoints.length);
        this.xScale = scaleDataTo1(tmp);
        this.x = tmp;

        tmp = Arrays.copyOf(yPoints, yPoints.length);
        this.yScale = scaleDataTo1(tmp);
        this.y = tmp;

        tmp = Arrays.copyOf(dXPoints, dXPoints.length);
        scaleDataTo1(tmp, xScale);
        this.dx = tmp;

        tmp = Arrays.copyOf(dYPoints, dYPoints.length);
        scaleDataTo1(tmp, yScale);
        this.dy = tmp;

        this.gev = new GeneralizedExtremeValue(x, y, dx, dy);
    }

    protected final float scaleDataTo1(float[] a) {
        float max = MiscMath.findMax(a);
        scaleDataTo1(a, max);
        return max;
    }
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
     * create the weight array used by the chi square sum method for the given
     * weight method.
     *
     * @param wdc
     * @param yModel
     * @return
     */
    float[] calcWeights(WEIGHTS_DURING_CHISQSUM wdc, float[] yModel) {

        float[] w = new float[yModel.length];

        if (wdc != null) {
            if (wdc.ordinal() == WEIGHTS_DURING_CHISQSUM.INVERSE_Y.ordinal()) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = 1.0f/(yScale*y[i]);
                }
            } else if (wdc.ordinal() == WEIGHTS_DURING_CHISQSUM.MODEL_Y.ordinal()) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = yModel[i]*yScale;
                }
            } else {
                boolean hasValidValues = hasValidValues(dy);
                if (hasValidValues) {
                    for (int i = 0; i < w.length; i++) {
                        w[i] = dy[i]*yScale;
                    }
                } else {
                    throw new IllegalStateException("dy has invalid values");
                }
            }
        } else {
            // defaults to using errors
            boolean hasValidValues = hasValidValues(dy);
            if (hasValidValues) {
                for (int i = 0; i < w.length; i++) {
                    w[i] = dy[i]*yScale;
                }
            } else {
                throw new IllegalStateException("dy has invalid values");
            }
        }

        return w;
    }

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
     * calculate the chi square of the yModel using the method for weights.
     * The best results for the fit are usually from using errors for the
     * fit.
     *
     * @param yModel
     * @param wdc
     * @return
     */
    public float calculateChiSquareSum(float[] yModel, WEIGHTS_DURING_CHISQSUM wdc) {

        if (yModel == null) {
            return Float.POSITIVE_INFINITY;
        }

        float[] w = calcWeights(wdc, yModel);

        float chiSum = 0.f;

        for (int i = 0; i < yModel.length; i++) {

            float z = yScale*(yModel[i] - y[i])/w[i];
            z *= z;

            chiSum += z;
        }

        return chiSum;
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
            float z = yScale*dy[i];
            z *= z;
            sum += z;
        }
        return sum;
    }

    /**
     * Fit the curve using a simple constraint for mu and a narrow range of possible values
     * for k and sigma.
     *
     * Mu is usually found to be the maximum of the y array of the data because the
     * code is intended to be used on GEV distributions, more specifically, those of
     * EV Type I and Type II with k > 0.
     *
     *
     * The initial k range is 0.001f to 10.f
     * and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     *
     * The method uses one iteration of a downhill simplex method for fitting.
     *
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {

        int yMaxIndex = MiscMath.findYMaxIndex(y);

        if (yMaxIndex == -1) {
            // all y's were zero in y
            return null;
        }
        float yNorm = y[yMaxIndex];
        float mu = x[yMaxIndex];

        GEVYFit bestFit = fitCurveKGreaterThanZero(weightMethod, mu, yNorm);

        return bestFit;
    }

    /**
     * fit a GEV curve to the data by using the given range of k and sigma.
     *
     * Internally, the method uses a downhill simplex algorithm for fitting.
     *
     * @param weightMethod the method for determing the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param mu
     * @param yNorm
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax, float mu, float yNorm)
        throws FailedToConvergeException, IOException {

        float yErrSquareSum = calcYErrSquareSum();

        GEVYFit bestFit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

        if (debug) {

            String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum,
                bestFit.getChiSqStatistic());

            System.out.println(str);
        }

        return bestFit;
    }

    /**
     * fit a GEV curve to the data by using the given range of k and sigma.
     *
     * Internally, the method uses a downhill simplex algorithm for fitting.
     *
     * @param weightMethod the method for determing the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @param k
     * @param sigma
     * @param mu
     * @param yNorm
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod,
        float k, float sigma, float mu, float yNorm) throws FailedToConvergeException, IOException {

        GEVYFit bestFit = calculateChiSqSumAndCurve(k, sigma, mu, weightMethod, yNorm);

        return bestFit;
    }

    /**
     * Fit the curve using a simple constraint for mu and a narrow range of possible values
     * for k and sigma.
     *
     * Mu is usually found to be the maximum of the y array of the data because the
     * code is intended to be used on GEV distributions, more specifically, those of
     * EV Type I and Type II with k > 0.
     *
     * The initial k range is 0.001f to 10.f
     * and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     *
     * The method uses one iteration of a downhill simplex method for fitting.
     *
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @param mu
     * @param yNorm
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod, float mu, float yNorm)
        throws FailedToConvergeException, IOException {

        // start with ranges useful for EV Type I and II for the two-point correlation functions
        float deltaX = x[1] - x[0];

        float kMin = kMinDefault;
        float kMax = kMaxDefault;
        float sigmaMin = sigmaMinDefault;
        float sigmaMax = sigmaMaxDefault;

        if (debug) {
            String str = String.format("kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f", kMin, kMax, sigmaMin, sigmaMax);
            System.out.println(str);
        }

        float yErrSquareSum = calcYErrSquareSum();

        GEVYFit bestFit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

        if (debug) {

            String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum,
                bestFit.getChiSqStatistic());

            System.out.println(str);
        }

        return bestFit;
    }

    /**
     * Fit the curve using a simple constraint for mu and a wide range of possible values
     * for k and sigma.
     *
     * Mu is usually found to be the maximum of the y array of the data because the
     * code is intended to be used on GEV distributions, more specifically, those of
     * EV Type I and Type II with k > 0.
     *
     * This method starts with a very large default range of values for k and sigma and
     * uses range reduction techniques to reduce the range while iteratively improving the
     * fit and the range.
     *
     * The initial k range is 0.00001f to 100.f
     * and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     *
     * The maximum number of iterations of range reduction + fitting is
     *    2 if the first chi square statistic is < 10
     *    10 if the first chi square statistic is larger than 10
     *
     * For the range reduction, a grid search through parameter space is
     * used and the outer boundaries of the k and sigma range of the best solution are kept.
     * For the fitting, a downhill simplex method is used.
     *
     * The iterations are performed until the chisquared statistic is < 2 or
     * the number of iterations is less than a maximum number.
     *
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZeroWithLargeKAndSigmaRange(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {

        int yMaxIndex = MiscMath.findYMaxIndex(y);

        if (yMaxIndex == -1) {
            // all y's were zero in y
            return null;
        }
        float mu = x[yMaxIndex];
        float yNorm = y[yMaxIndex];

        GEVYFit bestFit = fitCurveKGreaterThanZeroWithRangeReduction(weightMethod, mu, yNorm);

        return bestFit;
    }

    /**
     * Fit the curve using the provided index for mu and a wide range of possible values
     * for k and sigma.
     *
     * Mu is usually found to be the maximum of the y array of the data because the
     * code is intended to be used on GEV distributions, more specifically, those of
     * EV Type I and Type II with k > 0.  For that reason, the provided yConstIndex is
     * used as mu for the fit.
     *
     * This method starts with a very large default range of values for k and sigma and
     * uses range reduction techniques to reduce the range while iteratively improving the
     * fit and the range.
     *
     * The initial k range is 0.00001f to 100.f
     * and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     *
     * The maximum number of iterations of range reduction + fitting is
     *    2 if the first chi square statistic is < 10
     *    10 if the first chi square statistic is larger than 10
     *
     * For the range reduction, a grid search through parameter space is
     * used and the outer boundaries of the k and sigma range of the best solution are kept.
     * For the fitting, a downhill simplex method is used.
     *
     * The iterations are performed until the chisquared statistic is < 2 or
     * the number of iterations is less than a maximum number.
     *
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZeroWithRangeReduction(WEIGHTS_DURING_CHISQSUM weightMethod, float mu, float yNorm)
        throws FailedToConvergeException, IOException {

        // start with ranges useful for EV Type I and II for the two-point correlation functions
        float deltaX = x[1] - x[0];

        float kMin = 0.00001f;
        float kMax = 0.001f;
        float sigmaMin = 0.025f;
        float sigmaMax = 20.0f*sigmaMin;

        if (debug) {
            String str = String.format("kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f mu=%.7f", kMin, kMax, sigmaMin, sigmaMax, mu);
            System.out.println(str);
        }

        float yErrSquareSum = calcYErrSquareSum();

        float[] smallerMinMaxRanges = findSmallerRange(weightMethod, kMin, kMax,
            sigmaMin, sigmaMax, mu, yErrSquareSum, yNorm);

        kMin = smallerMinMaxRanges[0];
        kMax = smallerMinMaxRanges[1];
        sigmaMin = smallerMinMaxRanges[2];
        sigmaMax = smallerMinMaxRanges[3];

        if (debug) {
            String str2 = String.format("reduced=> kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f",
                kMin, kMax, sigmaMin, sigmaMax);
            System.out.println(str2);
        }

        GEVYFit bestFit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

        if (debug) {
            String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(),
                yErrSquareSum, bestFit.getChiSqStatistic());
            System.out.println(str);
        }

        boolean statisticIsHigh = (bestFit.getChiSqStatistic() > 10);

        int nIter = 0;
        int nIterMax = 10;

        int nIter2 = 0;
        int nIterMax2 = 2;

        while ( (statisticIsHigh || (nIter2 < nIterMax2)) && (nIter < nIterMax)) {

            float[] minMaxRanges = findSmallerRange(weightMethod, kMin, kMax, sigmaMin, sigmaMax, mu,
                yErrSquareSum, yNorm);

            kMin = minMaxRanges[0];
            kMax = minMaxRanges[1];
            sigmaMin = minMaxRanges[2];
            sigmaMax = minMaxRanges[3];

            if (debug) {
                String str = String.format("reduced=> kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f mu=%.7f",
                    kMin, kMax, sigmaMin, sigmaMax, mu);
                System.out.println(str);
                System.out.println("switch to fine grid of ranges + downhill simplex for each");
            }

            int kPower = MiscMath.findPowerOf10(bestFit.getK());
            float kFactorStep = ((kPower < 0) && (nIter > 0)) ? 0.01f : 0.1f;

            GEVYFit yfit = fitCurveKGreaterThanZeroUsingSmallSteps(weightMethod,
                kMin, kMax, sigmaMin, sigmaMax, mu,
                kFactorStep, 5.0f, yErrSquareSum, yNorm);

            if (yfit != null) {

                if (debug) {
                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                        yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum, yfit.getChiSqStatistic());
                    System.out.println(str);
                }

                if (yfit.getChiSqSum() < bestFit.getChiSqSum()) {

                    bestFit = yfit;

                    if (debug) {
                        String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                            bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum,
                            bestFit.getChiSqStatistic());
                        System.out.println(str);
                    }

                    if (bestFit.getChiSqStatistic() > 2) {
                        System.out.println("statistic is too high!");
                    } else {
                        break;
                    }
                }
            }
            nIter++;
            nIter2++;
            statisticIsHigh = (bestFit.getChiSqStatistic() > 10);
        }

        if (debug) {

            String label = String.format("k=%.1e s=%.1e m=%.1e chisq=%.1f yerrsq=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum);

            plotFit(bestFit, label);

            System.out.println("fitted gev has chisqsum=" + bestFit.getChiSqSum() + " while curve yerrsqsum=" + yErrSquareSum);
        }

        return bestFit;
    }

    /**
     * Fit curve within range of given parameters.  Internally, uses the Neder-Meade downhill
     * simplex method and no additional iterations combined with range reduction such as
     * is used by the fitCurveKGreaterThanZero method without min and max arguments.
     *
     * The mu parameter is fit for x[0] through x[yConstIndex] only.
     *
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZeroAndMu(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax) throws FailedToConvergeException, IOException {

        // using mu range useful for EV Type I and II for the two-point correlation functions

        if (debug) {
            String str = String.format("*kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f", kMin, kMax, sigmaMin, sigmaMax);
            System.out.println(str);
        }

        float yErrSquareSum = calcYErrSquareSum();

        GEVYFit bestFit = null;

        int yMaxIndex = MiscMath.findYMaxIndex(y);

        if (yMaxIndex == -1) {
            // all y's were zero in y
            return null;
        }
        float mu = x[yMaxIndex];
        float yPeak = y[yMaxIndex];

        int end = (yMaxIndex + 5);
        if (end > (y.length - 1)) {
            end = y.length - 1;
        }

        for (int i = 0; i <= end; i++) {

            float yNorm = y[i];

            for (int ii = 0; ii < 2; ii++) {

                float dmu;
                if (i < (y.length - 1)) {
                    dmu = (x[i] - x[i + 1])/2.0f;
                } else {
                    dmu = (x[0] - x[1])/2.0f;
                }
                mu = x[i] + ii*dmu;

                GEVYFit yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

                if (debug && (yfit != null)) {
                    String label = String.format("k=%.1e s=%.1e m=%.1e chisq=%.1f yerrsq=%.1f",
                        yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum);
                    System.out.println(label);
                }

                if ( (bestFit == null) || (yfit.getChiSqSum() < bestFit.getChiSqSum()) ) {
                    bestFit = yfit;
                }
            }
        }

        if (bestFit != null) {

            String label = String.format("k=%.1e s=%.1e m=%.1e chisq=%.1f yerrsq=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum);

            plotFit(bestFit, label);

            if (debug) {
                System.out.println("fitted gev has chisqsum=" + bestFit.getChiSqSum() + " while curve yerrsqsum=" + yErrSquareSum);
            }
        }

        return bestFit;
    }

    /**
     * fit the curve by randomly drawing  k and sigma parameters and keeping the
     * best fit.  The algorithm stops after a maximum number of iterations.
     * In general, this algorithm tends not to find the best fit for the GEV curve.
     *
     * @param weightMethod
     * @param yErrSquareSum
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    protected GEVYFit fitCurveKGreaterThanZeroUsingRandomParameters(WEIGHTS_DURING_CHISQSUM weightMethod, float mu,
        float yErrSquareSum, float yNorm) throws FailedToConvergeException, IOException, NoSuchAlgorithmException {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( System.currentTimeMillis() );

        GEVYFit bestFit = null;

        boolean go = true;
        int nIter = 0;
        int nMaxIter = 1000;
        while (go && (nIter < nMaxIter)) {
            float[] parameters = generateRandomParameters(sr, mu);

            GEVYFit yfit = calculateChiSqSumAndCurve(parameters[0], parameters[1], mu, weightMethod, yNorm);

            if ((bestFit == null) || (yfit.getChiSqSum() < bestFit.getChiSqSum())) {
                bestFit = yfit;

                //if ((yfit.getChiSqSum() < yErrSquareSum) && (yfit.getChiSqStatistic() < 10)){
                //    go = false;
                //}
            }

            nIter++;
        }

        return bestFit;
    }

    /**
     * fit the curve using the best fit within the given range of parameters.
     *
     * It uses a recursive nDimension grid search of the ranges and returns
     * the best fit.  The recursive algorithm completes either when a better
     * solution is not found after having descended into a finer grid element
     * range or when a maximum number of iterations have occurred.
     *
     * @param weightMethod
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param mu
     * @param yErrSquareSum
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    protected GEVYFit fitCurveKGreaterThanZeroUsingGrid(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax, float mu,
        float yErrSquareSum, float yNorm) throws FailedToConvergeException, IOException {

        float nDimension = 10;

        return fitCurveKGreaterThanZeroUsingGrid(weightMethod, kMin, kMax, sigmaMin, sigmaMax, mu,
            nDimension, yErrSquareSum, yNorm);
    }

    /**
     * fit the curve using the best fit within the given range of parameters.
     *
     * It uses a recursive nDimension grid search of the ranges and returns
     * the best fit.  The recursive algorithm completes either when a better
     * solution is not found after having descended into a finer grid element
     * range or when a maximum number of iterations have occurred.
     *
     * @param weightMethod
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param mu
     * @param yErrSquareSum
     * @param nDimension the number of divisions along one axis
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    protected GEVYFit fitCurveKGreaterThanZeroUsingGrid(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax, float mu,
        float nDimension, float yErrSquareSum, float yNorm) throws FailedToConvergeException, IOException {

        GEVYFit bestFit = null;

        float[] parameters = new float[2];

        float[] minMax = new float[4];

        boolean go = true;
        int nIter = 0;

        //NOTE:  this method doesn't always find the minima on the right side of the bisection
        //   of an area, but higher dimension helps that (==> implies finds closest center of
        //   minimum biased by "descent"...).
        // This might be a good method to start with a very large range of k and sigma,
        //   run once with nDimension = 10 and then use a faster method
        //   such as the downhill simplex afterwards

        int nMaxIter = 1000;
        while (go && (nIter < nMaxIter)) {

            GEVYFit currentBestFit = null;
            int currentBestFitNumber = 0;

            for (int i = 0; i < nDimension*nDimension; i++) {

                getSectionParametersFromGrid(kMin, kMax, sigmaMin, sigmaMax, i, parameters, nDimension);

                GEVYFit yfit = calculateChiSqSumAndCurve(parameters[0], parameters[1], mu, weightMethod, yNorm);

/*
                float[] tmp = new float[4];
                getSectionBoundariesFromGrid(kMin, kMax, sigmaMin, sigmaMax, i, tmp, nDimension);
                String dstr = String.format("[%d] next krange=[%.7f : %.7f] srange=[%.7f : %.7f] chiSqSum=%.1f k=%.7f s=%.7f",
                    i, tmp[0], tmp[1], tmp[2], tmp[3], yfit.getChiSqSum(), parameters[0], parameters[1]);
                System.out.println(dstr);
*/

                if ( (currentBestFit == null) || (yfit.getChiSqSum() < currentBestFit.getChiSqSum()) ) {
                    currentBestFit = yfit;
                    currentBestFitNumber = i;
                }
            }

            if ((bestFit == null) || currentBestFit.getChiSqSum() < bestFit.getChiSqSum()) {

                bestFit = currentBestFit;

                /*String dstr = String.format("===> choosing [%d] with chiSqSum=%.1f yerrsqsum=%.1f",
                    currentBestFitNumber, currentBestFit.getChiSqSum(), currentBestFit.yDataErrSq);
                System.out.println(dstr);*/

                getSectionBoundariesFromGrid(kMin, kMax, sigmaMin, sigmaMax, currentBestFitNumber, minMax, nDimension);

                // with complexity of GEV, seems to be a problem finding the correct side of the smallest subdivisions, so
                //   will take a less efficient approach of defining next bisection range as
                //   the kMin and current k and same for sigma
                kMin = minMax[0];
                kMax = minMax[1];
                sigmaMin = minMax[2];
                sigmaMax = minMax[3];

            } else {
                if (bestFit != null) {
                    bestFit.setChiSqStatistic(calculateChiSquareStatistic(bestFit.getYFit(), weightMethod));
                }
                return bestFit;
            }

            nIter++;
        }

        if (bestFit != null) {
            bestFit.setChiSqStatistic(calculateChiSquareStatistic(bestFit.getYFit(), weightMethod));
        }
        return bestFit;
    }

    protected float[] findSmallerRange(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax, float mu,
        float yErrSquareSum, float yNorm) throws FailedToConvergeException, IOException {

        float[] minMaxRanges = fitCurveKGreaterThanZeroUsingGridOuter(weightMethod, kMin, kMax, sigmaMin, sigmaMax, mu,
            10, yErrSquareSum, yNorm);

        return minMaxRanges;
    }

    /**
     * method used by findSmallerRange to reduce the k and sigma min and max ranges.
     *
     * It uses a recursive nDimension grid search of the ranges, and keeps the
     * outer range as the solution for each iteration.  It uses a small number
     * of iterations in an attempt to reduce the ranges, not solve for the best fit.
     *
     * The reduction instead of solution is necessary for some GEV curves over
     * a large range of values because the solution f(k, sigma, mu) is not unique,
     * and hence the local minimum found first using a rough division of parameter space
     * may not lead to the best fit.
     *
     * @param weightMethod
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param mu
     * @param nDimension
     * @param yErrSquareSum
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    protected float[] fitCurveKGreaterThanZeroUsingGridOuter(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax, float mu,
        int nDimension, float yErrSquareSum, float yNorm) throws FailedToConvergeException, IOException {

        float origKMaxHalf = kMax/2.0f;

        GEVYFit bestFit = null;

        float[] parameters = new float[2];

        float[] minMax = new float[4];

        boolean go = true;
        int nLoIter = 0;
        int nHiIter = 0;

        int nMaxIter = 3;
        while (go && (nLoIter < nMaxIter) && (nHiIter < 10) && (kMax > origKMaxHalf)) {

            GEVYFit currentBestFit = null;
            int currentBestFitNumber = 0;

            for (int i = 0; i < nDimension*nDimension; i++) {

                getSectionParametersFromGrid(kMin, kMax, sigmaMin, sigmaMax, i, parameters, nDimension);

                GEVYFit yfit = calculateChiSqSumAndCurve(parameters[0], parameters[1], mu, weightMethod, yNorm);

/*
                float[] tmp = new float[4];
                getSectionBoundariesFromGrid(kMin, kMax, sigmaMin, sigmaMax, i, tmp, nDimension);
                String dstr = String.format("[%d] next krange=[%.7f : %.7f] srange=[%.7f : %.7f] chiSqSum=%.1f k=%.7f s=%.7f",
                    i, tmp[0], tmp[1], tmp[2], tmp[3], yfit.getChiSqSum(), parameters[0], parameters[1]);
                System.out.println(dstr);
*/

                if ( (currentBestFit == null) || (yfit.getChiSqSum() < currentBestFit.getChiSqSum()) ) {
                    currentBestFit = yfit;
                    currentBestFitNumber = i;
                }
            }

            if ((bestFit != null) && !Float.isInfinite(currentBestFit.getChiSqSum()) &&
                currentBestFit.getChiSqSum() == bestFit.getChiSqSum()) {

                //TODO:  handle this case

                String dstr = String.format("  <==> same chisqsum=%.1f yerrsqsum=%.1f",
                    currentBestFit.getChiSqSum(), currentBestFit.getYDataErrSq());
                System.out.println(dstr);


                nLoIter++;

            } else if ((bestFit == null) || currentBestFit.getChiSqSum() < bestFit.getChiSqSum()) {

                bestFit = currentBestFit;
/*
                String dstr = String.format("===> choosing [%d] with chiSqSum=%.1f yerrsqsum=%.1f",
                    currentBestFitNumber, currentBestFit.getChiSqSum(), currentBestFit.yDataErrSq);
                System.out.println(dstr);
*/

                getSectionBoundariesFromGrid(kMin, kMax, sigmaMin, sigmaMax, currentBestFitNumber, minMax, nDimension);

                boolean lowDiv = (currentBestFitNumber > ((nDimension*nDimension)/2.0f));

                if (lowDiv) {
                    kMax = bestFit.getK();
                    sigmaMax = bestFit.getSigma();
                    nLoIter++;
                } else {
                    // solution range is too large.  try to reduce it
                    kMin = 1.1f*kMin;
                    kMax = 0.9f*kMax;
                    sigmaMin = 1.1f*sigmaMin;
                    sigmaMax = 0.9f*sigmaMax;
                    nHiIter++;
                }

            } else {
                return new float[]{kMin, kMax, sigmaMin, sigmaMax};
            }

        }

        return new float[]{kMin, kMax, sigmaMin, sigmaMax};
    }

    /**
     generate randomly values for k and sigma between a default of

     float kMin = 0.00001f;
     float kMax = 0.001f;
     float sigmaMin = 0.025f;
     float sMax = 20.0f*sigmaMin;

     * @param sr
     * @param xmax
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
     * divides the ranges into small steps and uses the simplex method to fit the
     * curve within each of those small steps.  returns the best fit.
     * Note that for the GEV curve the best fit may not be unique, that is
     * there may be another set of k and sigma which result in the same curve.
     *
     * @param weightMethod
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param kStepDelta
     * @param sigmaStepFactor
     * @param yConstIndex
     * @param yErrSquareSum
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    protected GEVYFit fitCurveKGreaterThanZeroUsingSmallSteps(WEIGHTS_DURING_CHISQSUM weightMethod,
        float kMin, float kMax, float sigmaMin, float sigmaMax, float mu,
        float kStepDelta, float sigmaStepFactor,
        float yErrSquareSum, float yNorm) throws FailedToConvergeException, IOException {

        //String str = String.format("**kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f", kMin, kMax, sigmaMin, sigmaMax);
        //System.out.println(str);

        GEVYFit bestFit = null;

        float kMin2 = kMin;
        float kMax2 = kMin2 + kStepDelta;
        while (kMax2 <= kMax) {
            float sigmaMin2 = sigmaMin;
            float sigmaMax2 = sigmaMin2 * sigmaStepFactor;
            while (sigmaMax2 <= sigmaMax) {
                GEVYFit fit = fitCurve(kMin2, kMax2, sigmaMin2, sigmaMax2, mu, yErrSquareSum, weightMethod, yNorm);

                /*String dstr = String.format("krange=[%.7f : %.7f] srange=[%.7f : %.7f] chiSqSum=%.1f k=%.7f s=%.7f",
                    kMin2, kMax2, sigmaMin2, sigmaMax2, fit.getChiSqSum(), fit.getK(), fit.getSigma());
                System.out.println(dstr);*/

                if ((bestFit == null) || (fit.getK() > 0) && (fit.getChiSqSum() < bestFit.getChiSqSum())) {

                    bestFit = fit;
                }
                sigmaMin2 = sigmaMax2;
                sigmaMax2 = sigmaMin2 * sigmaStepFactor;
            }
            kMin2 = kMax2;
            kMax2 = kMin2 + kStepDelta;
            //System.gc();
        }

        if (bestFit != null) {
            bestFit.setChiSqStatistic(calculateChiSquareStatistic(bestFit.getYFit(), weightMethod));
        }
        return bestFit;
    }

    protected void plotFit(GEVYFit yfit, String label) throws IOException {

        float xIntervalHalf = (yfit.getX()[1] - yfit.getX()[0]) / 2.0f;
        float xmin = yfit.getX()[0] - xIntervalHalf;
        float xmax = yfit.getX()[ yfit.getX().length - 1];
        float ymin = 0.0f;
        float ymax = MiscMath.findMax(y);

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);

        try {
            plotter.addPlot(x, y, dx, dy, yfit.getX(), yfit.getYFit(), label);
            plotter.writeFile2();
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }

        if (debug) {
            // print the x and y array points for debugging
            System.out.println("x=");
            for (int i = 0; i < x.length; i++) {
                if (i > 0) {
                    System.out.print(", ");
                }
                System.out.print(x[i] + "f");
            }
            System.out.println("\ny=");
            for (int i = 0; i < y.length; i++) {
                if (i > 0) {
                    System.out.print(", ");
                }
                System.out.print(y[i] + "f");
            }
            System.out.println("");
        }
    }

    public GEVYFit fitCurve(float kMin, float kMax, float sigmaMin, float sigmaMax, float mu, float yErrSquareSum,
        WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {

        return fitCurveForKGreaterThanZeroWithDownhillSimplex(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod);
    }

    public GEVYFit fitCurve(float kMin, float kMax, float sigmaMin, float sigmaMax, float mu, float yErrSquareSum,
        WEIGHTS_DURING_CHISQSUM weightMethod, float yNorm) throws FailedToConvergeException, IOException {

        //return fitCurveKGreaterThanZeroUsingGrid(weightMethod, kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, yNorm);
        return fitCurveForKGreaterThanZeroWithDownhillSimplex(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum,
            weightMethod, yNorm);
    }

    public GEVYFit calculateChiSqSumAndCurve(float k, float sigma, float mu, WEIGHTS_DURING_CHISQSUM wdc) {

        float yNorm = MiscMath.findMax(y);

        return calculateChiSqSumAndCurve(k, sigma, mu, wdc, yNorm);
    }

    public GEVYFit calculateChiSqSumAndCurve(float k, float sigma, float mu, WEIGHTS_DURING_CHISQSUM wdc, float yNorm) {

        float[] yGEV = gev.generateNormalizedCurve(k, sigma, mu, yNorm);

        float chisqsum = calculateChiSquareSum(yGEV, wdc);

        GEVYFit yfit = new GEVYFit();
        yfit.setChiSqSum(chisqsum);
        yfit.setK(k);
        yfit.setSigma(sigma);
        yfit.setMu(mu);
        yfit.setYFit(yGEV);
        yfit.setX(x);
        yfit.setXScale(xScale);
        yfit.setYScale(yScale);
        yfit.setYDataErrSq( calcYErrSquareSum() );

        return yfit;
    }
    /**
     * use downhill simplex method to find the best fitting parameters for the curve.
     *
     * At this point, the termination of the fit is by number of iterations, but
     * the code could be altered to use a termination test.
     *
     * Note that we don't fit for mu because the GEV curves in the two-point correlation
     * tend to be Extreme Value Type I of Type II so the peak is usually mu.
     * That is fit with fewer attempts in the wrapper to this method.
     *
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param mu
     * @param yErrSquareSum
     * @param weightMethod
     * @return
     */
    private GEVYFit fitCurveForKGreaterThanZeroWithDownhillSimplex(float kMin, float kMax, float sigmaMin, float sigmaMax,
        float mu, float yErrSquareSum, WEIGHTS_DURING_CHISQSUM weightMethod) {

        int yMaxIndex = MiscMath.findYMaxIndex(y);

        if (yMaxIndex == -1) {
            // all y's were zero in y
            return null;
        }
        float yNorm = y[yMaxIndex];

        return fitCurveForKGreaterThanZeroWithDownhillSimplex(kMin, kMax, sigmaMin, sigmaMax,
            mu, yErrSquareSum, weightMethod, yNorm);
    }

    /**
     * use downhill simplex method to find the best fitting parameters for the curve.
     *
     * At this point, the termination of the fit is by number of iterations, but
     * the code could be altered to use a termination test.
     *
     * Note that we don't fit for mu because the GEV curves in the two-point correlation
     * tend to be Extreme Value Type I of Type II so the peak is usually mu.
     * That is fit with fewer attempts in the wrapper to this method.
     *
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param mu
     * @param yErrSquareSum
     * @param weightMethod
     * @return
     */
    private GEVYFit fitCurveForKGreaterThanZeroWithDownhillSimplex(float kMin, float kMax, float sigmaMin, float sigmaMax,
        float mu, float yErrSquareSum, WEIGHTS_DURING_CHISQSUM weightMethod, float yNorm) {

        int nK = 100;
        float deltaK = (kMax - kMin)/(float)nK;

        int nSigma = 100;
        float deltaSigma = (sigmaMax - sigmaMin)/(float)nSigma;

        float k = kMin;
        float sigma = sigmaMin;

        // start with simplex for 3 points (fitting 2 parameters)
        GEVYFit[] yfits = new GEVYFit[3];
        yfits[0] = calculateChiSqSumAndCurve(k, sigma, mu, weightMethod, yNorm);

        yfits[1] = calculateChiSqSumAndCurve(k + deltaK, sigma + deltaSigma, mu, weightMethod, yNorm);

        yfits[2] = calculateChiSqSumAndCurve(k + 2*deltaK, sigma + 2*deltaSigma, mu, weightMethod, yNorm);

        /*  random choice starts do not work as well for GEV type I and type II
        SecureRandom sr = new SecureRandom();
        sr.setSeed(System.currentTimeMillis());
        yfits[0] = calculateChiSqSumAndCurve(sr.nextFloat()*(kMax - kMin)/2.f, sr.nextFloat()*(sigmaMax - sigmaMin)/2.f, mu, weightMethod);
        yfits[0].yDataErrSq = yErrSquareSum;

        yfits[1] = calculateChiSqSumAndCurve(sr.nextFloat()*(kMax - kMin)/2.f, sr.nextFloat()*(sigmaMax - sigmaMin)/2.f, mu, weightMethod);
        yfits[1].yDataErrSq = yErrSquareSum;

        yfits[2] = calculateChiSqSumAndCurve(sr.nextFloat()*(kMax - kMin)/2.f, sr.nextFloat()*(sigmaMax - sigmaMin)/2.f, mu, weightMethod);
        yfits[2].yDataErrSq = yErrSquareSum;
        */

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.9f; // 0 < beta < 1

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;

        while (go && (nIter < nMaxIter)) {

            sortFromMinToMax(yfits, 0, 2);

            int minSimplexIndex = 0;
            int maxSimplexIndex = 2;
            int midSimplexIndex = 1;

            // determine center for all points excepting the worse fit
            float ksum = 0.0f;
            float ssum = 0.0f;
            for (int i = 0; i < yfits.length - 1; i++) {
                ksum += yfits[i].getK();
                ssum += yfits[i].getSigma();
            }
            k = ksum / (yfits.length - 1);
            sigma = ssum / (yfits.length - 1);

            // "Reflection"
            float kReflect     = k     - (alpha * (yfits[maxSimplexIndex].getK()     - k));
            float sigmaReflect = sigma - (alpha * (yfits[maxSimplexIndex].getSigma() - sigma));
            GEVYFit yfitReflected = calculateChiSqSumAndCurve(kReflect, sigmaReflect, mu, weightMethod, yNorm);
            yfitReflected.setYDataErrSq(yErrSquareSum);

            if ((yfitReflected.getChiSqSum() < yfits[minSimplexIndex].getChiSqSum())
                && ((kReflect >= kMin) && (kReflect <= kMax) && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax))
            ) {

                // "Expansion"
                float kExpansion =     kReflect     - (gamma * (k     - kReflect));
                float sigmaExpansion = sigmaReflect - (gamma * (sigma - sigmaReflect));
                GEVYFit yfitExpansion = calculateChiSqSumAndCurve(kExpansion, sigmaExpansion, mu, weightMethod, yNorm);
                yfitExpansion.setYDataErrSq(yErrSquareSum);

                if ((yfitExpansion.getChiSqSum() < yfits[minSimplexIndex].getChiSqSum())
                    && ((kExpansion >= kMin) && (kExpansion <= kMax) && (sigmaExpansion >= sigmaMin) && (sigmaExpansion <= sigmaMax))
                ) {

                    yfits[maxSimplexIndex] = yfitExpansion;

                } else {

                    yfits[maxSimplexIndex] = yfitReflected;
                }

            } else if ((yfitReflected.getChiSqSum() > yfits[midSimplexIndex].getChiSqSum())
                && ((kReflect >= kMin) && (kReflect <= kMax) && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax))
            ) {

                if ((yfitReflected.getChiSqSum() <= yfits[maxSimplexIndex].getChiSqSum())
                    && ((kReflect >= kMin) && (kReflect <= kMax) && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax))
                ) {

                    yfits[maxSimplexIndex] = yfitReflected;
                }

                // "Contraction"
                float kContraction =     (beta * yfits[maxSimplexIndex].getK())     + (1 - beta)*k;
                float sigmaContraction = (beta * yfits[maxSimplexIndex].getSigma()) + (1 - beta)*sigma;
                GEVYFit yfitContraction = calculateChiSqSumAndCurve(kContraction, sigmaContraction, mu, weightMethod, yNorm);
                yfitContraction.setYDataErrSq(yErrSquareSum);

                if (yfitContraction.getChiSqSum() > yfits[maxSimplexIndex].getChiSqSum()
                    && ((kContraction >= kMin) && (kContraction <= kMax) && (sigmaContraction >= sigmaMin)
                    && (sigmaContraction <= sigmaMax))
                ) {

                    float ktmp = (yfits[midSimplexIndex].getK() + yfits[minSimplexIndex].getK())/2;
                    float stmp = (yfits[midSimplexIndex].getSigma() + yfits[minSimplexIndex].getSigma())/2;

                    yfits[midSimplexIndex] = calculateChiSqSumAndCurve(ktmp, stmp, mu, weightMethod, yNorm);
                    yfits[midSimplexIndex].setYDataErrSq(yErrSquareSum);

                } else {

                    yfits[maxSimplexIndex] = yfitContraction;
                }

            } else
                if ((kReflect >= kMin) && (kReflect <= kMax) && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax))
            {

                yfits[maxSimplexIndex] = yfitReflected;
            }

            nIter++;

            //if (yfits[minSimplexIndex].getChiSqSum() < yErrSquareSum) {
            //    go = false;
            //}

            if ((k > kMax) || (sigma > sigmaMax)) {
                go = false;
            }
        }

        yfits[0].kSolutionResolution = yfits[0].getK() - yfits[1].getK();
        yfits[0].sigmaSolutionResolution = yfits[0].getSigma() - yfits[1].getSigma();
        if (yfits[0].getYFit() != null) {
            yfits[0].setChiSqStatistic(calculateChiSquareStatistic(yfits[0].getYFit(), weightMethod));
        }

        return yfits[0];
    }

    /**
     * sort the array yfits by ascending chisquare sum using
     * the quick sort algorithm.
     *
     * @param yfits
     * @param p first index of the array yfits to be sorted
     * @param r the last index of the array yfits to be sorted, inclusive
     */
    void sortFromMinToMax(GEVYFit[] yfits, int p, int r) {
        int q = -1;

        if (p < r) {

            q = partition(yfits, p, r);

            sortFromMinToMax(yfits, p, q - 1);

            sortFromMinToMax(yfits, q + 1, r);
        }
    }
    /**
     * the partition function of the GEVYFit quick sort method.
     *
     * @param yfits
     * @param p
     * @param r
     * @return
     */
    int partition(GEVYFit[] yfits, int p, int r) {

        float x = yfits[r].getChiSqSum();

        int i = p - 1;

        for (int j = p; j < r ; j++ ) {
            if (yfits[j].getChiSqSum() <= x) {

                i++;

                GEVYFit swap = yfits[i];
                yfits[i] = yfits[j];
                yfits[j] = swap;
            }
        }

        GEVYFit swap = yfits[i + 1];
        yfits[i + 1] = yfits[r];
        yfits[r] = swap;

        return i + 1;
    }

    /**
     * populate parametersOut with the parameters k and sigma from the grid for position.
     * this is used with the recursive grid search.
     *
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param position
     * @param parametersOut
     * @param nDimension
     */
    void getSectionParametersFromGrid(float kMin, float kMax, float sigmaMin, float sigmaMax, int position,
        float[] parametersOut, float nDimension) {

        /*    col 0
         *     ||
         *     \/
         *    *   |    |    *
         *      0 |  1 | 2     <=== row 0
         *    ---------------
         *        |    |
         *  k   3 |  4 | 5
         *    ---------------
         *        |    |
         *      6 |  7 | 8
         *        |    |
         *    *   |    |    *
         *          s
         */

        float kInterval = (kMax - kMin)/nDimension;
        float sInterval = (sigmaMax - sigmaMin)/nDimension;

        int row = (int)(position/nDimension);// rounds down

        int col = (position % (int)nDimension);

        parametersOut[0] = kMax - (row * kInterval) - (kInterval/2.0f);

        parametersOut[1] = sigmaMin + (col * sInterval) + (sInterval/2.0f);
    }

    /**
     * populate minMaxOut with the grid's finer kMin, kMax, sigmaMin, and sigmaMax for the position
     * in the larger grid whose ranges are kMin, kMax, sigmaMin, and sigmaMax.
     * this is used with the recursive grid search.
     *
     * @param kMin
     * @param kMax
     * @param sigmaMin
     * @param sigmaMax
     * @param position
     * @param minMaxOut
     * @param nDimension
     */
    void getSectionBoundariesFromGrid(float kMin, float kMax, float sigmaMin, float sigmaMax, int position,
        float[] minMaxOut, float nDimension) {

        /*    col 0
         *     ||
         *     \/
         *    *   |    |    *
         *      0 |  1 | 2     <=== row 0
         *    ---------------
         *        |    |
         *  k   3 |  4 | 5
         *    ---------------
         *        |    |
         *      6 |  7 | 8
         *        |    |
         *    *   |    |    *
         *          s
         */

        float kInterval = (kMax - kMin)/nDimension;
        float sInterval = (sigmaMax - sigmaMin)/nDimension;

        int row = (int)(position/nDimension);// rounds down

        int col = (position % (int)nDimension);

        minMaxOut[0] = kMax - ((row + 1)*kInterval);
        minMaxOut[1] = kMax - (row*kInterval);

        minMaxOut[2] = (sigmaMin + (col)*sInterval);
        minMaxOut[3] = (sigmaMin + (col + 1)*sInterval);
    }

    protected boolean equals(float[] bestFitMinMax, float[] minMax) {
        if ((bestFitMinMax == null) || (minMax == null)) {
            return false;
        }
        for (int i = 0; i < bestFitMinMax.length; i++) {
            if (bestFitMinMax[i] != minMax[i]) {
                return false;
            }
        }
        return true;
    }

    public float getYScale() {
        return yScale;
    }
}

