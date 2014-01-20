package algorithms.curves;

import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
  <pre>
  A chi square minimization routine for the Generalized Extreme Value function.
  It has been tailored for use with the two-point correlation algorithm,
  that is for k > 0.
 
  If one knows the range of parameter space for k and sigma, it is better to
  use methods which calculate the fit with that better knowledge.
 
  mu is  the location parameter
  sigma is the scale parameter and is > 0
  k is the shape parameter
 
                           (   (      ( x-mu))-(1/k))
                           (-1*(1 + k*(-----))      )
                  1        (   (      (sigma))      )   (      ( x-mu))(-1-(1/k))
  y = y_const * ----- * exp                           * (1 + k*(-----))
                sigma                                   (      (sigma))
 </pre> 
 
 * @author nichole
 */
public class GEVChiSquareMinimization extends AbstractCurveFitter {

    public static final float kMinDefault = 0.00001f;
    public static final float kMaxDefault = 0.001f;
    public static final float kMinDefault2 = 0.001f;
    public static final float kMaxDefault2 = 0.1f;
    public static final float kMinDefault3 = 0.1f;
    public static final float kMaxDefault3 = 2.0f;
    public static final float sigmaMinDefault = 0.025f;
    public static final float sigmaMinDefault2 = 0.1f;
    public static final float sigmaMaxDefault = 0.5f;

    public static final int downhillSimplexStartDivisionsDefault = 1000;
    public static final int downhillSimplexStartDivisionsDefault2 = 3;
    protected int downhillSimplexStartDivisions = downhillSimplexStartDivisionsDefault;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public enum WEIGHTS_DURING_CHISQSUM {
        ERRORS, INVERSE_Y, MODEL_Y
    }

    public GEVChiSquareMinimization(float[] xPoints, float[] yPoints,
        float[] dXPoints, float[] dYPoints) {

        super(xPoints, yPoints, dXPoints, dYPoints);
    }

    /**
      <pre>
      Fit the curve using a simple constraint for mu and a narrow range of possible values
      for k and sigma.
     
      Mu is usually found to be the maximum of the y array of the data because the
      code is intended to be used on GEV distributions, more specifically, those of
      EV Type I and Type II with k > 0.
     
     
      The initial k range is 0.001f to 10.f
      and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     
      The method uses one iteration of a downhill simplex method for fitting.    
     </pre>
     
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
     <pre>
      Fit the curve using a simple constraint for mu and a narrow range of possible values
      for k and sigma.
     
      Mu is usually found to be the maximum of the y array of the data because the
      code is intended to be used on GEV distributions, more specifically, those of
      EV Type I and Type II with k > 0.
     
      The initial k range is 0.001f to 10.f
      and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     
      The method uses one iteration of a downhill simplex method for fitting.
     </pre>
     
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
        float kMin = kMinDefault;
        float kMax = kMaxDefault;
        float sigmaMin = sigmaMinDefault;
        float sigmaMax = sigmaMaxDefault;

        if (debug) {
            String str = String.format("kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f", kMin, kMax, sigmaMin, sigmaMax);
            log.fine(str);
        }

        float yErrSquareSum = calcYErrSquareSum();

        GEVYFit bestFit;

        bestFit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

        if (debug && (bestFit != null)) {

            String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum,
                bestFit.getChiSqStatistic());

            log.fine(str);
        }

        return bestFit;
    }

    /**
     <pre>
      Fit the curve using a simple constraint for mu and a narrow range of possible values
      for k and sigma.
     
      Mu is usually found to be the maximum of the y array of the data because the
      code is intended to be used on GEV distributions, more specifically, those of
      EV Type I and Type II with k > 0.
     
      The initial k range is 0.001f to 10.f
      and the initial sigma range is kMin * x[0] to kMax * x[x.length - 1];
     
      The method uses one iteration of a downhill simplex method for fitting.
     </pre>
     
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @param mu
     * @param yNorm
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZeroWithExtendedRanges(WEIGHTS_DURING_CHISQSUM weightMethod, float mu, float yNorm)
        throws FailedToConvergeException, IOException {

        // start with ranges useful for EV Type I and II for the two-point correlation functions
        float kMin = kMinDefault;
        float kMax = kMaxDefault;
        float sigmaMin = sigmaMinDefault;
        float sigmaMax = sigmaMaxDefault;

        float chiSqStatisticLimit = 0f;//1000f;

        if (debug) {
            String str = String.format("kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f", kMin, kMax, sigmaMin, sigmaMax);
            log.fine(str);
        }

        float yErrSquareSum = calcYErrSquareSum();

        GEVYFit bestFit;

        bestFit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

        if (debug && (bestFit != null)) {

            String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), yErrSquareSum,
                bestFit.getChiSqStatistic());

            log.fine(str);
        }

        if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

            GEVYFit yfit;

            yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

            if (debug && (yfit != null)) {

                String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum,
                    yfit.getChiSqStatistic());

                log.fine(str);
            }

            if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                bestFit = yfit;
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault2;
                kMax = kMaxDefault2;

                yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum,
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault3;
                kMax = kMaxDefault3;

                yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum,
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault3;
                kMax = kMaxDefault3;// might need kMinDefault3*10

                sigmaMin = sigmaMinDefault2;

                yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum,
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault3;
                kMax = 1.0f;

                sigmaMin = sigmaMinDefault;
                sigmaMax = 0.1f;

                yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yNorm);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f yerrsq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum,
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }
        }

        return bestFit;
    }

     /**
     <pre>
      Fit curve within range of given parameters.  Internally, uses the Neder-Meade downhill
      simplex method and ranges of k, sigma and mu.  If the best fitting chi-square statistic
      is larger than 100, it repeats the process using different start conditions for
      the simplex method and keeps the best result.
     </pre>
     
     * @param weightMethod the method for determining the tolerance of fit to a point.  the errors weight is
     *    the best choice if errors are available.
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZeroAndMu(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {

        GEVYFit bestFit = null;

        float yNorm = MiscMath.findMax(y);

        float yErrSquareSum = calcYErrSquareSum();


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

        float nDmu = 2;
        for (int i = 0; i <= end; i++) {

            //float yNorm = y[i];

            for (int ii = 0; ii < nDmu; ii++) {

                float dmu;
                if (i < (y.length - 1)) {
                    dmu = (x[i+1] - x[i])/nDmu;
                } else {
                    dmu = (x[1] - x[0])/nDmu;
                }
                mu = x[i] + ii*dmu;

                //GEVYFit yfit = fitCurveKGreaterThanZero(weightMethod, mu, yNorm);
                GEVYFit yfit = fitCurveKGreaterThanZeroWithExtendedRanges(weightMethod, mu, yNorm);

                if (debug && (yfit != null)) {

                    String label = String.format("k=%.2e s=%.2e m=%.2e chisq=%.2f yerrsq=%.2f",
                        yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum);

                    plotFit(yfit, label);

                    log.fine(label);
                }

                if ( (bestFit == null) || (yfit.getChiSqSum() < bestFit.getChiSqSum()) ) {
                    bestFit = yfit;
                }
            }
        }

        return bestFit;
    }

    /**
     <pre>
     Fit curve within range of given parameters.  Internally, uses the Neder-Meade downhill
     simplex method and no additional iterations combined with range reduction such as
      is used by the fitCurveKGreaterThanZero method without min and max arguments.
     
      The mu parameter is fit for x[0] through x[yConstIndex] only.
     </pre>
     
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
            log.fine(str);
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

        float nDmu = 2;
        for (int i = 0; i <= end; i++) {

            //float yNorm = y[i];

            for (int ii = 0; ii < nDmu; ii++) {

                float dmu;
                if (i < (y.length - 1)) {
                    dmu = (x[i+1] - x[i])/nDmu;
                } else {
                    dmu = (x[1] - x[0])/nDmu;
                }
                mu = x[i] + ii*dmu;

                GEVYFit yfit = fitCurve(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, weightMethod, yPeak);

                if (debug && (yfit != null)) {
                    String label = String.format("k=%.1e s=%.1e m=%.1e chisq=%.1f yerrsq=%.1f",
                        yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum);
                    log.fine(label);
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
                log.fine(label);
            }
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
            plotter.addPlot(x, y, xe, ye, yfit.getX(), yfit.getYFit(), label);
            plotter.writeFile2();
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }

        /*
        if (debug) {
            // print the x and y array points for debugging
            log.info("x=");
            for (int i = 0; i < x.length; i++) {
                if (i > 0) {
                    System.out.print(", ");
                }
                System.out.print(x[i] + "f");
            }
            log.info("\ny=");
            for (int i = 0; i < y.length; i++) {
                if (i > 0) {
                    System.out.print(", ");
                }
                System.out.print(y[i] + "f");
            }
            log.info("");
        }*/
    }

    public GEVYFit fitCurve(float kMin, float kMax, float sigmaMin, float sigmaMax, float mu, float yErrSquareSum,
        WEIGHTS_DURING_CHISQSUM weightMethod, float yNorm) throws FailedToConvergeException, IOException {

        //return fitCurveKGreaterThanZeroUsingGrid(weightMethod, kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum, yNorm);
        return fitCurveForKGreaterThanZeroWithDownhillSimplex(kMin, kMax, sigmaMin, sigmaMax, mu, yErrSquareSum,
            weightMethod, yNorm);
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
        float mu, float yErrSquareSum, WEIGHTS_DURING_CHISQSUM weightMethod, float yNorm) {

        int nK = this.downhillSimplexStartDivisions;
        float deltaK = (kMax - kMin)/(float)nK;

        int nSigma = this.downhillSimplexStartDivisions;
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
        yfits[0].muSolutionResolution = yfits[0].getMu() - yfits[1].getMu();
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

        float xxp = yfits[r].getChiSqSum();

        int i = p - 1;

        for (int j = p; j < r ; j++ ) {
            if (yfits[j].getChiSqSum() <= xxp) {

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

}

