package algorithms.curves;

import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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

    public static final int downhillSimplexStartDivisionsDefault = 1000;
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
      Fit the GEV curve iterating over a set of k, sigma, and mu parameters
      and fitting each set using a downhill simplex algorithm.
      The best fit is returned.   
     </pre>
     
     * @param weightMethod the method for determining the tolerance of fit to a 
     * point.  the errors weight is the best choice if errors are available.
     * @return yFit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    @Override
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) 
        throws FailedToConvergeException, IOException {
        
        String filePath = ResourceFinder.findFileInTestResources("sim_curve_params_01.txt");
        
        File f = new File(filePath);
        
        BufferedReader reader = null;
        FileReader fr = null;
        
        List<Float> kParams = new ArrayList<Float>();
        List<Float> sParams = new ArrayList<Float>();
        List<Float> mParams = new ArrayList<Float>();
        
        try {
            fr = new FileReader(f);
            reader = new BufferedReader(fr);
            String line = reader.readLine();
            line = reader.readLine();
            while (line != null) {
                
                String[] params = line.split("\\s+");
                
                kParams.add(Float.valueOf(params[0]));
                sParams.add(Float.valueOf(params[1]));
                mParams.add(Float.valueOf(params[2]));
                                
                line = reader.readLine();
            }
            
        } finally {
            if (fr != null) {
                fr.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
        
        int n = kParams.size();
        int bestFitIndex = -1;
        GEVYFit bestFit = null;
                
        for (int i = 0; i < n - 1; i++) {
            
            float kMin = kParams.get(i).floatValue();
            float kMax = kParams.get(i + 1).floatValue();
            float sMin = sParams.get(i).floatValue();
            float sMax = sParams.get(i + 1).floatValue();
            float mMin = mParams.get(i).floatValue();
            float mMax = mParams.get(i + 1).floatValue();
            
            if (sMax < sMin) {
                sMax = sMin;
            }
            if (mMax < mMin) {
                mMax = mMin;
            }
                        
            GEVYFit yFit = fitCurveKGreaterThanZeroAndMu(weightMethod,
                kMin, kMax, sMin, sMax, mMin, mMax);
            
            if (yFit != null) {
                if (bestFit == null) {
                    bestFit = yFit;
                    bestFitIndex = i;
                } else if (yFit.chiSqStatistic < bestFit.chiSqStatistic) {
                    bestFit = yFit;
                    bestFitIndex = i;
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
     * @param muMin
     * @param muMax
     * @return yfit the best fitting curve
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZeroAndMu(
        WEIGHTS_DURING_CHISQSUM weightMethod, float kMin, float kMax, 
        float sigmaMin, float sigmaMax, float muMin, float muMax) 
        throws FailedToConvergeException, IOException {

        // using mu range useful for EV Type I and II for the two-point correlation functions

        if (debug) {
            String str = String.format(
                "*kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f muMin=%.7f muMax=%.7f", 
                kMin, kMax, sigmaMin, sigmaMax, muMin, muMax);
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

        GEVYFit yfit = fitCurveForKGreaterThanZeroWithDownhillSimplex(kMin, 
            kMax, sigmaMin, sigmaMax, muMin, muMax, yErrSquareSum, 
            weightMethod, yPeak);
        
        if (debug && (yfit != null)) {
            String label = String.format("k=%.1e s=%.1e m=%.1e chisq=%.1f yerrsq=%.1f",
                yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(), yErrSquareSum);
            log.fine(label);
        }

        if ((bestFit == null) || (yfit.getChiSqSum() < bestFit.getChiSqSum())) {
            bestFit = yfit;
        }

        return bestFit;
    }

    protected void plotFit(GEVYFit yfit, String label) throws IOException {
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
            xmin, xmax, ymin, ymax);

        try {
            plotter.addPlot(x, y, xe, ye, yfit.getX(), yfit.getYFit(), label);
            plotter.writeFile2();
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }
    }

    public GEVYFit calculateChiSqSumAndCurve(float k, float sigma, float mu, 
        WEIGHTS_DURING_CHISQSUM wdc, float yNorm) {

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
    private GEVYFit fitCurveForKGreaterThanZeroWithDownhillSimplex(float kMin, 
        float kMax, float sigmaMin, float sigmaMax, float muMin, float muMax, 
        float yErrSquareSum, WEIGHTS_DURING_CHISQSUM weightMethod, 
        float yNorm) {

        int nK = this.downhillSimplexStartDivisions;
        float deltaK = (kMax - kMin)/(float)nK;

        int nSigma = this.downhillSimplexStartDivisions;
        float deltaSigma = (sigmaMax - sigmaMin)/(float)nSigma;
        
        int nMu = this.downhillSimplexStartDivisions;
        float deltaMu = (muMax - muMin)/(float)nMu;

        float k = kMin;
        float sigma = sigmaMin;
        float mu = muMin;

        // start with simplex for 3 points (fitting 2 parameters)
        GEVYFit[] yfits = new GEVYFit[3];
        yfits[0] = calculateChiSqSumAndCurve(k, sigma, mu, weightMethod, yNorm);

        yfits[1] = calculateChiSqSumAndCurve(k + deltaK, sigma + deltaSigma, 
            mu + deltaMu, weightMethod, yNorm);

        yfits[2] = calculateChiSqSumAndCurve(k + 2*deltaK, sigma + 2*deltaSigma, 
            mu + 2*deltaMu, weightMethod, yNorm);

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
            float muReflect = mu - (alpha * (yfits[maxSimplexIndex].getMu() - mu));
            GEVYFit yfitReflected = calculateChiSqSumAndCurve(kReflect, 
                sigmaReflect, muReflect, weightMethod, yNorm);
            yfitReflected.setYDataErrSq(yErrSquareSum);

            if ((yfitReflected.getChiSqSum() < yfits[minSimplexIndex].getChiSqSum())
                && 
                ((kReflect >= kMin) && (kReflect <= kMax) 
                && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax)
                && (muReflect >= muMin) && (muReflect <= muMax))
            ) {

                // "Expansion"
                float kExpansion =     kReflect     - (gamma * (k     - kReflect));
                float sigmaExpansion = sigmaReflect - (gamma * (sigma - sigmaReflect));
                float muExpansion = muReflect - (gamma * (mu - muReflect));
                GEVYFit yfitExpansion = calculateChiSqSumAndCurve(kExpansion, 
                    sigmaExpansion, muExpansion, weightMethod, yNorm);
                yfitExpansion.setYDataErrSq(yErrSquareSum);

                if ((yfitExpansion.getChiSqSum() < yfits[minSimplexIndex].getChiSqSum())
                    && (
                    (kExpansion >= kMin) && (kExpansion <= kMax) 
                    && (sigmaExpansion >= sigmaMin) && (sigmaExpansion <= sigmaMax)
                    && (muExpansion >= muMin) && (muExpansion <= muMax)
                    )
                ) {

                    yfits[maxSimplexIndex] = yfitExpansion;

                } else {

                    yfits[maxSimplexIndex] = yfitReflected;
                }

            } else if ((yfitReflected.getChiSqSum() > yfits[midSimplexIndex].getChiSqSum())
                && (
                (kReflect >= kMin) && (kReflect <= kMax) 
                && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax)
                && (muReflect >= muMin) && (muReflect <= muMax)
                )
            ) {

                if ((yfitReflected.getChiSqSum() <= yfits[maxSimplexIndex].getChiSqSum())
                    && (
                    (kReflect >= kMin) && (kReflect <= kMax) 
                    && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax)
                    && (muReflect >= muMin) && (muReflect <= muMax)
                    )
                ) {

                    yfits[maxSimplexIndex] = yfitReflected;
                }

                // "Contraction"
                float kContraction =     (beta * yfits[maxSimplexIndex].getK())     + (1 - beta)*k;
                float sigmaContraction = (beta * yfits[maxSimplexIndex].getSigma()) + (1 - beta)*sigma;
                float muContraction = (beta * yfits[maxSimplexIndex].getMu()) + (1 - beta)*mu;
                GEVYFit yfitContraction = calculateChiSqSumAndCurve(
                    kContraction, sigmaContraction, muContraction, weightMethod, yNorm);
                yfitContraction.setYDataErrSq(yErrSquareSum);

                if (yfitContraction.getChiSqSum() > yfits[maxSimplexIndex].getChiSqSum()
                    && (
                    (kContraction >= kMin) && (kContraction <= kMax) 
                    && (sigmaContraction >= sigmaMin) && (sigmaContraction <= sigmaMax)
                    && (muContraction >= muMin) && (muContraction <= muMax)
                    )
                ) {

                    float ktmp = (yfits[midSimplexIndex].getK() + yfits[minSimplexIndex].getK())/2;
                    float stmp = (yfits[midSimplexIndex].getSigma() + yfits[minSimplexIndex].getSigma())/2;
                    float mtmp = (yfits[midSimplexIndex].getMu() 
                        + yfits[minSimplexIndex].getMu())/2;

                    yfits[midSimplexIndex] = calculateChiSqSumAndCurve(
                        ktmp, stmp, mtmp, weightMethod, yNorm);
                    yfits[midSimplexIndex].setYDataErrSq(yErrSquareSum);

                } else {

                    yfits[maxSimplexIndex] = yfitContraction;
                }

            } else if (
                (kReflect >= kMin) && (kReflect <= kMax) 
                && (sigmaReflect >= sigmaMin) && (sigmaReflect <= sigmaMax)
                && (muReflect >= muMin) && (muReflect <= muMax)
                )
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
            yfits[0].setChiSqStatistic(
                calculateChiSquareStatistic(yfits[0].getYFit(), weightMethod));
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

        if (p < r) {

            int q = partition(yfits, p, r);

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
}
