package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.XY;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.curves.FailedToConvergeException;
import algorithms.curves.GEVChiSquareMinimization;
import algorithms.curves.GEVYFit;
import algorithms.curves.GeneralizedExtremeValue;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Class to estimate a background density for a set of points in which
 * the background surface density will be used by the calling program to find
 * clusters in the data.
 *
 * It calculates the two-point density function of rectangular voids, creates a
 * histogram from the distribution, fits a Generalized Extreme Value
 * distribution to the histogram, integrates the area under the curve to
 * estimate the background surface density. For denser background densities, a
 * slightly different integration limit is used.
 *
 * This assumes that the data contains background points and clustered points
 * whose 2 distributions are different from one another, but homogeneous within
 * their own. It assumes that the distributions do not need to be fit well, but
 * that only the rough range of the background surface density needs to be
 * learned.
 *
 <pre>
 * More specifically:
 * -- The location of the 'background' points in two dimensional space are likely
 *    Poisson, that is their locations in a fixed interval of space are
 *    independent of one another and occurred randomly.
 * -- The areas between the smallest voids in such a distribution are well fit by
 *    Generalized Extreme Value distributions. Extreme value distributions are used
 *    to describe the maximum or minimum of values drawn from a sample distribution
 *    that is essentially exponential.
 *    The fits improve as N, the number of data points, increases.
 * -- The GEV curve contains 3 independent fitting parameters and the curve is
 *    an exponential combined with a polynomial, so it's resulting fitted
 *    parameters are not unique, but the curve is useful for characterizing the
 *    background point distribution and then integrating under the curve.
 * -- The points within a cluster are distributed radially.
 * </pre>
 *
 *
 * Usage:
 *    TwoPointVoidStats stats = new TwoPointVoidStats(indexer);
 *    stats.calc();
 *
 * For a more detailed fit to the background, one can use TwoPointVoidStats:
 *     stats = new TwoPointVoidStats(indexer);
 *     stats.setUseCompleteSampling(true);
 *     stats.calc();
 *
 * The runtime for the default incomplete sampling method is N^2 - N
 * while the complete sampling method is N!/(2!(N-2)! * N!/(2!(N-2) which
 * is effectively N^4.  Note that for datasets smaller than 100 points,
 * the complete sampling is used even if the default was requested.
 *
 * The value returned by the code is used by TwoPointCorrelation as an error of
 * the background density which is then used to calculate a threshold to find
 * clusters.
 *
 * If debugging is turned on, a plot is generated and the path is printed to
 * standard out, and statements are printed to standard out.
 *
 * @author nichole
 */
public class TwoPointVoidStats extends AbstractPointBackgroundStats {

    // TODO:  REPLACE these with enums. if !useCompleteSampling && !useLeastCompleteSampling
    //            need a way to read structure of data to decide between
    //            range search and completeSamplingAlgWtihIncompleteSampling

    public enum Sampling {
        COMPLETE, LEAST_COMPLETE, SEMI_COMPLETE, SEMI_COMPLETE_RANGE_SEARCH /*default*/, SEMI_COMPLETE_RANGE_SEARCH_2
    }

    public enum State {
        POINTS_LOADED, DENSITIES_CALCULATED, HISTOGRAM_CREATED, HISTOGRAM_FITTED, STATS_FINALIZED
    }

    protected Sampling sampling = null;

    protected State state = null;

    protected float[] allTwoPointSurfaceDensities = null;
    protected float[] allTwoPointSurfaceDensitiesErrors = null;
    protected int nTwoPointSurfaceDensities = 0;
    protected int[] point1 = null;
    protected int[] point2 = null;
    protected StringArrayLite twoPointIdentities = null;

    protected int defaultNBins = 40;

    //for debugging, hold on to intermediate data
    protected HistogramHolder statsHistogram = null;
    protected GEVYFit bestFit = null;

    protected boolean adjustMuForDensity = true;

    protected boolean doLogPerformanceMetrics = false;

    //
    public HistogramHolder getStatsHistogram() {
        return statsHistogram;
    }
    public GEVYFit getBestFit() {
        return bestFit;
    }

    protected float[] gevRangeFittingParameters = null;

    /**
     * constructor for class
     *
     * @param indexedSortedPoints indexed points sorted by Y
     */
    public TwoPointVoidStats(DoubleAxisIndexer indexedSortedPoints) {
        super(indexedSortedPoints);

        state = State.POINTS_LOADED;
    }

    public TwoPointVoidStats(String persistedIndexerFilePath) throws IOException {
        super(persistedIndexerFilePath);

        state = State.POINTS_LOADED;
    }

    /**
     * Use the least complete sampling.  This is fine for datasets with many many outliers,
     * especially compared to the number within clusters.  It's the fastest, but
     * should not be used for datasets in which there are few to no-outliers as there
     * will be very few two-point voids that get sampled.
     * The runtime is O(N*log_2(N)) as it uses divide and conquer.
     */
    public void setUseLeastCompleteSampling() {
        this.sampling = Sampling.LEAST_COMPLETE;
    }

    /**
     * Use complete sampling.
     * The runtime for the complete sampling is N!/(2!(N-2)! * N!/(2!(N-2) which
     * is effectively N^4.
     *
     * Note: complete sampling is always used when there are less than 100 points.
     *
     */
    public void setUseCompleteSampling() {
        this.sampling = Sampling.COMPLETE;
    }

    /**
     * Use semi-complete range sampling.  This is good to use when the number of background
     * points, that is the number of outliers, is expected to be moderately dense or dense.
     * It performs poorly for datasets in which all points are expected to be in
     * clusters.
     * The runtime is O(n^2 - n).
     */
    public void setUseSemiCompleteRangeSampling() {
        this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH;
    }

    /**
     * Use the algorithm for Sampling.COMPLETE, but sample the 2nd index over a larger range to reduce the
     * runtime.  This algorithm is very good to use when Sampling.COMPLETE is needed,
     * but there are a large enough number of points to sample the background very
     * well even when we decrease the sampling by a large amount.  This works well
     * for datasets which have few to no background points, that is few to no outliers.
     *
     * The runtime is expected to be
     *
     */
    public void setUseSemiCompleteSampling() {
        this.sampling = Sampling.SEMI_COMPLETE;
    }

    public Sampling getSampling() {
        return Sampling.valueOf(sampling.name());
    }

    protected void logPerformanceMetrics() {
        this.doLogPerformanceMetrics = true;
    }

    protected void printPerformanceMetrics(long startTimeMillis, long stopTimeMillis, String methodName, String bigOh) {

        long diffSec = (stopTimeMillis - startTimeMillis)/1000;

        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        MemoryUsage heapUsage = mbean.getHeapMemoryUsage();
        MemoryUsage nonHeapUsage = mbean.getNonHeapMemoryUsage();

        String str = String.format("%35s:  N=%9d  %s  RT(sec)=%8d  instance estimates(bytes)=%9d   heapUsed(bytes)=%9d   memoryPoolsSum(bytes)=%9d",
            methodName,
            indexer.nXY, bigOh, diffSec, approximateMemoryUsed(),
            heapUsage.getUsed(), nonHeapUsage.getUsed() );

        Logger.getLogger(this.getClass().getSimpleName()).info(str);
    }

    // TODO:  replace with estimation using reflection one day
    public long approximateMemoryUsed() {

        int n = (allTwoPointSurfaceDensities == null) ? 0 : allTwoPointSurfaceDensities.length;

        long sumBytes = 4*16 + (3*8);

        if (statsHistogram != null) {
            sumBytes += statsHistogram.approximateMemoryUsed();
        }
        if (bestFit != null) {
            sumBytes += bestFit.approximateMemoryUsed();
        }
        if (twoPointIdentities != null) {
            sumBytes += twoPointIdentities.approximateMemoryUsed();
        }

        sumBytes += 2*(8 + 16);//enum reference + class

        long sumBits = (4*n*32) + 32 + 32 + 2 + (Sampling.values().length * (32))
            + (State.values().length * (32));

        sumBytes += (sumBits/8);

        // amount of padding needed to make it a round 8 bytes
        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    public void setGEVRangeParameters(float kMin, float kMax, float sigmaMin, float sigmaMax) {
        this.gevRangeFittingParameters = new float[]{kMin, kMax, sigmaMin, sigmaMax};
    }

    /**
     * calculate the 2-point void densities and then calculate the
     * statistics of those points.
     *
     * More specifically:
     * @see #findVoids()
     *
     * The sampled background is then fit with a GEV curve, and the peak centroid
     * is the background density estimate.
     *
     * @throws TwoPointVoidStatsException
     */
    public void calc() throws TwoPointVoidStatsException {

        calculateTwoPointVoidDensities();

        twoPointIdentities = null;

        calculateStats();
    }

    @Override
    protected void calculateStats() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        if (indexer.nXY > 999) {
            statsHistogram = createHistogramWithHigherPeakResolution();
        } else {
            statsHistogram = createHistogram();
        }

        state = State.HISTOGRAM_CREATED;

        //int yMaxBin = Histogram.findMax(histogram.getYHist());

        int yMaxBin = -1;
        float ymax = Float.MIN_VALUE;
        for (int i = 0; i < (statsHistogram.getYHist().length/2.); i++) {
            if (statsHistogram.getYHist()[i] > ymax) {
                ymax = statsHistogram.getYHist()[i];
                yMaxBin = i;
            }
        }

        calculateStatsForBackground(statsHistogram, yMaxBin);
    }

    /**
     * Sample the two-point voids to calculate their density.
     *
     * More specifically:
     * @see #findVoids()
     *
     * @throws TwoPointVoidStatsException
     */
    protected void calculateTwoPointVoidDensities() throws TwoPointVoidStatsException {

        if (state.ordinal() >= State.DENSITIES_CALCULATED.ordinal()) {
            return;
        }

        long startTimeMillis = System.currentTimeMillis();

        allTwoPointSurfaceDensities = new float[100];
        point1 = new int[100];
        point2 = new int[100];
        twoPointIdentities = new StringArrayLite(100);

        findVoids();

        if (nTwoPointSurfaceDensities == 0) {
            throw new TwoPointVoidStatsException("No pairs were found isolated within an area");
        }

        // condense arrays
        allTwoPointSurfaceDensities = Arrays.copyOf(allTwoPointSurfaceDensities, nTwoPointSurfaceDensities);
        point1 = Arrays.copyOf(point1, nTwoPointSurfaceDensities);
        point2 = Arrays.copyOf(point2, nTwoPointSurfaceDensities);

        allTwoPointSurfaceDensitiesErrors =
            calulateTwoPointDensityErrors(allTwoPointSurfaceDensities, point1, point2,
            indexer.getX(), indexer.getY(), indexer.getXErrors(), indexer.getYErrors());

        state = State.DENSITIES_CALCULATED;

        if (doLogPerformanceMetrics) {

            long stopTimeMillis = System.currentTimeMillis();

            String str = "";
            if (sampling.ordinal() == Sampling.LEAST_COMPLETE.ordinal()) {
                str = "O(n lg(n)) with n=";
            } else if (sampling.ordinal() == Sampling.COMPLETE.ordinal()) {
                str = "O(n^4) with n=";
            } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH.ordinal()) {
                str = "O(n^2 - n) with n=";
            } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE.ordinal()) {
                str = "not yet est with n=";
            }
            if (allTwoPointSurfaceDensities == null) {
                str = str + "0";
            } else {
                str = str + this.allTwoPointSurfaceDensities.length;
            }
            printPerformanceMetrics(startTimeMillis, stopTimeMillis, "calculateSurfaceDensities", str);
        }
    }

    /**
     * Calculate the 2 point density errors following the chain rule
     *
     *                                | df |^2               | df |^2         df   df
     *      (sigma_f)^2 =  (sigma_x)^2|----|   +  (sigma_y)^2|----|    +  2 * -- * -- * cov_ab
     *                                | dx |                 | dy |           dx   dy
     *
     *      For uncorrelated variables the covariance terms are zero.
     *
     * For two-point density:
     *                 N             dF      -N
     *      F(x,y) = ------   ==>   ---- =  -----
     *                 X             dx      X^2
     *
     *                                  |   -N   |^2
     *     (sigma_f)^2 =  (sigma_x)^2 * |--------|
     *                                  |   X^2  |
     *
     * @param densities - the two point densities
     * @param point1Indexes indexes to xp, yp of one of the 2 points in the two-point densities
     * @param point2Indexes indexes to xp, yp of the other of the 2 points in the two-point densities
     * @param xp x array of points
     * @param yp y array of points
     * @param xpe errors in x
     * @param ype errors in y
     * @return errors for densities
     */
    float[] calulateTwoPointDensityErrors(float[] densities, int[] point1Indexes, int[] point2Indexes,
        float[] xp, float[] yp, float[] xpe, float[] ype) {

        if (densities == null) {
            throw new IllegalArgumentException("densities cannot be null");
        }
        if (point1Indexes == null || point2Indexes == null) {
            throw new IllegalArgumentException("pointIndexes cannot be null");
        }
        if (xp == null || yp == null) {
            throw new IllegalArgumentException("xp and yp cannot be null");
        }
        if (xpe == null || ype == null) {
            throw new IllegalArgumentException("errors are required");
        }

        float[] densityErrors = new float[densities.length];

        for (int i = 0; i < densities.length; i++) {

            if (Float.isInfinite(densities[i])) {
                continue;
            }

            int pt1 = point1[i];
            int pt2 = point2[i];

            double xe1 = xpe[pt1] ; // xp[pt1];
            double xe2 = xpe[pt2] ; // xp[pt2];
            double ye1 = ype[pt1] ; // yp[pt1];
            double ye2 = ype[pt2] ; // yp[pt2];

            double xsigmasq = ( xe1 * xe1 ) + ( xe2 * xe2 );
            double ysigmasq = ( ye1 * ye1 ) + ( ye2 * ye2 );
            double xysigmasq = xsigmasq + ysigmasq;

            double linearDensity = allTwoPointSurfaceDensities[i];
            double linearDensitySq = Math.pow(linearDensity, 2);

            double err = Math.sqrt(xysigmasq * linearDensitySq);

            densityErrors[i] = (float)err;
        }
        return densityErrors;
    }

    protected HistogramHolder createHistogram() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        int nBins = (indexer.getNumberOfPoints() < 100) ? defaultNBins/2 : defaultNBins;

        //HistogramHolder histogram = Histogram.createHistogramForSkewedData(nBins, allTwoPointSurfaceDensities,
        //    allTwoPointSurfaceDensitiesErrors, false);

        HistogramHolder histogram = Histogram.createHistogramForSkewedData(
            nBins, allTwoPointSurfaceDensities, allTwoPointSurfaceDensitiesErrors, true);

        plotPairSeparations();

        return histogram;
    }

    protected HistogramHolder createHistogramWithHigherPeakResolution() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        int nBins = (indexer.getNumberOfPoints() < 100) ? defaultNBins/2 : defaultNBins;

        float[] minMax = MiscMath.calculateOuterRoundedMinAndMax(allTwoPointSurfaceDensities);

        HistogramHolder histogram = Histogram.createHistogramForSkewedDataForPeakResolution2(
            nBins, allTwoPointSurfaceDensities, allTwoPointSurfaceDensitiesErrors,
            minMax[0], minMax[1], 0);

        //HistogramHolder histogram = Histogram.createHistogramForSkewedDataForPeakResolution(
        //    nBins, allTwoPointSurfaceDensities, allTwoPointSurfaceDensitiesErrors,
        //    minMax[0], minMax[1]);

        plotPairSeparations();

        return histogram;
    }

    /**
     * internal method to calculate the statistics from the histogram using a
     * fit to the background distribution and a rough integration under the fit
     * to approximate a usable limit of the background surface density
     * distribution.
     *
     * @param histogram
     * @param yMaxBin
     * @throws TwoPointVoidStatsException
     */
    protected void calculateStatsForBackground(HistogramHolder histogram, int yMaxBin) throws TwoPointVoidStatsException {

        GEVYFit yfit = fitBackgroundHistogram(histogram, yMaxBin);

        state = State.HISTOGRAM_FITTED;

        finalizeStats(histogram, yfit);
    }

    protected GEVYFit fitBackgroundHistogram(HistogramHolder histogram, int yMaxBin) throws TwoPointVoidStatsException {

        try {

            GEVChiSquareMinimization chiSqMin = new GEVChiSquareMinimization(
                histogram.getXHist(), histogram.getYHistFloat(), histogram.getXErrors(), histogram.getYErrors());

            chiSqMin.setDebug(debug);

            GEVYFit yfit = null;

            if (gevRangeFittingParameters != null) {

                yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS,
                    gevRangeFittingParameters[0], gevRangeFittingParameters[1],
                    gevRangeFittingParameters[2], gevRangeFittingParameters[3]
                );

            } else {
                yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
            }

            if (yfit == null) {
                if (debug) {
                    plotPairSeparations();
                }
                throw new TwoPointVoidStatsException("histogram of linear densities was not fittable");
            } else {
                if (debug) {
                    System.out.println(yfit.toString());
                }
            }

            boolean didCalculateFitStats = yfit.calculateStatsAsStepFunction();
            if (!didCalculateFitStats) {
                if (debug) {
                    plotPairSeparations();
                }
                throw new TwoPointVoidStatsException("histogram of linear densities was not fittable");
            }

            /*
             * We are assigning the threshold based upon one of 3 different cases:
             *   (1) a sparsely populated background + clusters.  In this case, the distribution is
             *       almost all clusters.  The histogram width will be wider and there
             *       will be little to no tail.
             *       If there are many points, the distribution should approach Gaussian.
             *       The histogram shape is GEV w/ k > 0.
             *
             *   (2) a moderately populated background + clusters.
             *       The distribution will be GEV with k > 0 with small peaks in the tail.
             *
             *   (3) a densely populate background + clusters.
             *       The distribution will be a
             *
             */

            /*boolean increaseThreshhold = (chsqdiverr < 200.f) && (chsqstatdiverr < 100.f);

            if (increaseThreshhold) {
                // for higher background surface densities, should fit further down the distribution
                GEVYFit yfit2 = chiSqMin.fitCurveKGreaterThanZero(
                    GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, (histogram.yHist.length / 2));
                yfit = yfit2;
            }*/

            if (yfit == null) {
                throw new TwoPointVoidStatsException("did not find a solution");
            }

            return yfit;

        } catch (FailedToConvergeException e) {
            throw new TwoPointVoidStatsException(e);
        } catch (IOException e2) {
            throw new TwoPointVoidStatsException(e2);
        }
    }

    protected void finalizeStats(HistogramHolder histogram, GEVYFit yfit) throws TwoPointVoidStatsException {

        this.statsHistogram = histogram;

        this.bestFit = yfit;

        try {

            if (debug) {
                if (bestFit != null) {
                    System.out.println(bestFit.toString());
                }
                float xHalfInterval = (histogram.getXHist()[1] - histogram.getXHist()[0]) / 2.0f;
                float xmin = 0;
                float xmax = histogram.getXHist()[histogram.getXHist().length - 1] + xHalfInterval;
                float ymin = 0;
                float ymax = MiscMath.findMax(histogram.getYHistFloat());

                float[] xf = bestFit.getOriginalScaleX();
                float[] yf = bestFit.getOriginalScaleYFit();
                PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
                plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                    histogram.getXErrors(), histogram.getYErrors(), xf, yf, "");
                plotter.writeFile3();
            }

            String limitStr = "";

            boolean isSmallNumberHist = false;

            // for histograms with small number of points and small peak, need to use
            //   a pattern observed for just those and that is that the minimum before
            //   yPeak is the answer if lower than the first point.
            int yPeakIndex = MiscMath.findYMaxIndex(histogram.getYHist());
            if ((yPeakIndex > 0) && (histogram.getYHist()[yPeakIndex] < 100)) {
                float min = Float.MAX_VALUE;
                int minIndex = -1;
                for (int i = yPeakIndex; i > -1; i--) {
                    float a = histogram.getYHistFloat()[i];
                    if ((a < min) && !Float.isInfinite(a)) {
                        min = a;
                        minIndex = i;
                    } else {
                        break;
                    }
                }
                if (minIndex > 0) {
                    isSmallNumberHist = true;
                    this.backgroundSurfaceDensity = histogram.getXHist()[minIndex];
                }
            }

            float limit, limitError;

            if (!isSmallNumberHist) {

                // centroid of area defined by the top portion of the fit where y >= ypeak/2
                float[] areaAndXYTopCentroid = calculateCentroidOfTop(bestFit.getOriginalScaleX(), bestFit.getOriginalScaleYFit(), 0.5f);

                limitStr = "top centroid";

                limit = (areaAndXYTopCentroid != null) ? areaAndXYTopCentroid[1] : bestFit.getXPeak();

                //if (limitIndex == -1) {
                //    System.out.println("WARNING:  solution was not found");
                //    return;
                //}

                if (limit < 0) {
                    float min = Float.MAX_VALUE;
                    int minIndex = -1;
                    for (int i = yPeakIndex; i < histogram.getXHist().length; i++) {
                        float a = histogram.getYHistFloat()[i];
                        if ((a < min) && !Float.isInfinite(a)) {
                            min = a;
                            minIndex = i;
                        } else {
                            break;
                        }
                    }
                    if (minIndex > 0) {
                        isSmallNumberHist = true;
                        this.backgroundSurfaceDensity = histogram.getXHist()[minIndex];
                    }
                } else {

                    this.backgroundSurfaceDensity = limit;
                }
            }

            if (debug) {
                System.out.println(bestFit.toString());
            }

            int limitIndex = bestFit.getXPeakIndex();
            limitError = histogram.getXErrors()[limitIndex];


            /* Errors add in quadrature:
             *
             *    * errors in histogram
             *    * errors in GEV fitting
             *
             * Errors in histogram:
             *
             *    error in Y is sqrt(Y) and that is already in standard units.
             *    error in X is resolvability, which is bin size = (xHist[1]-xHist[0])/2.
             *
             *                                | df |^2               | df |^2         df   df
             *      (sigma_f)^2 =  (sigma_x)^2|----|   +  (sigma_y)^2|----|    +  2 * -- * -- * cov_ab
             *                                | dx |                 | dy |           dx   dy
             *
             *      For uncorrelated variables the covariance terms are zero.
             *
             * Errors in fitting:
             *
             *    It looks like the error in the fit could very roughly be math.sqrt(chisq)/nPoints, that is close to
             *        the chi square statistic which uses degrees of freedom instead of nPoints.
             *
             *    But, more formally, errors are based on fitting 2 variable parameters and 1 fixed parameter and on x and y.
             *    The 2 parameters are GEV.sigma and GEV.k (we use a fixed mu usually)
             *
             *                                    | df_fit |^2               | df_fit |^2
             *      (sigma_f_fit)^2 =  (sigma_x)^2|--------|   +  (sigma_y)^2|--------|
             *                                    |   dx   |                 |   dy   |
             *
             *                             | df_fit |^2                       |  df_fit  |^2
             *          +  (sigma_GEV.k)^2 |--------|   +  (sigma_GEV.sigma)^2|----------|
             *                             | dGEV.k |                         |dGEV.sigma|
             *
             *                                          (   (      ( x-mu))-(1/k))
             *                                          (-1*(1 + k*(-----))      )
             *                                 1        (   (      (sigma))      )  (      ( x-mu))(-1-(1/k))
             *      where f_fit = y_const * ----- * exp                           * (1 + k*(-----))
             *                               sigma                                  (      (sigma))
             *
             *      The deriv of f_fit w.r.t. sigma, k, and x is a very very long derivative...
             *
             *      If need that formally one day, it's tedious, but possible to work out with the chain rule.
             *      For now, will instead use chisq to approximate the standard deviation from the model.
             */

            // Error from one value is dependent upon having counted all in histogram.  we take the
            //  error for our single bin as an average over all because of that.
            float errorInEstimateFromHistogram = limitError;

            float gevTotalMeanFittingError =
                GeneralizedExtremeValue.calculateTotalMeanFittingError(bestFit.getX(), bestFit.getYFit());

            float errorInFitting = gevTotalMeanFittingError / histogram.getYHist().length;

            // add errors in quadrature
            this.backgroundSurfaceDensityError = (float) (Math.sqrt(
                errorInEstimateFromHistogram * errorInEstimateFromHistogram + errorInFitting * errorInFitting));

            if (debug) {
                System.out.println("estimating background minimum as location of " + limitStr + " of background profile counts = "
                + this.backgroundSurfaceDensity
                + " w/ x error in histogram bin =" + errorInEstimateFromHistogram
                + " gev fitting error for one point =" + errorInFitting);
            }

            if (debug && (bestFit.getChiSqSum() > bestFit.getYDataErrSq())) {
                System.out.println("WARNING:  chisq is larger than errors: "
                    + bestFit.getChiSqSum() + " (errsqsum=" + bestFit.getYDataErrSq() + ")");
            }

            state = State.STATS_FINALIZED;

        } catch (IOException e) {
            throw new TwoPointVoidStatsException(e);
        }
    }

    /**
     * find the rectangular 2 point densities using the method already set.
     * If the method has not yet been set, the following rules are used to
     * determine sampling:
     *
     *     if (indexer.nXY < 100) {
     *         sampling = Sampling.COMPLETE;
     *     } else {
     *         if (data are somewhat evenly distributed, in other words, expecting outliers) {
     *             sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH;
     *         } else {
     *             sampling = Sampling.SEMI_COMPLETE;
     *         }
     *     }
     *
     */
    protected void findVoids() {

        float[] x = indexer.getX();
        float[] y = indexer.getY();

        if (sampling == null) {

            // we need to decide between SEMI_COMPLETE and SEMI_COMPLETE_RANGE_SEARCH.


            // ***** NOTE ********
            // This method is not implemented yet, so until then, will solve first SEMI_COMPLETE
            //    and if many points are outliers or all points are in one cluster, will resolve it
            //    using SEMI_COMPLETE_RANGE_SEARCH.
            // This is unfortunately, logic built into the invokee, TwoPointCorrelation for now.



            /*boolean pointsAreSemiEvenlyDistributed = PointsUtil.pointsAreSemiEvenlyDistributed(indexer.getX(), indexer.getY());
            if (pointsAreSemiEvenlyDistributed) {
                // the rough range search works well with semi-even background points
                sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH;
            } else {
                // use the algorithm for complete sampling because it handles
                //    "no background points" and non-circular clusters, but use
                //    it with much less sampling to allow faster completion
                sampling = Sampling.SEMI_COMPLETE;
            }*/

            if (indexer.nXY > 9999) {
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH_2;
            } else if (indexer.nXY > 999) {
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH;
            } else {
                this.sampling = Sampling.SEMI_COMPLETE;
            }
        }

        if (debug) {
            System.out.println("findVoid sampling=" + sampling.name());
        }

        if (sampling.ordinal() == Sampling.LEAST_COMPLETE.ordinal()) {

            // uses divide and conquer
            findVoids(x, y, 0, x.length-1, 0, y.length-1);

        } else if (sampling.ordinal() == Sampling.COMPLETE.ordinal()) {

            int incr = 1;

            findVoidsUsingDoubleIndexes(incr);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH.ordinal()) {

            // uses a very rough range search
            findVoidsRoughRangeSearch(x, y, 2);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_2.ordinal()) {

            // uses a very rough range search
            findVoidsRoughRangeSearch(x, y, 10);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE.ordinal()) {

            // for datasets not evenly distributed, and with clusters >> outliers
            //  For this case, complete sampling works, but we want something faster.

            // N!/(2!( (N/8)-2)! * N!/(2!( (N/8) -2) ??  TODO: estimate this
            int incr = 8;

            findVoidsUsingDoubleIndexes(incr);
        }
    }

    /**
     * A divide and conquer approach to finding the rectangular areas
     * containing only two points.  it's a recursion with pattern
     * T(n) = 4T(n/2) + n  so the runtime is O(n^2).  It does not completely
     * sample every pair of points.
     *
     * The range search within void findVoids() is preferred even though slower
     * than this divide and conquer because the range search has a more
     * complete solution, that is a higher number of pairs bounding rectangular
     * voids are learned from the range search.
     * @param x
     * @param y
     * @param xIndexLo
     * @param xIndexHi
     * @param yIndexLo
     * @param yIndexHi
     */
    protected void findVoids(float[] x, float[] y, int xIndexLo, int xIndexHi,
        int yIndexLo, int yIndexHi) {
                                                                                 // cost     number of times
        if ((xIndexLo < xIndexHi) && (yIndexLo < yIndexHi)) {                    //

            int xIndexMid = (xIndexLo + xIndexHi)/2;                             //

            int yIndexMid = (yIndexLo + yIndexHi)/2;                             //

            findVoids(x, y, xIndexLo, xIndexMid, yIndexLo, yIndexMid);           // c4           N/2
            findVoids(x, y, xIndexLo, xIndexMid, yIndexMid + 1, yIndexHi);       // c5           N/2

            findVoids(x, y, xIndexMid + 1, xIndexHi, yIndexLo, yIndexMid);       // c6           N/2
            findVoids(x, y, xIndexMid + 1, xIndexHi, yIndexMid + 1, yIndexHi);   // c7           N/2

            processIndexedRegion(x, y, xIndexLo, xIndexHi, yIndexLo, yIndexHi);  //              N/2
        }
    }

    protected void findVoidsUsingDoubleIndexes(int incr) {
        // N!/(2!(N-2)! * N!/(2!(N-2)!

        float[] x = indexer.getX();
        float[] y = indexer.getY();

        for (int i = 0; i < indexer.nXY; i++) {
            if (debug) {
                System.out.println("findVoids i=" + i + "/" + indexer.nXY);
            }
            for (int ii = (i + 1); ii < indexer.nXY; ii+=incr) {
                for (int j = 0; j < indexer.nXY; j++) {
                    for (int jj = (j + 1); jj < indexer.nXY; jj+=incr) {
                        processIndexedRegion(x, y, i, ii, j, jj);
                    }
                }
            }
        }
    }

    /**
     *
     * @param x
     * @param y
     * @param nDiv used to form the number of intervals = nPoints/nDiv
     *        which should usually be 2 or so
     */
    protected void findVoidsRoughRangeSearch(float[] x, float[] y, int nDiv) {

        /* Divide y interval in half and execute the same size intervals in x over the full range

         Exaample, for array w/ 8 elements:
         y    x
         0:7  0:7

         0:3  4:7
         0:3  0:3
         4:7  0:3
         4:7  4:7

         0:1  0:1
         0:1  2:3
         0:1  4:5
         0:1  6:7

         2:3  0:1
         2:3  2:3
         2:3  4:5
         2:3  6:7
         ...

         start w/ smallest intervals of y and march across range.
         * increase by 2 units and march across range
         *
         * length is 8
         *   2's
         *   4's
         *   6's
         *   8's
         */
        int n = indexer.getNumberOfPoints();

        int nYIntervals = n / nDiv;                                          // cost     number of times

        for (int k = 0; k < nYIntervals; k++) {                             //               nDiv
            int binSz = (k + 1) * 2;                                        // c10
            binSz = (int)((k + 1) * 1.5f);

            int yLo = 0;
            while ((yLo + binSz) < n) {                                     //              (N/2)-1, (N/2)-2

                int nXIntervals = n / binSz;

                for (int j = 0; j < nXIntervals; j++) {                     //              N/2, N/4
                    int startX = j * binSz;
                    int endX = startX + binSz - 1;
//System.out.println("processIndexedRegion: " + startX + ":" + endX + ":" + yLo + ":" + (yLo + binSz));
                    processIndexedRegion(x, y, startX, endX, yLo, yLo + binSz);
                }
                yLo += 2;
            }
        }
    }

    /**
     * this returns the intersection of points within the area defined by the x
     * range and yrange where those are given in the reference frame of the
     * indexes sorted by each axis.
     *
     * @param x
     * @param y
     * @param xIndexLo
     * @param xIndexHi
     * @param yIndexLo
     * @param yIndexHi
     */
    private void processIndexedRegion(float[] x, float[] y, int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi) {

        boolean useCompleteSampling = (sampling.ordinal() == Sampling.COMPLETE.ordinal());

        int[] regionIndexes = indexer.findIntersectingRegionIndexesIfOnlyTwo(
            xIndexLo, xIndexHi, yIndexLo, yIndexHi, useCompleteSampling);

        int nPointsInRegion = regionIndexes.length;

        // an area which has 2 points in it and has the smallest nPoints/area, that is the
        // largest area is the value we are looking for.
        // we won't be finding the furthest pair because that 'rectangular area' would presumably
        // contain other points too.
        //    *an exception is made for useCompleteSampling if more than one same y is present
        if ( (nPointsInRegion < 2) || ((nPointsInRegion != 2) && !useCompleteSampling)) {
            return;
        }

        if (nPointsInRegion == 2) {

            processIndexedPair(x, y, regionIndexes[0], regionIndexes[1]);

        } else if (nPointsInRegion == 3) {

            // find the point that doesn't have the same y as the others.
            int uniqueYIndex = -1;
            for (int i = 0; i < regionIndexes.length; i++) {
                boolean foundSameYAsI = false;
                float ypi = y[ regionIndexes[i] ];
                for (int j = 0; j < regionIndexes.length; j++) {
                    if (i == j) {
                        continue;
                    }
                    float ypj = y[ regionIndexes[j] ];
                    if (ypj == ypi) {
                        foundSameYAsI = true;
                        break;
                    }
                }
                if (!foundSameYAsI) {
                    // we found the unique y
                    uniqueYIndex = i;
                    break;
                }
            }
            if (uniqueYIndex == -1) {
                StringBuilder err = new StringBuilder("ERROR: intersecting region contained more than 2 points, but unique y wasn't found");
                for (int i = 0; i < regionIndexes.length; i++) {
                    err.append("\n  (").append( x[ regionIndexes[i] ] ).append(", ").append( y[ regionIndexes[i] ] ).append(")");
                }
                throw new IllegalStateException(err.toString());
            }
            int regionIndex0 = regionIndexes[uniqueYIndex];
            for (int i = 0; i < regionIndexes.length; i++) {
                if (i != uniqueYIndex) {
                    int regionIndex1 = regionIndexes[i];
                    processIndexedPair(x, y, regionIndex0, regionIndex1);
                }
            }
        } else if (nPointsInRegion == 4) {
            // this is a rectangle, so write all combinations
            for (int i = 0; i < regionIndexes.length; i++) {
                int regionIndex0 = regionIndexes[i];
                for (int j = (i + 1); j < regionIndexes.length; j++) {
                    int regionIndex1 = regionIndexes[j];
                    processIndexedPair(x, y, regionIndex0, regionIndex1);
                }
            }
        }
    }

    /**
     * process the pair of points by calculating the linear density and storing it
     * in the instance arrays if not already stored.
     *
     * @param x
     * @param y
     * @param regionIndex0
     * @param regionIndex1
     */
    protected void processIndexedPair(float[] x, float[] y, int regionIndex0, int regionIndex1) {

        // Note: using 1-D instead of 2-D rectangles seems to be a better choice because it is not
        // dependent on rotation of the reference frame

        float d = (float) Math.sqrt(LinesAndAngles.distSquared( x[regionIndex0], y[regionIndex0],
            x[regionIndex1], y[regionIndex1]));

        if (d == 0) {
            return;
        }

        //float xc = (x[regionIndexes[0]] + x[regionIndexes[1]]) / 2.f;
        //float yc = (y[regionIndexes[0]] + y[regionIndexes[1]]) / 2.f;

        float linearDensity = 2.f / d;

        // expand arrays by 100 if needed
        if ((nTwoPointSurfaceDensities + 2) > allTwoPointSurfaceDensities.length) {
            allTwoPointSurfaceDensities = Arrays.copyOf(allTwoPointSurfaceDensities, nTwoPointSurfaceDensities + 100);
            point1 = Arrays.copyOf(point1, nTwoPointSurfaceDensities + 100);
            point2 = Arrays.copyOf(point2, nTwoPointSurfaceDensities + 100);
        }

        if (twoPointIdentities.storeIfDoesNotContain(regionIndex0, regionIndex1)) {
            allTwoPointSurfaceDensities[nTwoPointSurfaceDensities] = linearDensity;
            point1[nTwoPointSurfaceDensities] = regionIndex0;
            point2[nTwoPointSurfaceDensities] = regionIndex1;
            nTwoPointSurfaceDensities++;
        }
    }

    static float[] calculateCentroidOfTop(float[] xfit, float[] yfit, float frac) {

        XY xy = LinesAndAngles.createPolygonOfTopFWFractionMax(xfit, yfit, frac);

        return LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(xy.getX(), xy.getY());
    }

    protected float[] calcAreaAndCentroid(float[] x, float[] y, int[] regionIndexes) {

        if (regionIndexes.length == 1) {
            return null;
        } else if (regionIndexes.length == 2) {

            float xp0 = x[regionIndexes[0]];
            float yp0 = y[regionIndexes[0]];

            float xp1 = x[regionIndexes[1]];
            float yp1 = y[regionIndexes[1]];

            float xc = (xp0 + xp1) / 2.f;
            float yc = (yp0 + yp1) / 2.f;

            double dist = Math.sqrt(Math.pow((xp0 - xp1), 2) + Math.pow((yp0 - yp1), 2));
            double height = (dist > 0) ? 1 : Math.pow(1, Math.getExponent(dist));
            float area = (float) (height * dist);

            return new float[]{area, xc, yc};
        }

        float[] xPolygon = new float[regionIndexes.length];
        float[] yPolygon = new float[regionIndexes.length];
        for (int i = 0; i < regionIndexes.length; i++) {
            int index = regionIndexes[i];
            xPolygon[i] = x[index];
            yPolygon[i] = y[index];
        }

        if (xPolygon.length > 2) {

            try {
                GrahamScan scan = new GrahamScan();
                scan.computeHull(xPolygon, yPolygon);

                float[] xHull = scan.getXHull();
                float[] yHull = scan.getYHull();

                return LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(xHull, yHull);

            } catch (GrahamScanTooFewPointsException e) {
                return null;
            }

        } else {

            return null;
        }
    }

    public String persistTwoPointBackground() throws IOException {
        return serializeTwoPointDensities();
    }

    protected void serializeTwoPointBackground(ObjectOutputStream oos) throws IOException {

        oos.writeInt(nTwoPointSurfaceDensities);

        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            oos.writeFloat(allTwoPointSurfaceDensities[i]);
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            oos.writeInt(point1[i]);
            oos.writeInt(point2[i]);
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            oos.writeFloat(allTwoPointSurfaceDensitiesErrors[i]);
        }

        oos.flush();
    }

    protected String serializeTwoPointDensities() throws IOException {

        return serializeTwoPointBackground("stats_2pt_voids_");
    }

    public boolean readTwoPointBackground(String persistedFileName) throws IOException {

        boolean didDeserialize = deserializeTwoPointBackground(persistedFileName);;

        if (nTwoPointSurfaceDensities == 0) {
            throw new IOException("No pairs were found isolated within an area");
        } else {
            if (didDeserialize) {
                state = State.DENSITIES_CALCULATED;
            }
        }

        return didDeserialize;
    }

    protected void deserializeTwoPointBackground(ObjectInputStream ois) throws IOException {

        this.nTwoPointSurfaceDensities = ois.readInt();

        this.allTwoPointSurfaceDensities = new float[nTwoPointSurfaceDensities];
        this.allTwoPointSurfaceDensitiesErrors = new float[nTwoPointSurfaceDensities];
        this.point1 = new int[nTwoPointSurfaceDensities];
        this.point2 = new int[nTwoPointSurfaceDensities];

        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            this.allTwoPointSurfaceDensities[i] = ois.readFloat();
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            this.point1[i] = ois.readInt();
            this.point2[i] = ois.readInt();
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            this.allTwoPointSurfaceDensitiesErrors[i] = ois.readFloat();
        }
    }

    void plotFit(PolygonAndPointPlotter plotter) {

        if (statsHistogram == null) {
            return;
        }
        if (bestFit == null) {
            return;
        }
        try {
            float xHalfInterval = (statsHistogram.getXHist()[1] - statsHistogram.getXHist()[0]) / 2.0f;
            float xmin = statsHistogram.getXHist()[0] - xHalfInterval;
            float xmax = statsHistogram.getXHist()[statsHistogram.getXHist().length - 1];
            float ymin = MiscMath.findMin(statsHistogram.getYHistFloat());
            float ymax = MiscMath.findMax(statsHistogram.getYHistFloat());

            float[] xf = bestFit.getOriginalScaleX();
            float[] yf = bestFit.getOriginalScaleYFit();

            plotter.addPlot(statsHistogram.getXHist(), statsHistogram.getYHistFloat(),
                statsHistogram.getXErrors(), statsHistogram.getYErrors(), xf, yf, "");
            plotter.writeFile();

        } catch (IOException e) {
        }

    }
    void plot(TwoPointVoidStatsPlotter plotter, float xmin, float xmax, float ymin, float ymax) {

        try {
            if (statsHistogram != null) {

                plotPairSeparations(plotter, xmin, xmax, ymin, ymax);

                plotter.writeFile2();
            }

        } catch (IOException e) {
        }
    }
    protected void plotPairSeparations() {

        try {

            TwoPointVoidStatsPlotter plotter = new TwoPointVoidStatsPlotter();

            float[] mm = indexer.findXYMinMax();

            float xmin = MiscMath.roundDownByLargestPower(mm[0]);
            float xmax = MiscMath.roundUpByLargestPower(mm[1]);
            float ymin = MiscMath.roundDownByLargestPower(mm[2]);
            float ymax = MiscMath.roundUpByLargestPower(mm[3]);

            plotPairSeparations(plotter, xmin, xmax, ymin, ymax);

            plotter.writeFile();

        } catch (IOException e) {}
    }

    protected void plotPairSeparations(TwoPointVoidStatsPlotter plotter, float xmin,
        float xmax, float ymin, float ymax) {

        int[] t1 = Arrays.copyOf(point1, nTwoPointSurfaceDensities);
        int[] t2 = Arrays.copyOf(point2, nTwoPointSurfaceDensities);

        plotter.addTwoPointPlot(indexer.getX(), indexer.getY(), t1, t2,
            xmin, xmax, ymin, ymax);

        if (this.statsHistogram != null && allTwoPointSurfaceDensities != null) {

            float min = statsHistogram.getXHist()[0];
            float max = statsHistogram.getXHist()[statsHistogram.getXHist().length - 1] +
                (  (statsHistogram.getXHist()[1] - statsHistogram.getXHist()[0])/2.f);

            float[] tmp = new float[allTwoPointSurfaceDensities.length];
            int count = 0;
            for (int i = 0; i < tmp.length; i++) {
                if ((allTwoPointSurfaceDensities[i] >= min) && (allTwoPointSurfaceDensities[i] <= max)) {
                    tmp[count] = allTwoPointSurfaceDensities[i];
                    count++;
                }
            }
            tmp = Arrays.copyOf(tmp, count);
            plotter.addHistogram(tmp, max, statsHistogram.getXHist().length);
        }
    }
}
