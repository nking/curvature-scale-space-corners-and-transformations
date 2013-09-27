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
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
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
 * The value returned by the code is used by TwoPointCorrelation as an estimate
 * of the background density.  That density is "noise" which TwoPointCorrelation
 * looks for clusters in the signal higher it (densities > threshold*noise).
 *
 * The runtimes for the code are still in progress, but roughly approximated
 * in 2 stages:  (1) calculating and fitting the background voids
 * @see #findVoids@see(), and (2) finding the groups within the data using
 * a threshold from the background density.
 *
 * If debugging is turned on, a plot is generated and the path is printed to
 * standard out, and statements are printed to standard out.
 *
 * @author nichole
 */
public class TwoPointVoidStats extends AbstractPointBackgroundStats {

    public enum Sampling {
        COMPLETE, LEAST_COMPLETE, SEMI_COMPLETE, SEMI_COMPLETE_RANGE_SEARCH,
        SEMI_COMPLETE_RANGE_SEARCH_2, SEMI_COMPLETE_RANGE_SEARCH_3, SEMI_COMPLETE_RANGE_SEARCH_4
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

    //for debugging, hold on to intermediate data:  histogram and bestFit
    protected HistogramHolder statsHistogram = null;
    protected GEVYFit bestFit = null;

    protected boolean adjustMuForDensity = true;

    protected boolean doLogPerformanceMetrics = false;

    protected Logger log = Logger.getLogger(this.getClass().getName());

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
     * The runtime for the complete sampling is N!/(2!(N-2)! * N!/(2!(N-2) from
     * permutations of iterating over x and y, and the result is effectively N^4.
     *
     * Complete sampling takes quite awhile, so other sampling methods are usually
     * used by default if this is not set by the user, with the exception of
     * datasets with less than 100 points.  For those, complete sampling is
     * used by default.
     */
    public void setUseCompleteSampling() {
        this.sampling = Sampling.COMPLETE;
    }

    /**
     * Use semi-complete range sampling.  This is a sampling method which scans
     * over a range of values with varying bin sizes.  Its good to use for
     * datasets in which there are outliers, that is data outside of "clusters".
     *
     * Note, if using this method, you may want to consider instead letting the code find it by default instead,
     * as the code will vary the bin size and its factor of change
     * depending upon data size.
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
     * The runtime is not yet estimated, but first look suggests
     * that it is N!/(2!( (N/8)-2)! * N!/(2!( (N/8) -2)
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

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;

        /*
         * enums:  one has 7 items
         *         one has 5 items
         */
        long sumBits = 7*nbits;
        long tmpSumBytes = (sumBits/8) + overheadBytes;
        long padding = (tmpSumBytes % 8);
        long sumBytes = tmpSumBytes + padding;
        sumBits = 5*nbits;
        tmpSumBytes = (sumBits/8) + overheadBytes;
        padding = (tmpSumBytes % 8);
        sumBytes += (tmpSumBytes + padding);

        // a reference to each of the enums
        sumBytes += (2*intBytes);

        // 4 references
        sumBytes += (4*intBytes);

        // 4 variables at stack word size each
        sumBytes += (4*intBytes);


        int n = (allTwoPointSurfaceDensities == null) ? 0 : allTwoPointSurfaceDensities.length;

        if (statsHistogram != null) {
            sumBytes += statsHistogram.approximateMemoryUsed();
        }
        if (bestFit != null) {
            sumBytes += bestFit.approximateMemoryUsed();
        }
        if (twoPointIdentities != null) {
            sumBytes += twoPointIdentities.approximateMemoryUsed();
        }

        if (allTwoPointSurfaceDensities != null) {
            sumBytes += (arrayBytes + (allTwoPointSurfaceDensities.length*intBytes));
        }

        if (allTwoPointSurfaceDensitiesErrors != null) {
            sumBytes += (arrayBytes + (allTwoPointSurfaceDensitiesErrors.length*intBytes));
        }

        if (point1 != null) {
            sumBytes += (arrayBytes + (point1.length*intBytes));
        }

        if (point2 != null) {
            sumBytes += (arrayBytes + (point2.length*intBytes));
        }

        sumBytes += overheadBytes;

        // amount of padding needed to make it a round 8 bytes
        padding = (sumBytes % 8);

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
                str = "(n*bf/nDiv)^1.8) or O(n^2 - n) with n=";
            } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_2.ordinal()) {
                str = "? n=";
            } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_3.ordinal()) {
                str = "? n=";
            } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_4.ordinal()) {
                str = "(n*bf/nDiv)^2.2) n=";
            } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE.ordinal()) {
                str = "(n/8)^(2.25) n=";
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
                    log.info(yfit.toString());
                }
            }

            boolean didCalculateFitStats = yfit.calculateStatsAsStepFunction();
            if (!didCalculateFitStats) {
                if (debug) {
                    plotPairSeparations();
                }
                throw new TwoPointVoidStatsException("histogram of linear densities was not fittable");
            }

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
                    log.info(bestFit.toString());
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
            // TODO:  when the histograms are improved for small number datasets, perhaps this won't be necessary
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
                log.info(bestFit.toString());
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
                log.info("estimating background minimum as location of " + limitStr + " of background profile counts = "
                + this.backgroundSurfaceDensity
                + " w/ x error in histogram bin =" + errorInEstimateFromHistogram
                + " gev fitting error for one point =" + errorInFitting);
            }

            if (debug && (bestFit.getChiSqSum() > bestFit.getYDataErrSq())) {
                log.info("WARNING:  chisq is larger than errors: "
                    + bestFit.getChiSqSum() + " (errsqsum=" + bestFit.getYDataErrSq() + ")");
            }

            state = State.STATS_FINALIZED;

        } catch (IOException e) {
            throw new TwoPointVoidStatsException(e);
        }
    }

    /**
     * find the rectangular 2 point densities using the method already set.
     *
     * If the method has not yet been set, the following information is used to
     * determine sampling:
     *
     * When no sampling has been set, the method attempts to assign least amount
     * of sampling needed to build a decent histogram.
     * (One day it should quickly look at the overall structure of the data
     * to learn if outliers are likely, and hence choose between a default
     * sampling pattern similar to Sampling.COMPLETE for
     * no outliers, or sampling pattern similar to Sampling.SEMI_COMPLETE_RANGE_SEARCH
     * when outliers are expected, but it does not at this time.)
     * Until then, the number of data points are used to decide sampling.
     * If SEMI_COMPLETE was used, the ability to improve the sampling is done
     * in the invoker @see TwoPointCorrelation#bruteForceCalculateGroups() which uses
     * @see TwoPointCorrelation#temporaryWorkaroundForSampling().
     * There if SEMI_COMPLETE was used and there are outliers, the background fit
     * is redone using SEMI_COMPLETE_RANGE_SEARCH
     *
     * The rules when sampling has not been set are as follows (along with the
     * "refit" logic just described above):

            if (indexer.nXY > 10000) {
                // RT (n^1.8) to (n^2.2) roughly...
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH_4;
            } else if (indexer.nXY > 8000) {
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH_3;
            } else if (indexer.nXY > 999) {
                // RT (n^1.8) roughly...
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH;
            } else if (indexer.nXY < 100) {
                // RT O(n^4)
                this.sampling = Sampling.COMPLETE;
            } else {
                // RT O(n lg(n))
                this.sampling = Sampling.SEMI_COMPLETE;
            }

     */
    protected void findVoids() {

        float[] x = indexer.getX();
        float[] y = indexer.getY();

        if (sampling == null) {

            // ***** NOTE ********
            // This block attempts to assign least amount of sampling needed to build a decent
            // histogram.  One day it should quickly look at the overall structure of the data
            // to learn if outliers are likely, and hence choose
            // between a default sampling pattern similar to Sampling.COMPLETE for
            // no outliers, or sampling pattern similar to Sampling.SEMI_COMPLETE_RANGE_SEARCH
            // when outliers are expected.
            // Until then, the number of data points are used below and then
            // one style of sampling is used and the other is attempted if it
            // looks like that was not a good choice.  Note that the "try again" behavior
            // is pushed up to @see TwoPointCorrelation#bruteForceCalculateGroups() use
            // of @see TwoPointCorrelation#temporaryWorkaroundForSampling() for now.
            // There if SEMI_COMPLETE was used and there are outliers, the background fit
            // is redone using SEMI_COMPLETE_RANGE_SEARCH

            if (indexer.nXY > 10000) {
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH_4;
            } else if (indexer.nXY > 8000) {
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH_3;
            } else if (indexer.nXY > 999) {
                this.sampling = Sampling.SEMI_COMPLETE_RANGE_SEARCH;
            } else if (indexer.nXY < 100) {
                this.sampling = Sampling.COMPLETE;
            } else {
                this.sampling = Sampling.SEMI_COMPLETE;
            }
        }

        if (debug) {
            log.info("findVoid sampling=" + sampling.name() + " for " + indexer.nXY + " points");
        }

        if (sampling.ordinal() == Sampling.LEAST_COMPLETE.ordinal()) {

            // uses divide and conquer
            findVoids(0, x.length - 1, 0, y.length - 1);

        } else if (sampling.ordinal() == Sampling.COMPLETE.ordinal()) {

            int incr = 1;

            findVoidsUsingDoubleIndexes(incr);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH.ordinal()) {

            // uses a very rough range search
            findVoidsRoughRangeSearch(0, indexer.nXY - 1, 0, indexer.nXY - 1, 2, 1.5f);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_2.ordinal()) {

            // uses a very rough range search
            int nDiv = 10;
            float bFactor = 4;
            findVoidsRoughRangeSearch(0, indexer.nXY - 1, 0, indexer.nXY - 1, nDiv, bFactor);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_3.ordinal()) {

            // for area of data so large that randomly chosen patches are neccessary
            //   to reduce sample to decrease runtime
            findVoidsRandomSamples(10, 10);

        } else if (sampling.ordinal() == Sampling.SEMI_COMPLETE_RANGE_SEARCH_4.ordinal()) {

            // for area of data so large that randomly chosen patches are neccessary
            //   to reduce sample to decrease runtime
            //findVoidsRandomSamples(20, 10);
            findVoidsRoughRangeSearch(0, indexer.nXY - 1, 0, indexer.nXY - 1, 10, 4f);

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
     * @param xIndexLo
     * @param xIndexHi
     * @param yIndexLo
     * @param yIndexHi
     */
    protected void findVoids(int xIndexLo, int xIndexHi,
        int yIndexLo, int yIndexHi) {
                                                                                 // cost     number of times
        if ((xIndexLo < xIndexHi) && (yIndexLo < yIndexHi)) {                    //

            int xIndexMid = (xIndexLo + xIndexHi)/2;                             //

            int yIndexMid = (yIndexLo + yIndexHi)/2;                             //

            findVoids(xIndexLo, xIndexMid, yIndexLo, yIndexMid);           // c4           N/2
            findVoids(xIndexLo, xIndexMid, yIndexMid + 1, yIndexHi);       // c5           N/2

            findVoids(xIndexMid + 1, xIndexHi, yIndexLo, yIndexMid);       // c6           N/2
            findVoids(xIndexMid + 1, xIndexHi, yIndexMid + 1, yIndexHi);   // c7           N/2

            processIndexedRegion(xIndexLo, xIndexHi, yIndexLo, yIndexHi);  //              N/2
        }
    }

    protected void findVoidsUsingDoubleIndexes(int incr) {
        // N!/(2!(N-2)! * N!/(2!(N-2)!

        findVoidsUsingDoubleIndexes(0, indexer.nXY - 1, 0, indexer.nXY - 1, incr);
    }

    protected void findVoidsUsingDoubleIndexes(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, int incr) {
        // N!/(2!(N-2)! * N!/(2!(N-2)!

        for (int i = xIndexLo; i < xIndexHi; i++) {
            if (debug) {
                log.info("findVoids i=" + i + "/" + indexer.nXY);
            }
            for (int ii = (i + 1); ii < indexer.nXY; ii+=incr) {
                for (int j = yIndexLo; j < yIndexHi; j++) {
                    for (int jj = (j + 1); jj < indexer.nXY; jj+=incr) {
                        processIndexedRegion(i, ii, j, jj);
                    }
                }
            }
        }
    }

    /**
     * @param xIndexLo index given w.r.t. indexer.sortedXIndexes
     * @param xIndexHi index given w.r.t. indexer.sortedXIndexes
     * @param yIndexLo index given w.r.t. indexer.sortedYIndexes
     * @param yIndexHi index given w.r.t. indexer.sortedYIndexes
     * @param nDiv used to form the number of intervals = nPoints/nDiv
     *        which should usually be 2 or so
     * @param bFactor
     */
    protected void findVoidsRoughRangeSearch(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, int nDiv, float bFactor) {

        //Divide y interval in half and execute the same size intervals in x over the full range
        int nYIntervals = (yIndexHi - yIndexLo) / nDiv;                     // cost     number of times

        for (int k = 0; k < nYIntervals; k++) {                             //               nDiv
            int binSz = (int)((k + 1) * bFactor);                           // c10

            int yLo = yIndexLo;
            while ((yLo + binSz) < yIndexHi) {                              //               nDiv

                int nXIntervals = (xIndexHi - xIndexLo)/ binSz;

                for (int j = 0; j < nXIntervals; j++) {                     //              nDiv/bfactor
                    int startX = xIndexLo + (j * binSz);
                    int endX = startX + binSz;
//System.out.println("processIndexedRegion: " + startX + ":" + endX + ":" + yLo + ":" + (yLo + binSz));
                    processIndexedRegion(startX, endX, yLo, yLo + binSz);
                }
                yLo += binSz;
            }
        }
    }

    /**
     *
     * runtime estimation  is O(n^2) + O(n)
     *     where if evenly distributed, n = N/nSamples, 10*( (N/10)^2 + (N/10) ) ?  not yet checked...
     *
     * @param nSamples the number of samples to take
     * @param nDivisionsPerSide
     */
    protected void findVoidsRandomSamples(int nSamples, int nDivisionsPerSide) {

        try {
            int n = indexer.getNumberOfPoints();

            SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(System.nanoTime());

            // find nSamples non-overlapping regions to sample the area

            /*
             * divide into a grid of nDivisionsPerSide by nDivisionsPerSide and choose nSamples randomly from them
             *
             * Note that the division is in index space rather than value space to
             * use the indexer.
             *
             */
            int binSize = indexer.nXY/nDivisionsPerSide;

            /*    col 0
             *     ||
             *     \/
             *      0 |  1 | 2     <=== row 0
             *    ---------------
             *        |    |
             *      3 |  4 | 5
             *    ---------------
             *        |    |
             *      6 |  7 | 8
             */

            int nSSq = nDivisionsPerSide*nDivisionsPerSide;

            // choices are 0 through nSamples-1
            boolean[] selected = new boolean[nSSq];

            // q=0,1,2,3
            int quadNumber = 0;

            for (int i = 0; i < nSamples; i++) {

                quadNumber++;
                if (quadNumber > 4) {
                    quadNumber = 0;
                }

                // to create more evenly distributed random sampling,
                //     will try for each quadrant
                //   3*nDiv:4*nDiv
                //
                //   2*nDiv:3*nDiv
                //
                //   nDiv:2*nDiv
                //
                //   0:nDiv         nDiv:2*nDiv       2*nDiv:3*nDiv       3*nDiv:4*nDiv
                //
                //   0,1   1,1
                //   0,0   1,0

                /*int bin = sr.nextInt(nSSq);
                while (selected[bin]) {
                    bin = sr.nextInt(nSSq);
                }
                selected[bin] = true;
                int row = (bin/nDivisionsPerSide);
                int col = (bin % nDivisionsPerSide);
                */

                int row = 0;
                int col = 0;
                int bin = 0;

                boolean draw = true;
                while (draw) {
                    int dCol = sr.nextInt(nDivisionsPerSide/2);
                    int dRow = sr.nextInt(nDivisionsPerSide/2);
                    switch(quadNumber) {
                        case 0:
                            col = dCol;
                            row = dRow;
                            break;
                        case 1:
                            col = 2*dCol;
                            row = dRow;
                            break;
                        case 2:
                            col = dCol;
                            row = 2*dRow;
                            break;
                        default:
                            col = 2*dCol;
                            row = 2*dRow;
                            break;
                    }
                    bin = col + (row*nDivisionsPerSide);
                    draw = (selected[bin]);
                }
                selected[bin] = true;

                int startX = col*binSize;
                int endX = startX + binSize;
                int yLo = row*binSize;
                int yHi = yLo + binSize;

                if (debug) {
                    /*log.info("[" + quadNumber + "] " + " bin =" + bin + " nDiv=" + nDivisionsPerSide
                        + " " + String.format("  %4d : %4d", col, row)
                        + " " + String.format("  [X %.4f : %.4f] [Y %.4f : %.4f]",
                        indexer.x[indexer.sortedXIndexes[startX]], indexer.x[indexer.sortedXIndexes[endX]],
                        indexer.y[yLo], indexer.y[yHi])  );*/

                    log.info("processIndexedRegion: " + startX + ":" + endX + ":" + yLo + ":" + yHi);
                }

                findVoidsRoughRangeSearch(startX, endX, yLo, yHi, 2, 2);
            }

        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }

    /**
     * this returns the intersection of points within the area defined by the x
     * range and y range where those are given in the reference frame of the
     * indexes sorted by each axis.
     *
     * @param x
     * @param y
     * @param xIndexLo index given w.r.t. indexer.sortedXIndexes
     * @param xIndexHi index given w.r.t. indexer.sortedXIndexes
     * @param yIndexLo index given w.r.t. indexer.sortedYIndexes
     * @param yIndexHi index given w.r.t. indexer.sortedYIndexes
     */
    private void processIndexedRegion(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi) {

        boolean useCompleteSampling = (sampling.ordinal() == Sampling.COMPLETE.ordinal());

        // returns indexes w.r.t. indexer.x and indexer.y
        int[] regionIndexes = indexer.findIntersectingRegionIndexesIfOnlyTwo(
            xIndexLo, xIndexHi, yIndexLo, yIndexHi, useCompleteSampling);

        int nPointsInRegion = regionIndexes.length;

        // an area which has 2 points in it and has the smallest nPoints/area, that is the
        // largest area is the value we are looking for.
        // we won't be finding the furthest pair because that 'rectangular area' would presumably
        // contain other points too.
        //    *an exception is made for useCompleteSampling if more than one same y is present
        if ( (nPointsInRegion < 2) || ((nPointsInRegion != 2) && !useCompleteSampling) ) {
            return;
        }

        if (nPointsInRegion == 2) {

            processIndexedPair(regionIndexes[0], regionIndexes[1]);

        } else if (nPointsInRegion == 3) {

            // find the point that doesn't have the same y as the others.
            int uniqueYIndex = -1;
            for (int i = 0; i < regionIndexes.length; i++) {
                boolean foundSameYAsI = false;
                float ypi = indexer.y[ regionIndexes[i] ];
                for (int j = 0; j < regionIndexes.length; j++) {
                    if (i == j) {
                        continue;
                    }
                    float ypj = indexer.y[ regionIndexes[j] ];
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
                    err.append("\n  (").append( indexer.x[ regionIndexes[i] ] )
                        .append(", ").append( indexer.y[ regionIndexes[i] ] ).append(")");
                }
                throw new IllegalStateException(err.toString());
            }
            int regionIndex0 = regionIndexes[uniqueYIndex];
            for (int i = 0; i < regionIndexes.length; i++) {
                if (i != uniqueYIndex) {
                    int regionIndex1 = regionIndexes[i];
                    processIndexedPair(regionIndex0, regionIndex1);
                }
            }
        } else if (nPointsInRegion == 4) {
            // this is a rectangle, so write all combinations
            for (int i = 0; i < regionIndexes.length; i++) {
                int regionIndex0 = regionIndexes[i];
                for (int j = (i + 1); j < regionIndexes.length; j++) {
                    int regionIndex1 = regionIndexes[j];
                    processIndexedPair(regionIndex0, regionIndex1);
                }
            }
        }
    }

    /**
     * process the pair of points by calculating the linear density and storing it
     * in the instance arrays if not already stored.
     *
     * @param regionIndex0 index w.r.t. indexer.x and indexer.y
     * @param regionIndex1 index w.r.t. indexer.x and indexer.y
     */
    protected void processIndexedPair(int regionIndex0, int regionIndex1) {

        // Note: using 1-D instead of 2-D rectangles seems to be a better choice because it is not
        // dependent on rotation of the reference frame

        float d = (float) Math.sqrt(LinesAndAngles.distSquared(
            indexer.x[regionIndex0], indexer.y[regionIndex0],
            indexer.x[regionIndex1], indexer.y[regionIndex1]));

        if (d == 0) {
            return;
        }

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
