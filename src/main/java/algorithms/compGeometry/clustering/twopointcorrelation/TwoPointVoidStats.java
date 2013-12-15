package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.XY;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.curves.FailedToConvergeException;
import algorithms.curves.GEVChiSquareMinimization;
import algorithms.curves.GEVYFit;
import algorithms.curves.GeneralizedExtremeValue;
import algorithms.curves.ICurveFitter;
import algorithms.curves.NonQuadraticConjugateGradientSolver;
import algorithms.misc.DoubleAxisIndexerStats;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.misc.Statistic;
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
 * the background density will be used by the calling program to find
 * clusters, that is groups of points associated by proximity, in the data.
 *
 * It calculates the two-point density function of rectangular voids, creates a
 * histogram from the distribution, fits a Generalized Extreme Value
 * distribution to the histogram and interprets that based upon information
 * about the number density of points in both dimensions.
 *
 * Essentially one must learn before sampling and analysis, whether the majority of points
   are in groups or are background points ahead of time.
   The only way to know that ahead of time is by rough cell counts (in 2-dimensions) before
   making the 1-dimensional histogram of counts.
   -- If a significant number of cells are empty, this is easily seen as SPARSE_BACKGROUND sampling
      and then the background density is the lowest bin's x value in a well formed histogram
   -- If all cells have values within stdev of avg it's an even distribution of points.
      this is seen as needing COMPLETE sampling
      and then the background density is the peak bin's x value in a well formed histogram
   -- Else, we should avoid empty cells and avoid the cells with counts > avg + stdev
      and then use COMPLETE sampling of that subset of cells
      and then the background density is the peak bin's x value in a well formed histogram
      -- Note that the cells should not be too large in number because that would be
         biasing the histogram to remove the densities that represent the distance between
         groups.
      -- Note that for this case, we really want to be able to remove the groups after first approx
         and then do a complete sampling on the data without those and determine the
         background density as the peak's x bin as the answer from that.
         the problem is that one needs to fit a profile to the groups
         in order to subtract them and to believe that that profile is correct and differentiable
         from a uniform background distribution (which are details of the 2-point correlation function...)

         TODO: following up on extreme differences for distributions in the later before
 *       improving the sampling to match the above description
 *
 <pre>
 * More specifically for the true background points:
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
 *    background point distribution and analyzing the distribution for the most
 *    frequently occurring densities (the peak) and the smallest densities for
 *    sparse background data sets.
 * </pre>
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
 * of the background density.  That density is used as "noise" which TwoPointCorrelation
 * looks for clusters in the signal higher than it by a threshhold factor that's
 * usually 2 or 3 (densities > threshold*noise).
 *
 * The runtimes for the code are still in progress, but roughly approximated
 * in 2 stages:
 * (1) calculating and fitting the background voids @see #findVoids@see(),
 * and (2) finding the groups within the data using a threshold from the background density.
 *
 * If debugging is turned on, plots are generated and those file paths are printed to
 * standard out, and statements are printed to standard out.
 *
 * @author nichole
 */
public class TwoPointVoidStats extends AbstractPointBackgroundStats {

    public enum State {
        POINTS_LOADED, DENSITIES_CALCULATED, HISTOGRAM_CREATED, HISTOGRAM_FITTED, STATS_FINALIZED
    }

    protected VoidSampling sampling = null;

    protected Boolean interpretForSparseBackground = null;
    
    boolean automateTheFindMethodChoice = false;

    protected State state = null;

    protected IVoidFinder voidFinder = null;

    protected int defaultNBins = 40;

    //for debugging, hold on to intermediate data:  histogram and bestFit
    protected HistogramHolder statsHistogram = null;
    protected GEVYFit bestFit = null;

    protected boolean doLogPerformanceMetrics = false;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    // uses conjugate gradient method for non quadratic functions if = true, else downhill simplex.
    protected boolean useDefaultFitting = false;

    /**
     * the factor to use when comparing a value to average += standardDeviation*sigmaFactor.
     * It should be the same as TwoPointCorrelation.sigmaFactor.
     * It's value is usually between 2 and 3.
     */
    protected float sigmaFactor = 2.5f;

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

    public void setUseDownhillSimplexHistogramFitting() {
        this.useDefaultFitting = false;
    }

    /**
     * set the type of sampling to be used on the dataset to calculate the 2-point
     * densities.  Note that this is expected to be invoked only from
     * TwoPointCorrelation.
     * @param sampling
     */
    protected void setUseCompleteSampling() {
        this.sampling = VoidSampling.COMPLETE;
    }

    public VoidSampling getSampling() {
        return VoidSampling.valueOf(sampling.name());
    }

    public Boolean getInterpretForSparseBackground() {
        return interpretForSparseBackground;
    }
    /**
     * set the interpretation of the density histogram to the method
     * used for sparse backgrounds.  Note that this is expected to
     * only be called from TwoPointCorrelation.
     */
    protected void setInterpretForSparseBackgroundToTrue() {
        this.interpretForSparseBackground = Boolean.TRUE;
    }

    /**
     * this method is not recommended, but is to allow the code to choose between
     * the default method of finding clusters based upon a density above the background
     * density or the alternative method for finding clusters when there are no background
     * points between the groups.
     * The reason that it is not recommended is that the 2 types of datasets that the methods
     * would be used on are so different
     * that you likely know which you have before use of this code.
     */
    protected void automateTheFindMethodChoice() {
        this.automateTheFindMethodChoice = true;
    }
    /**
     * set the value of the variable 'sigmaFactor'.
     * @param stDevFactor
     */
    protected void setStandardDeviationFactor(float stDevFactor) {
        this.sigmaFactor = stDevFactor;
    }

    public void releaseLargeVariables() {
        if (voidFinder != null) {
            voidFinder.releaseLargeVariables();
        }
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
            indexer.getNXY(), bigOh, diffSec, approximateMemoryUsed(),
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

        if (statsHistogram != null) {
            sumBytes += statsHistogram.approximateMemoryUsed();
        }
        if (bestFit != null) {
            sumBytes += bestFit.approximateMemoryUsed();
        }
        if (voidFinder != null) {
            sumBytes += voidFinder.approximateMemoryUsed();
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

        calculateStats();
    }

    @Override
    protected void calculateStats() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        long startTimeMillis = System.currentTimeMillis();

        // may need to release more memory.  if so, release point1 and point2
        //long memAvail = Util.getAvailableHeapMemory();
        //log.fine("memory available = " + memAvail);

        if (voidFinder != null) {
            log.fine("nXY=" + indexer.getNXY() + " nD=" + voidFinder.getNumberOfTwoPointDensities());
        }

        /*
        if (indexer.getNXY() > 999) {
            statsHistogram = createHistogramWithHigherPeakResolution();
        } else {*/
            statsHistogram = createHistogram();
        //}

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

        if (doLogPerformanceMetrics) {

            long stopTimeMillis = System.currentTimeMillis();

            printPerformanceMetrics(startTimeMillis, stopTimeMillis,
                "calculateStats", Integer.toString(statsHistogram.getXHist().length) );
        }

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

        findVoids();

        state = State.DENSITIES_CALCULATED;

    }

    protected void findVoids() throws TwoPointVoidStatsException {

        /* If a significant number of cells are empty ==> SPARSE_BACKGROUND sampling.
           background density is the lowest bin's x value in a well formed histogram.
          If all cells have values within stdev of avg ==> COMPLETE sampling
           background density is the peak bin's x value in a well formed histogram.
          Else:
             attempt to sample only the cells that are not empty nor above avg + stdev
             and use COMPLETE sampling on them.
             background density is the peak bin's x value in a well formed histogram
             (though this later could be improved)
         */

        long startTimeMillis = System.currentTimeMillis();

        int nXY = indexer.getNXY();

        // for reduced sampling of large sets, need these in scope:
        int nCellsPerDimension = (int)Math.sqrt(indexer.nXY/1000);

        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();

        if (automateTheFindMethodChoice) {
            
            int nCellsPerDimensionForStats = (int)Math.sqrt(indexer.nXY/300.);
            if (nCellsPerDimensionForStats < 9) {
                nCellsPerDimensionForStats = 9;
            }
            Statistic statistic = stats.calculateCellDensities(nCellsPerDimensionForStats, indexer);

            float fractionEmpty = stats.fractionOfCellsWithoutPoints(statistic);

            if (fractionEmpty > 0.3f) {

                interpretForSparseBackground = Boolean.TRUE;
                
            }
            
System.out.println("nCellsPerDim=" + nCellsPerDimensionForStats + " fractionEmpty=" + fractionEmpty + " indexer.nXY=" + indexer.getNumberOfPoints());
System.out.println("fractionNotAvg=" + stats.fractionOfCellsOutSideOfAvgTolerance(statistic, sigmaFactor));
            
            sampling = VoidSampling.COMPLETE;
        }
        
        if (sampling == null) {

            sampling = VoidSampling.COMPLETE;
            
        }
        if (interpretForSparseBackground == null) {
            
            interpretForSparseBackground = Boolean.FALSE;
        }


        int nSampled = -1;

        if (debug) {
            log.info("findVoid sampling=" + sampling.name() + " for " + nXY + " points");
        }

        if (sampling.ordinal() == VoidSampling.LEAST_COMPLETE.ordinal()) {

            voidFinder = new DivideAndConquerVoidFinder();

        } else if (sampling.ordinal() == VoidSampling.COMPLETE.ordinal()) {

            voidFinder = new CompleteSamplingVoidFinder();

        } else if (sampling.ordinal() == VoidSampling.COMPLETE_ON_RANDOM_SUBSET.ordinal()) {

            // xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi
            int[] xyMinMaxCell = stats.chooseARandomCell(nCellsPerDimension, indexer);

            voidFinder = new SubsetSamplingVoidFinder();

            ((SubsetSamplingVoidFinder) voidFinder).setXSortedIdxLo(xyMinMaxCell[0]);
            ((SubsetSamplingVoidFinder) voidFinder).setXSortedIdxHi(xyMinMaxCell[1]);
            ((SubsetSamplingVoidFinder) voidFinder).setYSortedIdxLo(xyMinMaxCell[2]);
            ((SubsetSamplingVoidFinder) voidFinder).setYSortedIdxHi(xyMinMaxCell[3]);

            nSampled = (xyMinMaxCell[1] - xyMinMaxCell[0]) * (xyMinMaxCell[3] - xyMinMaxCell[2]);

        } else {

            throw new IllegalStateException("Did not configure a finder for " + sampling.toString() + "?");
        }

        voidFinder.setSampling(sampling);

        voidFinder.findVoids(indexer);


        if (doLogPerformanceMetrics) {

            long stopTimeMillis = System.currentTimeMillis();

            String str = "";
            if (sampling.ordinal() == VoidSampling.LEAST_COMPLETE.ordinal()) {
                str = "O(n lg(n)) with n=";
            } else if (sampling.ordinal() == VoidSampling.COMPLETE.ordinal()) {
                str = "O(n^2) with n=";
            } else if (sampling.ordinal() == VoidSampling.COMPLETE_ON_RANDOM_SUBSET.ordinal()) {
                str = "O(n1*n2) = O(" + nSampled + ")";
            }

            printPerformanceMetrics(startTimeMillis, stopTimeMillis, "calculateBackgroundVia2PtVoidFit-->calculateTwoPointVoidDensities", str);
        }
    }

    protected HistogramHolder createHistogram() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        //int nBins = (indexer.getNumberOfPoints() < 100) ? defaultNBins/2 : defaultNBins;

        //HistogramHolder histogram = Histogram.createHistogramForSkewedData(
        //    nBins, voidFinder.getTwoPointDensities(), voidFinder.getTwoPointDensityErrors(), true);

        HistogramHolder histogram = Histogram.calculateSturgesHistogramRemoveZeroTail(
            voidFinder.getTwoPointDensities(), voidFinder.getTwoPointDensityErrors());

        plotPairSeparations();

        return histogram;
    }

    protected HistogramHolder createHistogramWithHigherPeakResolution() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        if (voidFinder == null) {
            return null;
        }

        int nBins = (indexer.getNumberOfPoints() < 100) ? defaultNBins/2 : defaultNBins;

        float[] minMax = MiscMath.calculateOuterRoundedMinAndMax(voidFinder.getTwoPointDensities());

        HistogramHolder histogram = Histogram.createHistogramForSkewedDataForPeakResolution2(
            nBins, voidFinder.getTwoPointDensities(), voidFinder.getTwoPointDensityErrors(),
            minMax[0], minMax[1], 0);

        plotPairSeparations();

        return histogram;
    }

    /**
     * internal method to calculate the statistics from the histogram using a
     * fit to the background distribution and a rough integration under the fit
     * to approximate a usable limit of the background density
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

            GEVYFit yfit = null;

            ICurveFitter chiSqMin = null;

            if (useDefaultFitting) {

                chiSqMin = new NonQuadraticConjugateGradientSolver(
                    histogram.getXHist(), histogram.getYHistFloat(), histogram.getXErrors(), histogram.getYErrors());

            } else {

                chiSqMin = new GEVChiSquareMinimization(
                    histogram.getXHist(), histogram.getYHistFloat(), histogram.getXErrors(), histogram.getYErrors());
            }

            chiSqMin.setDebug(debug);

            if (gevRangeFittingParameters != null) {

                if (useDefaultFitting) {

                    yfit = ((NonQuadraticConjugateGradientSolver)chiSqMin)
                        .fitCurveParametersAllAtOnce(gevRangeFittingParameters[0], gevRangeFittingParameters[1],
                            gevRangeFittingParameters[2], gevRangeFittingParameters[3], 0.001f, 0.3f);
                } else {

                    yfit =
                        ((GEVChiSquareMinimization)chiSqMin).fitCurveKGreaterThanZeroAndMu(
                            GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS,
                        gevRangeFittingParameters[0], gevRangeFittingParameters[1],
                        gevRangeFittingParameters[2], gevRangeFittingParameters[3]
                    );
                }

            } else {

                if (useDefaultFitting) {

                    yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

                } else {

                    yfit = ((GEVChiSquareMinimization)chiSqMin).fitCurveKGreaterThanZeroAndMu(
                        GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);
                }
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

            return yfit;

        } catch (FailedToConvergeException e) {
            throw new TwoPointVoidStatsException(e);
        } catch (IOException e2) {
            throw new TwoPointVoidStatsException(e2);
        }
    }

    protected void finalizeStats(HistogramHolder histogram, GEVYFit yfit) throws TwoPointVoidStatsException {

        if (yfit == null) {
            // this should never happen from calculateStatsForBackground
            throw new TwoPointVoidStatsException("yfit cannot be null");
        }
        this.statsHistogram = histogram;

        this.bestFit = yfit;

        try {

            if (debug) {
                log.info(bestFit.toString());

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

            /* for interpretForSparseBackground:
                   background density is the lowest bin's x value in a well formed histogram.
               for all other sampling:
                   background density is the peak bin's x value in a well formed histogram.
            */

            float limit, limitError;

            if (interpretForSparseBackground == null) {
                interpretForSparseBackground = Boolean.FALSE;
            }

            if (interpretForSparseBackground) {

                int yPeakIndex = MiscMath.findYMaxIndex(histogram.getYHist());

                limit = histogram.getXHist()[0];

            } else {

                //limitStr = "top centroid";

                // centroid of area defined by the top portion of the fit where y >= ypeak/2
                float[] areaAndXYTopCentroid = calculateCentroidOfTop(bestFit.getOriginalScaleX(), bestFit.getOriginalScaleYFit(), 0.5f);

                float[] areaAndXYTopCentroid2 = calculateCentroidOfTop(histogram.getXHist(), histogram.getYHistFloat(), 0.5f);

                float xpeak = (areaAndXYTopCentroid[0] > 0) ? areaAndXYTopCentroid[1] : areaAndXYTopCentroid2[1];
                
System.out.println("GEV top centroid=" + areaAndXYTopCentroid[1] + " data top centroid=" + areaAndXYTopCentroid2[1] + " using xpeak=" + xpeak);

                limit = (xpeak > 0) ? xpeak : bestFit.getXPeak();
            }

            this.backgroundSurfaceDensity = limit;

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
                String interp = (interpretForSparseBackground) ? "sparse" : "complete";
                log.info("estimating background from sampling= " + sampling.toString()
                + " and interpretation=" + interp
                + this.backgroundSurfaceDensity
                + " w/ x error in histogram bin =" + errorInEstimateFromHistogram
                + " gev fitting error for one point =" + errorInFitting
                + " error in the histogram x bin at y peak = " + limitError);
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

    public int getNumberOfDensityPoints() {
        return (voidFinder != null) ? voidFinder.getNumberOfTwoPointDensities() : 0;
    }

    public String persistTwoPointBackground() throws IOException {
        return serializeTwoPointDensities();
    }

    protected void serializeTwoPointBackground(ObjectOutputStream oos) throws IOException {

        if (voidFinder == null) {
            throw new IllegalStateException("no voidFinder to persist");
        }

        oos.writeInt(voidFinder.getNumberOfTwoPointDensities());

        for (int i = 0; i < voidFinder.getNumberOfTwoPointDensities(); i++) {
            oos.writeFloat(voidFinder.getTwoPointDensities()[i]);
        }
        for (int i = 0; i < voidFinder.getNumberOfTwoPointDensities(); i++) {
            oos.writeInt(voidFinder.getPoint1()[i]);
            oos.writeInt(voidFinder.getPoint2()[i]);
        }
        for (int i = 0; i < voidFinder.getNumberOfTwoPointDensities(); i++) {
            oos.writeFloat(voidFinder.getTwoPointDensityErrors()[i]);
        }

        oos.flush();
    }

    protected String serializeTwoPointDensities() throws IOException {

        return serializeTwoPointBackground("stats_2pt_voids_");
    }

    public boolean readTwoPointBackground(String persistedFileName) throws IOException {

        boolean didDeserialize = deserializeTwoPointBackground(persistedFileName);;

        if (voidFinder.getNumberOfTwoPointDensities() == 0) {
            throw new IOException("No pairs were found isolated within an area");
        } else {
            if (didDeserialize) {
                state = State.DENSITIES_CALCULATED;
            }
        }

        return didDeserialize;
    }

    protected void deserializeTwoPointBackground(ObjectInputStream ois) throws IOException {

        voidFinder = new VoidReader(ois);

    }

    void plotFit(PolygonAndPointPlotter plotter) {

        if (statsHistogram == null) {
            return;
        }
        if (bestFit == null) {
            return;
        }
        try {

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

        } catch (IOException e) {
            Logger.getLogger(SerializerUtil.class.getName()).severe(e.getMessage());
        }
    }

    protected void plotPairSeparations(TwoPointVoidStatsPlotter plotter, float xmin,
        float xmax, float ymin, float ymax) {

        if (voidFinder == null) {
            return;
        }

        if (voidFinder.getPoint1() != null) {
            int[] t1 = Arrays.copyOf(voidFinder.getPoint1(), voidFinder.getNumberOfTwoPointDensities());
            int[] t2 = Arrays.copyOf(voidFinder.getPoint2(), voidFinder.getNumberOfTwoPointDensities());

            plotter.addTwoPointPlot(indexer.getX(), indexer.getY(), t1, t2,
                xmin, xmax, ymin, ymax);
        }

        if (this.statsHistogram != null && (voidFinder.getNumberOfTwoPointDensities() > 0)) {

            float min = statsHistogram.getXHist()[0];
            float max = statsHistogram.getXHist()[statsHistogram.getXHist().length - 1] +
                (  (statsHistogram.getXHist()[1] - statsHistogram.getXHist()[0])/2.f);

            float[] tmp = new float[voidFinder.getNumberOfTwoPointDensities()];
            int count = 0;
            for (int i = 0; i < tmp.length; i++) {
                if ((voidFinder.getTwoPointDensities()[i] >= min) && (voidFinder.getTwoPointDensities()[i] <= max)) {
                    tmp[count] = voidFinder.getTwoPointDensities()[i];
                    count++;
                }
            }
            tmp = Arrays.copyOf(tmp, count);
            plotter.addHistogram(tmp, max, statsHistogram.getXHist().length);
        }
    }
}
