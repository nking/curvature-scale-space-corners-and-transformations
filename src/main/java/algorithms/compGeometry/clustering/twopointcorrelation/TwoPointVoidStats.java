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
 * clusters in the data.
 *
 * It calculates the two-point density function of rectangular voids, creates a
 * histogram from the distribution, fits a Generalized Extreme Value
 * distribution to the histogram, integrates the area under the curve to
 * estimate the background density. For denser background densities, a
 * slightly different integration limit is used.
 *
 * This assumes that the data contains background points and clustered points
 * whose 2 distributions are different from one another, but homogeneous within
 * their own. It assumes that the distributions do not need to be fit well, but
 * that only the rough range of the background density needs to be
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

    public enum State {
        POINTS_LOADED, DENSITIES_CALCULATED, HISTOGRAM_CREATED, HISTOGRAM_FITTED, STATS_FINALIZED
    }

    protected VoidSampling sampling = null;

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
        this.sampling = VoidSampling.LEAST_COMPLETE;
    }
    
    public void setUseDownhillSimplexHistogramFitting() {
        this.useDefaultFitting = false;
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
        this.sampling = VoidSampling.COMPLETE;
    }

    public VoidSampling getSampling() {
        return VoidSampling.valueOf(sampling.name());
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
        
        if (indexer.getNXY() > 999) {
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

    protected HistogramHolder createHistogram() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        int nBins = (indexer.getNumberOfPoints() < 100) ? defaultNBins/2 : defaultNBins;

        //HistogramHolder histogram = Histogram.createHistogramForSkewedData(nBins, allTwoPointSurfaceDensities,
        //    allTwoPointSurfaceDensitiesErrors, false);
        
        HistogramHolder histogram = Histogram.createHistogramForSkewedData(
            nBins, voidFinder.getTwoPointDensities(), voidFinder.getTwoPointDensityErrors(), true);

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
            
            //TODO:  improve calcs for spatially separated groups without background points between them.
            //maybe take use findVoids of non empty cells only
            DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
            float fracEmpty = stats.fractionOfCellsWithoutPoints(6, indexer);
            log.finest("fracEmpty cells = " + fracEmpty); 

            if ((fracEmpty > 0.5) && (yPeakIndex > 0)) {
                isSmallNumberHist = true;
                this.backgroundSurfaceDensity = histogram.getXHist()[1];
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

    protected void findVoids() throws TwoPointVoidStatsException {

        long startTimeMillis = System.currentTimeMillis();
        
        int nXY = indexer.getNXY();
        
        // for reduced sampling of large sets, need these in scope:
        int nCellsPerDimension = (int)Math.sqrt(indexer.nXY/1000);

        DoubleAxisIndexerStats stats = new DoubleAxisIndexerStats();
                
        int nSampled = -1;
        
        if (nXY >= 4000) {
            
            // fits to void histograms of datasets where nXY > 1000 are the better fits.
            // for larger datasets, can take a sub-sample the void densities to reduce the runtime if can quickly see that the data is evenly sampled.
            // 2 cells holding 1000 points each is a dataset of size (2*sqrt(1000))^2 = 4000 at least.            
            
            boolean doesNotHaveLargeGaps = stats.doesNotHaveLargeGaps(nCellsPerDimension, indexer);
            
            log.info(" dataset spatial distribution doesNotHaveLargeGaps=" + doesNotHaveLargeGaps);
            
            if (doesNotHaveLargeGaps) {
                
                // reduce the sampling to a cell holding 1000 points 
                // and use the sampling that would be applied for 1000 points
                sampling = VoidSampling.SEMI_COMPLETE_SUBSET;               
            }
        }
        
        if (sampling == null) {

            this.sampling = VoidSampling.COMPLETE;
        }        

        if (debug) {
            log.info("findVoid sampling=" + sampling.name() + " for " + nXY + " points");
        }

        if (sampling.ordinal() == VoidSampling.LEAST_COMPLETE.ordinal()) {

            voidFinder = new DivideAndConquerVoidFinder();
            
            voidFinder.setSampling(sampling);
            
            voidFinder.findVoids(indexer);

        } else if (sampling.ordinal() == VoidSampling.COMPLETE.ordinal()) {

            voidFinder = new CompleteSamplingVoidFinder();
            
            voidFinder.setSampling(sampling);
            
            voidFinder.findVoids(indexer);

        } else if (sampling.ordinal() == VoidSampling.SEMI_COMPLETE_SUBSET.ordinal()) {

            // xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi
            int[] xyMinMaxCell = stats.chooseARandomCell(nCellsPerDimension, indexer);
           
            voidFinder = new SubsetSamplingVoidFinder();

            voidFinder.setSampling(sampling);

            ((SubsetSamplingVoidFinder) voidFinder).setXSortedIdxLo(xyMinMaxCell[0]);
            ((SubsetSamplingVoidFinder) voidFinder).setXSortedIdxHi(xyMinMaxCell[1]);
            ((SubsetSamplingVoidFinder) voidFinder).setYSortedIdxLo(xyMinMaxCell[2]);
            ((SubsetSamplingVoidFinder) voidFinder).setYSortedIdxHi(xyMinMaxCell[3]);

            voidFinder.findVoids(indexer);

            nSampled = (xyMinMaxCell[1] - xyMinMaxCell[0]) * (xyMinMaxCell[3] - xyMinMaxCell[2]);

        }
        
        if (doLogPerformanceMetrics) {

            long stopTimeMillis = System.currentTimeMillis();

            String str = "";
            if (sampling.ordinal() == VoidSampling.LEAST_COMPLETE.ordinal()) {
                str = "O(n lg(n)) with n=";
            } else if (sampling.ordinal() == VoidSampling.COMPLETE.ordinal()) {
                str = "O(n^2) with n=";
            } else if (sampling.ordinal() == VoidSampling.SEMI_COMPLETE_SUBSET.ordinal()) {
                str = "O(n1*n2) = O(" + nSampled + ")";
            }
            
            printPerformanceMetrics(startTimeMillis, stopTimeMillis, "calculateBackgroundVia2PtVoidFit-->calculateTwoPointVoidDensities", str);
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
