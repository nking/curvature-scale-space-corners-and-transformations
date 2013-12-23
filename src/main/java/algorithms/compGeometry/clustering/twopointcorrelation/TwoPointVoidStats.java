package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.curves.DerivGEV;
import algorithms.curves.FailedToConvergeException;
import algorithms.curves.GEVChiSquareMinimization;
import algorithms.curves.GEVYFit;
import algorithms.curves.GeneralizedExtremeValue;
import algorithms.curves.ICurveFitter;
import algorithms.curves.NonQuadraticConjugateGradientSolver;
import algorithms.misc.AxisIndexerStats;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.misc.Statistic;
import algorithms.util.ArrayPair;
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
   are in groups and have background points outside of the groups ahead of time or whether
   there are no background points outside of groups.
   The easiest way to approximate that ahead of time is by rough cell counts (in 2-dimensions) before
   making the 1-dimensional histogram of counts.
   -- (1) If a significant number of cells are empty, this seen as needing SPARSE_BACKGROUND interpretation
      and then the background density is the lowest bin's x value in a well formed histogram GEV fit
   -- (2) Else  the background density is the peak of the histogram or GEV fit in a well formed histogram
   
   The code automatically determines which of method (1) and (2) to use.
  
  If the user has better knowledge of which should be applied, can set that with:
     setInterpretForSparseBackgroundToTrue() or setInterpretForSparseBackgroundToFalse()
     
 <pre>
 * More specifically for the true background points:
 * -- The location of the 'background' points in two dimensional space are likely
 *    Poisson, that is their locations in a fixed interval of space are
 *    independent of one another and occurred randomly.
 * -- The areas between voids in such a distribution are well fit by
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
 * If debugging is turned on, plots are generated and those file paths are printed to
 * standard out, and statements are printed to standard out.
 *
 * @author nichole
 */
public class TwoPointVoidStats extends AbstractPointBackgroundStats {

    public enum State {
        POINTS_LOADED, DENSITIES_CALCULATED, HISTOGRAM_CREATED, HISTOGRAM_FITTED, STATS_FINALIZED
    }

    // null is signficant, so don't set a default unless change the code where check for null
    protected VoidSampling sampling = null;

    // null is signficant, so don't set a default unless change the code where check for null
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
     * the factor to use when comparing a density to backgroundDensity*sigmaFactor
     * or when comparing cell counts to average += standardDeviation*sigmaFactor.
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
    public TwoPointVoidStats(AxisIndexer indexedSortedPoints) {

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
        return sampling;
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
        if (automateTheFindMethodChoice) {
            throw new IllegalStateException("cannot have both 'automate' and 'set to sparse interpretation'");
        }
        this.interpretForSparseBackground = Boolean.TRUE;
    }
    protected void setInterpretForSparseBackgroundToFalse() {
        if (automateTheFindMethodChoice) {
            throw new IllegalStateException("cannot have both 'automate' and 'unset to sparse interpretation'");
        }
        this.interpretForSparseBackground = Boolean.FALSE;
    }

    protected void automateTheFindMethodChoice() {
        if (interpretForSparseBackground != null) {
            throw new IllegalStateException("cannot have both 'automate' and 'set or unset sparse interpretation'");
        }
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
     * statistics of those points to get the background density and error.
     *
     * More specifically:
     * @see #findVoids()
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

        statsHistogram = createHistogram();
        
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

        long startTimeMillis = System.currentTimeMillis();

        int nXY = indexer.getNXY();

        // for reduced sampling of large sets, need these in scope:
        int nCellsPerDimension = (int)Math.sqrt(indexer.nXY/1000);

        AxisIndexerStats stats = new AxisIndexerStats();

        if (automateTheFindMethodChoice) {
            
            int nCellsPerDimensionForStats = (int)Math.sqrt(indexer.nXY/300.);
            if (nCellsPerDimensionForStats < 9) {
                nCellsPerDimensionForStats = 9;
            }
            Statistic statistic = stats.calculateCellDensities(nCellsPerDimensionForStats, indexer);

            float fractionEmpty = stats.fractionOfCellsWithoutPoints(statistic);

            if (fractionEmpty > 0.3f) {

                interpretForSparseBackground = Boolean.TRUE;
                
                sampling = VoidSampling.COMPLETE;
                
            } else if (nXY > 10000) {
                
                sampling = VoidSampling.COMPLETE_ON_SUBSET;
                
                if (nXY > 50000) {
                    
                    nCellsPerDimension = 4;
                    
                } else if (nXY > 20000) {
                    
                    nCellsPerDimension = 3;
                    
                } else {
                    
                    nCellsPerDimension = 2;
                }
                
            } else {
                
                sampling = VoidSampling.COMPLETE;
                
            }
                                        
            log.fine("nCellsPerDim=" + nCellsPerDimensionForStats + " fractionEmpty=" + fractionEmpty 
                + " indexer.nXY=" + indexer.getNumberOfPoints() + "  avg of cells=" + statistic.getAverage()
                + "  cellXSize=" + statistic.getXSz() + " cell counts=" + Arrays.toString(statistic.getItems()));
            log.fine("fractionNotAvg=" + stats.fractionOfCellsOutSideOfAvgTolerance(statistic, sigmaFactor));
        }
        
//sampling = VoidSampling.COMPLETE_ON_SUBSET;

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

        } else if (sampling.ordinal() == VoidSampling.COMPLETE_ON_SUBSET.ordinal()) {

            // indexLo, int indexHi 
            // this chooses along the diagonal to keep x and y range the same to avoid needing to search along y more
            int[] xyMinMaxCell = stats.chooseARandomDiagonalCell(nCellsPerDimension, indexer);

            voidFinder = new SubsetSamplingVoidFinder();

            ((SubsetSamplingVoidFinder) voidFinder).setSortedIdxLo(xyMinMaxCell[0]);
            ((SubsetSamplingVoidFinder) voidFinder).setSortedIdxHi(xyMinMaxCell[1]);

            nSampled = (xyMinMaxCell[1] - xyMinMaxCell[0]) * (xyMinMaxCell[1] - xyMinMaxCell[0]);

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
            } else if (sampling.ordinal() == VoidSampling.COMPLETE_ON_SUBSET.ordinal()) {
                str = "O(n1*n2) = O(" + nSampled + ")";
            }

            printPerformanceMetrics(startTimeMillis, stopTimeMillis, "calculateBackgroundVia2PtVoidFit-->calculateTwoPointVoidDensities", str);
        }
    }

    protected HistogramHolder createHistogram() throws TwoPointVoidStatsException {

        if (state.ordinal() < State.DENSITIES_CALCULATED.ordinal()) {
            calculateTwoPointVoidDensities();
        }

        HistogramHolder histogram = Histogram.defaultHistogramCreator(
            voidFinder.getTwoPointDensities(), voidFinder.getTwoPointDensityErrors());

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
                   background density is the peak in a well formed histogram.
            */

            float limit, limitError;

            if (interpretForSparseBackground == null) {
                interpretForSparseBackground = Boolean.FALSE;
            }

            if (interpretForSparseBackground) {

                int yPeakIndex = MiscMath.findYMaxIndex(histogram.getYHist());

                limit = histogram.getXHist()[0];

                int limitIndex = bestFit.getXPeakIndex();
                
                limitError = histogram.getXErrors()[limitIndex];
                
                log.info("interpreting the density for a sparse background");
                
            } else {

                //limitStr = "top centroid";

                // centroid of area defined by the top portion of the fit or histogram where y >= ypeak/2
                float[] areaAndXYTopCentroid = calculateCentroidOfTop(bestFit.getOriginalScaleX(), bestFit.getOriginalScaleYFit(), 0.5f);

                float[] areaAndXYTopCentroid2 = calculateCentroidOfTop(histogram.getXHist(), histogram.getYHistFloat(), 0.5f);
                
                float xpeak = (areaAndXYTopCentroid[0] > 0) ? areaAndXYTopCentroid[1] : areaAndXYTopCentroid2[1];
                
                log.info("GEV top centroid=" + areaAndXYTopCentroid[1] + " data top centroid=" + areaAndXYTopCentroid2[1] + " using xpeak=" + xpeak);

                log.info("interpreting the density for a non-sparse background");
                
                if (xpeak > 0) {
                    
                    limit = xpeak;
                    
                    float areaAndXYTopCentroidError = (float)GeneralizedExtremeValue.calculateWidthFittingError(bestFit, 0.5f);
                    
                    limitError = areaAndXYTopCentroidError;
                    
                } else {
                    
                    limit = bestFit.getXPeak();
                    
                    float areaAndXYTopCentroid2Error = Histogram.calculateHistogramWidthYLimitError(
                        histogram.getXHist(), histogram.getYHistFloat(), histogram.getXErrors(), histogram.getYErrors(),
                        0.5f);
                    
                    limitError = areaAndXYTopCentroid2Error;
                }                                
            }

            this.backgroundDensity = limit;

            if (debug) {
                log.info(bestFit.toString());
            }            

            this.backgroundDensityError = limitError;
            
            if (debug) {

                // as comparison, log that roughly derived from chi square sum
                float gevTotalMeanFittingEmpiricalError = (float)Math.sqrt(
                    GeneralizedExtremeValue.calculateChiSq(bestFit.getX(), bestFit.getYFit()));

                // empirical error for one bin:
                float empiricalErrorInFitting = gevTotalMeanFittingEmpiricalError / histogram.getYHist().length;
                
                String interp = (interpretForSparseBackground) ? " sparse " : " complete ";
                log.info("\nestimating background from sampling= " + sampling.toString()
                    + " and interpretation=" + interp 
                    + "\ndens=" + this.backgroundDensity
                    + "\nw/ centroid error in area/y =" + backgroundDensityError
                    + "\ngev empirically estimated fitting error for one histogram bin =" + empiricalErrorInFitting
                );
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

        ArrayPair xy = LinesAndAngles.createPolygonOfTopFWFractionMax(xfit, yfit, frac);

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
