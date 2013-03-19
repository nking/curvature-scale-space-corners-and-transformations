package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
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
import java.util.Arrays;
import java.util.Locale;

/**
 * Class to estimate a background surface density for a set of points in which
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

    protected boolean useCompleteSampling = false;
    protected boolean allowTuningForThresholdDensity = true;

    protected float[] allTwoPointSurfaceDensities = null;
    protected float[] allTwoPointSurfaceDensitiesErrors = null;
    protected int nTwoPointSurfaceDensities = 0;
    protected int[] point1 = null;
    protected int[] point2 = null;
    protected StringArrayLite twoPointIdentities = null;

    protected int defaultNBins = 20;

    //for debugging, hold on to intermediate data
    protected HistogramHolder statsHistogram = null;
    protected GEVYFit bestFit = null;

    protected boolean adjustMuForDensity = true;

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
    }

    public TwoPointVoidStats(String persistedIndexerFileName) throws IOException {
        super(persistedIndexerFileName);
    }

    /**
     * Use complete sampling, else the default sampling.
     * the runtime for the complete sampling is N!/(2!(N-2)! * N!/(2!(N-2) which
     * is effectively N^4.  The runtime for the default sampling is N^2 - N.
     *
     * Complete sampling should be used and is used when there are less than 100 points.
     *
     * @param doUseCompleteSampling
     */
    public void setUseCompleteSampling(boolean doUseCompleteSampling) {

        if (indexer.nXY < 100) {

            useCompleteSampling = true;

        } else {

            useCompleteSampling = doUseCompleteSampling;
        }
    }

    public void setAllowTuningForThresholdDensity(boolean allowTuning) {
        this.allowTuningForThresholdDensity = allowTuning;
    }

    public void setGEVRangeParameters(float kMin, float kMax, float sigmaMin, float sigmaMax) {
        this.gevRangeFittingParameters = new float[]{kMin, kMax, sigmaMin, sigmaMax};
    }

    /**
     * calculate the 2-point void surface densities and then calculate the
     * statistics of those points.
     *
     * @throws TwoPointVoidStatsException
     */
    public void calc() throws TwoPointVoidStatsException {

        calculateSurfaceDensities();

        twoPointIdentities = null;

        calculateStats();
    }

    protected void calculateSurfaceDensities() throws TwoPointVoidStatsException {

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

    protected void calculateStats() throws TwoPointVoidStatsException {

        //HistogramHolder histogram = Histogram.createHistogramForSkewedData(20, allTwoPointSurfaceDensities,
        //    allTwoPointSurfaceDensitiesErrors);

        int nBins = (indexer.getNumberOfPoints() < 100) ? 10 : defaultNBins;

        HistogramHolder histogram = Histogram.createHistogramForSkewedData(
            nBins, allTwoPointSurfaceDensities, allTwoPointSurfaceDensitiesErrors);

        statsHistogram = histogram;

        plotPairSeparations();

        //int yMaxBin = Histogram.findMax(histogram.getYHist());

        int yMaxBin = -1;
        float ymax = Float.MIN_VALUE;
        for (int i = 0; i < (histogram.getYHist().length/2.); i++) {
            if (histogram.getYHist()[i] > ymax) {
                ymax = histogram.getYHist()[i];
                yMaxBin = i;
            }
        }

        calculateStatsForBackground(histogram, yMaxBin);
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

        statsHistogram = histogram;

        int muBin = yMaxBin;

        try {

            GEVChiSquareMinimization chiSqMin = new GEVChiSquareMinimization(
                histogram.getXHist(), histogram.getYHistFloat(), histogram.getXErrors(), histogram.getYErrors());

            chiSqMin.setDebug(debug);

            GEVYFit yfit = null;

            float kMin;
            float kMax;
            float sigmaMin;
            float sigmaMax;

            if (gevRangeFittingParameters != null) {
                kMin = gevRangeFittingParameters[0];
                kMax = gevRangeFittingParameters[1];
                sigmaMin = gevRangeFittingParameters[2];
                sigmaMax = gevRangeFittingParameters[3];
            } else {
                kMin = GEVChiSquareMinimization.kMinDefault;
                kMax = GEVChiSquareMinimization.kMaxDefault;
                sigmaMin = GEVChiSquareMinimization.sigmaMinDefault;
                sigmaMax = GEVChiSquareMinimization.sigmaMaxDefault;
            }

            yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(
                GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS,
                kMin, kMax, sigmaMin, sigmaMax
            );

            if (yfit == null) {
                if (debug) {
                    plotPairSeparations();
                }
                throw new TwoPointVoidStatsException("histogram of linear densities was not fittable");
            }

            boolean didCalculateFitStats = yfit.calculateStatsAsStepFunction();
            if (!didCalculateFitStats) {
                if (debug) {
                    plotPairSeparations();
                }
                throw new TwoPointVoidStatsException("histogram of linear densities was not fittable");
            }

            float chsqdiverr = yfit.getChiSqSum() / yfit.getYDataErrSq();
            float chsqstatdiverr = yfit.getChiSqStatistic() / yfit.getYDataErrSq();

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

            bestFit = yfit;

            if (debug) {
                if (bestFit != null) {
                    System.out.println(bestFit.toString());
                }
                float xHalfInterval = (histogram.getXHist()[1] - histogram.getXHist()[0]) / 2.0f;
                float xmin = 0;
                float xmax = histogram.getXHist()[histogram.getXHist().length - 1] + xHalfInterval;
                float ymin = 0;
                float ymax = MiscMath.findMax(histogram.getYHistFloat());

                float[] xf = yfit.getOriginalScaleX();
                float[] yf = yfit.getOriginalScaleYFit();
                PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
                plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                    histogram.getXErrors(), histogram.getYErrors(), xf, yf, "");
                plotter.writeFile3();
            }


            boolean increaseThreshold = false;
            boolean increaseThreshold2 = false;
            boolean decreaseThreshold = false;

            float medDivMean = yfit.getXMedian()/yfit.getXMean();
            float chistat = yfit.getChiSqStatistic();

            if (allowTuningForThresholdDensity) {
                if (chistat > 100000) {
                    increaseThreshold2 = true;
                } else if ((chistat > 2) && (medDivMean) > 15.0f) {
                    increaseThreshold2 = true;
                } else if ((chistat > 2) && (yfit.getOriginalScaleYFit()[yfit.getXPeakIndex()] < 100)) {
                    decreaseThreshold = true;
                } else if ((chistat > 2) && (medDivMean < 4.0f)) {
                    increaseThreshold2 = true;
                }
            }

            // Note, the following portion could be altered or improved for different uses.  One could fit more than one
            //   component, etc.
            // x% works well for sparsely populated backgrounds (number of background points < points in clusters)
            // and moderately populated backgrounds (number of background points approx equal to points in clusters)

            String limitStr = "10%";
            float limit, limitError;
            limit = yfit.getX10Percent();
            int limitIndex = yfit.getX10PercentIndex();
            limitError = histogram.getXErrors()[limitIndex];

            if (increaseThreshold) {
                limitStr = "peak+0.1";
                float xPeak = yfit.getXPeak();
                int xPeakIndex = yfit.getXPeakIndex();
                if ((xPeakIndex + 1) < yfit.getX().length) {
                    float delta = (yfit.getOriginalScaleX()[xPeakIndex + 1] - xPeak)/10.0f;
                    limit = xPeak + delta;
                    limitError = histogram.getXErrors()[xPeakIndex];
                } else {
                    limit = xPeak;
                    limitError = histogram.getXErrors()[xPeakIndex];
                }
            } else if (increaseThreshold2) {
                limitStr = "peak+0.25";
                float xPeak = yfit.getXPeak();
                int xPeakIndex = yfit.getXPeakIndex();
                if ((xPeakIndex + 1) < yfit.getX().length) {
                    float delta = (yfit.getOriginalScaleX()[xPeakIndex + 1] - xPeak)/4.0f;
                    limit = xPeak + delta;
                    limitError = histogram.getXErrors()[xPeakIndex];
                } else {
                    limit = xPeak;
                    limitError = histogram.getXErrors()[xPeakIndex];
                }
            } else if (decreaseThreshold) {
                limitStr = "0.5*10%";
                limit = 0.5f*yfit.getX10Percent();
            }

            //if (limitIndex == -1) {
            //    System.out.println("WARNING:  solution was not found");
            //    return;
            //}

            this.backgroundSurfaceDensity = limit;

            if (debug) {
                System.out.println(yfit.toString());
            }

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
                GeneralizedExtremeValue.calculateTotalMeanFittingError(yfit.getX(), yfit.getYFit());

            float errorInFitting = gevTotalMeanFittingError / histogram.getYHist().length;

            // add errors in quadrature
            this.backgroundSurfaceDensityError = (float) (Math.sqrt(
                errorInEstimateFromHistogram * errorInEstimateFromHistogram + errorInFitting * errorInFitting));

            if (debug) {
                System.out.println("estimating background minimum as location of " + limitStr + " of background profile counts = "
                + this.backgroundSurfaceDensity
                + " w/ x error in histogram bin =" + errorInEstimateFromHistogram
                + " gev fitting error for one point =" + errorInFitting + " combined error ===> "
                + this.backgroundSurfaceDensityError);
            }

            if (debug && (yfit.getChiSqSum() > chiSqMin.calcYErrSquareSum())) {
                System.out.println("WARNING:  chisq is larger than errors: "
                    + yfit.getChiSqSum() + " (errsqsum=" + chiSqMin.calcYErrSquareSum() + ")");
            }

        } catch (FailedToConvergeException e) {
            throw new TwoPointVoidStatsException(e);
        } catch (IOException e2) {
            throw new TwoPointVoidStatsException(e2);
        }
    }

    /**
     * find the rectangular 2 point surface densities using either complete
     * sampling or the default partial sampling via a range division
     * (default).
     *
     * The runtime for the complete sampling uses N!/(2!(N-2)! * N!/(2!(N-2) iterations
     *      N    number of iterations
     *      10       2   E3
     *    1000       2.5 E11
     * 1000000       2.5 E23
     *
     * The default method uses N^2 - N for number of iterations
     *      N    number of iterations
     *      10       2.8 E1
     *    1000       3.7 E5
     * 1000000       3.7 E11
     *
     * Roughly, the complete method runtime is proportional to N^4 and the incomplete
     * range sampling is proportional to N^2.
     * Simulations produce similar number density distributions within errors for
     * the smaller sample, so it's fine to use the default method.
     */
    protected void findVoids() {

        if (useCompleteSampling) {

            float[] x = indexer.getX();
            float[] y = indexer.getY();
                                                                                // cost     number of times
            for (int i = 0; i < indexer.nXY; i++) {
                System.out.println("findVoids i=" + i + "/" + indexer.nXY);
                for (int ii = (i + 1); ii < indexer.nXY; ii++) {
                    for (int j = 0; j < indexer.nXY; j++) {
                        for (int jj = (j + 1); jj < indexer.nXY; jj++) {        //           N!/(2!(N-2)! * N!/(2!(N-2)!
/*
System.out.println(" xsi=[" + i + ":" + ii + "] "
 + " ysi=[" + j + ":" + jj + "]  ==> "
 + " xi=[" + indexer.sortedXIndexes[i] + ":" + indexer.sortedXIndexes[ii] + "]"
 + " yi=[" + indexer.sortedYIndexes[j] + ":" + indexer.sortedYIndexes[jj] + "]"
 );
*/
                            processIndexedRegion(x, y, i, ii, j, jj);
                        }
                    }
                }
            }

        } else {

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
            float[] x = indexer.getX();
            float[] y = indexer.getY();
            int n = indexer.getNumberOfPoints();

            int nYIntervals = n / 2;                                            // cost     number of times

            for (int k = 0; k < nYIntervals; k++) {                             //               2

                int binSz = (k + 1) * 2;                                        // c10

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

        int[] regionIndexes = indexer.findIntersectingRegionIndexesIfOnlyTwo(xIndexLo, xIndexHi, yIndexLo, yIndexHi, useCompleteSampling);

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
                StringBuffer err = new StringBuffer("ERROR: intersecting region contained more than 2 points, but unique y wasn't found");
                for (int i = 0; i < regionIndexes.length; i++) {
                    float ypi = y[ regionIndexes[i] ];
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

        oos.flush();
    }

    protected String serializeTwoPointDensities() throws IOException {

        return serializeTwoPointBackground("stats_2pt_voids_");
    }

    public boolean readTwoPointBackground(String persistedFileName) throws IOException {

        boolean didDeserialize = deserializeTwoPointBackground(persistedFileName);;

        if (nTwoPointSurfaceDensities == 0) {
            throw new IOException("No pairs were found isolated within an area");
        }

        return didDeserialize;
    }

    protected void deserializeTwoPointBackground(ObjectInputStream ois) throws IOException {

        nTwoPointSurfaceDensities = ois.readInt();

        allTwoPointSurfaceDensities = new float[nTwoPointSurfaceDensities];
        point1 = new int[nTwoPointSurfaceDensities];
        point2 = new int[nTwoPointSurfaceDensities];

        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            allTwoPointSurfaceDensities[i] = ois.readFloat();
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            point1[i] = ois.readInt();
            point2[i] = ois.readInt();
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

        if (this.statsHistogram != null) {

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
