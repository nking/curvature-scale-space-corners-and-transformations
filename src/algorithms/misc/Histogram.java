package algorithms.misc;

import algorithms.util.Errors;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class Histogram {

    /**
     * create a histogram from the data that has little or no adjustment
     * for min and max.
     *
     * @param a
     * @param nBins
     * @param xHist
     * @param yHist
     */
    public static void createHistogram(float[] a, int nBins, float[] xHist, int[] yHist) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }

        float[] minMax = MiscMath.calculateOuterRoundedMinAndMax(a);

        createHistogram(a, nBins, minMax[0], minMax[1], xHist, yHist);
    }

    protected static float calculateBinWidth(float minValue, float maxValue, int nBins) {

        float xInterval = (maxValue - minValue)/nBins;

        // expand interval if necessary to make sure the last point is in the last bin
        if ((int) ((maxValue - minValue)/xInterval) != (nBins - 1)) {
            float t = (maxValue + minValue)/2.0f;
            int powDelta = MiscMath.findPowerOf10(t);
            float pow10 = (float)Math.pow(10, powDelta);
            xInterval = (maxValue - minValue + pow10)/nBins;
        }

        return xInterval;
    }

    /**
     * find the largest bin in which the width as determined by
     * calculateBinWidth is >= the matchBinWidth.  the return number
     * of bins will be no smaller than half of nBins given, which
     * means that the subsequent binWidth may not match the intended.
     *
     * @param minValue
     * @param maxValue
     * @param matchBinWidth
     * @param nBins
     * @return
     */
    protected static int calculateNBinsToMatchBinWidth(float minValue, float maxValue,
        float matchBinWidth, int nBins) {

        // the calculated binwidth has to be >= to matchBinWidth
        int n = nBins;

        float binWidth = calculateBinWidth(minValue, maxValue, n);

        while ( (binWidth < matchBinWidth) && (n >= ((nBins/2) + 1))) {
            n--;
            binWidth = calculateBinWidth(minValue, maxValue, n);
        }

        return n;
    }

    public static void createHistogram(float[] a, int nBins,
        float aMin, float aMax, float[] xHist, int[] yHist, float binWidth) {

        if (xHist == null || xHist.length != nBins) {
            throw new IllegalArgumentException("xHist has to be of size nBins and initialized");
        }
        if (yHist == null || yHist.length != nBins) {
            throw new IllegalArgumentException("yHist has to be of size nBins and initialized");
        }

        for (int i = 0; i < nBins; i++) {
            xHist[i] = aMin + i*binWidth + (binWidth/2.f);
        }

        for (int i = 0; i < a.length; i++) {
            int bin = (int) ((a[i] - aMin)/binWidth);
            if (bin < nBins) {
                yHist[bin]++;
            }
        }
    }

    public static void createHistogram(float[] a, int nBins,
        float aMin, float aMax, float[] xHist, int[] yHist) {

        if (xHist == null || xHist.length != nBins) {
            throw new IllegalArgumentException("xHist has to be of size nBins and initialized");
        }
        if (yHist == null || yHist.length != nBins) {
            throw new IllegalArgumentException("yHist has to be of size nBins and initialized");
        }

        Arrays.fill(yHist, 0);

        float xInterval = calculateBinWidth(aMin, aMax, nBins);

        createHistogram(a, nBins, aMin, aMax, xHist, yHist, xInterval);
    }

    /**
     * create a histogram adjusted to have a range
     * representing most of the data and a number of bins which attempts to have
     * counts per bin above 5.  The number of bins returned will be between
     * (preferredNumberOfBins)/2 and preferredNumberOfBins, inclusively.
     *
     * @param preferredNumberOfBins
     * @param values
     * @param valueErrors if null, the root mean square of values are used for the
     * point errors in order to calculate the histogram errors.
     * @param useAMinimumFromPeakResolution decide a minimum number of bins based upon
     *   having n bins across the peak profile. more specifically, n bins across
     *   the region x=0 to x=max(values) to x=max(values)/2, where n
     *   is min(10, 0.1*values.length)
     * @return
     */
    public static HistogramHolder createHistogramForSkewedData(int preferredNumberOfBins,
        float[] values, float[] valueErrors, boolean useAMinimumFromPeakResolution) {

        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        float[] minMax = MiscMath.calculateOuterRoundedMinAndMax(values);

        return createHistogramForSkewedData(preferredNumberOfBins, values, valueErrors,
            minMax[0], minMax[1], useAMinimumFromPeakResolution);
    }

    /**
     * create a histogram adjusted to have a range
     * representing most of the data and a number of bins which attempts to have
     * counts per bin above 5.  The number of bins returned will be between
     * (preferredNumberOfBins)/2 and preferredNumberOfBins, inclusively.
     *
     * @param preferredNumberOfBins
     * @param values
     * @param valueErrors if null, the root mean square of values are used for the
     * point errors in order to calculate the histogram errors.
     * @param minValue
     * @param maxValue maximum value
     * @param useAMinimumFromPeakResolution decide a minimum number of bins based upon
     *   having n bins across the peak profile. more specifically, n bins across
     *   the region x=0 to x=max(values) to x=max(values)/2, where n
     *   is min(10, 0.1*values.length)
     * @return
     */
    public static HistogramHolder createHistogramForSkewedData(int preferredNumberOfBins,
        float[] values, float[] valueErrors, float minValue, float maxValue, boolean useAMinimumFromPeakResolution) {

        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        if (useAMinimumFromPeakResolution) {
            return createHistogramForSkewedDataForPeakResolution(preferredNumberOfBins, values, valueErrors, minValue, maxValue);
        }

        float binWidth = calculateSturgesBinWidth(values, 0.25f);

        //float binWidth = calculateFreedmanDraconisBinWidth(values, 0.85f, maxValue);

        int nBins = (int)((maxValue - minValue)/binWidth);
        if (nBins < 5) {
            nBins = preferredNumberOfBins;
        }

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];

        // if there are many bins with counts < 5, should  reduce the number of bins, but not to less than half preferred
        boolean go = true;

        int lowerLimit = (int)(preferredNumberOfBins/2.5f);
        if (lowerLimit < 5) {
            lowerLimit = preferredNumberOfBins;
        }

        while (go) {
            int lessThan5Counts = 0;
            for (int i = 0; i < yHist.length; i++) {
                if (yHist[i] < 5) {
                    lessThan5Counts ++;
                }
            }
            int chk = (int)(0.25*nBins);
            if (lessThan5Counts >= chk) {
                int tNBins = (int)(nBins/1.5f);
                if (tNBins >= lowerLimit) {

                    nBins = tNBins;

                    xHist = new float[nBins];
                    yHist = new int[nBins];
                    Histogram.createHistogram(values, nBins, minValue, maxValue, xHist, yHist, binWidth);

                } else {
                    go = false;
                }
            } else {
                go = false;
            }
        }

        // change bin size so that preferredNumberOfBins is used.
        float factor = nBins/preferredNumberOfBins;
        if (factor > 1) {
            nBins = preferredNumberOfBins;
            binWidth *= factor;
        }
        xHist = new float[nBins];
        yHist = new int[nBins];
        Histogram.createHistogram(values, nBins, minValue, maxValue, xHist, yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        if (valueErrors == null) {
            valueErrors = Errors.populateYErrorsBySqrt(yHistFloat);
        }
        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);

        return histogram;
    }

    /**
     * create a histogram adjusted to have a range
     * representing most of the data and a number of bins which attempts to have
     * counts per bin above 5.  The number of bins returned will be between
     * (preferredNumberOfBins)/2 and preferredNumberOfBins, inclusively.
     *
     * decide a minimum number of bins based upon
     *   having n bins across the peak profile. more specifically, n bins across
     *   the region x=0 to x=max(values) to x=max(values)/2, where n
     *   is min(10, 0.1*values.length)
     *
     * @param preferredNumberOfBins
     * @param values
     * @param valueErrors if null, the root mean square of values are used for the
     * point errors in order to calculate the histogram errors.
     * @param minValue
     * @param maxValue maximum value
     * @return
     */
    public static HistogramHolder createHistogramForSkewedDataForPeakResolution(int preferredNumberOfBins,
        float[] values, float[] valueErrors, float minValue, float maxValue) {

        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        // -- start with a histogram based upon sturges or draconis.
        // -- then reduce the bin width if necessary such that there are 10 or
        //    so bins across the width of half max of the profile.
        // -- reduce nbins until the last bins hold values > 5

        float binWidth = calculateSturgesBinWidth(values, 0.25f);

        //float binWidth = calculateFreedmanDraconisBinWidth(values, 0.85f, maxValue);

        int nBins = (int)((maxValue - minValue)/binWidth);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];

        float xMax = maxValue;

        Histogram.createHistogram(values, nBins, minValue, maxValue, xHist, yHist, binWidth);

        // calculate a minimum number of bins based upon resolving the peak if useAMinimumFromPeakResolution=true
        float peakResBinWidth = calculateBinWidthForPeakResolution(xHist, yHist);

        if (binWidth > peakResBinWidth) {
            binWidth = peakResBinWidth;
        }

        nBins = (int)((xMax - minValue)/binWidth);

        int minNumberOfBins = preferredNumberOfBins/2;

        // if there are many bins with counts < 5, should  reduce the number of bins, but not to less than half preferred
        boolean go = true;

        int yPeakIndex = MiscMath.findYMaxIndex(yHist);

        while (go) {
            int lessThan5Counts = 0;
            int indexStartingConsecutiveLowCounts = -1;
            int prevIndexStartingConsecutiveLowCounts = yHist.length;
            for (int i = (yHist.length - 1); i > -1; i--) {
                if (yHist[i] < 5) {

                    lessThan5Counts++;

                    // avoid using points with x < x at ypeak
                    if ( (i > yPeakIndex) &&
                        ((i == (prevIndexStartingConsecutiveLowCounts - 1)) ||
                        (i == (prevIndexStartingConsecutiveLowCounts - 2)))) {

                        indexStartingConsecutiveLowCounts = i;
                        prevIndexStartingConsecutiveLowCounts = indexStartingConsecutiveLowCounts;
                    } else {
                        prevIndexStartingConsecutiveLowCounts = i;
                    }
                }
            }
            int chk = (int)(0.1*nBins);
            chk = 4;

            if (indexStartingConsecutiveLowCounts > -1) {

                if ((indexStartingConsecutiveLowCounts + 1) < yHist.length) {
                    indexStartingConsecutiveLowCounts++;
                }

                // reduce xMax
                xMax = xHist[indexStartingConsecutiveLowCounts];

                // only the histograms with nBins > min should be repeated to remove empty bins at tail
                nBins = (int)((xMax - minValue)/binWidth);
                if (nBins < minNumberOfBins) {
                    nBins = minNumberOfBins;
                    binWidth = (xMax - minValue)/nBins;
                }

                xHist = new float[nBins];
                yHist = new int[nBins];
                Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);

                go = false;

            } else if (lessThan5Counts >= chk) {

                int tNBins = (int)(nBins/1.5f);
                if (tNBins >= minNumberOfBins) {

                    // reduce the number of bins
                    nBins = tNBins;

                    xMax = (nBins*binWidth) + minValue;

                    xHist = new float[nBins];
                    yHist = new int[nBins];
                    Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);

                } else {
                    go = false;
                }
            } else {
                go = false;
            }
        }

        // change bin size so that preferredNumberOfBins or minNumberOfBins is used.
        float factor = (float)nBins/preferredNumberOfBins;
        if (factor > 1) {
            nBins = preferredNumberOfBins;
            binWidth *= factor;
        }
        xHist = new float[nBins];
        yHist = new int[nBins];
        Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        if (valueErrors == null) {
            valueErrors = Errors.populateYErrorsBySqrt(yHistFloat);
        }
        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);

        return histogram;
    }

    protected static float calculateSturgesBinWidth(float[] values, float binWidthFactor) {

        if (values == null ) {
            throw new IllegalArgumentException("values cannot be null and must be the same length");
        }
        /*
        Sturges:
                               stdev of the n data values
           bin width = 3.49 * ----------------------------
                                        n^(1/3)
        */
        float[] meanStdev = MiscMath.findMeanAndStDev(values);

        float binWidth = (float) (binWidthFactor*3.49f* meanStdev[1]/Math.pow(values.length, (1./3.)));

        return binWidth;
    }

    public static float calculateFreedmanDraconisBinWidth(float[] values, float binWidthFactor) {

        if (values == null) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        /*
        Freedman-Diaconis:
           IQR is interquartile range
           Q1 = 1/4 of data
           Q3 = 3/4 of data
           IQR = Q3 - Q1
           fdBinWidth = 2*( xHist[q3Index] - xHist[q1Index])
               / (float)(Math.pow(values.length, (1.0/3.0)));
        */
        // use Freedman-Diaconis choice to see match with this distr

        float[] tmp = Arrays.copyOf(values, values.length);
        Arrays.sort(tmp);
        float q1 = tmp[ values.length/4];
        float q2 = tmp[ values.length/2];
        float q3 = tmp[ 3*values.length/4];

        float binWidth = (float) (binWidthFactor*2.0f*(q3-q1)/Math.pow(values.length, (1./3.)));

        return binWidth;
    }

    public static float calculateFreedmanDraconisBinWidth(float[] values, float binWidthFactor, float maxValue) {

        if (values == null) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        /*
        Freedman-Diaconis:
           IQR is interquartile range
           Q1 = 1/4 of data
           Q3 = 3/4 of data
           IQR = Q3 - Q1
           fdBinWidth = 2*( xHist[q3Index] - xHist[q1Index])
               / (float)(Math.pow(values.length, (1.0/3.0)));
        */
        // use Freedman-Diaconis choice to see match with this distr

        float[] tmp = Arrays.copyOf(values, values.length);
        Arrays.sort(tmp);
        for (int i = 0; i < tmp.length; i++) {
            if (tmp[i] > maxValue) {
                tmp = Arrays.copyOf(tmp, i);
                break;
            }
        }

        float q1 = tmp[ values.length/4];
        float q2 = tmp[ values.length/2];
        float q3 = tmp[ 3*values.length/4];

        float binWidth = (float) (binWidthFactor*2.0f*(q3-q1)/Math.pow(values.length, (1./3.)));

        return binWidth;
    }

    /**
     * Populate the arrays xHistErrorsOutput and yHistErrorsOutput with errors
     * for the x and y bins calculated from valueErrors.
     *
     *
     * Trying to adjust errors to use the fact that a bin with a value of zero is
     * actually determined from all bins, so the error in a bin with a value of
     * zero is a positive number less than infinity and is due to all errors of
     * points that go into making the histogram.
     *
     * For each bin:
     *   (sigma)^2 = (ave sigma_from_all)^2 + (sigma from all points in bin, in function that reduces it by the number of points)
     *
     * The total of all bins should be equal to the sigma_from_all;
     *
     * N = 3 bins:
     *    [0] (sigma_bin0)^2 = (sigma_from_all/nBins)^2 + (sigma_bin0 * Fnc(bin0))^2
     *    [1] (sigma_bin1)^2 = (sigma_from_all/nBins)^2 + (sigma_bin1 * Fnc(bin1))^2
     *    [2] (sigma_bin2)^2 = (sigma_from_all/nBins)^2 + (sigma_bin2 * Fnc(bin2))^2
     * Where (sigma_bin0)^2 is the sum of all points that went into bin0 added in quadrature.
     *
     * The total of all 3 should equal the total from all points added in quadrature and no more.
     *
     *    (sigma_from_all_points)^2 = nBins*(sigma_from_all/nBins)^2 + (sigma_bin0 * Fnc(bin0))^2
     *                                + (sigma_bin1 * Fnc(bin1))^2 + (sigma_bin2 * Fnc(bin2))^2
     *
     *    (sigma_from_all_points)^2 * (1 - (1/nBins)) = (sigma_bin0 * Fnc(bin0))^2 + (sigma_bin1 * Fnc(bin1))^2 + (sigma_bin2 * Fnc(bin2))^2
     *
     * Not easy to see how to solve, so returning to first principles:
     *     Each bin should have a contribution from all points and then a contribution from its own.
     *     In order for the total to not exceed the sum in quadrature of all points, a bin's sigma squared should
     *     be the ave from all + it's own times (1-1/N) roughly, which is what the equation above suggests.
     *
     * @param xHist
     * @param yHist
     * @param values
     * @param valueErrors errors that are on the same scale as the values, that is, these are
     *   NOT percent errors
     * @param xHistErrorsOutput
     * @param yHistErrorsOutput
     */
    public static void calulateHistogramBinErrors(float[] xHist, int[] yHist,
        float[] values, float[] valueErrors, float[] xHistErrorsOutput, float[] yHistErrorsOutput) {

        float xInterval = xHist[1] - xHist[0];
        float xmin = xHist[0] - (xInterval/2.0f);

        float[] sumErrorPerBin = new float[xHist.length];

        float sumErrorAllPoints = 0;
        //float sumPercentErrorAllPoints = 0;

        for (int i = 0; i < valueErrors.length; i++) {

            int bin = (int) ((values[i] - xmin)/xInterval);

            // in units of y
            float a = valueErrors[i];
            a *= a;

            if (bin < xHist.length) {

                sumErrorPerBin[bin] += a;

                sumErrorAllPoints += a;

                //float b = valueErrors[i]/values[i];
                //sumPercentErrorAllPoints += (b*b);
            }
        }

        // the x value was determined from all points, so the error should be taken
        //   as the average error of all points
        float aveErrorOfAllPoints = (float)Math.sqrt(sumErrorAllPoints)/yHist.length;

        // the percent errors, divided over all bins, which is what was done to learn the binwdith
        //float avePercentErrorOfAllPoints = (float)Math.sqrt(sumPercentErrorAllPoints)/yHist.length;

        float sumAlt = 0;
        float sumOrigSquared = 0;
        for (int i = 0; i < yHist.length; i++) {

            float sumErrorBinSquared = sumErrorPerBin[i];

            // compare to aveErrorOfAllPoints... should be similar
            //float be = (float) Math.sqrt(sumErrorBinSquared);

            float c = aveErrorOfAllPoints;

            xHistErrorsOutput[i] = c;

            //float a = (float)Math.sqrt(sumErrorBinSquared * (1.0f - (1.0f/yHist.length));
            float ai = sumErrorBinSquared;
            float af = sumErrorBinSquared/yHist.length;
            float a = (float)Math.sqrt(ai - af);

            float yBinError = (yHist[i] == 0) ? 0 : a/yHist[i];

            yHistErrorsOutput[i] = yBinError;

            sumAlt += yBinError;


            float yBinErrorOrig = (yHist[i] == 0) ? 0 : (float)Math.sqrt(sumErrorPerBin[i])/yHist[i];

            sumOrigSquared += yBinErrorOrig;
        }

        float contribFromAllToEachBin = (sumOrigSquared - sumAlt)/yHist.length;

        for (int i = 0; i < yHist.length; i++) {

            yHistErrorsOutput[i] += contribFromAllToEachBin;
        }
    }

    public static float calculateHistogramAreaError(float[] xHist, int[] yHist,
        float[] xErrors, float[] yErrors) {

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
         *      if f = XY, and X and Y are not correlated, we have:
         *          sigma^2  =  xError^2*(Y^2)  +  yError^2*(X^2) + xError^2*yError^2 <==== covariance term correct?
         */
        float binSize = xHist[1] - xHist[0];
        float sum = 0.0f;
        for (int i = 0; i < xHist.length; i++) {

            float xe = xErrors[i];
            float yeSquare = yErrors[i]*yErrors[i];

            float sigmaDfDxSquared = (xe*xe)*(yHist[i] * yHist[i]);
            float sigmaDfDySquared = yeSquare*(binSize * binSize);

            float covXY = (xe * xe) * yeSquare;

            float s = sigmaDfDxSquared + sigmaDfDySquared + covXY;

            sum += s;
        }

        sum = (float) Math.sqrt(sum);

        return sum;
    }

    /**
     *   *  |  *
     *      |
     * -----|.....  <---- angle is w.r.t y=0, x=xc.  increases in CCW order
     *  *   |
     *      |  *
     *
     * @param xc
     * @param yc
     * @param radius
     * @param angleInDegreesFromYEQ0XGT0  angle in degrees, CCW from point y=0, x=xc
     * @return
     */
    static float[] calculateXAndYFromXcYcAndRadiusCCW(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        return calculateXAndYFromXcYcAndRadius(xc, yc, radius, 360 - angleInDegreesFromYEQ0XGT0);
    }
    public static float[] calculateXAndYFromXcYcAndRadius(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        double dx = radius * Math.cos(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));
        double dy = radius * Math.sin(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));

        float x = (float) (xc + dx);
        float y = (float) (yc - dy);
        return new float[]{x, y};
    }

    /**
     * decide a minimum number of bins based upon
     *   having n bins across the peak profile. more specifically, n bins across
     *   the region x=0 to x=max(values) to x=max(values)/2, where n
     *   is min(10, 0.1*values.length)
     *
     * @param xHist
     * @param yHist
     * @return
     */
    public static float calculateBinWidthForPeakResolution(float[] xHist, int[] yHist) {

        int nBins;
        if ( (int)(0.1*yHist.length) < 10 ) {
            if ( (int)(0.1*yHist.length) < 3 ) {
                nBins = 3;
            } else {
                nBins = (int)(0.1*yHist.length);
            }
        } else {
            nBins = 10;
        }

        // find value max(values)/2 beyond peak
        int maxValueIndex = MiscMath.findYMaxIndex(yHist);

        float halfMaxValue = yHist[maxValueIndex]/2.0f;
        int halfMaxValueIndex = -1;

        for (int i = maxValueIndex; i < yHist.length; i++) {

            float value = yHist[i];

            if (value > halfMaxValue) {

                halfMaxValueIndex = i;

            } else {
                break;
            }
        }

        float binWidth = (xHist[halfMaxValueIndex])/nBins;

        return binWidth;
    }
}
