package algorithms.misc;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 *  TODO:  improve this class...
 *
 * @author nichole
 */
public class Histogram {

    protected static Logger log = Logger.getLogger(Histogram.class.getName());

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

        Arrays.fill(yHist, 0);

        for (int i = 0; i < nBins; i++) {
            xHist[i] = aMin + i*binWidth + (binWidth/2.f);
        }

        for (int i = 0; i < a.length; i++) {
            int bin = (int) ((a[i] - aMin)/binWidth);
            if ((bin > -1) && (bin < nBins)) {
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

        /*if (useAMinimumFromPeakResolution) {
            return createHistogramForSkewedDataForPeakResolution(preferredNumberOfBins, values, valueErrors, minValue, maxValue);
            //return createHistogramForSkewedDataForPeakResolution2(preferredNumberOfBins, values, valueErrors, minValue, maxValue, 1);
        }*/

        return createHistogramForSkewedDataIgnorePeakResolution(values, valueErrors, minValue, maxValue);
    }

    /**
     * this is the histogram used for the small number histograms, that is for values.length < 1000
     *
     * @param preferredNumberOfBins
     * @param values
     * @param valueErrors
     * @param minValue min x histogram value
     * @param maxValue max x histogram value
     * @return
     */
    private static HistogramHolder createHistogramForSkewedDataIgnorePeakResolution(
        float[] values, float[] valueErrors, float minValue, float maxValue) {

        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        /*
str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/
// Sturges   good for very low numbers
        int nIntervalsSturges = (int)Math.ceil( Math.log(values.length)/Math.log(2) );

        int nBins = Math.max(nIntervalsSturges, 10);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        float minx = minValue;
        float maxx = maxValue;

        float binWidth = calculateBinWidth(minx, maxx, nBins);

        Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);

        boolean go = true;

        int nIter = 0;
        
        boolean didFlatten20 = false;
        boolean didFlatten20To1 = false;
        boolean didFlatten1 = false;

        while (go && (nIter < 100) && (minx < maxx)) {

            int yPeakIndex = MiscMath.findYMaxIndex(yHist);
            
            float yHistMax = yHist[yPeakIndex];
            int minCounts = (yHistMax < 50) ? 3 : 20;
            
            // correct for the over-correcting that reduces the histogram to all zeros and a single peak
            if (yHistMax <= 20.0f && !didFlatten20 && !didFlatten1) {
                minx = minValue;
                maxx = maxValue;
                binWidth = calculateBinWidth(minx, maxx, nBins);
                Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
                if (yHistMax > 1.0f) {
                    didFlatten20 = true;
                    minCounts = 5;
                } else {
                    didFlatten1 = true;
                    minCounts = 10;
                }
            } else if (didFlatten20) {
                if (yHistMax == 1.0f) {
                    minx = minValue;
                    maxx = maxValue;
                    binWidth = calculateBinWidth(minx, maxx, nBins);
                    Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
                    didFlatten20 = false;
                    didFlatten20To1 = true;
                    minCounts = 10;
                } else {
                    minCounts = 5;
                }
            } else if (didFlatten1) {
                minCounts = 10;
            } else if (didFlatten20To1) {
                minCounts = 10;
            }

            int startsWithLessThanMinCounts = 0;

            int endsWithLessThanMinCounts = 0;
            int indexStartingConsecutiveEndLowCounts = -1;
            int prevIndexStartingConsecutiveEndLowCounts = yHist.length;

            for (int i = (yHist.length - 1); i > -1; i--) {

                if (yHist[i] < minCounts) {

                    endsWithLessThanMinCounts++;
                    // avoid using points with x < x at ypeak
                    if ( (i > yPeakIndex) && ((i == (prevIndexStartingConsecutiveEndLowCounts - 1)) ||
                    (i == (prevIndexStartingConsecutiveEndLowCounts - 2)))) {
                        
                        indexStartingConsecutiveEndLowCounts = i;
                        
                        prevIndexStartingConsecutiveEndLowCounts = indexStartingConsecutiveEndLowCounts;
                    
                    } else {
                        prevIndexStartingConsecutiveEndLowCounts = i;
                    }

                    // count the consecutive low counts at start of histogram
                    if (i > 0) {
                        // if there is a low count in i-1 too, increment this one, else reset it to zero (no skips)
                        if (yHist[i-1] < minCounts) {
                            startsWithLessThanMinCounts++;
                        } else {
                            startsWithLessThanMinCounts = 0;
                        }
                    }
                }
            }

            if (startsWithLessThanMinCounts == 0 && endsWithLessThanMinCounts == 0) {
                go = false;
            } else {
                if (endsWithLessThanMinCounts > 0) {
                    maxx -= 1.01*binWidth;
                }
                if (startsWithLessThanMinCounts > 0) {
                    minx += 1.01*binWidth;
                }

                binWidth = calculateBinWidth(minx, maxx, nBins);

                Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);

                nIter++;
            }
        }

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
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

        int minNumberOfBins = (preferredNumberOfBins > 30) ? 20 : ((preferredNumberOfBins > 20) ? 10 : preferredNumberOfBins/2);

        if (binWidth > peakResBinWidth) {
            binWidth = peakResBinWidth;
            nBins = (int)((xMax - minValue)/binWidth);
        }

        // if there are many bins with counts < 5, should  reduce the number of bins, but not to less than half preferred.
        // also, if there are many empty bins at start of histogram, increase the bin size.
        boolean go = true;

/*
algorithms.util.PolygonAndPointPlotter plotter = null;
try {
plotter = new algorithms.util.PolygonAndPointPlotter();
} catch (Exception e) {}
String str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/

        while (go) {

            int yPeakIndex = MiscMath.findYMaxIndex(yHist);

            int startsWithLessThan5Counts = 0;
            int indexStartingConsecutiveStartLowCounts = -1;

            int endsWithLessThan5Counts = 0;
            int indexStartingConsecutiveEndLowCounts = -1;
            int prevIndexStartingConsecutiveEndLowCounts = yHist.length;
            for (int i = (yHist.length - 1); i > -1; i--) {

                if (yHist[i] < 5) {

                    // count the consecutive low counts at the end of the histogram:
                    endsWithLessThan5Counts++;
                    // avoid using points with x < x at ypeak
                    if ( (i > yPeakIndex) && ((i == (prevIndexStartingConsecutiveEndLowCounts - 1)) ||
                    (i == (prevIndexStartingConsecutiveEndLowCounts - 2)))) {
                        indexStartingConsecutiveEndLowCounts = i;
                        prevIndexStartingConsecutiveEndLowCounts = indexStartingConsecutiveEndLowCounts;
                    } else {
                        prevIndexStartingConsecutiveEndLowCounts = i;
                    }

                    // count the consecutive low counts at start of histogram
                    if (i > 0) {
                        // if there is a low count in i-1 too, increment this one, else reset it to zero (no skips)
                        if (yHist[i-1] < 5) {
                            startsWithLessThan5Counts++;
                            indexStartingConsecutiveStartLowCounts = i-1;
                        } else {
                            startsWithLessThan5Counts = 0;
                            indexStartingConsecutiveStartLowCounts = -1;
                        }
                    }
                }
            }

            /*
             * Handle cases:
             *
             * (0) zeros at start  && zeros at end
             *     -- reduce xmax
             *        --> reduce nBins
             *
             * (1) zeros at start  only
             *        --> reduce nBins
             *
             * (2) zeros at end only
             *     -- reduce xmax
             *        --> reduce nBins
             *
             * Note, don't want cases to undo each others work upon iteration, so chose
             * to reduce the number of bins (rather than increase and decrease the bin width).
             */

            int chk = (int)(0.1*nBins);
            chk = 4;

            if ((indexStartingConsecutiveStartLowCounts == 0) && (indexStartingConsecutiveEndLowCounts > -1)) {
                // case 0

                // reduce xMax
                float tmpXMax = xHist[indexStartingConsecutiveEndLowCounts];;

                int tmpNBins = (int)((tmpXMax - minValue)/binWidth);

                if (tmpNBins >= minNumberOfBins) {

                    xMax = tmpXMax;

                    nBins = tmpNBins;

                    //binWidth = calculateBinWidth(minValue, xMax, nBins);
                } else {
                    go = false;
                }

            } else if (indexStartingConsecutiveStartLowCounts == 0) {
                // case 1

                if (startsWithLessThan5Counts > 2) {

                    // reduce nBins
                    float firstNonZero = xHist[startsWithLessThan5Counts];
                    float minBinWidth = firstNonZero - minValue;

                    int tmpNBins = (int)((xMax - minValue)/minBinWidth);

                    if (tmpNBins >= minNumberOfBins) {

                        nBins = tmpNBins;

                        binWidth = minBinWidth;

                    } else {
                        go = false;
                    }

                } else {
                    go = false;
                }

            } else if (indexStartingConsecutiveEndLowCounts > -1) {
                // case 2

                // reduce xMax
                float tmpXMax = xHist[indexStartingConsecutiveEndLowCounts];;

                int tmpNBins = (int)((tmpXMax - minValue)/binWidth);

                if (tmpNBins >= minNumberOfBins) {

                    xMax = tmpXMax;

                    nBins = tmpNBins;

                } else {
                    go = false;
                }

            } else {
                go = false;
            }

            if (go) {

                xHist = new float[nBins];
                yHist = new int[nBins];
                Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);
/*
str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/
            }

        }

        xHist = new float[nBins];
        yHist = new int[nBins];
        Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
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
/*
str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/

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
    public static HistogramHolder createHistogramForSkewedDataForPeakResolution2(int preferredNumberOfBins,
        float[] values, float[] valueErrors, float minValue, float maxValue, int nIter) {

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

        int minNumberOfBins = (preferredNumberOfBins > 30) ? 30 : ((preferredNumberOfBins > 20) ? 20 : preferredNumberOfBins/2);

        if (binWidth > peakResBinWidth) {

            float tmpBinWidth = peakResBinWidth;

            int tmpNBins = (int)((xMax - minValue)/tmpBinWidth);

            if ((nIter == 1) || (tmpNBins < (2*nBins) )) {
                nBins = tmpNBins;
                binWidth = tmpBinWidth;

            } else if (tmpNBins < minNumberOfBins) {
                nBins = minNumberOfBins;
                binWidth = (maxValue - minValue)/nBins;

            } else {
                nBins = preferredNumberOfBins;
                binWidth = (xMax - minValue)/nBins;
            }
        }

        // if there are many bins with counts < 5, should  reduce the number of bins, but not to less than half preferred.
        // also, if there are many empty bins at start of histogram, increase the bin size.
/*
algorithms.util.PolygonAndPointPlotter plotter = null;
try {
plotter = new algorithms.util.PolygonAndPointPlotter();
} catch (Exception e) {}
String str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/

        boolean retry = (nIter < 200) ? true : false;

        int yPeakIndex = MiscMath.findYMaxIndex(yHist);

        int startsWithLessThan5Counts = 0;
        int indexStartingConsecutiveStartLowCounts = -1;

        int endsWithLessThan5Counts = 0;
        int indexStartingConsecutiveEndLowCounts = -1;
        int prevIndexStartingConsecutiveEndLowCounts = yHist.length;
        for (int i = (yHist.length - 1); i > -1; i--) {

            if (yHist[i] < 5) {

                // count the consecutive low counts at the end of the histogram:
                endsWithLessThan5Counts++;
                // avoid using points with x < x at ypeak
                if ( (i > yPeakIndex) && ((i == (prevIndexStartingConsecutiveEndLowCounts - 1)) ||
                (i == (prevIndexStartingConsecutiveEndLowCounts - 2)))) {
                    indexStartingConsecutiveEndLowCounts = i;
                    prevIndexStartingConsecutiveEndLowCounts = indexStartingConsecutiveEndLowCounts;
                } else {
                    prevIndexStartingConsecutiveEndLowCounts = i;
                }

                // count the consecutive low counts at start of histogram
                if (i > 0) {
                    // if there is a low count in i-1 too, increment this one, else reset it to zero (no skips)
                    if (yHist[i-1] < 5) {
                        startsWithLessThan5Counts++;
                        indexStartingConsecutiveStartLowCounts = i-1;
                    } else {
                        startsWithLessThan5Counts = 0;
                        indexStartingConsecutiveStartLowCounts = -1;
                    }
                }
            }
        }

        /*
         * Handle cases:
         *
         * (0) zeros at start  && zeros at end
         *     -- reduce xmax
         *        --> reduce nBins
         *
         * (1) zeros at start  only
         *        --> reduce nBins
         *
         * (2) zeros at end only
         *     -- reduce xmax
         *        --> reduce nBins
         *
         * Note, don't want cases to undo each others work upon iteration, so chose
         * to reduce the number of bins (rather than increase and decrease the bin width).
         */

        int chk = (int)(0.1*nBins);
        chk = 4;

        if ((indexStartingConsecutiveStartLowCounts == 0) && (indexStartingConsecutiveEndLowCounts > -1)) {
            // case 0

            // reduce xMax
            xMax = xHist[indexStartingConsecutiveEndLowCounts];;

            int tmpNBins = (int)((xMax - minValue)/binWidth);

            if (tmpNBins >= minNumberOfBins) {

                nBins = tmpNBins;

            } else {

                nBins = minNumberOfBins;

                binWidth = (xMax - minValue)/nBins;

                retry = false;
            }

        } else if (indexStartingConsecutiveStartLowCounts == 0) {
            // case 1

            if (startsWithLessThan5Counts > 2) {

                // reduce nBins
                float firstNonZero = xHist[startsWithLessThan5Counts];
                float minBinWidth = firstNonZero - minValue;

                int tmpNBins = (int)((xMax - minValue)/minBinWidth);

                if (tmpNBins >= minNumberOfBins) {

                    nBins = tmpNBins;

                    binWidth = minBinWidth;

                } else {

                    nBins = minNumberOfBins;

                    binWidth = (xMax - minValue)/nBins;

                    retry = false;
                }

            } else {
                retry = false;
            }

        } else if (indexStartingConsecutiveEndLowCounts > -1) {
            // case 2

            // reduce xMax
            xMax = xHist[indexStartingConsecutiveEndLowCounts];;

            int tmpNBins = (int)((xMax - minValue)/binWidth);

            if (tmpNBins >= minNumberOfBins) {

                nBins = tmpNBins;

                binWidth = minNumberOfBins;

            } else {
                nBins = minNumberOfBins;

                binWidth = (xMax - minValue)/nBins;

                retry = false;
          }

        } else {
            retry = false;
        }

        if (retry) {

            float[] tmpValues = new float[values.length];
            float[] tmpValueErrors = new float[values.length];
            int count = 0;
            for (int i = 0; i < tmpValues.length; i++) {
                if (values[i] <= xMax) {
                    tmpValues[count] = values[i];
                    tmpValueErrors[count] = valueErrors[i];
                    count++;
                }
            }
            tmpValues = Arrays.copyOf(tmpValues, count);
            tmpValueErrors = Arrays.copyOf(tmpValueErrors, count);

            nIter++;

            return createHistogramForSkewedDataForPeakResolution2(preferredNumberOfBins,
                tmpValues, tmpValueErrors, minValue, xMax, nIter);

        } else {

            xHist = new float[nBins];
            yHist = new int[nBins];
            Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);
/*
str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/

            xHist = new float[nBins];
            yHist = new int[nBins];
            Histogram.createHistogram(values, nBins, minValue, xMax, xHist, yHist, binWidth);  
            
            float[] yHistFloat = new float[yHist.length];
            for (int i = 0; i < yHist.length; i++) {
                yHistFloat[i] = (float) yHist[i];
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
/*str = String.format("n=%10d w=%.5f", nBins, binWidth);
plotter.addPlot(xHist, yHist, null, null, str);
try {
plotter.writeFile();
} catch (Exception e) {}*/

           
            return histogram;
        }
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

        float binWidth = (float) (binWidthFactor * 3.49f * meanStdev[1]/Math.pow(values.length, (1./3.)));

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

            if ((bin > -1) && (bin < xHist.length)) {

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

    public static HistogramHolder calculateSturgesHistogramWithoutTrimming(
        float[] values, float[] valueErrors) {
    
        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        int nIntervalsSturges = (int)Math.ceil( Math.log(values.length)/Math.log(2) );

        int nBins = Math.max(nIntervalsSturges, 10);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        float minx = MiscMath.findMin(values);
        float maxx = MiscMath.findMax(values);

        float binWidth = calculateBinWidth(minx, maxx, nBins);

        Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
        
        
        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
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
    
    public static HistogramHolder calculateSturgesHistogramRemoveZeroTail(
        float[] values, float[] valueErrors) {
    
        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException("values and valueErrors cannot be null and must be the same length");
        }

        int nIntervalsSturges = (int)Math.ceil( Math.log(values.length)/Math.log(2) );
        
        //int nItervalsRice = (int)(2*Math.pow(values.length, 0.3333));
        
        int nBins = 25;
        
        if (values.length > 10000) {
            nBins = 40;
        }

        nBins = Math.max(nIntervalsSturges, nBins);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
       
        float minx = MiscMath.findMin(values);
        float maxx = MiscMath.findMax(values);

        float binWidth = calculateBinWidth(minx, maxx, nBins);

        Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
        
        float maxy = MiscMath.findMax(yHist);

        int minCountsLimit = (int)Math.max(5, 0.03f*maxy);
        int countsBelowMinAtTail = 0;
        int lastLowCountIdx = yHist.length - 1;
        for (int i = (yHist.length - 1); i > -1; i--) {
            if (yHist[i] < minCountsLimit) {
                countsBelowMinAtTail++;
                lastLowCountIdx = i;
            } else {
                break;
            }
        }
        if (countsBelowMinAtTail > 0) {
            
            maxx = xHist[lastLowCountIdx];
            
            // keep nbins the same?            
            binWidth = calculateBinWidth(minx, maxx, nBins);

            Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
        }
        
        // if there are a large number of points, we'd like to increase the resolution of the peak if needed
        if (true || values.length > 10000) {
            int nLeftOfPeak = MiscMath.findYMaxIndex(yHist);
            int nIter = 0;
            while (nIter < 30 && nLeftOfPeak < 3 && (yHist[nLeftOfPeak] > 100)) {
                binWidth *= 0.8f;
                Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
                nLeftOfPeak = MiscMath.findYMaxIndex(yHist);
                nIter++;
            }
        }
          
        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
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
}
