package com.climbwithyourfeet.clustering.util;

import java.util.Arrays;

/**
 * smaller version of Histogram class to be packaged in a library jar
 * @author nichole
 */
public class Histogram {
    
    /**
     *
     * @param xMin
     * @param xMax
     * @param nBins
     * @param values
     * @param valueErrors
     * @return
     */
    public static HistogramHolder createSimpleHistogram(
        final float xMin, final float xMax, int nBins, 
        float[] values, float[] valueErrors) {

        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
                
        float binWidth = calculateBinWidth(xMin, xMax, nBins);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, xMin, xMax,
            xHist, yHist, binWidth);

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
     *
     * @param minValue
     * @param maxValue
     * @param nBins
     * @return
     */
    protected static float calculateBinWidth(float minValue, float maxValue, int nBins) {

        float xInterval = (maxValue - minValue)/(float)nBins;

        // expand interval if necessary to make sure the last point is in the last bin
        if ((int) ((maxValue - minValue)/xInterval) != (nBins - 1)) {
            float t = (maxValue + minValue)/2.0f;
            int powDelta = MiscMath.findPowerOf10(t);
            float pow10 = (float)Math.pow(10, powDelta);
            xInterval = (maxValue - minValue + pow10)/(float)nBins;
        }

        return xInterval;
    }
    
    /**
     *
     * @param a
     * @param nBins
     * @param aMin
     * @param aMax
     * @param xHist
     * @param yHist
     * @param binWidth
     */
    public static void createHistogram(float[] a, int nBins,
        float aMin, float aMax, float[] xHist, int[] yHist, float binWidth) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (xHist == null || xHist.length != nBins) {
            throw new IllegalArgumentException("xHist has to be of size nBins and initialized");
        }
        if (yHist == null || yHist.length != nBins) {
            throw new IllegalArgumentException("yHist has to be of size nBins and initialized");
        }

        Arrays.fill(yHist, 0);

        for (int i = 0; i < nBins; i++) {
            xHist[i] = aMin + (float)i*binWidth + (binWidth/2.f);
        }

        for (int i = 0; i < a.length; i++) {
            int bin = (int) ((a[i] - aMin)/binWidth);
            if ((bin > -1) && (bin < nBins)) {
                yHist[bin]++;
            }
        }
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
        float[] values, float[] valueErrors, float[] xHistErrorsOutput, 
        float[] yHistErrorsOutput) {

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

            // estimating it as binWidth/2. any point has x value = bin center +- binWidth/2
            xHistErrorsOutput[i] = xInterval/2.0f; //c;

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
    
     /**
     * make assumption that errors are dominated by shot noise and so the noise is sqrt(y[i]).
     *
     * @param y
     * @return
     */
    public static float[] populateYErrorsBySqrt(float[] y) {

        float[] dy = new float[y.length];

        float maxError = 0.f;

        for (int i = 0; i < dy.length; i++) {
            if (y[i] > 0) {
                dy[i] =(float)(Math.sqrt(y[i]));
                if (dy[i] > maxError) {
                    maxError = dy[i];
                }
            }
        }
        for (int i = 0; i < dy.length; i++) {
            if (y[i] == 0) {
                dy[i] = maxError;
            }
        }

        return dy;
    }

}
