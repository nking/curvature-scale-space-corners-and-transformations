package algorithms.imageProcessing;

import algorithms.util.Errors;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class ImageStatisticsHelper {
    
    /**
     * calculates the mean of values and returns it.
     * @param img
     * @return [meanR, meanG, meanB]
     */
    public static int getMean(final GreyscaleImage img) {
        
        long sum = 0;
        
        for (int i = 0; i < img.getNPixels(); i++) {
            sum += img.getValue(i);
        }

        return (int)(sum/img.getNPixels());
    }
    
    /**
     * calculates the median of values.
     * @param img
     * @return [meanR, meanG, meanB]
     */
    public static int getMedian(final GreyscaleImage img) {
        
        int[] values = new int[img.getNPixels()];
        for (int i = 0; i < img.getNPixels(); ++i) {
            values[i] = img.getValue(i);
        }
        
        return getMedian(values); 
    }
    
    public static int getMean(int[] a) {
        long sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i];
        }
        return (int)(sum/a.length);
    }
    
    public static float getMean(float[] a) {
        double sum = 0.;
        for (int i = 0; i < a.length; i++) {
            sum += a[i];
        }
        return (float)(sum/a.length);
    }
    
    public static int getMedian(int[] a) {
        int[] c = Arrays.copyOf(a, a.length);
        Arrays.sort(c);
        return c[c.length/2];
    }
    
    public static float getMedian(float[] a) {
        float[] c = Arrays.copyOf(a, a.length);
        Arrays.sort(c);
        return c[c.length/2];
    }
    
    /**
     * returns the Q1, Q2, Q3 and Q4 of the data a
     * 
     * @param a
     * @return 
     */
    public static float[] getQuartiles(float[] a) {
        
        if (a.length < 3) {
            throw new IllegalArgumentException("a.length must be at least 3");
        }
        
        float[] c = Arrays.copyOf(a, a.length);
        
        Arrays.sort(c);
        
        /*
                      median
             min        .         max
               .        .         .
               .   |    .    |    .
                q1   q2   q3   q4
        */
        
        int medianIdx = c.length >> 1;
        
        int q12Idx = (medianIdx - 1) >> 1;
        
        int q34Idx = (c.length + (medianIdx + 1))/2;
                
        return new float[]{c[q12Idx], c[medianIdx], c[q34Idx], c[c.length - 1]};
    }
    
    public static int[] getQuartiles(int[] a) {
        
        int[] c = Arrays.copyOf(a, a.length);
        
        Arrays.sort(c);
        
        /*
                      median
             min        .         max
               .        .         .
               .   |    .    |    .
                q1   q2   q3   q4
        */
        
        int medianIdx = c.length >> 1;
        
        int q12Idx = (medianIdx - 1) >> 1;
        
        int q34Idx = (c.length + (medianIdx + 1))/2;
                
        return new int[]{c[q12Idx], c[medianIdx], c[q34Idx], c[c.length - 1]};
    }
   
    /**
     * examine the statistics of pixels in a border of width borderWidth
     * around the borders of the image and return the statistics.
     * 
     * @param input
     * @param borderWidth
     * @param useSturges
     * @return 
     */
    public static ImageStatistics examineImageBorders(final GreyscaleImage input, 
        int borderWidth, boolean useSturges) {
                       
        float[] values = new float[input.getNPixels()];
        
        int count = 0;
        
        /**
         * | |
         * | |
         * | |
         * | |
         */
        for (int i = 0; i < borderWidth; i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                values[count] = input.getValue(i, j);                
                count++;
            }
        }
        
        /**
         * | |        | |
         * | |        | |
         * | |        | |
         * | |        | |
         */
        for (int i = (input.getWidth() - borderWidth); i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                values[count] = input.getValue(i, j);                
                count++;
            }
        }
        
        /**
         *   _________
         * | |________| |
         * | |        | |
         * | |        | |
         * | |        | |
         */
        for (int i = borderWidth; i < (input.getWidth() - borderWidth); i++) {
            for (int j = 0; j < borderWidth; j++) {
                values[count] = input.getValue(i, j);                
                count++;
            }
        }
        
        /**
         *   _________
         * | |________| |
         * | |        | |
         * | |________| |
         * | |________| |
         */
        for (int i = borderWidth; i < (input.getWidth() - borderWidth); i++) {
            for (int j = (input.getHeight() - borderWidth); j < input.getHeight(); j++) {
                values[count] = input.getValue(i, j);                
                count++;
            }
        }
        
        values = Arrays.copyOf(values, count);
        
        return examine(values, useSturges);
    }
    
    /**
     * examine the statistics of pixels in a border of width borderWidth
     * around the borders of the image and return the statistics.
     * 
     * @param input
     * @param useSturges
     * @return 
     */
    public static ImageStatistics examineImage(final GreyscaleImage img, 
        boolean useSturges) {
                       
        float[] values = new float[img.getNPixels()];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            values[i] = img.getValue(i);
        }
        
        return examine(values, useSturges);
    }
    
    /**
     * examine the statistics of pixels in a border of width borderWidth
     * around the borders of the image and return the statistics.
     * 
     * @param pixValues
     * @param useSturges
     * @return 
     */
    public static ImageStatistics examine(float[] pixValues, boolean useSturges) {
        
        ImageStatistics stats = new ImageStatistics();
          
        stats.setMedian(getMedian(pixValues));
        
        stats.setMean(getMean(pixValues));
        
        stats.setMin(MiscMath.findMin(pixValues));
        
        float xMax = MiscMath.findMax(pixValues);
        
        stats.setMax(xMax);
        
        float[] simulatedErrors = Errors.populateYErrorsBySqrt(pixValues);

        HistogramHolder hist = useSturges ?
            Histogram.calculateSturgesHistogram(0.0f, 256.0f, pixValues, 
                simulatedErrors)
            : Histogram.createSimpleHistogram(0.0f, 256.0f,
                10, pixValues, simulatedErrors);
        
        // think we probably want to remove the highest intensity bin, so
        // can think         
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        
        float mode = hist.getXHist()[yMaxIdx];
        
        stats.setMode(mode);
        
        stats.setHistogram(hist);
        
        stats.setQuartiles(ImageStatisticsHelper.getQuartiles(pixValues));
        
        stats.setHistogramAreaFraction(hist.getHistArea(xMax, 2));
        
        return stats;
    }
    
    /**
     * examine the statistics of pixels in a border of width borderWidth
     * around the borders of the image and return the statistics.
     * 
     * @param pixValues
     * @param useSturges
     * @return 
     */
    public static ImageStatistics examine(int[] pixValues, boolean useSturges) {
        
        ImageStatistics stats = new ImageStatistics();
          
        stats.setMedian(getMedian(pixValues));
        
        stats.setMean(getMean(pixValues));
        
        stats.setMin(MiscMath.findMin(pixValues));
        
        float xMax = MiscMath.findMax(pixValues);
        
        stats.setMax(xMax);
        
        float[] pixValuesF = new float[pixValues.length];
        for (int i = 0; i < pixValuesF.length; ++i) {
            pixValuesF[i] = pixValues[i];
        }
        
        float[] simulatedErrors = Errors.populateYErrorsBySqrt(pixValuesF);

        HistogramHolder hist = useSturges ?
            Histogram.calculateSturgesHistogram(0.0f, 256.0f, pixValuesF, 
                simulatedErrors)
            : Histogram.createSimpleHistogram(0.0f, 256.0f,
                10, pixValuesF, simulatedErrors);
        
        // think we probably want to remove the highest intensity bin, so
        // can think         
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        
        float mode = hist.getXHist()[yMaxIdx];
        
        stats.setMode(mode);
        
        stats.setHistogram(hist);
        
        stats.setQuartiles(ImageStatisticsHelper.getQuartiles(pixValuesF));
        
        stats.setHistogramAreaFraction(hist.getHistArea(xMax, 2));
        
        return stats;
    }
    
    /**
     * examine a width and height of pixels around the border of the image in
     * order to look for a low level intensity of the image, that is an effective
     * bias level due to the ambient lighting that can be subtracted from 
     * other pixels.  Note that if there are real zeros in the border histograms,
     * no 'bias' level should be subtracted from each pixel, but the histogram
     * is still useful for finding a lower threshold.
     * 
     * @param input 
     * @param useSturges 
     * @return  
     */
    public static ImageStatistics examineImageBorders(final GreyscaleImage input,
        boolean useSturges) {
                        
        if (input.getWidth() < 5) {
            return null;
        }
        
        // if <= 256x256, use whole image
        if ((input.getWidth() * input.getHeight()) < 65537) {
            return examineImage(input, useSturges);
        }
        
        int width = 10;
        
        if (input.getWidth() < 20) {
            width = 1;
        } else if (input.getWidth() < 50) {
            width = 5;
        } else if (input.getWidth() < 1000) {
            width = 10;
        } else {
            // choose 5 percent of image width or a default of 100 pixels?
            width = 100;
        }
        
        return examineImageBorders(input, width, useSturges);
    }
     
    public static int countPixels(final GreyscaleImage img, int lowValue, 
        int highValue) {
        
        int c = 0;
        
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                int v = img.getValue(col, row);
                if ((v >= lowValue) && (v <= highValue)) {
                    c++;
                }
            }
        }
        
        return c;
    }

}
