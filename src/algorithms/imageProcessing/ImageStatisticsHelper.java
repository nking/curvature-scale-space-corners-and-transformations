package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.util.Errors;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import algorithms.imageProcessing.transform.Transformer;

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
    
    /**
     * returns the Q1, Q2, Q3 and Q4 of the data a
     * 
     * @param a
     * @return 
     */
    public static int[] getQuartiles(List<Integer> a) {
        
        int n = a.size();
        
        if (n < 3) {
            throw new IllegalArgumentException("a.length must be at least 3");
        }
                
        Collections.sort(a);
        
        /*
                      median
             min        .         max
               .        .         .
               .   |    .    |    .
                q1   q2   q3   q4
        */
        
        int medianIdx = n >> 1;
        
        int q12Idx = (medianIdx - 1) >> 1;
        
        int q34Idx = (n + (medianIdx + 1))/2;
                
        return new int[]{a.get(q12Idx).intValue(), a.get(medianIdx).intValue(), 
            a.get(q34Idx).intValue(), a.get(n - 1).intValue()};
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

    public static float[] calculateAvgCIEXY(ImageExt input, Set<PairInt> points) {
        
        if (points.isEmpty()) {
            return new float[]{Float.NaN, Float.NaN};
        }
            
        float cieXSum = 0;
        float cieYSum = 0;
        for (PairInt p : points) {
            cieXSum += input.getCIEX(p.getX(), p.getY());
            cieYSum += input.getCIEY(p.getX(), p.getY());
        }
        
        cieXSum /= (float)points.size();
        cieYSum /= (float)points.size();
        
        return new float[]{cieXSum, cieYSum};
    }

    static float[] calculateAvgRGB(ImageExt input, Set<PairInt> points) {
        
        if (points.isEmpty()) {
            return new float[]{Float.NaN, Float.NaN, Float.NaN};
        }
            
        float rSum = 0;
        float gSum = 0;
        float bSum = 0;
        for (PairInt p : points) {
            rSum += input.getR(p.getX(), p.getY());
            gSum += input.getG(p.getX(), p.getY());
            bSum += input.getB(p.getX(), p.getY());
        }
        
        rSum /= (float)points.size();
        gSum /= (float)points.size();
        bSum /= (float)points.size();
        
        return new float[]{rSum, gSum, bSum};
    }
    
    /**
     * a method to determine the intersection of transformed image 1 with
     * image 2 and then examine the distribution of stats's points in 
     * 4 quadrants of the intersection to return whether stats are present in
     * all quadrants.  A caveat of the method is that not all of the 
     * intersection necessarily has image details which could be matched, for 
     * example, clear sky does not have corners using the methods here.
     * @param parameters
     * @param stats
     * @param img1Width
     * @param img1Height
     * @param img2Width
     * @param img2Height
     * @return 
     */
    public static boolean statsCoverIntersection(TransformationParameters 
        parameters, List<FeatureComparisonStat> stats,
        int img1Width, int img1Height, int img2Width, int img2Height) {
        
        if (parameters == null || stats.size() < 4) {
            return false;
        }
        int[] counts = getQuadrantCountsForIntersection(parameters, stats,
            img1Width, img1Height, img2Width, img2Height);
        
        int nq = 0;
        for (int c : counts) {
            if (c > 0) {
                nq++;
            }
        }
        
        return (nq == 4);
    }
    
    /**
     * a method to determine the intersection of transformed image 1 with
     * image 2 and then examine the distribution of stats's points in 
     * 4 quadrants of the intersection to return whether stats are present in
     * all quadrants.  A caveat of the method is that not all of the 
     * intersection necessarily has image details which could be matched, for 
     * example, clear sky does not have corners using the methods here.
     * @param parameters
     * @param stats
     * @param img1Width
     * @param img1Height
     * @param img2Width
     * @param img2Height
     * @return 
     */
    public static int[] getQuadrantCountsForIntersection(TransformationParameters 
        parameters, List<FeatureComparisonStat> stats,
        int img1Width, int img1Height, int img2Width, int img2Height) {
        
        if (parameters == null) {
            return new int[4];
        }
        
        /*
        calculate the intersection of the 2 images.
        divide the region into 4 parts (2 vertical and 2 horizontal) by noting
        the 4 boundary points for each and making a polygon for each.
        
        then use point in polygon tests to count the number of stats.point2's
        in each of the 4 regions.        
        */
        
        /*
       / \   ( tr )    ( tr )            (x2q2, y2q2)  d5   (x2q3, y2q3)
        |
        |                                 d2           d3             d4
        |
        0    ( tr )    ( tr )            (x2q1, y2q1)  d1   (x2q0, y2q0)
          0 -->
        */
        
        double[][] img2Intersection = getBoundsOfIntersectionInFrame2(
            parameters, img1Width, img1Height, img2Width, img2Height);
        
        float[] d1 = new float[]{
            (float)((img2Intersection[0][0] + img2Intersection[1][0])/2.f),
            (float)((img2Intersection[0][1] + img2Intersection[1][1])/2.f)};     
        float[] d2 = new float[]{
            (float)((img2Intersection[1][0] + img2Intersection[2][0])/2.f),
            (float)((img2Intersection[1][1] + img2Intersection[2][1])/2.f)};
        float[] d4 = new float[]{
            (float)((img2Intersection[0][0] + img2Intersection[3][0])/2.f),
            (float)((img2Intersection[0][1] + img2Intersection[3][1])/2.f)};
        float[] d5 = new float[]{
            (float)((img2Intersection[2][0] + img2Intersection[3][0])/2.f),
            (float)((img2Intersection[2][1] + img2Intersection[3][1])/2.f)};
        float[] d3 = new float[]{(d2[0] + d4[0])/2.f, (d1[1] + d5[1])/2.f};
        
        float[] xPoly0 = new float[5];
        float[] yPoly0 = new float[5];
        xPoly0[0] = (float)img2Intersection[0][0];
        yPoly0[0] = (float)img2Intersection[0][1];
        xPoly0[1] = d1[0];
        yPoly0[1] = d1[1];
        xPoly0[2] = d3[0];
        yPoly0[2] = d3[1];
        xPoly0[3] = d4[0];
        yPoly0[3] = d4[1];
        xPoly0[4] = xPoly0[0];
        yPoly0[4] = yPoly0[0];

        /*
       / \   ( tr )    ( tr )            (x2q2, y2q2)  d5   (x2q3, y2q3)
        |
        |                                 d2           d3             d4
        |
        0    ( tr )    ( tr )            (x2q1, y2q1)  d1   (x2q0, y2q0)
          0 -->
        */
        
        float[] xPoly1 = new float[5];
        float[] yPoly1 = new float[5];
        xPoly1[0] = d1[0];
        yPoly1[0] = d1[1];
        xPoly1[1] = (float)img2Intersection[1][0];
        yPoly1[1] = (float)img2Intersection[1][1];
        xPoly1[2] = d2[0];
        yPoly1[2] = d2[1];
        xPoly1[3] = d3[0];
        yPoly1[3] = d3[1];
        xPoly1[4] = xPoly1[0];
        yPoly1[4] = yPoly1[0];
        
        float[] xPoly2 = new float[5];
        float[] yPoly2 = new float[5];
        xPoly2[0] = d3[0];
        yPoly2[0] = d3[1];
        xPoly2[1] = d2[0];
        yPoly2[1] = d2[1];
        xPoly2[2] = (float)img2Intersection[2][0];
        yPoly2[2] = (float)img2Intersection[2][1];
        xPoly2[3] = d5[0];
        yPoly2[3] = d5[1];
        xPoly2[4] = xPoly2[0];
        yPoly2[4] = yPoly2[0];
        
        /*
       / \   ( tr )    ( tr )            (x2q2, y2q2)  d5   (x2q3, y2q3)
        |
        |                                 d2           d3             d4
        |
        0    ( tr )    ( tr )            (x2q1, y2q1)  d1   (x2q0, y2q0)
          0 -->
        */
        
        float[] xPoly3 = new float[5];
        float[] yPoly3 = new float[5];
        xPoly3[0] = d4[0];
        yPoly3[0] = d4[1];
        xPoly3[1] = d3[0];
        yPoly3[1] = d3[1];
        xPoly3[2] = d5[0];
        yPoly3[2] = d5[1];
        xPoly3[3] = (float)img2Intersection[3][0];
        yPoly3[3] = (float)img2Intersection[3][1];
        xPoly3[4] = xPoly3[0];
        yPoly3[4] = yPoly3[0];
        
        PointInPolygon poly = new PointInPolygon();
        
        int[] count = new int[4];
        for (FeatureComparisonStat stat : stats) {
            int x = stat.getImg2Point().getX() * stat.getBinFactor2();
            int y = stat.getImg2Point().getY() * stat.getBinFactor2();
            boolean isIn = poly.isInSimpleCurve(x, y, xPoly0, yPoly0, 5);
            if (isIn) {
                count[0]++;
            } else {
                isIn = poly.isInSimpleCurve(x, y, xPoly1, yPoly1, 5);
                if (isIn) {
                    count[1]++;
                } else {
                    isIn = poly.isInSimpleCurve(x, y, xPoly2, yPoly2, 5);
                    if (isIn) {
                        count[2]++;
                    } else {
                        isIn = poly.isInSimpleCurve(x, y, xPoly3, yPoly3, 5);
                        if (isIn) {
                            count[3]++;
                        }
                    }
                }
            }
        }
        
        return count;
    }
    
    /**
     transforms the bounds of img1 to the frame of img2 and calculates the
     intersection of them.  The result is the pixel coordinates of the
     4 corners of the intersection bounds in the order
     (x2q1, y2q1), (x2q2, y2q2), (x2q3, y2q3), (x2q4, y2q4).
        
       / \  (x2q3, y2q3)      (x2q4, y2q4)
        |
        |
        0   (x2q2, y2q2)      (x2q1, y2q1)
          0 -->
        
     * @param parameters
     * @param img1Width
     * @param img1Height
     * @param img2Width
     * @param img2Height
     * @return 
     */
    public static double[][] getBoundsOfIntersectionInFrame2(TransformationParameters 
        parameters, int img1Width, int img1Height, int img2Width, int img2Height) {
        
        //calculate the intersection of the 2 images
        
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();
        
        Transformer transformer = new Transformer();
        
        TransformationParameters revParams = tc.swapReferenceFrames(parameters);
        
        /*
        
       / \   ( tr )    ( tr )            (x2q3, y2q3)      (x2q4, y2q4)
        |
        |
        0    ( tr )    ( tr )            (x2q2, y2q2)      (x2q1, y2q1)
          0 -->
        
        */
        
        // determine intersection of img2 with img1 in img1 reference frame
        double[] q1Tr = transformer.applyTransformation(revParams, 
            img2Width - 1, 0);
        
        double[] q2Tr = transformer.applyTransformation(revParams, 
            0, 0);
        
        double[] q3Tr = transformer.applyTransformation(revParams, 
            0, img2Height - 1);
        
        double[] q4Tr = transformer.applyTransformation(revParams, 
            img2Width - 1, img2Height - 1);
        
        // if the transformed bounds are off image, reset the bounds to img1 bounds
        double[][] img1Intersection = new double[4][2];
        img1Intersection[0] = q1Tr;
        img1Intersection[1] = q2Tr;
        img1Intersection[2] = q3Tr;
        img1Intersection[3] = q4Tr;
        
        for (double[] xyTr : img1Intersection) {
            if (xyTr[0] < 0) {
                xyTr[0] = 0;
            } else if (xyTr[0] > (img1Width - 1)) {
                xyTr[0] = (img1Width - 1);
            }
            if (xyTr[1] < 0) {
                xyTr[1] = 0;
            } else if (xyTr[1] > (img1Height - 1)) {
                xyTr[1] = (img1Height - 1);
            }
        }
        
        // transform the img1 intersection into reference frame of img2
        double[] q1TrTr = transformer.applyTransformation(parameters, q1Tr[0], q1Tr[1]);
        
        double[] q2TrTr = transformer.applyTransformation(parameters, q2Tr[0], q2Tr[1]);
        
        double[] q3TrTr = transformer.applyTransformation(parameters, q3Tr[0], q3Tr[1]);
        
        double[] q4TrTr = transformer.applyTransformation(parameters, q4Tr[0], q4Tr[1]);
        
        double[][] img2Intersection = new double[4][2];
        img2Intersection[0] = q1TrTr;
        img2Intersection[1] = q2TrTr;
        img2Intersection[2] = q3TrTr;
        img2Intersection[3] = q4TrTr;
        
        for (double[] xyTr : img2Intersection) {
            if (xyTr[0] < 0) {
                xyTr[0] = 0;
            } else if (xyTr[0] > (img2Width - 1)) {
                xyTr[0] = (img2Width - 1);
            }
            if (xyTr[1] < 0) {
                xyTr[1] = 0;
            } else if (xyTr[1] > (img2Height - 1)) {
                xyTr[1] = (img2Height - 1);
            }
        }
        
        return img2Intersection;
    }

}
