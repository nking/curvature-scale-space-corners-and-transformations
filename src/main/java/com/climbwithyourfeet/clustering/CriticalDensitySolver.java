package com.climbwithyourfeet.clustering;

import algorithms.sorting.MultiArrayMergeSort;
import com.climbwithyourfeet.clustering.util.Histogram;
import com.climbwithyourfeet.clustering.util.HistogramHolder;
import com.climbwithyourfeet.clustering.util.MiscMath;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CriticalDensitySolver {
        
    private boolean debug = false;
    
    /**
     *
     */
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     *
     */
    public CriticalDensitySolver() {
    }
    
    /**
     *
     */
    public void setToDebug() {
        debug = true;
    }
    
    /**
     * using histograms of inverse sqrt of distance transform, find the center 
     * of the first peak and return it, else return 0
     * (0 as a critical density should result in an infinite critical separation
     * so no clusters).
     * @param distTrans
     * @param nPoints
     * @param width
     * @param height
     * @return 
     */
    float findCriticalDensity(int[][] distTrans, int nPoints, int width, 
        int height) {
        
        float critDens = findCriticalDensity(distTrans);

        return critDens;
        
    }
    
    /**
     * using histograms of 1/sqrt(distanceTransform[i][j]), find the center of 
     * the first peak and return it, else
     * return 0;
     * @param values
     * @return 
     */
    protected float findCriticalDensity(float[] values) {
        
        if (values == null || values.length < 10) {
            throw new IllegalArgumentException("values length must be 10 or more");
        }
        
        /*
        the goal of this method is to form a histogram comparable to the
        GeneralizedExtremeValue function to find the critical density
        (which is near the peak).
        
        TODO: This method needs improvements, especially for small numbers.
        */
        
        float[] vErrors = Histogram.populateYErrorsBySqrt(values);

        float xl = MiscMath.findMax(values);
        int nb = 40;
        if (nb > values.length) {
            nb = values.length/4;
            if (nb == 0) {
                nb = 1;
            }
        }
                
        HistogramHolder hist = Histogram.createSimpleHistogram(
            0, xl, nb, values, vErrors);
        
        if (debug) {
            String outFileSuffix = "_cluster_";
            hist.plotHistogram("clstr", outFileSuffix);
        }

        int len = hist.getXHist().length;
        
        int yFirstPeakIdx = Histogram.findFirstPeakIndex(hist);
                
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());

        //System.out.println("y1=" + yFirstPeakIdx + " ymx=" +
        //    yMaxIdx + " len=" + len);
        
        if (yMaxIdx > yFirstPeakIdx && (yMaxIdx > (len/2))) {
            
            nb = 8;
            hist = Histogram.createSimpleHistogram(
                0, xl, nb, values, vErrors);
            
            if (hist == null) {
                throw new IllegalStateException("error in algorithm");
            }

            if (debug) {
                String outFileSuffix = "_cluster_2_";
                hist.plotHistogram("clstr", outFileSuffix);
            }
        
            yFirstPeakIdx = Histogram.findFirstPeakIndex(hist);
            
            yMaxIdx = Histogram.findFirstMinimaFollowingPeak(hist, yFirstPeakIdx);

        } else {
                    
            // calculate the y quartiles above zero
            float[] quartiles = calcXQuartilesAboveZero(hist);

            if (debug) {
                System.out.println("quartiles=" + Arrays.toString(quartiles));
            }
            
            xl = Math.max(quartiles[0], quartiles[1]);

            nb /= 2;
            if (nb == 0) {
                nb = 1;
            }

            // make another histogram w/ x range being the 2nd quartile
            hist = Histogram.createSimpleHistogram(
                0, xl, nb, values, vErrors);
        
            yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        }
        
        if (yMaxIdx == -1) {
            return 0;
        }
        
        return 1.1f * hist.getXHist()[yMaxIdx];
    }

    private float findCriticalDensity(int[][] distTrans) {
        
        int w = distTrans.length;
        int h = distTrans[0].length;
        
        float[] values = new float[w * h];
        int count2 = 0;
        for (int i0 = 0; i0 < w; ++i0) {
            for (int j0 = 0; j0 < h; ++j0) {
                int v = distTrans[i0][j0];
                values[count2] = (float) (1. / Math.sqrt(v));
                count2++;
            }
        }

        return findCriticalDensity(values);
    }

    private float[] calcXQuartilesAboveZero(HistogramHolder hist) {

        // making a cumulative array of y.
        
        int n = hist.getXHist().length;
        double[] ys = new double[n];
        int count = 0;
        ys[0] = hist.getYHist()[0];
        for (int i = 1; i < n; ++i) {
            ys[i] = hist.getYHist()[i] + ys[i - 1];
        }
        
        double yTot = ys[n - 1];
        
        // where y cumulative is yTot/2
        int medianIdx = Arrays.binarySearch(ys, yTot/2);
        if (medianIdx < 0) {
            // idx = -*idx2 - 1
            medianIdx = -1*(medianIdx + 1);
        }
        if (medianIdx > (n - 1)) {
            medianIdx = n - 1;
        }
        
        // where y curmulative is yTot/4
        int q12Idx = Arrays.binarySearch(ys, yTot/2);
        if (q12Idx < 0) {
            // idx = -*idx2 - 1
            q12Idx = -1*(q12Idx + 1);
        }
        if (q12Idx > (n - 1)) {
            q12Idx = n - 1;
        }
        
        // where y curmulative is 3*yTot/4
        int q34Idx = Arrays.binarySearch(ys, 3.*yTot/4.);
        if (q34Idx < 0) {
            // idx = -*idx2 - 1
            q34Idx = -1*(q34Idx + 1);
        }
        if (q34Idx > (n - 1)) {
            q34Idx = n - 1;
        }
        
        float[] xs = hist.getXHist();        
        
        return new float[]{xs[q12Idx], xs[medianIdx], xs[q34Idx], xs[n - 1]};
    }

}
