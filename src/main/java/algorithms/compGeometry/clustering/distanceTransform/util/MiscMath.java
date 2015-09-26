package algorithms.compGeometry.clustering.distanceTransform.util;

import algorithms.misc.HistogramHolder;
import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MiscMath {
    
    public static float findMax(float[] a) {
        float max = Float.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
            }
        }
        return max;
    }
    
    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(float[] a) {
        if (a == null) {
            return -1;
        }
        float max = Float.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(int[] a) {
        if (a == null) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && (a[i] < Integer.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }
    
    public static double calculateArea(HistogramHolder hist, int idxStart, int idxStopIncl) {
        
        if (hist == null) {
            throw new IllegalArgumentException("hist cannot be null");
        }
        if (idxStart < 0 || (idxStart > (hist.getXHist().length - 1))) {
            throw new IllegalArgumentException("idxStart is out of bounds of hist");
        }
        if (idxStopIncl < 0 || (idxStopIncl > (hist.getXHist().length - 1))) {
            throw new IllegalArgumentException("idxStopIncl is out of bounds of hist");
        }
        
        double area = 0;
        
        for (int idx = idxStart; idx < idxStopIncl; ++idx) {
            float yTerm = hist.getYHistFloat()[idx + 1] + hist.getYHistFloat()[idx];
            float xLen = hist.getXHist()[idx + 1] - hist.getXHist()[idx];
            if (xLen < 0) {
                xLen *= -1;
            }
                        
            area += (yTerm * xLen)*0.5;
        }
        
        return area;
    }
    
    public static double[] calculateXYCentroids(Set<PairInt> points) {
        
        double xc = 0;
        double yc = 0;

        for (PairInt p : points) {
            xc += p.getX();
            yc += p.getY();
        }
        xc /= (double)points.size();
        yc /= (double)points.size();

        return new double[]{xc, yc};
    }
}
