package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 * @author nichole
 */
public class ColorHistogram {
    
    /**
     * histogram of 16 bins each of r' = r/(r+g+b) and g' = g/(r+g+b).
     * 
     * @param img
     * @param points 
     * @return histogram accessed as hist[0][bin] and hist[1][bin] where
     * 0 is for r' histogram and 1 is for the g' histogram.
     */
    public int[][] histogram(Image img, Set<PairInt> points) {
        
        int nBins = 16;
        
        int[][] hist = new int[2][];
        for (int i = 0; i < 2; ++i) {
            hist[i] = new int[nBins];
        }
       
        float binWidth = 1.f/(float)nBins;
                
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            float r = img.getR(x, y);
            float g = img.getG(x, y);
            float b = img.getB(x, y);
            float tot = r + g + b;
            
            float rPrime, gPrime;
            if (tot == 0) {
                rPrime = 0;
                gPrime = 0;
            } else {
                rPrime = r/tot;
                gPrime = g/tot;
            }
            
            int binNumber = Math.round(rPrime/binWidth);
            if (binNumber > (nBins - 1)) {
                binNumber = nBins - 1;
            }
            hist[0][binNumber]++;
            
            binNumber = Math.round(gPrime/binWidth);
            if (binNumber > (nBins - 1)) {
                binNumber = nBins - 1;
            }
            hist[1][binNumber]++;
        }
        
        return hist;
    }
    
    /**
     * histogram of 16 bins each of CIE LAB colors where the bins
     * of each color are taken to be in the range of min and max
     * possible values,
     * (0,0,0) to (28.51 3.28 2.15).
     * 
     * @param img
     * @param points 
     * @return histogram accessed as hist[0][bin] and hist[1][bin] and
     * hist[2][bin] where 0, 1, 2 are histograms for L, A, and B, respectively
     * using a range of values (0,0,0) to (28.51 3.28 2.15) and 16 bins.
     */
    public int[][] histogramCIELAB(ImageExt img, Set<PairInt> points) {
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        float maxL = 28.51f;
        float maxA = 3.28f;
        float maxB = 2.15f;
       
        float binWidthL = maxL/(float)nBins;
        float binWidthA = maxA/(float)nBins;
        float binWidthB = maxB/(float)nBins;
                
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            float[] lab = img.getCIELAB(x, y);
            
            int binNumberL = Math.round(lab[0]/binWidthL);
            if (binNumberL > (nBins - 1)) {
                binNumberL = nBins - 1;
            }
            hist[0][binNumberL]++;
            
            int binNumberA = Math.round(lab[1]/binWidthA);
            if (binNumberA > (nBins - 1)) {
                binNumberA = nBins - 1;
            }
            hist[1][binNumberA]++;
            
            int binNumberB = Math.round(lab[2]/binWidthB);
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
    /**
     * summation over colors of the differences between the
     * two histograms (note that the histograms are each normalized internally
     * by their total number of counts so that a difference can be calculated).
     * @param hist0
     * @param hist1
     * @return 
     */
    public float similarity(int[][] hist0, int[][] hist1) {
        
        if ((hist0.length != hist1.length)) {
            throw new IllegalArgumentException(
                "hist0 and hist1 must be same dimensions");
        }
        for (int i = 0; i < hist0.length; ++i) {
            if ((hist0[i].length != hist1[i].length)) {
                throw new IllegalArgumentException(
                    "hist0 and hist1 must be same dimensions");
            }
        }
        
        // the values need to be normalized by the total number of points
        // in the histogram, so calculat that first
        int[] n0 = new int[hist0.length];
        int[] n1 = new int[hist0.length];
        for (int i = 0; i < hist0.length; ++i) {
            for (int j = 0; j < hist0[i].length; ++j) {
                n0[i] += hist0[i][j];
                n1[i] += hist1[i][j];
            }
        }
            
        float sum = 0;
        for (int i = 0; i < hist0.length; ++i) {
            for (int j = 0; j < hist0[i].length; ++j) {
                float y0 = (float)hist0[i][j]/(float)n0[i];
                float y1 = (float)hist1[i][j]/(float)n1[i];
                float yDiff = y1 - y0;
                if (yDiff < 0) {
                    yDiff *= -1;
                }
                sum += yDiff;
            }
        }
        
        return sum;
    }
    
}
