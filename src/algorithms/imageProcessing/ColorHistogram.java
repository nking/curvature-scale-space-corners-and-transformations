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
     * summation over colors of the differences between the
     * two histograms (note that the histograms are each normalized internally
     * by their total number of counts so that a difference can be calculated).
     * @param hist0
     * @param hist1
     * @return 
     */
    public float similarity(int[][] hist0, int[][] hist1) {
        
        if ((hist0.length != hist1.length) || (hist0[0].length != hist1[0].length)) {
            throw new IllegalArgumentException(
                "hist0 and hist1 must be same dimensions");
        }
        
        float sum = 0;
        for (int i = 0; i < hist0.length; ++i) {
            int n0 = 0;
            int n1 = 0;
            for (int j = 0; j < hist0[i].length; ++j) {
                n0 += hist0[i][j];
                n1 += hist1[i][j];
            }
            for (int j = 0; j < hist0[i].length; ++j) {
                float y0 = (float)hist0[i][j]/(float)n0;
                float y1 = (float)hist1[i][j]/(float)n1;
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
