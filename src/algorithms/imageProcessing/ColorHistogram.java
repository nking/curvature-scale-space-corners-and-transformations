package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 * @author nichole
 */
public class ColorHistogram {
    
    private static float maxL = 28.51f;
    private static float maxA = 3.28f;
    private static float maxB = 2.15f;
        
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
    
    public void add2To1(int[][] hist1, int[][] hist2) {
        
        for (int i = 0; i < hist1.length; ++i) {
            for (int j = 0; j < hist1[i].length; ++j) {
                hist1[i][j] += hist2[i][j];
            }
        }
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
        
        float binWidthL = maxL/(float)nBins;
        float binWidthA = maxA/(float)nBins;
        float binWidthB = maxB/(float)nBins;
                
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            float[] lab = img.getCIELAB(x, y);
            
            int binNumberL = Math.abs(Math.round(lab[0]/binWidthL));
            if (binNumberL > (nBins - 1)) {
                binNumberL = nBins - 1;
            }
            hist[0][binNumberL]++;
            
            int binNumberA = Math.abs(Math.round(lab[1]/binWidthA));
            if (binNumberA > (nBins - 1)) {
                binNumberA = nBins - 1;
            }
            hist[1][binNumberA]++;
            
            int binNumberB = Math.abs(Math.round(lab[2]/binWidthB));
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
    /**
     * summation over color histograms of the property min(h_i, h_j)
     * @param hist0
     * @param hist1
     * @return 
     */
    public float intersection(int[][] hist0, int[][] hist1) {
        
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
        
        /*
        After histograms are normalized to same total number
        
        K(a,b) = 
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */
            
        float sum = 0;
        float sum0 = 0;
        float sum1 = 0;
        for (int i = 0; i < hist0.length; ++i) {
            for (int j = 0; j < hist0[i].length; ++j) {
                float y0 = (float)hist0[i][j]/(float)n0[i];
                float y1 = (float)hist1[i][j]/(float)n1[i];
                sum += Math.min(y0, y1);
                sum0 += y0;
                sum1 += y1;
            }
        }
        
        float sim = sum / ((float)Math.min(sum0, sum1));
        
        return sim;
    }

    public int[][] histogramRGB(ImageExt img, Set<PairInt> points) {
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        float binWidth = 256.f/(float)nBins;

        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            float r = img.getR(x, y);
            float g = img.getG(x, y);
            float b = img.getB(x, y);
            
            int binNumberR = Math.abs(Math.round(r/binWidth));
            if (binNumberR > (nBins - 1)) {
                binNumberR = nBins - 1;
            }
            hist[0][binNumberR]++;
            
            int binNumberG = Math.abs(Math.round(g/binWidth));
            if (binNumberG > (nBins - 1)) {
                binNumberG = nBins - 1;
            }
            hist[1][binNumberG]++;
            
            int binNumberB = Math.abs(Math.round(b/binWidth));
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
    public int[][] histogramHSV(ImageExt img, Set<PairInt> points) {
        
        // range of Color values of h, s, b are 0:1 for all
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        float binWidth = 1.f/(float)nBins;

        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            float h = img.getHue(x, y);
            float s = img.getSaturation(x, y);
            float b = img.getBrightness(x, y);
            
            int binNumberH = Math.abs(Math.round(h/binWidth));
            if (binNumberH > (nBins - 1)) {
                binNumberH = nBins - 1;
            }
            hist[0][binNumberH]++;
            
            int binNumberS = Math.abs(Math.round(s/binWidth));
            if (binNumberS > (nBins - 1)) {
                binNumberS = nBins - 1;
            }
            hist[1][binNumberS]++;
            
            int binNumberB = Math.abs(Math.round(b/binWidth));
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
}
