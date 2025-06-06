package algorithms.imageProcessing;

import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * 
 * NOTE: instead of this class, consider using a polar
 * theta of cie luv to better cover the distribution
 * of colors.
 * 
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
                
        int[][] hist = createWithDefaultSize();
        
        // 16
        int nBins = hist[0].length;
       
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
    
    public void add2To1(int[] hist1, int[] hist2) {
        for (int i = 0; i < hist1.length; ++i) {
            hist1[i] += hist2[i];
        }
    }
    
    /**
     * histogram of 16 bins each of CIE LAB colors where the bins
     * of each color are taken to be in the range of min and max
     * possible values
     * <pre>
     *    L    0 to 28.5
     *    A  -46.9  62.5
     *    B  -45.7  48.0
     * </pre>
     * 
     * @param img
     * @param points 
     * @return histogram accessed as hist[0][bin] and hist[1][bin] and
     * hist[2][bin] where 0, 1, 2 are histograms for L, A, and B, respectively
     * using a range of values (0,-50,-50) to (28.5 62.5 48.0) and 16 bins.
     */
    public int[][] histogramCIELAB(ImageExt img, Set<PairInt> points) {
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        //(0,-50,-50) to (28.5 62.5 48.0) 
        
        float binWidthL = 28.5f/(float)nBins;
        float binWidthA = (62.5f + 50.f)/(float)nBins;
        float binWidthB = (48.0f + 50.f)/(float)nBins;
                
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            float[] lab = img.getCIELAB(x, y);
            
            int binNumberL = (int)(lab[0]/binWidthL);
            if (binNumberL > (nBins - 1)) {
                binNumberL = nBins - 1;
            }
            hist[0][binNumberL]++;
            
            //(0,-50,-50) to (28.5 62.5 48.0)
            int binNumberA = (int)((lab[1] - -50.f)/binWidthA);
            if (binNumberA > (nBins - 1)) {
                binNumberA = nBins - 1;
            }
            hist[1][binNumberA]++;
            
            int binNumberB = (int)((lab[2] - -50.f)/binWidthB);
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
    /**
     * histogram of 16 bins each of CIE LAB colors where the bins
     * of each color are taken to be in the range of min and max
     * possible values
     * <pre>
     *    L    0 to 28.5
     *    A  -46.9  62.5
     *    B  -45.7  48.0
     * </pre>
     * 
     * @param img
     * @return histogram accessed as hist[0][bin] and hist[1][bin] and
     * hist[2][bin] where 0, 1, 2 are histograms for L, A, and B, respectively
     * using a range of values (0,-50,-50) to (28.5 62.5 48.0) and 16 bins.
     */
    public int[][] histogramCIELAB(ImageExt img, TIntList xPoints,
        TIntList yPoints) {
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        //(0,-50,-50) to (28.5 62.5 48.0) 
        
        float binWidthL = 28.5f/(float)nBins;
        float binWidthA = (62.5f + 50.f)/(float)nBins;
        float binWidthB = (48.0f + 50.f)/(float)nBins;
                
        for (int i = 0; i < xPoints.size(); ++i) {
            int x = xPoints.get(i);
            int y = yPoints.get(i);
            
            float[] lab = img.getCIELAB(x, y);
            
            int binNumberL = (int)(lab[0]/binWidthL);
            if (binNumberL > (nBins - 1)) {
                binNumberL = nBins - 1;
            }
            hist[0][binNumberL]++;
            
            //(0,-50,-50) to (28.5 62.5 48.0)
            int binNumberA = (int)((lab[1] - -50.f)/binWidthA);
            if (binNumberA > (nBins - 1)) {
                binNumberA = nBins - 1;
            }
            hist[1][binNumberA]++;
            
            int binNumberB = (int)((lab[2] - -50.f)/binWidthB);
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
    /**
     * histogram of 16 bins each of CIE LUV colors where the bins
     * of each color are taken to be in the range of min and max
     * possible values
     * <pre>
     * range of return values when using default standard illumination of
     * D65 daylight is:
     * L       0 to 104.5
     * u   -86.9 to 183.8
     * v  -141.4 to 112.3
     * </pre>
     * 
     * @param img
     * @param points 
     * @return histogram accessed as hist[0][bin] and hist[1][bin] and
     * hist[2][bin] where 0, 1, 2 are histograms for L, A, and B, respectively
     * using a range of values and 16 bins.
     */
    public int[][] histogramCIELUV(Image img, Set<PairInt> points) {
        
        int nBins = 16;
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
               
        float binWidthL = 104.5f/(float)nBins;
        float binWidthU = (183.8f + 86.9f)/(float)nBins;
        float binWidthV = (112.3f + 141.4f)/(float)nBins;
                
        for (PairInt p : points) {
            
            float[] lab = cieC.rgbToPolarCIELUV(
                img.getR(p), img.getG(p), img.getB(p));
            
            int binNumberL = (int)(lab[0]/binWidthL);
            if (binNumberL > (nBins - 1)) {
                binNumberL = nBins - 1;
            }
            hist[0][binNumberL]++;
            
            //(0,-50,-50) to (28.5 62.5 48.0)
            int binNumberU = (int)((lab[1] - -50.f)/binWidthU);
            if (binNumberU > (nBins - 1)) {
                binNumberU = nBins - 1;
            }
            hist[1][binNumberU]++;
            
            int binNumberV = (int)((lab[2] - -50.f)/binWidthV);
            if (binNumberV > (nBins - 1)) {
                binNumberV = nBins - 1;
            }
            hist[2][binNumberV]++;
        }
        
        return hist;
    }
    
    /**
     * histogram of 16 bins each of CIE LAB colors where the bins
     * of each color are taken to be in the range of min and max
     * possible values,
     * range of values (0,-50,-50) to (28.5 62.5 48.0).
     * 
     * @param img
     * @param points 
     * @return histogram accessed as hist[0][bin] and hist[1][bin] and
     * hist[2][bin] where 0, 1, 2 are histograms for L, A, and B, respectively
     * using a range of values range of values (0,-50,-50) to (28.5 62.5 48.0) and 16 bins.
     */
    public int[][] histogramCIELAB(ImageExt img, TIntSet points) {
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        //(0,-50,-50) to (28.5 62.5 48.0) 
        
        float binWidthL = 28.5f/(float)nBins;
        float binWidthA = (62.5f + 50.f)/(float)nBins;
        float binWidthB = (48.0f + 50.f)/(float)nBins;
                
        TIntIterator iter = points.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            
            float[] lab = img.getCIELAB(pixIdx);
            
            int binNumberL = (int)(lab[0]/binWidthL);
            if (binNumberL > (nBins - 1)) {
                binNumberL = nBins - 1;
            }
            hist[0][binNumberL]++;
            
            //(0,-50,-50) to (28.5 62.5 48.0)
            int binNumberA = (int)((lab[1] - -50.f)/binWidthA);
            if (binNumberA > (nBins - 1)) {
                binNumberA = nBins - 1;
            }
            hist[1][binNumberA]++;
            
            int binNumberB = (int)((lab[2] - -50.f)/binWidthB);
            if (binNumberB > (nBins - 1)) {
                binNumberB = nBins - 1;
            }
            hist[2][binNumberB]++;
        }
        
        return hist;
    }
    
    /**
     * histogram of CIE LUV C and H which are the polar angle of U,V and the 
     * magnitude.  The histogram is a one dimensional binning of 
     * divisions of C and H in combination.  8 divisions in C or H results
     * in 64 bins, for example.
     * 
     * The ranges of the bins are the ranges from use of the D65 standard
     * illuminant which are
     * <pre>
     * * L* 0 to 105
     * c  0 to 139
     * h  0 to 359
     * </pre>
     * 
     * @param img
     * @param points 
     * @return histogram of CIE LUV C and H which are the polar angle of 
     * U,V and the magnitude.  The histogram is a one dimensional binning 
     * of divisions of C and H in combination.  8 divisions in C or H 
     * results in 64 bins, for example.
     */
    public int[] histogramCIECH64(ImageExt img, Set<PairInt> points) {
        
        int nBins = 8;
        
        int nTot = nBins * nBins;

        int[] hist = new int[nTot];
        
        float binWidthC = 139.f/(float)nBins;
        float binWidthH = 359.f/(float)nBins;
        
        // idx = (cBin * nBins) + hBin
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (PairInt p : points) {
            
            int pixIdx = img.getInternalIndex(p);
            
            float[] lch = cieC.rgbToCIELCH(
                img.getR(pixIdx), img.getG(pixIdx), img.getR(pixIdx));
            
            int binNumberC = (int)(lch[1]/binWidthC);
            if (binNumberC > (nBins - 1)) {
                binNumberC = nBins - 1;
            }
            
            int binNumberH = (int)(lch[2]/binWidthH);
            if (binNumberH > (nBins - 1)) {
                binNumberH = nBins - 1;
            }
            
            int hIdx = (binNumberC * nBins) + binNumberH;
            
            assert(hIdx < nTot);
            
            hist[hIdx]++;
        }
        
        return hist;
    }
    
    /**
     * histogram of CIE LUV C and H which are the polar angle of U,V and the 
     * magnitude.  The histogram is a one dimensional binning of 
     * divisions of C and H in combination.  8 divisions in C or H results
     * in 64 bins, for example.
     * 
     * The ranges of the bins are the ranges from use of the D65 standard
     * illuminant which are
     * <pre>
     * * L* 0 to 105
     * c  0 to 139
     * h  0 to 359
     * </pre>
     * 
     * @param img
     * @return histogram of CIE LUV C and H which are the polar angle of 
     * U,V and the magnitude.  The histogram is a one dimensional binning 
     * of divisions of C and H in combination.  8 divisions in C or H 
     * results in 64 bins, for example.
     */
    public int[] histogramCIECH64(ImageExt img, TIntList xPoints,
        TIntList yPoints) {
        
        int nBins = 8;
        
        int nTot = nBins * nBins;

        int[] hist = new int[nTot];
        
        float binWidthC = 139.f/(float)nBins;
        float binWidthH = 359.f/(float)nBins;
        
        // idx = (cBin * nBins) + hBin
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < xPoints.size(); ++i) {
            
            int pixIdx = img.getInternalIndex(xPoints.get(i), yPoints.get(i));
            
            float[] lch = cieC.rgbToCIELCH(
                img.getR(pixIdx), img.getG(pixIdx), img.getR(pixIdx));
            
            int binNumberC = (int)(lch[1]/binWidthC);
            if (binNumberC > (nBins - 1)) {
                binNumberC = nBins - 1;
            }
            
            int binNumberH = (int)(lch[2]/binWidthH);
            if (binNumberH > (nBins - 1)) {
                binNumberH = nBins - 1;
            }
            
            int hIdx = (binNumberC * nBins) + binNumberH;
            
            assert(hIdx < nTot);
            
            hist[hIdx]++;
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
        // in the histogram, so calculate that first
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
                float y0 = 0;
                float y1 = 0;
                if (n0[i] > 0) {
                    y0 = (float)hist0[i][j]/(float)n0[i];
                }
                if (n1[i] > 0) {
                    y1 = (float)hist1[i][j]/(float)n1[i];
                }
                sum += Math.min(y0, y1);
                sum0 += y0;
                sum1 += y1;
            }
        }
        
        float d = Math.min(sum0, sum1);
        
        float sim = (d == 0.f) ? 0 : sum/d;
                
        return sim;
    }
    
    public int[][] createWithDefaultSize() {
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        return hist;
    }
    
    /**
     * note, hist0 and hist1 must be normalized before use here 
     * (unlike intersection which internally normalizes the histograms).
     * @param hist0
     * @param hist1
     * @return 
     */
    public float chiSquaredSum(int[][] hist0, int[][] hist1) {
     
        float[] x2 = new float[3];
        chiSquaredSum(hist0, hist1, x2);
        
        float sum = 0;
        for (int i = 0; i < hist0.length; ++i) {
            sum += x2[i];
        }
        
        return sum;
    }
    
    /**
     * note, hist0 and hist1 must be normalized before use here 
     * (unlike intersection which internally normalizes the histograms).
     * @param hist0
     * @param hist1
     * @param output 
     */
    public void chiSquaredSum(int[][] hist0, int[][] hist1, float[] output) {
        
        if ((hist0.length != hist1.length)) {
            throw new IllegalArgumentException(
                "hist0 and hist1 must be same dimensions");
        }
        if (output.length != hist0.length) {
            throw new IllegalArgumentException("output must be same length as "
                + " hist0 first dimension");
        }
        for (int i = 0; i < hist0.length; ++i) {
            if ((hist0[i].length != hist1[i].length)) {
                throw new IllegalArgumentException(
                    "hist0 and hist1 must be same dimensions");
            }
        }
        
        //1/2 times sum over all bins of : (h1 - h2)^2/(h1 + h2)
        for (int i = 0; i < hist0.length; ++i) {
            for (int k = 0; k < hist0[i].length; ++k) {
                float diff = hist0[i][k] - hist1[i][k];
                float add = hist0[i][k] + hist1[i][k];
                if (add > 0) {
                    output[i] += (diff * diff)/add;
                }
            }
            output[i] *= 0.5f;
        }        
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
    public int[][] histogramHSV(ImageExt img, TIntList xPoints,
        TIntList yPoints) {
        
        // range of Color values of h, s, b are 0:1 for all
        
        int nBins = 16;
        
        int[][] hist = new int[3][];
        for (int i = 0; i < 3; ++i) {
            hist[i] = new int[nBins];
        }
        
        float binWidth = 1.f/(float)nBins;

        for (int i = 0; i < xPoints.size(); ++i) {
            
            int x = xPoints.get(i);
            int y = yPoints.get(i);
            
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
    
    /**
     * 1D histogram of 16 bins of greyscale image.  max possible
     * value in image should be provided also.
     * 
     * @param img
     * @param points 
     * @return histogram accessed as hist[bin]
     */
    public int[] histogram1D(GreyscaleImage img, Set<PairInt> points,
        int maxV) {
             
        int nBins = 16;
        
        int[] hist = new int[nBins];
        
        float binWidth = (float)maxV/(float)nBins;
           
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            float v = img.getValue(x, y);
            
            int binNumber = Math.round(v/binWidth);
            if (binNumber > (nBins - 1)) {
                binNumber = nBins - 1;
            }
            hist[binNumber]++;           
        }
        
        return hist;
    }
    
    /**
     * 1D histogram of 16 bins of greyscale image.  max possible
     * value in image should be provided also.
     * 
     * @param values from 0 to maxV, inclusive 
     * @return histogram accessed as hist[bin]
     */
    public int[] histogram1D(TIntList values, int maxV) {
             
        int nBins = 16;
        
        int[] hist = new int[nBins];
        
        float binWidth = (float)maxV/(float)nBins;
           
        for (int i = 0; i < values.size(); ++i) {
            
            float v = values.get(i);
            
            int binNumber = Math.round(v/binWidth);
            if (binNumber > (nBins - 1)) {
                binNumber = nBins - 1;
            }
            hist[binNumber]++;           
        }
        
        return hist;
    }
    
    /**
     * summation over color histograms of the property min(h_i, h_j)
     * @param hist0
     * @param hist1
     * @return 
     */
    public float intersection(int[] hist0, int[] hist1) {
        
        if ((hist0.length != hist1.length)) {
            throw new IllegalArgumentException(
                "hist0 and hist1 must be same dimensions");
        }
                
        // the values need to be normalized by the total number of points
        // in the histogram, so calculate that first
        int n0 = 0;
        int n1 = 0;
        for (int i = 0; i < hist0.length; ++i) {
            n0 += hist0[i];
            n1 += hist1[i];
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
        for (int j = 0; j < hist0.length; ++j) {
            float y0 = 0;
            float y1 = 0;
            if (n0 > 0) {
                y0 = (float)hist0[j]/(float)n0;
            }
            if (n1 > 0) {
                y1 = (float)hist1[j]/(float)n1;
            }
            sum += Math.min(y0, y1);
            sum0 += y0;
            sum1 += y1;
        }
        
        float d = Math.min(sum0, sum1);
        
        float sim = (d == 0.f) ? 0 : sum/d;
                
        return sim;
    }

    /**
     * 
     * @param hist0
     * @param hist1
     * @param bins0Shift the number of bins to shift bins1 by during
     * comparison.
     * @return 
     */
    public float intersection(int[] hist0, int[] hist1, int bins0Shift) {
    
        if ((hist0.length != hist1.length)) {
            throw new IllegalArgumentException(
                "hist0 and hist1 must be same dimensions");
        }
                
        // the values need to be normalized by the total number of points
        // in the histogram, so calculate that first
        int n0 = 0;
        int n1 = 0;
        for (int i = 0; i < hist0.length; ++i) {
            n0 += hist0[i];
            n1 += hist1[i];
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
        for (int j = 0; j < hist0.length; ++j) {
            float y0 = 0;
            float y1 = 0;
            int idx0 = j + bins0Shift;
            if (idx0 > (hist0.length - 1)) {
                idx0 -= hist0.length;
            }
            if (n0 > 0) {
                y0 = (float)hist0[idx0]/(float)n0;
            }
            if (n1 > 0) {
                y1 = (float)hist1[j]/(float)n1;
            }
            sum += Math.min(y0, y1);
            sum0 += y0;
            sum1 += y1;
        }
        
        float d = Math.min(sum0, sum1);
        
        float sim = (d == 0.f) ? 0 : sum/d;
                
        return sim;
    }
    
    /**
     ptImg values for histogram bins:
     0:  red = 0 - 18
     1:  orange = 18 - 40
     2:  yellow = 41 - 60ish
     3:  green = 61 - 106
     4:  blue = 107 - 192
     5:  purple = 193 - 255
        
     * @param ptImg, the "H" image of LCH color space
     * @param labeledSets
     * @return 
     */
    public static List<OneDIntArray> createPTHistograms(GreyscaleImage ptImg, 
        List<TIntSet> labeledSets) {
        
        /*
        ptImg values for histogram bins:
         0:  red = 0 - 18
         1:  orange = 18 - 40
         2:  yellow = 41 - 60ish
         3:  green = 61 - 106
         4:  blue = 107 - 192
         5:  purple = 193 - 255
        */
        
        List<OneDIntArray> output = new ArrayList<OneDIntArray>();
        
        for (int i = 0; i < labeledSets.size(); ++i) {
            
            TIntSet set = labeledSets.get(i);
            TIntIterator iter = set.iterator();
            
            int[] hist = new int[6];
            output.add(new OneDIntArray(hist));
            
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                int v = ptImg.getValue(pixIdx);
                if (v < 19) {
                    hist[0]++;
                } else if (v < 41) {
                    hist[1]++;
                } else if (v < 61) {
                    hist[2]++;
                } else if (v < 107) {
                    hist[3]++;
                } else if (v < 193) {
                    hist[4]++;
                } else {
                    hist[5]++;
                }
            }
        }
        
        return output;
    }

    /**
     ptImg values for histogram bins:
     0:  red = 0 - 18
     1:  orange = 18 - 40
     2:  yellow = 41 - 60ish
     3:  green = 61 - 106
     4:  blue = 107 - 192
     5:  purple = 193 - 255
        
     * @param ptImg, the "H" image of LCH color space
     * @param labeledSet
     * @return 
     */
    public static int[] createPTHistogram(GreyscaleImage ptImg, 
        TIntSet labeledSet) {
        
        /*
        ptImg values for histogram bins:
         0:  red = 0 - 18
         1:  orange = 18 - 40
         2:  yellow = 41 - 60ish
         3:  green = 61 - 106
         4:  blue = 107 - 192
         5:  purple = 193 - 255
        */
        
        TIntIterator iter = labeledSet.iterator();

        int[] hist = new int[6];

        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int v = ptImg.getValue(pixIdx);
            if (v < 19) {
                hist[0]++;
            } else if (v < 41) {
                hist[1]++;
            } else if (v < 61) {
                hist[2]++;
            } else if (v < 107) {
                hist[3]++;
            } else if (v < 193) {
                hist[4]++;
            } else {
                hist[5]++;
            }
        }
        
        return hist;
    }

    /**
     ptImg values for histogram bins:
     0:  red = 0 - 18
     1:  orange = 18 - 40
     2:  yellow = 41 - 60ish
     3:  green = 61 - 106
     4:  blue = 107 - 192
     5:  purple = 193 - 255
        
     * @param ptImg, the "H" image of LCH color space
     * @param labeledSet
     * @return 
     */
    public static int[] createPTHistogram(GreyscaleImage ptImg, 
        Set<PairInt> labeledSet) {
        
        /*
        ptImg values for histogram bins:
         0:  red = 0 - 18
         1:  orange = 18 - 40
         2:  yellow = 41 - 60ish
         3:  green = 61 - 106
         4:  blue = 107 - 192
         5:  purple = 193 - 255
        */
        
        int[] hist = new int[6];

        for (PairInt p : labeledSet) {
            int v = ptImg.getValue(p);
            if (v < 19) {
                hist[0]++;
            } else if (v < 41) {
                hist[1]++;
            } else if (v < 61) {
                hist[2]++;
            } else if (v < 107) {
                hist[3]++;
            } else if (v < 193) {
                hist[4]++;
            } else {
                hist[5]++;
            }
        }
        
        return hist;
    }

}
