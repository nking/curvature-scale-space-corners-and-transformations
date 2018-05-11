package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import gnu.trove.set.TLongSet;
import java.util.Arrays;

/**
 algorithm based upon paper
 "Integral Histogram: A Fast Way to Extract Histograms in
   Cartesian Spaces" by Porikli 2005

 Algorithm uses a pattern similar to a summed area table
 to cumulatively add the content from one lower pixel in
 x and in y with a runtime complexity that is O(N_pixels),
 but a space complexity that is O(N_pixels * nBins) where
 nBins is the number of bins in a histogram (default is 16 for
 default value range of 0 to 255, inclusive).
 
 The logic is for greyscale images that are the polar
 * angles of the X,Y chromaticity
 * of a colorspace such that the brightness is not part
 * of the image measurement.  The image values are color encoding.
 * 
 * 
 * @author nichole
 */
public class PolarThetaIntegralHistograms {
    
    public PolarThetaIntegralHistograms() {
        
    }
    
    /**
     * runtime complexity is O(N_pixels).
     * 
     * @param img greyscale image where intensity is polar angles of the X,Y 
     * chromaticity
     * @param nBins
     * @return 
     */
    public int[][] createHistograms(GreyscaleImage img, int nBins) {
        
        TLongSet includePixels = null;
        
        return createHistograms(img, includePixels, nBins);
    }
    
    /**
     * runtime complexity is O(N_pixels).
     * 
     * @param img greyscale image where intensity is polar angles of the X,Y 
     * chromaticity
     * @param includePixels set of pixel coords to include, but if this is null
     * all pixels are included
     * @param nBins
     * @return 
     */
    public int[][] createHistograms(GreyscaleImage img, TLongSet includePixels,
        int nBins) {

        int w = img.getWidth();
        int h = img.getHeight();
        int nPix = img.getNPixels();
        
        int binWidth = 256/nBins;
        
        int[][] out = new int[nPix][];
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                
                int pixIdx = img.getInternalIndex(x, y);
                
                int v = img.getValue(pixIdx);
                 
                if (v < 0 || v > 255) {
                    throw new IllegalArgumentException("img values muse be "
                        + "between 1 and 255, inclusive");
                }
                
                int bin = v/binWidth;
                if (bin >= nBins) {
                    bin = nBins - 1;
                }
                
                boolean incl = (includePixels == null) || 
                    includePixels.contains(pixIdx);
                    
                if (pixIdx == 0) {
                    
                    out[pixIdx] = new int[nBins];
                    
                    if (incl) {
                        out[pixIdx][bin]++;
                    }
                    
                } else if (x > 0 && y > 0) {
                    
                    int pixIdxL = img.getInternalIndex(x - 1, y);
                    int pixIdxB = img.getInternalIndex(x, y - 1);
                    int pixIdxLB = img.getInternalIndex(x - 1, y - 1);
                    
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                                        
                    if (incl) {
                        out[pixIdx][bin]++;
                    }
                    
                    // add bin, add pixIdxB and subtract pixIdxLB
                    HOGUtil.add(out[pixIdx], out[pixIdxB]);
                    HOGUtil.subtract(out[pixIdx], out[pixIdxLB]);
                    
                } else if (x > 0) {
                    
                    int pixIdxL = img.getInternalIndex(x - 1, y);
                                   
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                    if (incl) {
                        out[pixIdx][bin]++;
                    }
                
                } else if (y > 0) {
                
                    int pixIdxB = img.getInternalIndex(x, y - 1);
                                    
                    out[pixIdx] = Arrays.copyOf(out[pixIdxB], nBins);
                    
                    if (incl) {
                        out[pixIdx][bin]++;
                    }
                }
            }
        }

        return out;
    }
    
    /**
     * 
     * runtime complexity is O(nBins)
     * 
     * @param startX
     * @param stopX
     * @param startY
     * @param stopY
     * @param output
     * @param outputN 
     */
    public void extractWindow(int[][] histograms, int startX, int stopX, 
        int startY, int stopY, int w, int h, 
        int output[], int[] outputN) {

        HOGUtil.extractWindow(histograms, startX, stopX, startY, stopY, w, h, 
            output, outputN);
    }
    
    /**
     * apply a windowed sum across the integral histogram image,
     * where the window size is N_PIX_PER_CELL_DIM.
     *  
     * @param histograms
     * @param w
     * @param h
     * @param N_PIX_PER_CELL_DIM
     */    
    public void applyWindowedSum(int[][] histograms, int w, int h, 
        int N_PIX_PER_CELL_DIM) {
        
        if (N_PIX_PER_CELL_DIM < 1) {
            throw new IllegalArgumentException("N_PIX_PER_CELL_DIM must be >= 1");
        }
        
        HOGUtil.applyWindowedSum(histograms, w, h, N_PIX_PER_CELL_DIM);
    }
    
}
