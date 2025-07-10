package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
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
 
 * Added to the algorithm is a bilinear splitting of values into 2 histograms
 * bins instead of 1.
 * 
 * @author nichole
 */
public class GradientIntegralHistograms {
    
    public GradientIntegralHistograms() {
        
    }
    
    // t, b0, b1, f0, f1;
    void calculateBinsAndFractions(int x, int y, GreyscaleImage theta, int nBins,
        double[] out) {
    
        float binWidth = 180.f/nBins;
        int w = theta.getWidth();
        int h = theta.getHeight();
        
        double b, f0, f1;
        double halfBinWidth = binWidth/2.;
        int pixIdx, t, b0, b1, c0, v, v0, v1;
        
        pixIdx = theta.getInternalIndex(x, y);
                
        t = theta.getValue(pixIdx);

        if (t < 0 || t > 180) {
            throw new IllegalArgumentException("theta values muse be "
                + "between 1 and 179, inclusive");
        }
        if (t == 180) {
            t = 0;
        }
        v = theta.getValue(pixIdx);
        assert(v >= 0.);

        b = (double)t/binWidth;
        b0 = (int)b;
        c0 = (int)(((double)b0 * binWidth) + halfBinWidth);

        f1 = Math.abs((double)(c0 - t)/binWidth);
        f0 = 1. - f1;
        assert(f1 >= 0.);
        assert(f0 >= 0.);

        if (t == c0) {
            // t is centered in bin
            b1 = b0;
        } else if (t < c0) {
            b1 = b0 - 1;
            if (b1 < 0) {
                b1 = nBins - 1;
            }
        } else {
            b1 = b0 + 1;
            if (b1 == nBins) {
                b1 = 0;
            }
        }
        out[0] = t;
        out[1] = b0;
        out[2] = b1;
        out[3] = f0;
        out[4] = f1;
    }
    
    /**
     * runtime complexity is O(N_pixels).
     * 
     * @param gradient
     * @param theta
     * @param nBins
     * @return 2D histogram with first dimension being the pixel index and 
     * 2nd dimension being the angle based histogram for the pixel.
     */
    public int[][] createHistograms(GreyscaleImage gradient,  
        GreyscaleImage theta, int nBins) {

        int w = gradient.getWidth();
        int h = gradient.getHeight();
        
        if (w != theta.getWidth() || h != theta.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }
        
        // if min < 0, creates a new copy and adds a bias level of min
        gradient = HOGUtil.applyBiasLevelToMakePositive(gradient);
        
        int nPix = gradient.getNPixels();
                
        /*
        NOTE: have changed to place the counts into the 2 bins that the angle is closest
        to, that is, a binlinear interpolation of the contribution to theta.
        */
        
        int[][] out = new int[nPix][];        
        double[] bilinearParams = new double[5];
        double f0, f1;
        int pixIdx, t, v, b0, b1;
                
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                
                pixIdx = gradient.getInternalIndex(x, y);
                
                v = gradient.getValue(pixIdx);
                assert(v >= 0.);
                
                // t, b0, b1, f0, f1;
                calculateBinsAndFractions(x, y, theta, nBins, bilinearParams);
            
                t = (int)bilinearParams[0];
                b0 = (int)bilinearParams[1];
                b1 = (int)bilinearParams[2];
                f0 = bilinearParams[3];
                f1 = bilinearParams[4];
                
                assert(Math.abs((f0 + f1) - 1) < 0.01);
                
                if (pixIdx == 0) {
                    //x==0 && y==0
                    out[pixIdx] = new int[nBins];
                    
                    out[pixIdx][b0] += (int)Math.round(f0 * v);
                    out[pixIdx][b1] += (int)Math.round(f1 * v);
                    
                } else if (x > 0 && y > 0) {
                    
                    int pixIdxL = gradient.getInternalIndex(x - 1, y);
                    int pixIdxB = gradient.getInternalIndex(x, y - 1);
                    int pixIdxLB = gradient.getInternalIndex(x - 1, y - 1);
                    
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                                        
                    out[pixIdx][b0] += (int)Math.round(f0 * v);
                    out[pixIdx][b1] += (int)Math.round(f1 * v);
                    
                    // add bin, add pixIdxB and subtract pixIdxLB
                    HOGUtil.add(out[pixIdx], out[pixIdxB]);
                    HOGUtil.subtract(out[pixIdx], out[pixIdxLB]);
                    
                } else if (x > 0) {
                    
                    int pixIdxL = gradient.getInternalIndex(x - 1, y);
                         
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                    
                    out[pixIdx][b0] += (int)Math.round(f0 * v);
                    out[pixIdx][b1] += (int)Math.round(f1 * v);
                    
                } else if (y > 0) {
                
                    int pixIdxB = gradient.getInternalIndex(x, y - 1);
                     
                    out[pixIdx] = Arrays.copyOf(out[pixIdxB], nBins);
                    
                    out[pixIdx][b0] += (int)Math.round(f0 * v);
                    out[pixIdx][b1] += (int)Math.round(f1 * v);
                }
            }
        }

        return out;
    }
    
    /**
     * extract the sum of histograms in the window inclusively defined as
     * (startX:stopX, startY:stopY).
     * 
     * runtime complexity is O(nBins)
     * 
     * @param histograms the 2D histogram integral image.  the first dimension 
     *    is the pixel index and the 2nd is the histogram bin, e.g.
     *    histograms[pixIdx][binIdx].
     * @param startX
     * @param stopX the last x pixel in the window, inclusive
     * @param startY
     * @param stopY the last y pixel in the window, inclusive
     * @param output
     * @param outputN an empty 1 dimensional array of size 1 to return the 
     * number of pixels in the cell
     */
    public void extractWindow(int[][] histograms, int startX, int stopX, 
        int startY, int stopY, int w, int h, 
        int[] output, int[] outputN) {

        HOGUtil.extractWindow(histograms, startX, stopX, startY, stopY, w, h, 
            output, outputN);
    }
    
    /**
     * extract the sum of histograms in the window inclusively defined as
     * (startX:stopX, startY:stopY).
     * 
     * runtime complexity is O(nBins)
     * 
     * @param histograms the 2D histogram integral image.  the first dimension 
     *    is the pixel index and the 2nd is the histogram bin, e.g.
     *    histograms[pixIdx][binIdx].
     * @param startX
     * @param stopX the last x pixel in the window, inclusive
     * @param startY
     * @param stopY the last y pixel in the window, inclusive
     * @param output
     * @param outputN an empty 1 dimensional array of size 1 to return the 
     * number of pixels in the cell
     */
    public void extractWindow(int[][] histograms, int startX, int stopX, 
        int startY, int stopY, int w, int h, 
        long[] output, int[] outputN) {

        HOGUtil.extractWindow(histograms, startX, stopX, startY, stopY, w, h, 
            output, outputN);
    }
    
    /**
     * apply a windowed sum across the gradient integral histogram image,
     * where the window size is N_PIX_PER_CELL_DIM.
     * The result is an image of histograms, where each histogram represents
     * the integration over the surrounding N_PIX_PER_CELL_DIM window.
     * For windows near the image edge, a factor is applied to bring the
     * counts up by factor (N_PIX_PER_CELL_DIM/n_pix_in_window).
     * The resulting 2D histograms are then made into an integral image again.
     *  
     * @param histograms
     * @param w
     * @param h
     * @param N_PIX_PER_CELL_DIM
     */    
    public void applyWindowedSum(int[][] histograms, int w, int h, 
        int N_PIX_PER_CELL_DIM) {
        
        int[][] img2 = applyWindowedSum0(histograms, w, h, N_PIX_PER_CELL_DIM);
        
        img2 = transformIntoIntegral2DHist(img2, w, h);
        
        for (int i = 0; i < histograms.length; ++i) {
            System.arraycopy(img2[i], 0, histograms[i], 0, img2[i].length);
        }
    }
    
    /**
     * apply a windowed sum across the gradient integral histogram image,
     * where the window size is N_PIX_PER_CELL_DIM.
     * The result is an image of histograms, where each histogram represents
     * the integration over the surrounding N_PIX_PER_CELL_DIM window.
     * The result is NOT a 2D integral histogram image, just a 2D histogram image;
     *  
     * @param histograms
     * @param w
     * @param h
     * @param N_PIX_PER_CELL_DIM
     */    
    int[][] applyWindowedSum0(int[][] histograms, int w, int h, 
        int N_PIX_PER_CELL_DIM) {
        
        return HOGUtil.applyWindowedSum0(histograms, w, h, N_PIX_PER_CELL_DIM);
    }
    
    int[][] transformIntoIntegral2DHist(int[][] hist, int imageWidth,
        int imageHeight) {
        
        return HOGUtil.transformIntoIntegral2DHist(hist, imageWidth, 
            imageHeight);
    }

}
