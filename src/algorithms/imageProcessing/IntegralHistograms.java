package algorithms.imageProcessing;

import algorithms.imageProcessing.features.HOGUtil;
import gnu.trove.set.TLongSet;
import java.util.Arrays;

/**
 algorithm based upon paper
 "Integral Histogram: A Fast Way to Extract Histograms in
   Cartesian Spaces" by Porikli 2005

 Algorithm uses a pattern similar to a summed area table
 to cumulatively add the content from one lower pixel in
 x and in y with a runtime complexity that is O(N_pixels),
 * but a space complexity that is O(N_pixels * nBins) where
 * nBins is the number of bins in a histogram (default is 16 for
 * default value range of 0 to 255, inclusive).
 
 * @author nichole
 */
public class IntegralHistograms {
        
    /**
     * for a default range in values of 0 to 255, inclusive, and a default
     * bin size of 16, calculate the integral histograms.
     * 
     * @return output two dimensional array with first dimension being 
     * the pixel index and the second being the histogram at that pixel.
     */
    public int[][] create(GreyscaleImage img) {
     
        return create(img, 0, 255, 16);
    }
    
    /**
     * for a default range in values of 0 to 255, inclusive, and a default
     * bin size of 16, calculate the integral histograms.
     * 
     * @return output two dimensional array with first dimension being 
     * the pixel index and the second being the histogram at that pixel.
     */
    public int[][] create(GreyscaleImage img, int nBins) {
             
        return create(img, 0, 255, nBins);
    }
    
    /**
     * NOT TESTED YET
     * 
     * for a default range in values of 0 to 255, inclusive, and a default
     * bin size of 16, calculate the integral histograms.
     * 
     * @param img
     * @param minValue
     * @param maxValue
     * @param nBins
     * @return output two dimensional array with first dimension being 
     * the pixel index and the second being the histogram at that pixel.
     */
    public int[][] create(GreyscaleImage img, int minValue, int maxValue, int nBins) {
        
        TLongSet includePixels = null;
        
        return create(img, includePixels, minValue, maxValue, nBins);
    }
    
    /**
     * NOT TESTED YET
     * 
     * for a default range in values of 0 to 255, inclusive, and a default
     * bin size of 16, calculate the integral histograms.
     * 
     * @param img
     * @param includePixels set of pixel coords to include, but if this is null
     * all pixels are included
     * @param minValue
     * @param maxValue
     * @param nBins
     * @return output two dimensional array with first dimension being 
     * the pixel index and the second being the histogram at that pixel.
     */
    public int[][] create(GreyscaleImage img, TLongSet includePixels,
        int minValue, int maxValue, int nBins) {

        //NOTE: because there is little change between the data in one pixel
        // and the next, it should be possible encode and compress this
        // data structure significantly.
        
        int w = img.getWidth();
        int h = img.getHeight();
        int nPix = img.getNPixels();
        int binWidth = (int)Math.ceil(
            ((float)maxValue - (float)minValue + 1)/(float)nBins);
        
        int[][] out = new int[nPix][];
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int pixIdx = img.getInternalIndex(x, y);
                int bin = (img.getValue(pixIdx) - minValue)/binWidth;
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
     * NOT YET TESTED
     * extract the sum of a window centered at (x,y) of x dimension d and y
     * dimension d and return that summed histogram 
     * and the number of pixels in the
     * aperture in the output variables, output and outputN.
     * NOTE GreyscaleImage, x, and y are in column major format
     * @param integralHistograms
     * @param width the width of the image that integeralHistograms
     * was made from.
     * @param height the height of the image that integeralHistograms
     * was made from.
     * @param x coordinate for x center of window
     * @param y coordinate for y center of window
     * @param d diameter of window in x and y
     * @param output one dimensional array of size nBins in which the
     * sum of the histograms in the window will be returned.
     * @param outputN the number of pixels in the window.  if the window
     * size extends beyond the borders of the "image", this number will be
     * smaller than d * d.
     */
    public void extractWindowFromIntegralHistograms(int[][] integralHistograms, 
        int width, int height, int x, int y, int d, int output[], int[] outputN) {
             
        if (outputN == null || outputN.length != 1) {
            throw new IllegalArgumentException(
                "outputN must be initialized and of size=1");
        }
        if (output == null || output.length != integralHistograms[0].length) {
            throw new IllegalArgumentException(
                "output must be initialized to size nBins which should be equal "
                    + "to size integralHistograms[0].length");
        }
        
        if (d < 0) {
            throw new IllegalArgumentException(
                "d must be a non-negative number");
        }
        
        int w = width;
        int h = height;
        
        if (x < 0 || y < 0 || (x > (w - 1)) || (y > (h - 1))) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "image. x=" + x + " y=" + y + " w=" + w + " h=" + h);
        }
        
        final int r = (d >> 1);
        
        int startX, stopX, startY, stopY;
        if ((r & 1) == 1) {
            startX = x - r;
            stopX = x + r;
            startY = y - r;
            stopY = y + r;
        } else {
            startX = x - r - 1;
            stopX = x + r;
            startY = y - r - 1;
            stopY = y + r;
        }
        if (startX < 0) {
            startX = 0;
        }
        if (startY < 0) {
            startY = 0;
        }
        if (startX >= width) {
            startX = width - 1;
        }
        if (startY >= height) {
            startY = height - 1;
        }
        if (stopX < 0) {
            stopX = 0;
        }
        if (stopY < 0) {
            stopY = 0;
        }
        if (stopX >= width) {
            stopX = width - 1;
        }
        if (stopY >= height) {
            stopY = height - 1;
        }
        
        HOGUtil.extractWindow(integralHistograms, startX, stopX, startY, stopY, 
            w, h, output, outputN);
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
    
}
