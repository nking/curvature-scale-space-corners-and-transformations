package algorithms.imageProcessing;

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
     * NOT TESTED YET
     * 
     * for a default range in values of 0 to 255, inclusive, and a default
     * bin size of 16, calculate the integral histograms.
     * 
     * @param img
     * @param minValue
     * @param maxValue
     * @param binWidth
     * @return output two dimensional array with first dimension being 
     * the pixel index and the second being the histogram at that pixel.
     */
    public int[][] create(GreyscaleImage img, int minValue, int maxValue,
        int binWidth) {

        //NOTE: because there is little change between the data in one pixel
        // and the next, it should be possible encode and compress this
        // data structure significantly.
        
        int w = img.getWidth();
        int h = img.getHeight();
        int nPix = img.getNPixels();
        int nBins = (maxValue - minValue + 1)/binWidth;
        
        int[][] out = new int[nPix][];
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int pixIdx = img.getInternalIndex(x, y);
                int bin = (img.getValue(pixIdx) - minValue)/binWidth;
                if (pixIdx == 0) {
                    out[pixIdx] = new int[nBins];
                    out[pixIdx][bin]++;
                } else if (x > 0 && y > 0) {
                    int pixIdxL = img.getInternalIndex(x - 1, y);
                    int pixIdxB = img.getInternalIndex(x, y - 1);
                    int pixIdxLB = img.getInternalIndex(x - 1, y - 1);
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                    out[pixIdx][bin]++;
                    // add bin, add pixIdxB and subtract pixIdxLB
                    add(out[pixIdx], out[pixIdxB]);
                    subtract(out[pixIdx], out[pixIdxLB]);
                } else if (x > 0) {
                    int pixIdxL = img.getInternalIndex(x - 1, y);
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                    out[pixIdx][bin]++;
                } else if (y > 0) {
                    int pixIdxB = img.getInternalIndex(x, y - 1);
                    out[pixIdx] = Arrays.copyOf(out[pixIdxB], nBins);
                    out[pixIdx][bin]++;
                }
            }
        }

        return out;
    }

    private void add(int[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    private void subtract(int[] subtractFrom, int[] subtract) {
        for (int i = 0; i < subtractFrom.length; ++i) {
            subtractFrom[i] += subtract[i];
        }
    }
    
}
