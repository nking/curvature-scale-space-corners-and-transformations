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
 * but a space complexity that is O(N_pixels * nBins) where
 * nBins is the number of bins in a histogram (default is 16 for
 * default value range of 0 to 255, inclusive).
 
 * @author nichole
 */
public class GradientIntegralHistograms {
    
    private final int[][] histograms;
    
    private final int w;
    private final int h;
    
    /**
     * assumes a theta range of 0 to 180 and creates a histogram image
     * with nThetaBins and using a frequency unit calculated as the
     * gradient value, that is a gradient of 124 adds a count of 124
     * to the histogram bin for theta value for its pixel.
     * 
     * Note that the values in gradient are expected to be in range 0 to 255,
     * inclusive.
     * 
     * runtime complexity is O(N_pixels).
     * 
     * @param gradient values in range 0 to 255, inclusive.
     * @param theta angle of gradient in range 0 to 180.
     * @param nThetaBins
     */
    public GradientIntegralHistograms(GreyscaleImage gradient,  GreyscaleImage theta,
        int nThetaBins) {
        
        this.histograms = createHistograms(gradient, theta, nThetaBins);
    
        this.w = gradient.getWidth();
        this.h = gradient.getHeight();
    }
        
    /**
     * runtime complexity is O(N_pixels).
     * 
     * @param gradient
     * @param theta
     * @param nBins
     * @return 
     */
    private int[][] createHistograms(GreyscaleImage gradient,  
        GreyscaleImage theta, int nBins) {

        int w = gradient.getWidth();
        int h = gradient.getHeight();
        
        if (w != theta.getWidth() || h != theta.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }
        
        int nPix = gradient.getNPixels();
        
        int binWidth = 180/nBins;
        
        int[][] out = new int[nPix][];
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                
                int pixIdx = gradient.getInternalIndex(x, y);
                
                int t = theta.getValue(pixIdx);
                 
                if (t < 0 || t > 180) {
                    throw new IllegalArgumentException("theta values muse be "
                        + "between 1 and 179, inclusive");
                }
                if (t == 180) {
                    t = 0;
                }
                
                int bin = t/binWidth;
                
                int v = gradient.getValue(pixIdx);
               
                if (pixIdx == 0) {
                    
                    out[pixIdx] = new int[nBins];
                    
                    out[pixIdx][bin] += v;
                    
                } else if (x > 0 && y > 0) {
                    
                    int pixIdxL = gradient.getInternalIndex(x - 1, y);
                    int pixIdxB = gradient.getInternalIndex(x, y - 1);
                    int pixIdxLB = gradient.getInternalIndex(x - 1, y - 1);
                    
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                                        
                    out[pixIdx][bin] += v;
                    
                    // add bin, add pixIdxB and subtract pixIdxLB
                    add(out[pixIdx], out[pixIdxB]);
                    subtract(out[pixIdx], out[pixIdxLB]);
                    
                } else if (x > 0) {
                    
                    int pixIdxL = gradient.getInternalIndex(x - 1, y);
                                   
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                    out[pixIdx][bin] += v;
                
                } else if (y > 0) {
                
                    int pixIdxB = gradient.getInternalIndex(x, y - 1);
                                    
                    out[pixIdx] = Arrays.copyOf(out[pixIdxB], nBins);
                    
                    out[pixIdx][bin] += v;
                }
            }
        }

        return out;
    }
    
    /**
     * NOT READY FOR USE. 
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
    public void extractWindow(int startX, int stopX, int startY, int stopY, 
        int output[], int[] outputN) {
        
        int dx, dy;
        if (startX > -1) {
            dx = stopX - startX;
        } else {
            dx = stopX;
        }
        if (startY > -1) {
            dy = stopY - startY;
        } else {
            dy = stopY;
        }
        
        if (dx > 0 && dy > 0) {
            if ((startX > 0) && (stopX < w) && (startY > 0) && (stopY < h)) {
                int nPix = (dx + 1) * (dy + 1);
                outputN[0] = nPix;
                add(output, histograms[getPixIdx(stopX, stopY)]);
                subtract(output, histograms[getPixIdx(startX - 1, stopY)]);
                subtract(output, histograms[getPixIdx(stopX, startY - 1)]);
                add(output, histograms[getPixIdx(startX - 1, startY - 1)]);
                return;
            }
        }
        
        if (dx == 0 && dy == 0) {
            if (stopX == 0) {
                if (stopY == 0) {
                    outputN[0] = 1;
                    add(output, histograms[getPixIdx(stopX, stopY)]);
                    return;
                }
                outputN[0] = 1;
                add(output, histograms[getPixIdx(stopX, stopY)]);
                subtract(output, histograms[getPixIdx(stopX, stopY - 1)]);
                return;
            } else if (stopY == 0) {
                //stopX == 0 && stopY == 0 has been handles in previous block
                outputN[0] = 1;
                add(output, histograms[getPixIdx(stopX, stopY)]);
                subtract(output, histograms[getPixIdx(stopX - 1, stopY)]);
                return;
            } else {
                // stopX > 0
                outputN[0] = 1;
                add(output, histograms[getPixIdx(stopX, stopY)]);
                subtract(output, histograms[getPixIdx(startX - 1, stopY)]);
                subtract(output, histograms[getPixIdx(stopX, startY - 1)]);
                add(output, histograms[getPixIdx(startX - 1, startY - 1)]);
                return;
            }
            //System.out.println(" --> startX=" + startX +
            //    " stopX=" + stopX + " startY=" + startY + " stopY=" + stopY);            
        }
        
        if (startX > 0 && startY > 0) {
            int nPix = (dx + 1) * (dy + 1);
            outputN[0] = nPix;
            add(output, histograms[getPixIdx(stopX, stopY)]);
            subtract(output, histograms[getPixIdx(startX, stopY)]);
            subtract(output, histograms[getPixIdx(stopX, startY)]);
            add(output, histograms[getPixIdx(startX, startY)]);
            return;
        } else if (startX > 0) {
            // startY is < 0
            int nPix = (dx + 1) * (stopY + 1);
            outputN[0] = nPix;
            add(output, histograms[getPixIdx(stopX, stopY)]);
            subtract(output, histograms[getPixIdx(startX - 1, stopY)]);
            return;
        } else if (startY > 0) {
            // startX < 0
            int nPix = (stopX + 1) * (dy + 1);
            outputN[0] = nPix;
            add(output, histograms[getPixIdx(stopX, stopY)]);
            subtract(output, histograms[getPixIdx(stopX, startY - 1)]);
            return;
        } else {
            // startX < 0 && startY < 0
            int nPix = (stopX + 1) * (stopY + 1);
            outputN[0] = nPix;
            add(output, histograms[getPixIdx(stopX, stopY)]);
            return;
        }
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
    /*public void extractWindowFromIntegralHistograms(int[][] integralHistograms, 
        int width, int height, int x, int y, int d, int output[], int[] outputN) {
     
        //TODO: rewrite using simpler pattern used in SummedAreaTable
        
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
        
        Arrays.fill(output, 0);
        
        // extract the summed histograms in area of dxd window centered on x,y
        if (r > 0) {
            if (x > r && x < (w-r) && (y > r) && (y < (h-r))) {
                outputN[0] = d * d;
                add(output, integralHistograms[getPixIdx(x+r, y+r, width)]);
                subtract(output, integralHistograms[getPixIdx(x-r, y+r, width)]);
                subtract(output, integralHistograms[getPixIdx(x+r, y-r, width)]);
                add(output, integralHistograms[getPixIdx(x-r, y-r, width)]);
                return;
            }
        }
                
        // handling borders separately
        
        int startX = x - r - 1;
        int stopX = x + r;
        int startY = y - r - 1;
        int stopY = y + r;
        
        if (stopX > (w - 1)) {
            stopX = w - 1;
        }
        if (stopY > (h - 1)) {
            stopY = h - 1;
        }
        
        //System.out.println("x=" + x + " y=" + y + " r=" + r
        //    + " startX=" + startX +
        //    " stopX=" + stopX + " startY=" + startY + " stopY=" + stopY);
       
        // when r == 0, bounds need another edit or immediate return
        //
        //2            2           2           2           2       *
        //1            1 *         1           1    *      1
        //0 *          0           0    *      0           0
        //  0  1  2      0  1  2     0  1  2     0  1  2     0  1  2
        //
        if (r == 0) {
            if (stopX == 0) {
                if (stopY == 0) {
                    outputN[0] = 1;
                    add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
                    return;
                }
                startY = stopY - 1;
                if (startY == 0) {
                    outputN[0] = 1;
                    add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
                    subtract(output, integralHistograms[getPixIdx(stopX, startY, width)]);
                    return;
                }
            } else {
                // stopX > 0
                startX = stopX - 1;
                if (stopY == 0) {
                    outputN[0] = 1;
                    add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
                    subtract(output, integralHistograms[getPixIdx(startX, stopY, width)]);
                    return;
                }
                startY = stopY - 1;
            }
            //System.out.println(" --> startX=" + startX +
            //    " stopX=" + stopX + " startY=" + startY + " stopY=" + stopY);            
        }
        
        if (startX >= 0 && startY >= 0) {
            int nPix = (r == 0) ? 1 : (stopX - startX) * (stopY - startY);
            outputN[0] = nPix;
            add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
            subtract(output, integralHistograms[getPixIdx(startX, stopY, width)]);
            subtract(output, integralHistograms[getPixIdx(stopX, startY, width)]);
            add(output, integralHistograms[getPixIdx(startX, startY, width)]);
            return;
        } else if (startX >= 0) {
            // startY is < 0
            int nPix = (r == 0) ? 1 : (stopX - startX) * (stopY + 1);
            outputN[0] = nPix;
            add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
            subtract(output, integralHistograms[getPixIdx(startX, stopY, width)]);
            return;
        } else if (startY >= 0) {
            // startX < 0
            int nPix = (r == 0) ? 1 : (stopX + 1) * (stopY - startY);
            outputN[0] = nPix;
            add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
            subtract(output, integralHistograms[getPixIdx(stopX, startY, width)]);
            return;
        } else {
            // startX < 0 && startY < 0
            int nPix = (r == 0) ? 1 : (stopX + 1) * (stopY + 1);
            outputN[0] = nPix;
            add(output, integralHistograms[getPixIdx(stopX, stopY, width)]);
            return;
        }
    }
    */
   
    protected void add(int[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    protected void subtract(int[] subtractFrom, int[] subtract) {
        for (int i = 0; i < subtractFrom.length; ++i) {
            subtractFrom[i] -= subtract[i];
        }
    }
    
    protected int getPixIdx(int col, int row) {
        return (row * w) + col;
    }
    
}
