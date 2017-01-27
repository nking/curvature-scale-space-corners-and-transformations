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
 
 * @author nichole
 */
public class GradientIntegralHistograms {
    
    public GradientIntegralHistograms() {
        
    }
    
    /**
     * runtime complexity is O(N_pixels).
     * 
     * @param gradient
     * @param theta
     * @param nBins
     * @return 
     */
    public int[][] createHistograms(GreyscaleImage gradient,  
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

        if (output.length != histograms[0].length) {
            throw new IllegalArgumentException("output.length must == nThetaBins");
        }
        
        if (stopX < startX || stopY < startY) {
            throw new IllegalArgumentException("stopX must be >= startX and "
                + "stopY >= startY");
        }
        
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
                System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                    output, 0, output.length);
                subtract(output, histograms[getPixIdx(startX - 1, stopY, w)]);
                subtract(output, histograms[getPixIdx(stopX, startY - 1, w)]);
                add(output, histograms[getPixIdx(startX - 1, startY - 1, w)]);
                return;
            }
        }
        
        if (dx == 0 && dy == 0) {
            if (stopX == 0) {
                if (stopY == 0) {
                    outputN[0] = 1;
                    System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                        output, 0, output.length);
                    return;
                }
                outputN[0] = 1;
                System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                    output, 0, output.length);
                subtract(output, histograms[getPixIdx(stopX, stopY - 1, w)]);
                return;
            } else if (stopY == 0) {
                //stopX == 0 && stopY == 0 has been handles in previous block
                outputN[0] = 1;
                System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                    output, 0, output.length);
                subtract(output, histograms[getPixIdx(stopX - 1, stopY, w)]);
                return;
            } else {
                // stopX > 0
                outputN[0] = 1;
                System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                    output, 0, output.length);
                subtract(output, histograms[getPixIdx(startX - 1, stopY, w)]);
                subtract(output, histograms[getPixIdx(stopX, startY - 1, w)]);
                add(output, histograms[getPixIdx(startX - 1, startY - 1, w)]);
                return;
            }
            //System.out.println(" --> startX=" + startX +
            //    " stopX=" + stopX + " startY=" + startY + " stopY=" + stopY);            
        }
        
        if (startX > 0 && startY > 0) {
            int nPix = (dx + 1) * (dy + 1);
            outputN[0] = nPix;
            System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                output, 0, output.length);
            subtract(output, histograms[getPixIdx(startX, stopY, w)]);
            subtract(output, histograms[getPixIdx(stopX, startY, w)]);
            add(output, histograms[getPixIdx(startX, startY, w)]);
            return;
        } else if (startX > 0) {
            // startY is < 0
            int nPix = (dx + 1) * (stopY + 1);
            outputN[0] = nPix;
            System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                output, 0, output.length);
            subtract(output, histograms[getPixIdx(startX - 1, stopY, w)]);
            return;
        } else if (startY > 0) {
            // startX < 0
            int nPix = (stopX + 1) * (dy + 1);
            outputN[0] = nPix;
            System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                output, 0, output.length);
            subtract(output, histograms[getPixIdx(stopX, startY - 1, w)]);
            return;
        } else {
            // startX < 0 && startY < 0
            int nPix = (stopX + 1) * (stopY + 1);
            outputN[0] = nPix;
            System.arraycopy(histograms[getPixIdx(stopX, stopY, w)], 0, 
                output, 0, output.length);
            return;
        }
    }
    
    /**
     * apply a windowed sum across the gradient integral histogram image,
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
        
        int nBins = histograms[0].length;
                
        int[][] img2 = new int[w * h][];
        
        int[] outN = new int[1];
        
        // a centered window sum
        int r = N_PIX_PER_CELL_DIM >> 1;
        int r0, r1;
        if (r == 0) {
            r0 = 0;
            r1 = 0;
        } else if ((N_PIX_PER_CELL_DIM & 1) == 1) {
            r0 = -r;
            r1 = r;
        } else {
            r0 = -r;
            r1 = r - 1;
        }
                        
        // extract the summed area of each dxd window centered on x,y
        for (int x = 0; x < w; ++x) {
            
            int x2 = x + r0;
            int x3 = x + r1;
            if (x3 < 0) {
                continue;
            } else if (x2 < 0) {
                x2 = 0;
            } else if (x2 >= w) {
                break;
            }
            if (x3 > (w - 1)) {
                x3 = w - 1;
            }

            for (int y = 0; y < h; ++y) {
                
                int y2 = y + r0;
                int y3 = y + r1;
                if (y3 < 0) {
                    continue;
                } else if (y2 < 0) {
                    y2 = 0;
                } else if (y2 >= h) {
                    break;
                }
                if (y3 > (h - 1)) {
                    y3 = h - 1;
                }
                                
                int pixIdx = (y * w) + x;

                int[] out = new int[nBins];

                extractWindow(histograms, x2, x3, y2, y3, w, h, out, outN);
                
                img2[pixIdx] = out;
            }
        }
        
        for (int i = 0; i < histograms.length; ++i) {
            System.arraycopy(img2[i], 0, histograms[i], 0, nBins);
        }
    }
    
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
    
    protected int getPixIdx(int col, int row, int w) {
        return (row * w) + col;
    }
    
}
