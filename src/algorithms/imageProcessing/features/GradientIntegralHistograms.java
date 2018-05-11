package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import algorithms.util.PixelHelper;
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
                                        
                    out[pixIdx][b0] += (f0 * v);
                    out[pixIdx][b1] += (f1 * v);
                    
                    // add bin, add pixIdxB and subtract pixIdxLB
                    HOGUtil.add(out[pixIdx], out[pixIdxB]);
                    HOGUtil.subtract(out[pixIdx], out[pixIdxLB]);
                    
                } else if (x > 0) {
                    
                    int pixIdxL = gradient.getInternalIndex(x - 1, y);
                         
                    out[pixIdx] = Arrays.copyOf(out[pixIdxL], nBins);
                    
                    out[pixIdx][b0] += (f0 * v);
                    out[pixIdx][b1] += (f1 * v);
                    
                } else if (y > 0) {
                
                    int pixIdxB = gradient.getInternalIndex(x, y - 1);
                     
                    out[pixIdx] = Arrays.copyOf(out[pixIdxB], nBins);
                    
                    out[pixIdx][b0] += (f0 * v);
                    out[pixIdx][b1] += (f1 * v);
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
        
        if (N_PIX_PER_CELL_DIM < 1) {
            throw new IllegalArgumentException("N_PIX_PER_CELL_DIM must be >= 1");
        }
        
        int nBins = histograms[0].length;
                
        int[][] img2 = new int[w * h][];
        
        int[] outN = new int[1];
        
        int windowSize = N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM;
        
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
        
        float factor;
                        
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
            if (x3 >= w) {
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
                if (y3 >= h) {
                    y3 = h - 1;
                }
                                
                int pixIdx = (y * w) + x;
                
                img2[pixIdx] = new int[nBins];

                extractWindow(histograms, x2, x3, y2, y3, w, h, img2[pixIdx], outN);
                
                if (outN[0] < windowSize) {
                    factor = (float)windowSize/(float)outN[0];
                    HOGUtil.mult(img2[pixIdx], factor);
                }             
            }
        }
        
        return img2;
    }
    
    int[][] transformIntoIntegral2DHist(int[][] hist, int imageWidth,
        int imageHeight) {
        
        int nPix = imageWidth * imageHeight;
        PixelHelper ph = new PixelHelper();
        
        int[][] out = new int[nPix][];
        for (int i = 0; i < nPix; ++i) {
            out[i] = Arrays.copyOf(hist[i], hist[i].length);
        }
        
        int[] tmp;
        int pixIdx;
                
        for (int x = 0; x < imageWidth; ++x) {
            for (int y = 0; y < imageHeight; ++y) {
                
                pixIdx = (int)ph.toPixelIndex(x, y, imageWidth);
                
                tmp = out[pixIdx];
                                
                if (x > 0 && y > 0) {
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x - 1, y, imageWidth)]);
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x, y - 1, imageWidth)]);
                    HOGUtil.subtract(tmp, out[(int)ph.toPixelIndex(x - 1, y - 1, 
                        imageWidth)]);
                    //int v = out.getValue(x, y)
                    //        + out.getValue(x - 1, y) 
                    //        + out.getValue(x, y - 1) 
                    //        - out.getValue(x - 1, y - 1);
                    //out.setValue(x, y, v);
                } else if (x > 0) {
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x - 1, y, imageWidth)]);
                    //int v = out.getValue(x, y) 
                    //        + out.getValue(x - 1, y);
                    //out.setValue(x, y, v);
                } else if (y > 0) {
                    HOGUtil.add(tmp, out[(int)ph.toPixelIndex(x, y - 1, imageWidth)]);
                    //int v = out.getValue(x, y)
                    //        + out.getValue(x, y - 1);
                    //out.setValue(x, y, v);
                }
            }
        }
        return out;
    }

}
