package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.IntegralHistograms;
import algorithms.util.OneDIntArray;
import algorithms.util.OneDLongArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 A class holding histograms of a greyscale image 
 * and methods to extract
 * features.
 * The algorithm is similar to those mentioned in HOGs.java except
 * that the histograms below have bins of color space and the added
 * units are "1" unit of the pixel color.
 
  @author nichole
*/
public class HGS {
   
    private static float eps = 0.000001f;
    
    private final int nBins;
    
    private final int N_PIX_PER_CELL_DIM;
    
    private final int N_CELLS_PER_BLOCK_DIM;
    
    // histogrm integral images with a windowed sum of N_PIX_PER_CELL_DIM
    private final int[][] gHists;
    
    private final int w;
    private final int h;
    
    private boolean debug = false;
    
    //TODO: calculate the limits in nPixels this can handle due to
    //   using integers instead of long for storage.
    //  8.4 million pix, roughly 2900 X 2900
    
    /**
     * constructor
     * @param img gradient image or greyscale
     */
    public HGS(GreyscaleImage img) {
        
        // binWidth of 16
        nBins = 16;
        N_PIX_PER_CELL_DIM = 4;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = img.getWidth();
        h = img.getHeight();
        
        gHists = init(img);
    }
    
    public HGS(GreyscaleImage img, int nCellsPerDim, int nPixPerCellDim,
        int nBins) {
        
        // binWidth of 16
        this.nBins = nBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = img.getWidth();
        h = img.getHeight();
        
        gHists = init(img);
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    private int[][] init(GreyscaleImage img) {
       
        IntegralHistograms gh = new IntegralHistograms();
        
        int[][] histograms = gh.create(img, nBins);

        //apply a windowed avg across the integral image
        gh.applyWindowedSum(histograms, w, h, N_PIX_PER_CELL_DIM);
        
        return histograms;  
    }

    /**     * 
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * This uses the block normalization of Dalal & Triggs.
     * 
     * @param x
     * @param y
     * @param outHist 
     */
    public void extractBlock(int x, int y, int[] outHist) {
                
        if (outHist.length != nBins) {
            throw new IllegalArgumentException("outHist.length != nBins");
        }

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normaliza each cell by that sum.
        
        int nH = N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM;

        double blockTotal = 0;        
                
        List<OneDIntArray> cells = new ArrayList<OneDIntArray>(nH);
        
        for (int cX = 0; cX < N_CELLS_PER_BLOCK_DIM; ++cX) {
            
            int cXOff = -(N_CELLS_PER_BLOCK_DIM/2) + cX;
        
            int x2 = x + (cXOff * N_PIX_PER_CELL_DIM);
            
            if ((x2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                break;
            } else if (x2 < 0) {
                x2 = 0;
            } else if (x2 >= w) {
                break;
            }
            
            for (int cY = 0; cY < N_CELLS_PER_BLOCK_DIM; ++cY) {
                    
                int cYOff = -(N_CELLS_PER_BLOCK_DIM/2) + cY;
                
                int y2 = y + (cYOff * N_PIX_PER_CELL_DIM);

                if ((y2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                    break;
                } else if (y2 < 0) {
                    y2 = 0;
                } else if (y2 >= h) {
                    break;
                }
                                
                int pixIdx = (y2 * w) + x2;
                
                int[] out = Arrays.copyOf(gHists[pixIdx], gHists[pixIdx].length);

                cells.add(new OneDIntArray(out));
            
                int t = sumCounts(out);
                
                blockTotal += (t * t);
            }
        }
        
        blockTotal /= (double)cells.size();
        blockTotal = Math.sqrt(blockTotal);
        double norm = 1./(blockTotal + eps);
                
        float maxBlock = 255.f;
            //(N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM) *
            //(N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM) * 255,f;
   
        norm *= maxBlock;
        
        Arrays.fill(outHist, 0, outHist.length, 0);
        
        for (int i = 0; i < cells.size(); ++i) {
            int[] a = cells.get(i).a;
            for (int j = 0; j < a.length; ++j) {
                //v /= Math.sqrt(blockTotal + 0.0001);
                a[j] = (int)Math.round(norm * a[j]);
            } 
            add(outHist, a);
        }
        
        /*        
        part of a block of 3 X 3 cells
        
           2        2        2        2
           1        1        1        1
           0        0        0        0
          -9 -8 -7 -6 -5 -4 -3 -2 -1  *  1  2  3  4  5  6  7  9
                                      *
        */        
    }
    
    /**
     * NOT READY FOR USE
     * 
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * This uses the block normalization of Dalal & Triggs.
     * 
     * @param x
     * @param y
     * @param outHist 
     */
    public void extractBlock(int x, int y, long[] outHist) {
                
        if (outHist.length != nBins) {
            throw new IllegalArgumentException("outHist.length != nBins");
        }

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normaliza each cell by that sum.
        
        int nH = N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM;

        double blockTotal = 0;        
                
        List<OneDLongArray> cells = new ArrayList<OneDLongArray>(nH);
        
        for (int cX = 0; cX < N_CELLS_PER_BLOCK_DIM; ++cX) {
            
            int cXOff = -(N_CELLS_PER_BLOCK_DIM/2) + cX;
        
            int x2 = x + (cXOff * N_PIX_PER_CELL_DIM);
            
            if ((x2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                break;
            } else if (x2 < 0) {
                x2 = 0;
            } else if (x2 >= w) {
                break;
            }
            
            for (int cY = 0; cY < N_CELLS_PER_BLOCK_DIM; ++cY) {
                    
                int cYOff = -(N_CELLS_PER_BLOCK_DIM/2) + cY;
                
                int y2 = y + (cYOff * N_PIX_PER_CELL_DIM);

                if ((y2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                    break;
                } else if (y2 < 0) {
                    y2 = 0;
                } else if (y2 >= h) {
                    break;
                }
                                
                int pixIdx = (y2 * w) + x2;
                
                long[] out = copyHist(pixIdx);

                cells.add(new OneDLongArray(out));
            
                long t = sumCounts(out);
                
                blockTotal += (t * t);
            }
        }
        
        blockTotal /= (double)cells.size();
        blockTotal = Math.sqrt(blockTotal);
        double norm = 1./(blockTotal + eps);
                
        float maxBlock = 255.f;
            //(N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM) *
            //(N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM) * 255,f;
   
        norm *= maxBlock;
        
        Arrays.fill(outHist, 0, outHist.length, 0);
        
        for (int i = 0; i < cells.size(); ++i) {
            long[] a = cells.get(i).a;
            for (int j = 0; j < a.length; ++j) {
                //v /= Math.sqrt(blockTotal + 0.0001);
                a[j] = (int)Math.round(norm * a[j]);
            } 
            add(outHist, a);
        }
        
        /*        
        part of a block of 3 X 3 cells
        
           2        2        2        2
           1        1        1        1
           0        0        0        0
          -9 -8 -7 -6 -5 -4 -3 -2 -1  *  1  2  3  4  5  6  7  9
                                      *
        */        
    }
    
    /**
     * NOT READY FOR USE
     * 
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * 
     * The normalization is the number of pixels visited during
     * construction.   The result is better for uses needing a signal level.
     * TODO: consider adding errors for this.
     * 
     * @param x
     * @param y
     * @param outHist 
     */
    public void extractBlock2(int x, int y, int[] outHist) {
                
        if (outHist.length != nBins) {
            throw new IllegalArgumentException("outHist.length != nBins");
        }

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        int nH = N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM;
                
        List<OneDIntArray> cells = new ArrayList<OneDIntArray>(nH);
        
        for (int cX = 0; cX < N_CELLS_PER_BLOCK_DIM; ++cX) {
            
            int cXOff = -(N_CELLS_PER_BLOCK_DIM/2) + cX;
        
            int x2 = x + (cXOff * N_PIX_PER_CELL_DIM);
            
            if ((x2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                break;
            } else if (x2 < 0) {
                x2 = 0;
            } else if (x2 >= w) {
                break;
            }
            
            for (int cY = 0; cY < N_CELLS_PER_BLOCK_DIM; ++cY) {
                    
                int cYOff = -(N_CELLS_PER_BLOCK_DIM/2) + cY;
                
                int y2 = y + (cYOff * N_PIX_PER_CELL_DIM);

                if ((y2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                    break;
                } else if (y2 < 0) {
                    y2 = 0;
                } else if (y2 >= h) {
                    break;
                }
                                
                int pixIdx = (y2 * w) + x2;
                
                int[] out = Arrays.copyOf(gHists[pixIdx], gHists[pixIdx].length);

                cells.add(new OneDIntArray(out));                            
            }
        }
        
        //TODO: might need to consider a factor to place the significant
        //    digits into integers
        
        // max value in a cell's histogram will be
        //   number of pixels * 255
        //   = N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM * 255
        //
        // the normalized counts would be divided by 
        //     N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM * 255
        //     but the result is smaller than 1 and we want to use integers,
        //     so can multiply the result by 255.
        
        double norm = 1./(double)cells.size();
            //* N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM);
               
        Arrays.fill(outHist, 0, outHist.length, 0);
        
        for (int i = 0; i < cells.size(); ++i) {
            int[] a = cells.get(i).a;
            add(outHist, a);
        }
        
        for (int i = 0; i < outHist.length; ++i) {
            outHist[i] = (int)Math.round(outHist[i] * norm);
        }
    }
    
    /**
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     * The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     * 
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     * 
     * @param histA
     * @param histB
     * @return 
     */
    public float intersection(int[] histA, int[] histB) {
       
        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }
        
        int nBins = histA.length;
        
        int binWidth = 256/nBins;
        
        /*
        histograms are already normalized
        
        K(a,b) = 
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */
            
        float sum = 0;
        float sumA = 0;
        float sumB = 0;
        for (int j = 0; j < nBins; ++j) {
            
            float yA = histA[j];
            float yB = histB[j];
            
            sum += Math.min(yA, yB);
            sumA += yA;
            sumB += yB;
            
            //System.out.println(" " + yA + " -- " + yB + " sum="+sum + ", " + sumA + "," + sumB);
        }
        
        float d = eps + Math.min(sumA, sumB);
        
        float sim = sum/d;
        
        return sim;
    }
    
    public float[] diff(int[] histA, int[] histB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        double sumDiff = 0;
        double err = 0;
                        
        for (int j = 0; j < nBins; ++j) {
            
            float yA = histA[j];
            float yB = histB[j];
            
            float maxValue = Math.max(yA, yB) + eps;

            float diff = Math.abs((yA - yB)/maxValue);
            
            //sumDiff += (diff * diff);
            sumDiff += diff;

            //      already squared
            err += (diff/maxValue);           
        }
        
        sumDiff /= (double)nBins;

        //sumDiff = Math.sqrt(sumDiff);

        err /= (double)nBins;
        err = Math.sqrt(err);

        return new float[]{(float)sumDiff, (float)err};
    }

    public int[] extractFeature(int xCenter, int yCenter, int detectorWidth,
        int detectorHeight) {
        
        int hw = detectorWidth/2;
        int hh = detectorHeight/2;

        if ((xCenter - hw) < 0 || (yCenter - hh) < 0 
            || (xCenter + hw) >= w || (yCenter + hh) >= h) {
            throw new IllegalArgumentException("out of bounds of "
                + "original image");
        }
        
        int hc = N_PIX_PER_CELL_DIM/2;
        
        /*        
                          xc,yc            
             |         |         |         |
        */
        int nX0 = (hw - hc)/N_PIX_PER_CELL_DIM;
        int startX = xCenter - (nX0 * N_PIX_PER_CELL_DIM);
        if (startX < hc) {
            nX0 = (xCenter - hc)/N_PIX_PER_CELL_DIM;
            startX = xCenter - (nX0 * N_PIX_PER_CELL_DIM);
        }
        int nX1 = (hw - hc)/N_PIX_PER_CELL_DIM;
        int stopX = xCenter + (nX1 * N_PIX_PER_CELL_DIM);
        if (stopX >= (this.w - hc)) {
            nX1 = (w - 1 - xCenter - hc)/N_PIX_PER_CELL_DIM;
            stopX = xCenter + (nX1 * N_PIX_PER_CELL_DIM);
        }
        int nY0 = (hh - hc)/N_PIX_PER_CELL_DIM;
        int startY = yCenter - (nY0 * N_PIX_PER_CELL_DIM);
        if (startY < hc) {
            nY0 = (yCenter - hc)/N_PIX_PER_CELL_DIM;
            startY = yCenter - (nY0 * N_PIX_PER_CELL_DIM);
        }
        int nY1 = (hh - hc)/N_PIX_PER_CELL_DIM;
        int stopY = yCenter + (nY1 * N_PIX_PER_CELL_DIM);
        if (stopY >= (this.h - hc)) {
            nY1 = (h - 1 - yCenter - hc)/N_PIX_PER_CELL_DIM;
            stopY = yCenter + (nY1 * N_PIX_PER_CELL_DIM);
        }
        
        //System.out.println(" startX=" + startX + " stopX=" + stopX
        //    + " startY=" + startY + " stopY=" + stopY
        //    + " HC=" + hc
        //);
        
        int nH = (nX0 + nX1 + 1) * (nY0 + nY1 + 1) * nBins;
        
        int[] tmp = new int[nBins];
        int[] out = new int[nH];
        
        int count = 0;
        double blockTotal = 0;
                
        // scan forward by 1 cell
        for (int x = startX; x <= stopX; x += N_PIX_PER_CELL_DIM) {
            for (int y = startY; y <= stopY; y += N_PIX_PER_CELL_DIM) {
                
                extractBlock(x, y, tmp);
                
                System.arraycopy(tmp, 0, out, count * nBins, nBins);
                
                double t = sumCounts(tmp);
                blockTotal += (t * t);               
                count++;                
            }
        }
        
        //System.out.println("NH=" + nH + " count=" + count + " blockTotal=" + blockTotal);
        
        // normalize over detector
        if (count > 0) {
            blockTotal = Math.sqrt(blockTotal/(double)count);
        }
        
        double norm = 1./(blockTotal + eps);

        float maxBlock = 255.f;
            //(N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM) *
            //(N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM);

        norm *= maxBlock;
        
        assert(!Double.isNaN(norm));

        for (int i = 0; i < out.length; ++i) {
            out[i] *= norm;
            assert(out[i] >= 0);
        }

        return out;
    }
    
    public float intersectionOfFeatures(int[] featureA, int[] featureB) {

        if ((featureA.length != featureB.length)) {
            throw new IllegalArgumentException(
                "featureA and featureB must be same dimensions");
        }
        
        int[] tmpA = new int[nBins];
        int[] tmpB = new int[nBins];
        
        float t;
        double sum = 0;
        for (int j = 0; j < featureA.length; j += nBins) {
            System.arraycopy(featureA, j, tmpA, 0, nBins);
            System.arraycopy(featureB, j, tmpB, 0, nBins);
            t = intersection(tmpA, tmpB);
            //System.out.println("    inter=" + t);
            sum += (t * t);
        }

        sum /= (double)(featureA.length/nBins);
        sum = Math.sqrt(sum);

        return (float)sum;
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * calculate the difference of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     *
     * Note that because the feature contains spatially ordered concatenation of
     * histograms, the registration of featureA and featureB to the same 
     * orientation must be done before this method (more specifically, before
     * extraction to features).
     *      *
     * @param featureA
     * @param featureB
     * @return
     */
    public float[] diffOfFeatures(int[] featureA, int[] featureB) {

        if ((featureA.length != featureB.length)) {
            throw new IllegalArgumentException(
                "featureA and featureB must be same dimensions");
        }
        
        int[] tmpA = new int[nBins];
        int[] tmpB = new int[nBins];
        
        float[] t;
        double sum = 0;
        double sumSqErr = 0;
        for (int j = 0; j < featureA.length; j += nBins) {
            System.arraycopy(featureA, j, tmpA, 0, nBins);
            System.arraycopy(featureB, j, tmpB, 0, nBins);
            t = diff(tmpA, tmpB);
            //System.out.println("    inter=" + t);
            sum += t[0];
            sumSqErr += (t[1] * t[1]);
        }

        sum /= (double)(featureA.length/nBins);
        //sum = Math.sqrt(sum);
        
        //TODO: check normalization by nBins here
        sumSqErr /= (double)(featureA.length/nBins);
        sumSqErr = Math.sqrt(sumSqErr);

        return new float[]{(float)sum, (float)sumSqErr};
    }
    
    private int sumCounts(int[] hist) {
        
        int sum = 0;
        for (int v : hist) {
            sum += v;
        }
        
        return sum;
    }
     
    private long sumCounts(long[] hist) {
        
        long sum = 0;
        for (long v : hist) {
            sum += v;
        }
        
        return sum;
    }
    
    private void add(int[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    private void add(long[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    private void add(long[] addTo, long[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    private long[] copyHist(int pixIdx) {
        int[] a = gHists[pixIdx];
        long[] out = new long[a.length];
        for (int i = 0; i < a.length; ++i) {
            out[i] = a[i];
        }
        return out;
    }
    
    public int getNumberOfBins() {
        return nBins;
    }

    public int getImageWidth() {
        return w;
    }

    public int getImageHeight() {
        return h;
    }

}
