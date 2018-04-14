package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.IntegralHistograms;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
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
                
        float maxBlock = (N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM) *
            (N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM) * 255,f;
   
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
    
    public float[] diff(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        if (orientationA < 0 || orientationA > 180 || orientationB < 0 ||
            orientationB > 180) {
            throw new IllegalArgumentException("orientations must be in range 0 to 180,"
                + "  inclusive,  or!=" + orientationA + " orB=" + orientationB);
        }
        if (orientationA == 180) {
            orientationA = 0;
        }
        if (orientationB == 180) {
            orientationB = 0;
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        int shiftA = (orientationA - 90)/binWidth;
        int shiftB = (orientationB - 90)/binWidth;
        
        double sumDiff = 0;
        double err = 0;
                        
        for (int j = 0; j < nBins; ++j) {
            
            int idxA = j + shiftA;
            if (idxA < 0) {
                idxA += nBins;
            } else if (idxA > (nBins - 1 )) {
                idxA -= nBins;
            }

            int idxB = j + shiftB;
            if (idxB < 0) {
                idxB += nBins;
            } else if (idxB > (nBins - 1 )) {
                idxB -= nBins;
            }

            float yA = histA[idxA];
            float yB = histB[idxB];
            
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

    private int sumCounts(int[] hist) {
        
        int sum = 0;
        for (int v : hist) {
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
    
    public int getNumberOfBins() {
        return nBins;
    }

}
