package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.OneDIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 NOT READY FOR USE.
 
 An implementation of Histograms of Oriented Gradients
 constructed from reading the following papers:
 <pre>
 "Histograms of Oriented Gradients for Human Detection" 
 by Dalal and Triggs, 2010 
 and
 "Distinctive Image Features from Scale-Invariant Keypoints"
 by Lowe, 2004
 </pre>

    The histograms of oriented gradients are constructed to compare
    the 2 part data of gradient angle and gradient magnitude
    between regions within different images as a feature descriptor.
    The feature is defined over several pixels surrounding the central
    keypoint in a pattern to improve the measurement.

   The GradientIntegralHistogram is used to store histograms of values that
   that are counts defined by gradient magnitudes in bins of
   orientation.
   
   As recommended by Dalal & Triggs, 9 bins are used for 180 degrees of gradient 
   angle range.
   
   The building of the integral image has a runtime complexity of O(N_pixels) 
   for the gradient and O(N_pixels) for the histogram integral image.
   
   Extraction of histogram data is 4 steps, just as in summed area tables, but 
   there is additionally the copy which is nBins steps.
  
   The extraction of data is at the "cell" level, which is recommended to be 
   6 X 6 pixels^2 by Dalal and Triggs.
   
   a block of cells is gathered for a point and that is an addition of the 
   N_cells X N_cells histograms.
   
   Before the addition, block level normalization is calculated for each cell.
   
   The block level normalization uses L2NormHys.
   They found best results using normalization of each cell histogram
   by the total over the block.
   The total number of values is summed over all cells within the
   block and a normalization factor for each cell is then computed using 
   that block total.  The individually normalized cells are then added
   over the block to create the block histogram.
       for each cell calculate cell_total_count.
       total_block_count = sum over cells ( cell_total_count )
       for each cell, normalization is 1/total_block_count
   then the block histogram is added over the same bins in each cell.
   To keep the block histogram as integer but normalized to same 
   max value could apply a further factor of 
   max possible value for a block being, 
   for example (2X2)*(6X6)*(255) = 36720.

   Note that a shift is needed for identifying the bin that is the
   canonical angle 0, that is a shift specific to a 
   dominant angle correction for the histogram.
   That shift will be applied during the addition stage to produce a
   canonicalized descriptor in reference frame w.r.t. rotation correction.
   (Note that for the use case here, the reference frame orientation will
   be supplied to the method. it's learned from the mser ellipse in one
   use case for example. so the application of a dominant orientation
   will be the same, but the calculation will not be performed for it.
   though that could be added as a method later...not necessary for
   current requirements).

   Comparison of block descriptors is then a histogram intersection,
   where normally 0 is no intersection, hence maximally different,
   and an intersection equal to the max value is maximally similar.
   (see the method ColorHistogram.intersection, but here, the normalization
   will already have been applied instead of determined in the method).
  
  Other details are in converting the intersection to a cost or score
  and specialized methods for that specific to this project will be
  present in this class.
  
  @author nichole
*/
public class HOGs {
   
    // 9 is default
    private final int nAngleBins;
    
    // 6 x 6 is recommended
    private final int N_PIX_PER_CELL_DIM;
    
    // 2x2 or 3x3 is recommended
    private final int N_CELLS_PER_BLOCK_DIM;
    
    private final GradientIntegralHistograms gHists;
    
    private final int w;
    private final int h;
    
    public HOGs(GreyscaleImage rgb) {
        
        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = 6;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = rgb.getWidth();
        h = rgb.getHeight();
        
        gHists = init(rgb);
    }
    
    public HOGs(GreyscaleImage gradientXY, GreyscaleImage theta) {
        
        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = 6;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = gradientXY.getWidth();
        h = gradientXY.getHeight();
        
        gHists = init(gradientXY, theta);
    }
    
    private GradientIntegralHistograms init(GreyscaleImage rgb) {
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private GradientIntegralHistograms init(GreyscaleImage gradientXY, 
        GreyscaleImage theta) {
        
        if (w != gradientXY.getWidth() || h != gradientXY.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }
        
        if (w != theta.getWidth() || h != theta.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms(gradientXY,
            theta, nAngleBins);

        return gh;        
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * 
     * @param x
     * @param y
     * @param outHist 
     */
    public void extractFeature(int x, int y, int[] outHist) {

        if (outHist.length != nAngleBins) {
            throw new IllegalArgumentException("outHist.length != nAngleBins");
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
        
        int cXOff = -(N_CELLS_PER_BLOCK_DIM/2);
        int[] outN = new int[1];
        
        List<OneDIntArray> cells = new ArrayList<OneDIntArray>(nH);
        
        for (int cX = 0; cX < N_CELLS_PER_BLOCK_DIM; ++cX) {
            
            int cYOff = -(N_CELLS_PER_BLOCK_DIM/2);
        
            int x2 = x + (cXOff * N_PIX_PER_CELL_DIM);
            
            if ((x2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                break;
            } else if (x2 < 0) {
                x2 = 0;
            } else if (x2 >= w) {
                break;
            }
            int x3 = x2 + N_PIX_PER_CELL_DIM - 1;
            if (x3 > (w - 1)) {
               x3 = w - 1;
            }
            
            for (int cY = 0; cY < N_CELLS_PER_BLOCK_DIM; ++cY) {
                    
                int y2 = y + (cYOff * N_PIX_PER_CELL_DIM);

                if ((y2 + N_PIX_PER_CELL_DIM - 1) < 0) {
                    break;
                } else if (y2 < 0) {
                    y2 = 0;
                } else if (y2 >= h) {
                    break;
                }
                int y3 = y2 + N_PIX_PER_CELL_DIM - 1;
                if (y3 > (h - 1)) {
                    y3 = h - 1;
                }

                int[] out = new int[nAngleBins];

                gHists.extractWindow(x2, x3, y2, y3, out, outN);

                cells.add(new OneDIntArray(out));
            
                int t = sumCounts(out);
                
                blockTotal += (t * t);
            }
            cYOff++;
        }
        
        double norm = 1./Math.sqrt(blockTotal + 0.0001);
                
        float maxBlock = (N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM) *
            (N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM) * 255,f;
   
        norm *= maxBlock;
        
        //NOTE: this may need to be revised:
        // one more term to the normalization needed in case the number of
        //  cells is fewer than nH (due to the point being near the image edge).
        //  will scale the results by the ratio.
        norm *= (double)nH/(double)cells.size();
        
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
        
        throw new UnsupportedOperationException("not yet implemented");
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
}
