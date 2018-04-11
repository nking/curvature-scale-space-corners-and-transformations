package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 CAVEAT: small amount of testing done, not yet throughly tested.

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
   That shift will be applied during the intersection stage to produce a
   canonicalized feature in a rotation corrected reference frame.
   (Note that for the use case here, the reference frame orientation will
   be supplied to the method. it's learned from the mser ellipse in one
   use case for example. so the application of a dominant orientation
   will be the same, but the calculation will not be performed for it.
   though that could be added as a method later...not necessary for
   current requirements).

   Comparison of block feature is then a histogram intersection,
   where normally 0 is no intersection, hence maximally different,
   and an intersection equal to the max value is maximally similar.
   (see the method ColorHistogram.intersection, but here, the normalization
   will already have been applied instead of determined in the method).

  Other details are in converting the intersection to a cost or score
  and specialized methods for that specific to this project will be
  present in this class.
   
  TODO: need to add an improved feature comparison method that considers the
  empty bins significant too.  the error calculations in Histogram.java include
  that

  @author nichole
*/
public class HOGs {

    private static float eps = 0.000001f;
    
    // 9 is default
    private final int nAngleBins;

    // 6 x 6 is recommended
    private final int N_PIX_PER_CELL_DIM;

    // 2x2 or 3x3 is recommended
    private final int N_CELLS_PER_BLOCK_DIM;

    // histogrm integral images with a windowed sum of N_PIX_PER_CELL_DIM
    private final int[][] gHists;

    private final int w;
    private final int h;

    private boolean debug = true;

    //TODO: calculate the limits in nPixels this can handle due to
    //   using integers instead of long for storage.
    // e.g. 8.4 million pix, roughly 2900 X 2900

    public HOGs(GreyscaleImage rgb) {

        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = 4;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = rgb.getWidth();
        h = rgb.getHeight();

        gHists = init(rgb);
    }

    public HOGs(GreyscaleImage rgb, int nCellsPerDim, int nPixPerCellDim) {

        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = rgb.getWidth();
        h = rgb.getHeight();

        gHists = init(rgb);
    }
    
    public HOGs(GreyscaleImage rgb, int nCellsPerDim, int nPixPerCellDim,
            int nAngleBins) {

        this.nAngleBins = nAngleBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = rgb.getWidth();
        h = rgb.getHeight();

        gHists = init(rgb);
    }

    public HOGs(GreyscaleImage gradientXY, GreyscaleImage theta) {

        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = 4;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = gradientXY.getWidth();
        h = gradientXY.getHeight();

        gHists = init(gradientXY, theta);
    }

    public void setToDebug() {
        debug = true;
    }

    private int[][] init(GreyscaleImage rgb) {

        ImageProcessor imageProcessor = new ImageProcessor();

        // instead of sobel, using 1st deriv
        GreyscaleImage[] gXgY = imageProcessor.createCentralDifferenceGradients(rgb);

        GreyscaleImage theta = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);

        GreyscaleImage gXY = 
            imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);

        if (debug) {
            algorithms.misc.MiscDebug.writeImage(gXgY[0], "_gX_");
            algorithms.misc.MiscDebug.writeImage(gXgY[1], "_gY_");
            algorithms.misc.MiscDebug.writeImage(gXY, "_gXY_");
            algorithms.misc.MiscDebug.writeImage(theta, "_theta_");
        }

        return init(gXY, theta);
    }

    private int[][] init(GreyscaleImage gradientXY, GreyscaleImage theta) {

        if (w != gradientXY.getWidth() || h != gradientXY.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }

        if (w != theta.getWidth() || h != theta.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }

        if (debug) {
            algorithms.misc.MiscDebug.writeImage(gradientXY, "_gXY_");
            algorithms.misc.MiscDebug.writeImage(theta, "_theta_");
        }

        GradientIntegralHistograms gh = new GradientIntegralHistograms();

        int[][] histograms = gh.createHistograms(gradientXY, theta, nAngleBins);

        //apply a windowed sum across the integral image.
        // result is that at each pixel is a histogram holding the sum of histograms
        //    from the surrounding N_PIX_PER_CELL_DIM window. 
        gh.applyWindowedSum(histograms, w, h, N_PIX_PER_CELL_DIM);

        return histograms;
    }

    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in 
     * block were set during
     * construction.
     * 
     * The feature is nAngleBins in length for 180 degrees
     * and the bin with the largest value
     * is the bin holding the angle perpendicular to the windowed point.
     * (for example: a horizontal line, the feature of a point on the
     * line has largest bin being the 90 degrees bin).
     *
     * @param x
     * @param y
     * @param outHist
     */
    public void extractBlock(int x, int y, int[] outHist) {

        if (outHist.length != nAngleBins) {
            throw new IllegalArgumentException("outHist.length != nAngleBins");
        }

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }

        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normalize each cell by that sum.

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

        if (!cells.isEmpty()) {
            blockTotal /= (double)cells.size();
            blockTotal = Math.sqrt(blockTotal);
        }
        
        double norm = 1./(blockTotal + eps);

        float maxBlock = 255.f * cells.size();
            //(N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM) *
            //(N_PIX_PER_CELL_DIM * N_PIX_PER_CELL_DIM);

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
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * Extract the blocks within the feature, add and normalize them and place
     * the result in the output array.
     * 
     * The feature is nAngleBins in length for 180 degrees
     * and the bin with the largest value
     * is the bin holding the angle perpendicular to the windowed point.
     * (for example: a horizontal line, the feature of a point on the
     * line has largest bin being the 90 degrees bin).
     *
     * Each block is extracted with "extractBlock" which is normalized by its 
     * own sum, then the sum of the blocks is used here to normalize the
     * feature.
     * 
     * @param xCenter
     * @param yCenter
     * @param detectorWidth
     * @param detectorHeight
     * @return outHist
     */
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
        
        int nH = (nX0 + nX1 + 1) * (nY0 + nY1 + 1) * nAngleBins;
        
        int[] tmp = new int[nAngleBins];
        int[] out = new int[nH];
        
        int count = 0;
        double blockTotal = 0;
                
        // scan forward by 1 cell
        for (int x = startX; x <= stopX; x += N_PIX_PER_CELL_DIM) {
            for (int y = startY; y <= stopY; y += N_PIX_PER_CELL_DIM) {
                
                extractBlock(x, y, tmp);
                
                System.arraycopy(tmp, 0, out, count * nAngleBins, nAngleBins);
                
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

        float maxBlock = 255.f * count;
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
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
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
     * Note also that you may want to try the rotation of oppossite direction.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return
     */
    public float intersection(int[] histA, int orientationA, int[] histB,
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

            sum += Math.min(yA, yB);
            sumA += yA;
            sumB += yB;

            //System.out.println(" " + yA + " -- " + yB + " sum="+sum + ", " + sumA + "," + sumB);
        }

        float d = eps + Math.min(sumA, sumB);
        float sim = sum/d;

        return sim;
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return difference, error
     */
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

        /*
        histograms are already normalized

        K(a,b) =
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */

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
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * calculate the intersection of histA and histB which have already
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
    public float intersectionOfFeatures(int[] featureA, int[] featureB) {

        if ((featureA.length != featureB.length)) {
            throw new IllegalArgumentException(
                "featureA and featureB must be same dimensions");
        }
        
        int[] tmpA = new int[nAngleBins];
        int[] tmpB = new int[nAngleBins];
        
        float t;
        double sum = 0;
        for (int j = 0; j < featureA.length; j += nAngleBins) {
            System.arraycopy(featureA, j, tmpA, 0, nAngleBins);
            System.arraycopy(featureB, j, tmpB, 0, nAngleBins);
            t = intersection(tmpA, 0, tmpB, 0);
            //System.out.println("    inter=" + t);
            sum += (t * t);
        }

        sum /= (double)(featureA.length/nAngleBins);
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
        
        int[] tmpA = new int[nAngleBins];
        int[] tmpB = new int[nAngleBins];
        
        float[] t;
        double sum = 0;
        double sumSqErr = 0;
        for (int j = 0; j < featureA.length; j += nAngleBins) {
            System.arraycopy(featureA, j, tmpA, 0, nAngleBins);
            System.arraycopy(featureB, j, tmpB, 0, nAngleBins);
            t = diff(tmpA, 0, tmpB, 0);
            //System.out.println("    inter=" + t);
            sum += t[0];
            sumSqErr += (t[1] * t[1]);
        }

        sum /= (double)(featureA.length/nAngleBins);
        //sum = Math.sqrt(sum);
        
        sumSqErr /= (double)(featureA.length/nAngleBins);
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
        return nAngleBins;
    }

    /**
     * NOT READY FOR USE.
     *
     * @param xy
     * @return
     */
    public int calculateDominantOrientation(Collection<PairInt> xy) {

        long[] combined = new long[nAngleBins];

        for (PairInt p : xy) {
            int pixIdx = (p.getY() * w) + p.getX();
            add(combined, gHists[pixIdx]);
        }

        TIntList maxIdxs = new TIntArrayList();
        long maxValue = Long.MIN_VALUE;
        for (int i = 0; i < combined.length; ++i) {
            long v = combined[i];
            if (v > maxValue) {
                maxIdxs.clear();
                maxIdxs.add(i);
                maxValue = v;
            } else if (v == maxValue) {
                maxIdxs.add(i);
            }
        }

        int binWidth = 180 / nAngleBins;

        //TODO: revise to fit max peak using neighboring bins
        if (maxIdxs.size() == 1) {
            return Math.round((maxIdxs.get(0) + 0.5f) * binWidth);
        }

        //TODO: revise this
        if (maxIdxs.size() == 2) {

            // min of average with and without a wraparound phase
            float ang0 = (maxIdxs.get(0) + 0.5f) * binWidth;
            float ang1 = (maxIdxs.get(1) + 0.5f) * binWidth;

            if (ang0 < ang1) {
                float diff0 = 360 + ang0 - ang1;
                float diff1 = ang1 - ang0;
                if (diff0 < diff1) {
                    return Math.round(0.5f * (360 + ang0 + ang1));
                }

                return Math.round(0.5f * (ang0 + ang1));
            }

            float diff0 = 360 + ang1 - ang0;
            float diff1 = ang0 - ang1;
            if (diff0 < diff1) {
                return Math.round(0.5f * (360 + ang1 + ang0));
            }

            return Math.round(0.5f * (ang0 + ang1));

        } else {

            // average of all indexes

            double[] angles = new double[maxIdxs.size()];
            for (int i = 0; i < maxIdxs.size(); ++i) {
                angles[i] = ((maxIdxs.get(i) + 0.5f) * binWidth);
            }

            double angleAvg =
                AngleUtil.calculateAverageWithQuadrantCorrections(
                    angles, false);

            return (int)Math.round(angleAvg);
        }
    }

    public TIntSet calculateDominantOrientations(Collection<PairInt> xy) {

        long[] combined = new long[nAngleBins];

        for (PairInt p : xy) {
            int pixIdx = (p.getY() * w) + p.getX();
            add(combined, gHists[pixIdx]);
        }

        int maxIdx = MiscMath.findYMaxIndex(combined);

        if (maxIdx == -1) {
            throw new IllegalArgumentException("histogram is full of "
                + " min value longs");
        }

        // if any bins have values within 80% of max, add to maxIdxs
        TIntList maxIdxs = new TIntArrayList();

        long max = combined[maxIdx];
        double limit = 0.8 * max;

        for (int i = 0; i < combined.length; ++i) {
            long v = combined[i];
            if (v >= limit) {
                maxIdxs.add(i);
            }
        }

        int binWidth = 180 / nAngleBins;

        TIntSet orientations = new TIntHashSet();

        for (int i = 0; i < maxIdxs.size(); ++i) {
            int idx = maxIdxs.get(i);
            int angle = Math.round((idx + 0.5f) * binWidth);
            orientations.add(angle);
        }

        return orientations;
    }

    public void _printHistograms() {
        for (int i = 0; i < gHists.length; ++i) {
            int[] gh0 = gHists[i];
            if (hasNonZeroes(gh0)) {
                System.out.format("%d) %s\n", i, Arrays.toString(gh0));
            }
        }
    }
    private boolean hasNonZeroes(int[] a) {
        for (int b : a) {
            if (b != 0) {
                return true;
            }
        }
        return false;
    }
    
    int[] extractHistogram(int x, int y) {
        int pixIdx = (y * w) + x;
        int[] out = Arrays.copyOf(gHists[pixIdx], gHists[pixIdx].length);
        return out;
    }
    
}
