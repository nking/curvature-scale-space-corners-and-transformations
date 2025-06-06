package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.util.AngleUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.Collection;

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
   
  TODO: improve distance function used in feature comparison method.  need to consider the
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

    private boolean debug = false;

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

    /**
     * create a HOG structure for patch comparisons. 
     * @param gradientXY
     * @param theta
     */
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

        // note, there may be negative values in the 2D integral histograms, 
        //    but when extracted to 1D for their single pixels or windows, 
        //    the extracted histogram should be all positive
        gh.applyWindowedSum(histograms, w, h, N_PIX_PER_CELL_DIM);
        
        //_printHistograms_xy(histograms);
        
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

        Arrays.fill(outHist, 0, outHist.length, 0);
                
        int r = N_CELLS_PER_BLOCK_DIM >> 1;
        int stopY = y + r;
        int stopX = x + r;
        int startX = x - r;
        int startY = y - r;
        if ((h & 1) == 0) {
            startX--;
            startY--;            
        }
        if (startX < 0) {
            startX = 0;
        }
        if (startY > 0) {
            startY = 0;
        }
        if (stopX >= w) {
            stopX = w - 1;
        }
        if (stopY >= h) {
            stopY = h - 1;
        }
        
        int[] outputN = new int[1];
        
        HOGUtil.extractWindow(gHists, startX, stopX, startY, stopY, w, h, 
            outHist, outputN);
        
        assert(isAllPositive(outHist));
            
        double blockTotal = sumCounts(outHist);
        blockTotal *= blockTotal;

        double norm;
        if (blockTotal > 0) {
            blockTotal /= (double)outputN[0];
            blockTotal = Math.sqrt(blockTotal);
            norm = 255./blockTotal;
        } else {
            norm = 255.;
        }
        
        for (int i = 0; i < outHist.length; ++i) {
            outHist[i] = (int)Math.round(norm * outHist[i]);
        }
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
    public void extractBlock(int x, int y, long[] outHist) {

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

        Arrays.fill(outHist, 0, outHist.length, 0);
                
        int r = N_CELLS_PER_BLOCK_DIM >> 1;
        int stopY = y + r;
        int stopX = x + r;
        int startX = x - r;
        int startY = y - r;
        if ((h & 1) == 0) {
            startX--;
            startY--;            
        }
        if (startX < 0) {
            startX = 0;
        }
        if (startY > 0) {
            startY = 0;
        }
        if (stopX >= w) {
            stopX = w - 1;
        }
        if (stopY >= h) {
            stopY = h - 1;
        }
        
        int[] outputN = new int[1];
        
        HOGUtil.extractWindow(gHists, startX, stopX, startY, stopY, w, h, 
            outHist, outputN);
            
        double blockTotal = sumCounts(outHist);
        blockTotal *= blockTotal;

        double norm;
        if (blockTotal > 0) {
            blockTotal /= (double)outputN[0];
            blockTotal = Math.sqrt(blockTotal);
            norm = 255./blockTotal;
        } else {
            norm = 255.;
        }

        for (int i = 0; i < outHist.length; ++i) {
            outHist[i] = (int)Math.round(norm * outHist[i]);
        }
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
        
        //assert(isAllPositive(gHists));
        
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
                
                //assert(isAllPositive(tmp));
                
                System.arraycopy(tmp, 0, out, count * nAngleBins, nAngleBins);
                
                double t = sumCounts(tmp);
                blockTotal += (t * t);               
                count++;                
            }
        }
        
        assert(blockTotal >= 0.);
        
        //System.out.println("NH=" + nH + " count=" + count + " blockTotal=" + blockTotal);
        
        // normalize over detector
        double norm;
        if (count > 0) {
            blockTotal = Math.sqrt(blockTotal/(double)count);
            norm = 255./blockTotal;
        } else {
            norm = 255.;
        }
        
        assert(!Double.isNaN(norm));

        for (int i = 0; i < out.length; ++i) {
            out[i] *= (int)Math.round(norm);
            //assert(out[i] >= 0);
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

        return HOGUtil.intersection(histA, orientationA, histB, orientationB);
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     * 
     * calculate the difference of histA and histB (normalized by the maximum value
       in their histograms).
     * A result of 0 is maximally similar and a result of 1 is maximally disssimilar.
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

        return HOGUtil.diff(histA, orientationA, histB, orientationB);
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
     * calculate the difference of featureA and featureB (normalized by the maximum value
       in their vectors).
     * A result of 0 is maximally similar and a result of 1 is maximally dissimilar.
     *
     * Note that because the feature contains spatially ordered concatenation of
     * histograms, the registration of featureA and featureB to the same 
     * orientation must be done before this method (more specifically, before
     * extraction to features).
     *
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
    
    public static long sumCounts(long[] hist) {

        long sum = 0;
        for (long v : hist) {
            sum += v;
        }

        return sum;
    }

    public static void add(int[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }

    public static void add(long[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    public static void add(long[] addTo, long[] addFrom) {
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
        
        //NOTE: adding an orientation for the center of points
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        PairInt xyCen = ch.calculateXYCentroids2(xy);
        if (xy.contains(xyCen)) {
            int pixIdx = (xyCen.getY() * w) + xyCen.getX();
            maxIdx = MiscMath.findYMaxIndex(gHists[pixIdx]);
            if (maxIdx > -1) {
                int angle = Math.round((maxIdx + 0.5f) * binWidth);
                orientations.add(angle);
            }
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
    public void _printHistograms_xy(int[][] a) {
        PixelHelper ph = new PixelHelper();
        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {
                int pixIdx = (int)ph.toPixelIndex(col, row, w);
                int[] gh0 = a[pixIdx];
                if (hasNonZeroes(gh0)) {
                    System.out.format("%d,%d) %s\n", col, row, Arrays.toString(gh0));
                }
            }
        }
    }
    public void _printHistograms_xy() {
        _printHistograms_xy(gHists);
    }
    private boolean hasNonZeroes(int[] a) {
        for (int b : a) {
            if (b != 0) {
                return true;
            }
        }
        return false;
    }

    public int getImageWidth() {
        return w;
    }
    public int getImageHeight() {
        return h;
    }

    private boolean isAllPositive(int[] a) {
        for (int i = 0; i < a.length; ++i) {
            if (a[i] < 0) {
                System.out.println("  " + Arrays.toString(a));
                return false;
            }
        }
        return true;
    }
    private boolean isAllPositive(int[][] h) {
        for (int i = 0; i < h.length; ++i) {
            for (int j = 0; j < h[i].length; ++j) {
                if (h[i][j] < 0) {
                    return false;
                }
            }
        }
        return true;
    }

}
