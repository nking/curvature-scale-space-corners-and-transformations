package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.OneDLongArray;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import algorithms.util.TwoDIntArray;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

/**
 * This is a version of HOGs which requires a CRegion list as input
 * and creates masked integral histograms for each.
 * It does not cover the entire given image, only the CRegions.
 * 

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
public class HOGRegions {

    private static float eps = 0.000001f;
    
    // 9 is default
    private final int nAngleBins;

    // 6 x 6 is recommended
    private final int N_PIX_PER_CELL_DIM;

    // 2x2 or 3x3 is recommended
    private final int N_CELLS_PER_BLOCK_DIM;

    /**
     * key = region index
     */
    private final TIntObjectMap<TwoDIntArray> regionIndexHist = 
        new TIntObjectHashMap<TwoDIntArray>();
    private final TIntObjectMap<OneDIntArray> regionIndexMinMaxXY =
        new TIntObjectHashMap<OneDIntArray>();
    private final TIntObjectMap<CRegion> regionIndexRegions;
    //In reference frame of subImage
    private final TIntObjectMap<TLongSet> regionCoords 
        = new TIntObjectHashMap<TLongSet>();

    private final int w;
    private final int h;

    private boolean debug = false;

    //TODO: calculate the limits in nPixels this can handle due to
    //   using integers instead of long for storage.
    // e.g. 8.4 million pix, roughly 2900 X 2900

    public HOGRegions(GreyscaleImage rgb, TIntObjectMap<CRegion> regionMap) {

        regionIndexRegions = regionMap;
        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = 4;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = rgb.getWidth();
        h = rgb.getHeight();
        init(rgb);
    }

    public HOGRegions(GreyscaleImage rgb, int nCellsPerDim, int nPixPerCellDim, 
        TIntObjectMap<CRegion> regionMap) {

        regionIndexRegions = regionMap;
        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = rgb.getWidth();
        h = rgb.getHeight();
        init(rgb);
    }
    
    public HOGRegions(GreyscaleImage rgb, int nCellsPerDim, int nPixPerCellDim,
        int nAngleBins, TIntObjectMap<CRegion> regionMap) {

        regionIndexRegions = regionMap;

        this.nAngleBins = nAngleBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = rgb.getWidth();
        h = rgb.getHeight();
        init(rgb);
    }

    public HOGRegions(GreyscaleImage gradientXY, GreyscaleImage theta, 
        TIntObjectMap<CRegion> regionMap) {

        regionIndexRegions = regionMap;
        nAngleBins = 9;
        N_PIX_PER_CELL_DIM = 4;
        N_CELLS_PER_BLOCK_DIM = 2;
        w = gradientXY.getWidth();
        h = gradientXY.getHeight();
        init(gradientXY, theta);
    }
    
    public HOGRegions(TIntObjectMap<CRegion> regionMap, int imageWidth,
        int imageHeight, int nCellsPerDim, int nPixPerCellDim, int nAngleBins) {

        regionIndexRegions = regionMap;
        this.nAngleBins = nAngleBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = imageWidth;
        h = imageHeight;
    }

    public void setToDebug() {
        debug = true;
    }

    private void init(GreyscaleImage rgb) {

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

        init(gXY, theta);
    }

    private void init(GreyscaleImage gradientXY, GreyscaleImage theta) {

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
                        
        TIntObjectIterator<CRegion> iter = regionIndexRegions.iterator();
        for (int i = 0; i < regionIndexRegions.size(); ++i) {
            iter.advance();
            
            int rIndex = iter.key();
            CRegion r = iter.value();
            
            int[] minMaxXY = new int[4];
        
            //In reference frame of subImage
            TLongSet rCoords = new TLongHashSet();
                
            int[][] histograms = createHistogram(gradientXY, theta, 
                r.extractCoords(), minMaxXY, rCoords);
            
            regionCoords.put(rIndex, rCoords);
            regionIndexHist.put(rIndex, new TwoDIntArray(histograms));
            regionIndexMinMaxXY.put(rIndex, new OneDIntArray(minMaxXY));
        }
    }
    
    int[][] createHistogram(GreyscaleImage gradientXY, 
        GreyscaleImage theta, Collection<PairInt> points,
        int[] outputMinMaxXY, TLongSet outputRefFramePixs) {

        if (w != gradientXY.getWidth() || h != gradientXY.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }

        if (w != theta.getWidth() || h != theta.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }

        PixelHelper ph = new PixelHelper();
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms();
        
        outputMinMaxXY[0] = Integer.MAX_VALUE;
        outputMinMaxXY[1] = Integer.MIN_VALUE;
        outputMinMaxXY[2] = Integer.MAX_VALUE;
        outputMinMaxXY[3] = Integer.MIN_VALUE;
        for (PairInt xy : points) {
            if (xy.getX() < outputMinMaxXY[0]) {
                outputMinMaxXY[0] = xy.getX();
            }
            if (xy.getX() > outputMinMaxXY[1]) {
                outputMinMaxXY[1] = xy.getX();
            }
            if (xy.getY() < outputMinMaxXY[2]) {
                outputMinMaxXY[2] = xy.getY();
            }
            if (xy.getY() > outputMinMaxXY[3]) {
                outputMinMaxXY[3] = xy.getY();
            }
        } 
        
        GreyscaleImage gXY = gradientXY.subImage2(outputMinMaxXY[0], 
            outputMinMaxXY[1], outputMinMaxXY[2], outputMinMaxXY[3]);

        GreyscaleImage th = theta.subImage2(outputMinMaxXY[0], 
            outputMinMaxXY[1], outputMinMaxXY[2], outputMinMaxXY[3]);

        int w2 = outputMinMaxXY[1] - outputMinMaxXY[0] + 1;
        int h2 = outputMinMaxXY[3] - outputMinMaxXY[2] + 1;
        int xOffset = outputMinMaxXY[0];
        int yOffset = outputMinMaxXY[2];

        //In reference frame of subImage
        for (PairInt xy : points) {
            long pixIdx = ph.toPixelIndex(xy.getX() - xOffset, 
                xy.getY() - yOffset, w2);
            if (pixIdx < 0) {
                int z = 0;
            } else if (pixIdx >= (w2 * h2)) {
                int z = 0;
            }
            outputRefFramePixs.add(pixIdx);
        }
        
        // mask out pixels not in the region
        for (int i2 = 0; i2 < w2; ++i2) {
            for (int j2 = 0; j2 < h2; ++j2) {
                long pixIdx = ph.toPixelIndex(i2, j2, w2);
                if (!outputRefFramePixs.contains(pixIdx)) {
                    gXY.setValue(i2, j2, 0);
                    th.setValue(i2, j2, 0);
                }
            }
        }

        int[][] histograms = gh.createHistograms(gXY, th, nAngleBins);

        //apply a windowed sum across the integral image.
        // result is that at each pixel is a histogram holding the sum of histograms
        //    from the surrounding N_PIX_PER_CELL_DIM window. 
        gh.applyWindowedSum(histograms, gXY.getWidth(), gXY.getHeight(), 
            N_PIX_PER_CELL_DIM);

        return histograms;
    }
    
    public void addARegion(GreyscaleImage gradientXY, GreyscaleImage theta,
        Canonicalizer.RegionPoints regionPoints) {

        if (w != gradientXY.getWidth() || h != gradientXY.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }

        if (w != theta.getWidth() || h != theta.getHeight()) {
            throw new IllegalArgumentException("gradient and theta must be same size");
        }
    
        int[] minMaxXY = new int[4];
        
        //In reference frame of subImage
        TLongSet rCoords = new TLongHashSet();
                
        int[][] histograms = createHistogram(gradientXY, theta, 
            regionPoints.points, minMaxXY, rCoords);
        
        int w2 = minMaxXY[1] - minMaxXY[0] + 1;
        int h2 = minMaxXY[3] - minMaxXY[2] + 1;
        //int xOffset = minMaxXY[0];
        //int yOffset = minMaxXY[2];
         
        Canonicalizer canonicalizer = new Canonicalizer();
        
        Set<CRegion> cRegions = canonicalizer.canonicalizeRegions4(
            regionPoints, w2, h2);
        
        for (CRegion cRegion : cRegions) {
        
            int rIndex = this.regionIndexRegions.size();
            cRegion.dataIdx = rIndex;
        
            regionCoords.put(rIndex, rCoords);
            regionIndexHist.put(rIndex, new TwoDIntArray(histograms));
            regionIndexMinMaxXY.put(rIndex, new OneDIntArray(minMaxXY));
            regionIndexRegions.put(rIndex, cRegion);
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
     * @param rIndex
     * @param x
     * @param y
     * @param outHist
     */
    public void extractBlock(int rIndex, int x, int y, int[] outHist) {

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

        int[][] hist = regionIndexHist.get(rIndex).a;
        int[] minMaxXY = regionIndexMinMaxXY.get(rIndex).a;
        TLongSet rCoords = regionCoords.get(rIndex);
        
        int width = minMaxXY[1] - minMaxXY[0] + 1;
        int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        List<OneDIntArray> cells = new ArrayList<OneDIntArray>(nH);

        for (int cX = 0; cX < N_CELLS_PER_BLOCK_DIM; ++cX) {

            int cXOff = -(N_CELLS_PER_BLOCK_DIM/2) + cX;

            int x2 = x + (cXOff * N_PIX_PER_CELL_DIM);

            int xSub = x2 - xOffset;
            
            if ((xSub + N_PIX_PER_CELL_DIM - 1) < 0) {
                break;
            } else if (xSub < 0) {
                xSub = 0;
            } else if (xSub >= width) {
                break;
            }

            for (int cY = 0; cY < N_CELLS_PER_BLOCK_DIM; ++cY) {

                int cYOff = -(N_CELLS_PER_BLOCK_DIM/2) + cY;

                int y2 = y + (cYOff * N_PIX_PER_CELL_DIM);

                int ySub = y2 - yOffset;
                
                if ((ySub + N_PIX_PER_CELL_DIM - 1) < 0) {
                    break;
                } else if (ySub < 0) {
                    ySub = 0;
                } else if (ySub >= height) {
                    break;
                }
                
                //if (!rCoords.contains(ySub * width + xSub)) {
                //    continue;
                //}

                int pixIdx = (ySub * width) + xSub;

                int[] out = Arrays.copyOf(hist[pixIdx], hist[pixIdx].length);

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

        return HOGUtil.diff(histA, orientationA, histB, orientationB);
    }
    
    private int sumCounts(int[] hist) {

        int sum = 0;
        for (int v : hist) {
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

    /**
     * NOT READY FOR USE.
     *
     * @param rIndex
     * @return
     */
    public int calculateDominantOrientation(int rIndex) {

        int[][] hist = regionIndexHist.get(rIndex).a;
        int[] minMaxXY = regionIndexMinMaxXY.get(rIndex).a;
        
        int width = minMaxXY[1] - minMaxXY[0] + 1;
        //int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        CRegion r = regionIndexRegions.get(rIndex);
        
        long[] combined = new long[nAngleBins];

        for (Entry<PairInt, PairInt> entry : r.offsetsToOrigCoords.entrySet()) {
            PairInt xy = entry.getValue();
            int pixIdx = ((xy.getY() - yOffset) * width) + (xy.getX() - xOffset);
            add(combined, hist[pixIdx]);
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

    public TIntSet calculateDominantOrientations(int rIndex) {

        int[][] hist = regionIndexHist.get(rIndex).a;
        int[] minMaxXY = regionIndexMinMaxXY.get(rIndex).a;
        
        int width = minMaxXY[1] - minMaxXY[0] + 1;
        //int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        CRegion r = regionIndexRegions.get(rIndex);
        
        Set<PairInt> xys = new HashSet<PairInt>();
        
        long[] combined = new long[nAngleBins];

        for (Entry<PairInt, PairInt> entry : r.offsetsToOrigCoords.entrySet()) {
            PairInt xy = entry.getValue();
            int pixIdx = ((xy.getY() - yOffset) * width) + (xy.getX() - xOffset);
            add(combined, hist[pixIdx]);
            xys.add(xy);
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
        PairInt xyCen = ch.calculateXYCentroids2(xys);
        if (xys.contains(xyCen)) {
            int pixIdx = ((xyCen.getY() - yOffset) * width) + (xyCen.getX() - xOffset);
            maxIdx = MiscMath.findYMaxIndex(hist[pixIdx]);
            if (maxIdx > -1) {
                int angle = Math.round((maxIdx + 0.5f) * binWidth);
                orientations.add(angle);
            }
        }
        
        return orientations;
    }

    public int getImageWidth() {
        return w;
    }
    
    public int getImageHeight() {
        return h;
    }

    public int getNumberOfBins() {
        return nAngleBins;
    }

}
