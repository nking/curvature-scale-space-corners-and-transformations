package algorithms.imageProcessing.segmentation;

import algorithms.QuickSort;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.Gaussian1DFirstDeriv;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a variant of kmeans whose goal is to make k super-pixels in the image
 * based upon greyscale similarity and x,y proximity.
 *
 * The code implements the algorithm of:
 * "SLIC Superpixels Compared to State-of-the-Art Superpixel Methods"
   by Achanta, Appu Shaji,Smith,  Lucchi, Fua, and Su ̈sstrunk,
 * but replaces the cielab colorspace with greyscale color space.
 * 
 * @author nichole
 */
public class SLICSuperPixelsGS {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected int[] labels = null;

    protected double[] distances = null;

    protected final int k;

    protected final int s;

    protected final int nXs;

    protected final int nYs;

    protected final float[][] seedDescriptors;

    protected final GreyscaleImage img;

    protected double[][] gradient = null;

    protected final double threshold;

    protected final double maxBrightness;

    /**
     * constructor for the super-pixel algorithm SLIC that uses the
     * color space of greyscale and creates approximately nClusters (a.k.a.
     * super-pixels).
     * The number of
     * clusters may be adjusted for even intervals over width and height.
     * @param img
     * @param nClusters
     * @param clrNorm a smaller number than 10 finds details such as textures
     * and a number closer to 40 produces very blocky segmentation mostly
     * lacking curves. (NOTE: these numbers need to be updated for GS, but
     * the directions of change are same).
     */
    public SLICSuperPixelsGS(GreyscaleImage img, int nClusters, double clrNorm) {

        if (clrNorm < 1 || (clrNorm > 40)) {
            throw new IllegalArgumentException("clrNorm should be in range "
                + "1 to 40");
        }

        maxBrightness = clrNorm;

        double sampling = Math.sqrt(( (float)img.getNPixels()/(float)nClusters));

        if (sampling < 1) {
            sampling = 1;
        }

        this.s = (int)Math.round(sampling);

        nXs = Math.round((float)img.getWidth()/(float)s);
        nYs = Math.round((float)img.getHeight()/(float)s);
        this.k = nXs * nYs;

        log.info("k = " + k + " s=" + this.s + " nXs=" + nXs + " nYs=" + nYs);

        this.img = img;

        // v, x, y
        seedDescriptors = new float[k][];
        for (int i = 0; i < k; ++i) {
            seedDescriptors[i] = new float[3];
        }

        // max error would be ( maxClr * maxClr * k) + 2*( s/2 * s/2 * k)
        double maxError = 2*(maxBrightness * maxBrightness * k);
        maxError = Math.sqrt(maxError);
        threshold = 0.01 * maxError;
    }

    /**
     * an optional method to set the gradient, else is
     * calculated internal to the class and discarded.
     * @param gradientImg
     */
    public void setGradient(GreyscaleImage gradientImg) {

        if (gradientImg.getNPixels() != img.getNPixels()
            || gradientImg.getWidth() != img.getWidth()
            || gradientImg.getHeight() != img.getHeight()) {
            throw new IllegalArgumentException(
                "gradientImg must be same size as img");
        }

        int width = gradientImg.getWidth();
        int height = gradientImg.getHeight();

        gradient = new double[width][];
        for (int i = 0; i < width; ++i) {
            gradient[i] = new double[height];
            for (int j = 0; j < height; ++j) {
                gradient[i][j] = gradientImg.getValue(i, j);
            }
        }
    }

    /**
     * constructor for the super-pixel algorithm SLIC that uses the default
     * color space of greyscale and creates approximately nClusters (a.k.a.
     * super-pixels).
     * The number of
     * clusters may be adjusted for even intervals over width and height.
     * Uses a default clrNorm = 10;
     * @param img
     * @param nClusters
     * 
     */
    public SLICSuperPixelsGS(GreyscaleImage img, int nClusters) {

        if  (nClusters > img.getNPixels()) {
            throw new IllegalArgumentException(
                "nClusters must be smaller than number of pixels in img");
        }

        maxBrightness = 10;

        //TOOD: after have an implementation as authors suggest,
        //  change to use deltaE instead of sqrt sum diffs of CIE lab
        //  and compare differences in results and runtime (many more flops...)

        double sampling = Math.sqrt(( (float)img.getNPixels()/(float)nClusters));

        if (sampling < 1) {
            sampling = 1;
        }

        this.s = (int)Math.round(sampling);

        nXs = Math.round((float)img.getWidth()/(float)s);
        nYs = Math.round((float)img.getHeight()/(float)s);
        this.k = nXs * nYs;

        log.info("k = " + k);

        this.img = img;

        // v, x, y
        seedDescriptors = new float[k][];
        for (int i = 0; i < k; ++i) {
            seedDescriptors[i] = new float[3];
        }

        // max error would be ( maxClr * maxClr * k) + 2*( s/2 * s/2 * k)
        double maxError = 2*(maxBrightness * maxBrightness * k);
        maxError = Math.sqrt(maxError);
        threshold = 0.01 * maxError;
    }

     public void calculate() {

        init();

        int nIterMax = 20;

        int nIter = 0;

        while (nIter < nIterMax) {

            assignPixelsNearSeeds();

            double l2Norm = adjustClusterCenters();

            log.fine("l2Norm=" + l2Norm + " nIter=" + nIter);

            if (l2Norm < threshold) {
                break;
            }

            nIter++;
        }

        assignTheUnassigned();

        //ImageSegmentation imageSegmentation = new ImageSegmentation();
        //imageSegmentation.replaceSinglePixelLabelsCIELAB(labels, img);

    }

    protected void init() {

        if (labels != null) {
            throw new IllegalStateException("variables have been initialized");
        }

        int nPix = img.getNPixels();

        labels = new int[nPix];
        Arrays.fill(labels, -1);

        distances = new double[nPix];
        Arrays.fill(distances, Double.MAX_VALUE);

        // init cluster centers, seedDescriptors for grid with cell size s
        populateSeedDescriptors();
    }

    private void populateSeedDescriptors() {

        if (gradient == null) {
            gradient = calcGradient();
        }

        /*
        sampled on a regular grid spaced S pixels apart.
        To produce roughly equally sized superpixels,
        the grid interval is S.

        ** The centers are moved to seed locations corresponding to the lowest
        gradient position in a 3 X 3 neighborhood. This is done to avoid
        centering a superpixel on an edge and to reduce the chance of
        seeding a superpixel with a noisy pixel.
        */

        int w = img.getWidth();
        int h = img.getHeight();

        // determine the centers of each s x s cell within search range of 3x3 in gradient

        final int sHalf = s/2;
        int dx, dy;
        if (s < 3) {
            dx = 0;
            dy = 0;
        } else {
            dx = 1;
            dy = 1;
        }

        // kCurrent = (iNy * nX) + iNx;
        for (int iNy = 0; iNy < nYs; ++iNy) {

            int y1 = sHalf + iNy*s;
            if ((dy == 0) && y1 > (h - 1)) {
                y1 = h - 1;
            }

            for (int iNx = 0; iNx < nXs; ++iNx) {

                int x1 = sHalf + iNx*s;
                if ((dx == 0) && x1 > (w - 1)) {
                    x1 = w - 1;
                }

                int kCurrent = (iNy * nXs) + iNx;

                // find smallest gradient within range (x1, y1) +-1 to set center for seed
                double minG = Double.MAX_VALUE;
                int minX2 = -1;
                int minY2 = -1;
                for (int x2 = (x1 - dx); x2 <= (x1 + dx); ++x2) {
                    if (x2 < 0 || (x2 > (w - 1))) {
                        continue;
                    }
                    for (int y2 = (y1 - dy); y2 <= (y1 + dy); ++y2) {
                        if (y2 < 0 || (y2 > (h - 1))) {
                            continue;
                        }
                        double g = gradient[x2][y2];
                        if (g < minG) {
                            minG = g;
                            minX2 = x2;
                            minY2 = y2;
                        }
                    }
                }

                // [v, x, y]
                seedDescriptors[kCurrent][0] = img.getValue(minX2, minY2);
                seedDescriptors[kCurrent][1] = minX2;
                seedDescriptors[kCurrent][2] = minY2;

                int pixIdx2 = img.getInternalIndex(minX2, minY2);
                labels[pixIdx2] = kCurrent;
                distances[pixIdx2] = 0;

                log.fine("seed " + kCurrent + " x=" + minX2 + " y=" + minY2);
            }
        }
    }

   /**
    * O(N)
    *
    * @return
    */
   private double[][] calcGradient() {

       int nPix = img.getNPixels();
       int width = img.getWidth();
       int height = img.getHeight();

       double[][] gradient = new double[width][];
       for (int i = 0; i < width; ++i) {
           gradient[i] = new double[height];
       }

       float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

       int h = (kernel.length - 1) >> 1;

       // TODO: edit to operate on one direction, then operate on that result
       //        in other direction

       for (int i = 0; i < nPix; ++i) {
           final int x1 = img.getCol(i);
           final int y1 = img.getRow(i);

           float vXSum = 0;
           float vYSum = 0;

           for (int g = 0; g < kernel.length; ++g) {
               float gg = kernel[g];
               if (gg == 0) {
                   continue;
               }

               int x2, y2;
               // calc for X gradient first
               int delta = g - h;
                x2 = x1 + delta;
                y2 = y1;
                // edge corrections.  use replication
                if (x2 < 0) {
                    x2 = -1 * x2 - 1;
                } else if (x2 >= width) {
                    int diff = x2 - width;
                    x2 = width - diff - 1;
                }
                vXSum += gg * img.getValue(x2, y2);

                // calc for y
                y2 = y1 + delta;
                x2 = x1;
                // edge corrections.  use replication
                if (y2 < 0) {
                    y2 = -1 * y2 - 1;
                } else if (y2 >= height) {
                    int diff = y2 - height;
                    y2 = height - diff - 1;
                }
                vYSum += gg * img.getValue(x2, y2);
           }

           double vC = Math.sqrt(vXSum * vXSum + vYSum * vYSum);

           // presumably, greyscale gradient is fine

           gradient[x1][y1] = vC;
       }

        return gradient;
    }

    private void assignPixelsNearSeeds() {

        int w = img.getWidth();
        int h = img.getHeight();

        for (int kCurrent = 0; kCurrent < k; ++kCurrent) {

            float[] desc1 = seedDescriptors[kCurrent];
            int x1 = (int) desc1[1];
            int y1 = (int) desc1[2];

            for (int x2 = (x1 - s); x2 <= (x1 + s); ++x2) {
                for (int y2 = (y1 - s); y2 <= (y1 + s); ++y2) {

                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    int pixIdx2 = img.getInternalIndex(x2, y2);

                    double dist = calcDist(desc1,
                        img.getValue(x2, y2), x2, y2);

                    if (dist < distances[pixIdx2]) {
                        distances[pixIdx2] = dist;
                        labels[pixIdx2] = kCurrent;
                    }
                }
            }
        }

    }

    private double adjustClusterCenters() {

        // L2 norm is the residuals added in quadratur

        double[][] meanDescriptors = new double[k][];
        // v, x, y
        for (int i = 0; i < k; ++i) {
            meanDescriptors[i] = new double[3];
        }

        // calculate new colors
        int[] count = new int[k];

        for (int i = 0; i < img.getNPixels(); ++i) {
            int label = labels[i];
            if (label == -1) {
                // this is still unassigned
                continue;
            }
            int x = img.getCol(i);
            int y = img.getRow(i);
            int v = img.getValue(i);

            meanDescriptors[label][0] += v;
            meanDescriptors[label][1] += x;
            meanDescriptors[label][2] += y;
            count[label]++;
        }

        // v, x, y
        double[] sqDiffSum = new double[3];

        double diff;
        for (int kCurrent = 0; kCurrent < k; ++kCurrent) {
            assert(count[kCurrent] > 0);
            float nc = count[kCurrent];
            for (int m = 0; m < 3; ++m) {
                meanDescriptors[kCurrent][m] /= nc;
                diff = meanDescriptors[kCurrent][m] - seedDescriptors[kCurrent][m];
                sqDiffSum[m] += (diff * diff);
            }
        }

        double l2Norm = 0;
        for (double sd : sqDiffSum) {
            l2Norm += sd;
        }
        l2Norm = Math.sqrt(l2Norm);

        // calc sqrt of sum of sq diffs with old centers and reset old centers
        for (int kCurrent = 0; kCurrent < k; ++kCurrent) {
            for (int m = 0; m < 3; ++m) {
                seedDescriptors[kCurrent][m] = (float)meanDescriptors[kCurrent][m];
            }
        }

        return l2Norm;
    }
    
    private double calcDist(float[] desc1, int v, int x2, int y2) {

        double dClrSq = 0;
        float diff = desc1[0] - v;
        dClrSq += (diff * diff);

        float diffX = desc1[1] - x2;
        float diffY = desc1[2] - y2;
        double dXYSq = diffX * diffX + diffY * diffY;

        double dComb = dClrSq + (dXYSq * maxBrightness * maxBrightness)/((float)s * s);

        dComb = Math.sqrt(dComb);

        return dComb;
    }

    private void assignTheUnassigned() {

        Map<PairInt, TIntSet> unassignedMap = new HashMap<PairInt, TIntSet>();

        int w = img.getWidth();
        int h = img.getHeight();
        int n = img.getNPixels();

        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                if (labels[pixIdx] == -1) {
                    addNeighobLabelsForPoint(unassignedMap, i, j, dxs, dys);
                }
            }
        }

        ArrayDeque<PairInt> queue0 = populateByNumberOfNeighbors(unassignedMap);

        ArrayDeque<PairInt> queue1 = new ArrayDeque<PairInt>();

        int nIter = 0;

        Set<PairInt> visited = new HashSet<PairInt>();

        while (!queue0.isEmpty() || !queue1.isEmpty()) {

            PairInt p;
            if (!queue1.isEmpty()) {
                p = queue1.poll();
            } else {
                p = queue0.poll();
            }

            if (visited.contains(p)) {
                continue;
            }
            visited.add(p);

            int x1 = p.getX();
            int y1 = p.getY();

            TIntSet adjLabels;
            if (nIter == 0) {
                adjLabels = unassignedMap.get(p);
                assert (adjLabels != null);
            } else {
                adjLabels = new TIntHashSet();
                addNeighobLabelsForPoint(adjLabels, x1, y1, dxs, dys);
            }

            double minD = Double.MAX_VALUE;
            int minLabel2 = -1;

            TIntIterator iter = adjLabels.iterator();
            while (iter.hasNext()) {
                int label2 = iter.next();
                double dist = calcDist(seedDescriptors[label2],
                    img.getValue(x1, y1), x1, y1);

                if (dist < minD) {
                    minD = dist;
                    minLabel2 = label2;
                }
            }

            int pixIdx1 = img.getInternalIndex(p.getX(), p.getY());
            labels[pixIdx1] = minLabel2;
            distances[pixIdx1] = minD;

            unassignedMap.remove(p);

            for (int m = 0; m < dxs.length; ++m) {
                int x2 = p.getX() + dxs[m];
                int y2 = p.getY() + dys[m];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int pixIdx2 = img.getInternalIndex(x2, y2);
                if (labels[pixIdx2] == -1) {
                    PairInt p2 = new PairInt(x2, y2);
                    queue1.add(p2);
                    //assert (!visited.contains(p2));
                }
            }
            nIter++;
        }

        assert(unassignedMap.isEmpty());
    }

    private ArrayDeque<PairInt> populateByNumberOfNeighbors(
        Map<PairInt, TIntSet> unassignedMap) {

        int n = unassignedMap.size();

        PairInt[] points = new PairInt[n];
        int[] nN = new int[n];

        int count = 0;
        for (Entry<PairInt, TIntSet> entry : unassignedMap.entrySet()) {
            points[count] = entry.getKey();
            nN[count] = entry.getValue().size();
            count++;
        }

        QuickSort.sortBy1stArg(nN, points);

        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();

        for (int i = (n - 1); i > -1; --i) {

            int nP = nN[i];
            if (nP == 0) {
                break;
            }

            queue.add(points[i]);
        }

        return queue;
    }

    public int[] getLabels() {
        return labels;
    }

    private void addNeighobLabelsForPoint(Map<PairInt, TIntSet> unassignedMap,
        int i, int j, int[] dxs, int[] dys) {

        int w = img.getWidth();
        int h = img.getHeight();

        PairInt p = new PairInt(i, j);

        TIntSet adjLabels = unassignedMap.get(p);
        if (adjLabels == null) {
            adjLabels = new TIntHashSet();
            unassignedMap.put(p, adjLabels);
        }

        addNeighobLabelsForPoint(adjLabels, i, j, dxs, dys);
    }

    private void addNeighobLabelsForPoint(TIntSet adjLabels,
        int i, int j, int[] dxs, int[] dys) {

        int w = img.getWidth();
        int h = img.getHeight();

        PairInt p = new PairInt(i, j);

        for (int m = 0; m < dxs.length; ++m) {
            int x2 = i + dxs[m];
            int y2 = j + dys[m];
            if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                continue;
            }
            int pixIdx2 = img.getInternalIndex(x2, y2);
            if (labels[pixIdx2] > -1) {
                adjLabels.add(labels[pixIdx2]);
            }
        }
    }
}
