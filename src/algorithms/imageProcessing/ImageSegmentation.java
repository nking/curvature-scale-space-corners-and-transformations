package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.NearestPoints1D;
import algorithms.compGeometry.NearestPointsInLists;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusColor;
import algorithms.compGeometry.clustering.KMeansPlusPlusFloat;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.features.BlobMedialAxes;
import algorithms.imageProcessing.features.BlobsAndPerimeters;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.features.IntensityClrFeatures;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MedianSmooth;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.KDTree;
import algorithms.search.KDTreeNode;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class holding several different image segmentation methods.  Note that
 * some other techniques involving contrast for example, are elsewhere.
 *
 * @author nichole
 */
public class ImageSegmentation {

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * applies KMeansPlusPlus algorithm to the values in input
     * (greyscale intensities) to create kBands of clustered pixels
     * (operates on input).
     * @param input
     * @param kBands
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    public void applyUsingKMPP(GreyscaleImage input, int kBands)
        throws IOException, NoSuchAlgorithmException {

        KMeansPlusPlus instance = new KMeansPlusPlus();
        instance.computeMeans(kBands, input);
        
        int[] seeds = instance.getCenters();
        
        int[] imgSeedIndexAssignments = instance.getImgPixelSeedIndexes();

        for (int pixIdx = 0; pixIdx < input.getNPixels(); ++pixIdx) {

            int seedIdx = imgSeedIndexAssignments[pixIdx];

            int seedValue = seeds[seedIdx];

            input.setValue(pixIdx, seedValue);
        }
    }
    
    public void applyUsingKMPP(Image input, int kBands) throws IOException, 
        NoSuchAlgorithmException {
        
        KMeansPlusPlusColor instance = new KMeansPlusPlusColor();
        instance.computeMeans(kBands, input);
        
        int[][] seeds = instance.getCenters();
        
        int[] imgSeedIndexAssignments = instance.getImgPixelSeedIndexes();
        
        for (int pixIdx = 0; pixIdx < input.getNPixels(); ++pixIdx) {

            int seedIdx = imgSeedIndexAssignments[pixIdx];

            int r = seeds[0][seedIdx];
            int g = seeds[1][seedIdx];
            int b = seeds[2][seedIdx];
            
            input.setRGB(pixIdx, r, g, b);
        }
    }

    /**
     * applies binary algorithm (simple thresholding) to the values in input
     * (greyscale intensities) to create pixels given above or below highest
     * frequency value.  Note that some images may need to be pre-processed
     * in order to use this one (for example, correct for illumination and
     * remove items like sky if main objects are not sky).
     * (operates on input).
     * @param input
     */
    public void applyBinaryUsingFrequency(GreyscaleImage input) {

        PairIntArray valueCounts = Histogram.createADescendingSortbyFrequencyArray(
            input);

        if (valueCounts == null || valueCounts.getN() == 0) {
            return;
        }

        //96, 30
        int divider = valueCounts.getX(0);
        int v0 = 255/4;
        int v1 = 3*v0;

        for (int i = 0; i < input.getNPixels(); ++i) {
            if (input.getValue(i) < divider) {
                input.setValue(i, v0);
            } else {
                input.setValue(i, v1);
            }
        }

    }

    /**
     * applies DTClustering algorithm to the values in input
     * (greyscale intensities) to create kBands of clustered pixels
     * (operates on input).
     * (This one is competitive with applyUsingPolarCIEXYAndFrequency
     * with lowFreqLimit 0.1f)
     * @param input
     * @param kBands
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    public void applyUsingDTClustering(GreyscaleImage input, int kBands)
        throws IOException, NoSuchAlgorithmException {

        PairIntArray valueCounts = Histogram.createADescendingSortbyFrequencyArray(
            input);

        // first, trying clustering by value and frequency,
        // second, compare that to clustering just by giving points as value,value
        //    which is effectivly a 1D clustering
        //Set<PairInt> points = new HashSet<PairInt>();
        Set<com.climbwithyourfeet.clustering.util.PairInt> points =
            new HashSet<com.climbwithyourfeet.clustering.util.PairInt>();

        for (int i = 0; i < valueCounts.getN(); ++i) {
            com.climbwithyourfeet.clustering.util.PairInt p = new
                com.climbwithyourfeet.clustering.util.PairInt(
                    valueCounts.getX(i), valueCounts.getY(i));

            points.add(p);
        }

        DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt> cFinder
            = new DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt>(
                points, input.getWidth(), input.getHeight());

        cFinder.calculateCriticalDensity();

        cFinder.findClusters();

        int n = cFinder.getNumberOfClusters();

        int[] centers = new int[n];
        for (int i = 0; i < n; ++i) {

            Set<com.climbwithyourfeet.clustering.util.PairInt> set = cFinder.getCluster(i);

            // find centeroid for x
            double xc = 0;
            for (com.climbwithyourfeet.clustering.util.PairInt p : set) {
                double x1 = p.getX();
                xc += x1;
            }

            centers[i] = (int)Math.round(xc/(double)set.size());
        }

        if (n > kBands) {
            n = kBands;
        }

        Arrays.sort(centers);
        int[] kCenters = new int[n];
        int count = 0;
        for (int i = (centers.length - 1); i > (centers.length - 1 - n); --i) {
            kCenters[count] = centers[i];
            count++;
        }

        assignToNearestCluster(input, kCenters);
    }

    /**
     * places points by their proximity to cluster centers
     * @param input
     * @param binCenters
     */
    public void assignToNearestCluster(GreyscaleImage input, int[] binCenters) {

        NearestPoints1D np = new NearestPoints1D(binCenters);

        for (int col = 0; col < input.getWidth(); col++) {

            for (int row = 0; row < input.getHeight(); row++) {

                int v = input.getValue(col, row);

                int vc = np.findClosestValue(v);

                input.setValue(col, row, vc);
            }
        }
    }
                
    private KDTreeNode findNearestNeighborBruteForce(int[] xs, int[] ys, 
        int x, int y) {
        
        int minDistSq = Integer.MAX_VALUE;
        int minDistIdx = -1;
        
        for (int i = 0; i < xs.length; ++i) {
            int diffX = xs[i] - x;
            int diffY = ys[i] - y;
            int distSq = diffX * diffX + diffY * diffY;
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistIdx = i;
            }
        }
        
        KDTreeNode node = new KDTreeNode();
        node.setX(xs[minDistIdx]);
        node.setY(ys[minDistIdx]);
        
        return node;
    }

    public List<Set<PairInt>> assignToNearestPolarCIECluster(
        Map<PairInt, Integer> polarCIEXYMap, int[] binCenters) {

        Arrays.sort(binCenters);

        int nc = binCenters.length;

        List<Set<PairInt>> groups = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < nc; ++i) {
            groups.add(new HashSet<PairInt>());
        }

        for (Entry<PairInt, Integer> entry : polarCIEXYMap.entrySet()) {

            int theta = entry.getValue().intValue();

            int idx = Arrays.binarySearch(binCenters, theta);
            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }
            if (idx > (nc - 1)) {
                idx = nc - 1;
            }

            int vc = binCenters[idx];

            if (idx == 0) {

                int bisectorBelowHalfLength = (360 - binCenters[nc - 1] + vc)/2;
                if ((vc - theta) > bisectorBelowHalfLength) {
                    idx = nc - 1;
                }

            } else {

                int bisectorBelow = ((binCenters[idx - 1] + vc) / 2);

                if (theta < bisectorBelow) {
                    idx = idx - 1;
                }
            }

            Set<PairInt> set = groups.get(Integer.valueOf(idx));
            set.add(entry.getKey());
        }

        for (int i = (groups.size() - 1); i > -1; --i) {
            if (groups.get(i).isEmpty()) {
                groups.remove(i);
            }
        }

        return groups;
    }

    /**
     * applies a blur of sigma=1 to image,
     * converts each pixel color to the polar angle of CIE XY Lab color space
     * with an origin of (0.35, 0.35) and uses a histogram binning of kColors=8,
     * then maps those bins to 0 to 255,
     * then replaces a pixel if 5,6 or its neighbors have the same color,
     * then applies histogram equalization to stretch the values to range
     * 0 to 255.
     * @param input
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistEq(ImageExt input) {

        int kColors = 8;

        return applyUsingCIEXYPolarThetaThenHistEq(input, kColors, true);
    }

    /**
     * converts each pixel color to the polar angle of CIE XY Lab color space
     * with an origin of (0.35, 0.35) and uses a histogram binning of kColors,
     * then maps those bins to 0 to 255,
     * then replaces a pixel if 5,6 or its neighbors have the same color,
     * then applies histogram equalization to stretch the values to range
     * 0 to 255.
     * runtime complexity is about O(N) + O(N * lg_2(N)) though the later
     * term has small constant multiples of it.
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @param useBlur if true, a blur of sigma=1 is applied to the image before
     * processing.
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistEq(ImageExt input,
        int kColors, boolean useBlur) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img;

        int minNeighborLimit;

        if (useBlur) {

            Image input2 = input.copyImage();

            //O(N) + O(N * lg_2(N)) where N=n_pixels and k=kernel.length
            imageProcessor.blur(input2, 1/*(float)Math.sqrt(2)/2.f*/);

            //O(N) + O(N * lg_2(N))
            img = applyUsingCIEXYPolarThetaThenHistogram(input2, kColors);

            minNeighborLimit = 6;

        } else {

            //O(N) + O(N * lg_2(N))
            img = applyUsingCIEXYPolarThetaThenHistogram(input, kColors);

            minNeighborLimit = 5;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // ----replace pixel, if 5,6 or more neighbors have same color -----
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        Map<Integer, Integer> freqMap = new HashMap<Integer, Integer>();

        int nChanged = 1;
        int nIterMax = 100;
        int nIter = 0;

        // change value of pixels largely surrounded by another value

        // runtime complexity is nIterations * O(N*8)
        // number of iterations appears to be small multiple of minNeighborLimit
        while (!useBlur && (nIter < nIterMax) && (nChanged > 0)) {

            log.fine("***nIter=" + nIter + " nChanged=" + nChanged +
                " minNeighbotLimit=" + minNeighborLimit);

            nChanged = 0;

            //O(N*8)
            for (int col = 0; col < w; col++) {
                for (int row = 0; row < h; row++) {

                    freqMap.clear();

                    Integer maxCountValue = null;
                    int maxCount = Integer.MIN_VALUE;

                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        int x = dxs[nIdx] + col;
                        int y = dys[nIdx] + row;

                        if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                            break;
                        }

                        Integer v = Integer.valueOf(img.getValue(x, y));

                        Integer c = freqMap.get(v);
                        if (c == null) {
                            c = Integer.valueOf(1);
                        } else {
                            c = Integer.valueOf(c.intValue() + 1);
                        }
                        freqMap.put(v, c);

                        if (c.intValue() > maxCount) {
                            maxCount = c.intValue();
                            maxCountValue = v;
                        }
                    }

                    if ((maxCount >= minNeighborLimit) &&
                        (img.getValue(col, row) != maxCountValue.intValue())) {

                        img.setValue(col, row, maxCountValue.intValue());
                        nChanged++;
                    }
                }
            }

            nIter++;
        }

        // rescale the image
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();

        return img;
    }

    /**
     * applies a blur of sigma=1 to image,
     * converts each pixel's color into CIE XY polar theta, then uses KMeansPlusPlus
     * to create kColors=8 bins points then remaps the points to use the
     * range 0 to 255, then replaces a pixel if it has 7 neighbors of same
     * color, then applies histogram equalization to rescale the range to be
     * between 0 and 255.
     * @param input
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenKMPPThenHistEq(ImageExt
        input) {

        int kColors = 8;

        return applyUsingCIEXYPolarThetaThenKMPPThenHistEq(input, kColors, true);
    }

    /**
     * converts each pixel's color into CIE XY polar theta, then uses KMeansPlusPlus
     * to create kColors bins points then remaps the points to use the
     * range 0 to 255, then replaces a pixel if it has 5,6 neighbors of same
     * color, then applies histogram equalization to rescale the range to be
     * between 0 and 255.
     *
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @param useBlur
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenKMPPThenHistEq(ImageExt input,
        int kColors, boolean useBlur) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img;

        int minNeighborLimit;

        if (useBlur) {

            Image input2 = input.copyImage();

            imageProcessor.blur(input2, 1/*(float)Math.sqrt(2)/2.f*/);

            img = applyUsingCIEXYPolarThetaThenKMPP(input2, kColors);

            minNeighborLimit = 6;

        } else {

            img = applyUsingCIEXYPolarThetaThenKMPP(input, kColors);

            minNeighborLimit = 5;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // ----replace pixel, if 5,6 or more neighbors have same color -----
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        Map<Integer, Integer> freqMap = new HashMap<Integer, Integer>();

        int nChanged = 1;
        int nIterMax = 100;
        int nIter = 0;

        while (!useBlur && (nIter < nIterMax) && (nChanged > 0)) {

            log.fine("nIter=" + nIter + " nChanged=" + nChanged);

            nChanged = 0;

            for (int col = 0; col < w; col++) {
                for (int row = 0; row < h; row++) {

                    freqMap.clear();

                    Integer maxCountValue = null;
                    int maxCount = Integer.MIN_VALUE;

                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        int x = dxs[nIdx] + col;
                        int y = dys[nIdx] + row;

                        if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                            break;
                        }

                        Integer v = Integer.valueOf(img.getValue(x, y));

                        Integer c = freqMap.get(v);
                        if (c == null) {
                            c = Integer.valueOf(1);
                        } else {
                            c = Integer.valueOf(c.intValue() + 1);
                        }
                        freqMap.put(v, c);

                        if (c.intValue() > maxCount) {
                            maxCount = c.intValue();
                            maxCountValue = v;
                        }
                    }

                    if ((maxCount >= minNeighborLimit) &&
                        (img.getValue(col, row) != maxCountValue.intValue())) {

                        img.setValue(col, row, maxCountValue.intValue());
                        nChanged++;
                    }
                }
            }

            nIter++;
        }

        // rescale the image
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();

        return img;
    }

    /**
     * converts each pixel's color into CIE XY polar theta,
     * then applies histogram mapping of kColors to remap the pixels to
     * values between 0 and 255.
     *
     * @param input
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistogram(Image input) {

        return applyUsingCIEXYPolarThetaThenHistogram(input, 253);
    }

    /**
     * converts each pixel's color into CIE XY polar theta,
     * then applies histogram mapping of kColors to remap the pixels to
     * values between 0 and 255.
     * runtime complexity is O(N) + O(N * lg_2(N)).
     * @param input
     * @param kColors the number of color bins to use for the image segmentation.
     * The minimum allowed value is 2 and the maximum allowed value is 253.
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistogram(Image input,
        int kColors) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        /*
        TODO: needs some improvements in color mapping.
           -- for regions that are very spotty, might consider using the
              intensity image to help define the region and take the highest
              number density color in that region and assign it to all.
        */

        int w = input.getWidth();
        int h = input.getHeight();

        float[] tmpColorBuffer = new float[2];

        GreyscaleImage output = new GreyscaleImage(w, h);

        Map<PairInt, Float> pixThetaMap = new HashMap<PairInt, Float>();

        CIEChromaticity cieC = new CIEChromaticity();

        float[] thetaValues = new float[input.getNPixels()];
        int thetaCount = 0;

        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {

                PairInt p = new PairInt(col, row);

                int r = input.getR(col, row);
                int g = input.getG(col, row);
                int b = input.getB(col, row);

                if ((r < 25) && (g < 25) && (b < 25)) {
                    continue;
                } else if ((r > 230) && (g > 230) && (b > 230)) {//might need to use 195 as lower limit
                    output.setValue(col, row, 255);
                    continue;
                }

                float[] cieXY = tmpColorBuffer;
                cieC.rgbToXYChromaticity(r, g, b, cieXY);

                if (cieC.isWhite(cieXY[0], cieXY[1])) {
                    output.setValue(col, row, 255);
                } else {

                    double thetaRadians = cieC.calculateXYTheta(cieXY[0], cieXY[1]);

                    thetaValues[thetaCount] = (float)thetaRadians;

                    pixThetaMap.put(p, Float.valueOf((float)thetaRadians));

                    thetaCount++;
                }
            }
        }

        thetaValues = Arrays.copyOf(thetaValues, thetaCount);

        // O(N * lg_2(N))
        createAndApplyHistMapping(output, pixThetaMap, thetaValues, kColors);

        return output;
    }

    /**
     * converts each pixel's color into CIE XY polar theta, then uses KMeansPlusPlus
     * to create kColors bins points then remaps the points to use the
     * range 0 to 255.
     *
     * runtime complexity is O(N) + O(N*lg_2(N))
     *
     * @param input
     * @param kColors the number of color bins to use for the image segmentation.
     * The minimum allowed value is 2 and the maximum allowed value is 253.
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenKMPP(Image input, int kColors) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        /*
        TODO: needs some improvements in color mapping.
           -- for regions that are very spotty, might consider using the
              intensity image to help define the region and take the highest
              number density color in that region and assign it to all.
        */

        int w = input.getWidth();
        int h = input.getHeight();

        float[] tmpColorBuffer = new float[2];

        GreyscaleImage output = new GreyscaleImage(w, h);

        Map<PairInt, Float> pixThetaMap = new HashMap<PairInt, Float>();

        CIEChromaticity cieC = new CIEChromaticity();

        float[] thetaValues = new float[input.getNPixels()];
        int thetaCount = 0;

        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {

                PairInt p = new PairInt(col, row);

                int r = input.getR(col, row);
                int g = input.getG(col, row);
                int b = input.getB(col, row);

                if ((r < 25) && (g < 25) && (b < 25)) {
                    continue;
                } else if ((r > 230) && (g > 230) && (b > 230)) {//might need to use 195 as lower limit
                    output.setValue(col, row, 255);
                    continue;
                }

                float[] cieXY = tmpColorBuffer;
                cieC.rgbToXYChromaticity(r, g, b, cieXY);

                if (cieC.isWhite(cieXY[0], cieXY[1])) {
                    output.setValue(col, row, 255);
                } else {

                    double thetaRadians = cieC.calculateXYTheta(cieXY[0], cieXY[1]);

                    thetaValues[thetaCount] = (float)thetaRadians;

                    pixThetaMap.put(p, Float.valueOf((float)thetaRadians));

                    thetaCount++;
                }
            }
        }

        thetaValues = Arrays.copyOf(thetaValues, thetaCount);

        createAndApplyKMPPMapping(output, pixThetaMap, thetaValues, kColors);

        return output;
    }

    private void createAndApplyKMPPMapping(GreyscaleImage output,
        Map<PairInt, Float> pixThetaMap, float[] thetaValues,
        final int kColors) {

        //TODO: assert kColors.  The invoker is reserving 2 bands for
        // B & W, so nBins should probably be (kColors - 2)...
        // correct this for the invoker when testing
        int nBins = kColors;

        KMeansPlusPlusFloat kmpp = new KMeansPlusPlusFloat();
        kmpp.computeMeans(nBins, thetaValues);

        float minValue = kmpp.getMinValue();
        float maxValue = kmpp.getMaxValue();

        float[] binCenters = kmpp.getCenters();

        Iterator<Map.Entry<PairInt, Float> > iter = pixThetaMap.entrySet().iterator();

        while (iter.hasNext()) {

            Map.Entry<PairInt, Float> entry = iter.next();

            PairInt p = entry.getKey();

            float theta = entry.getValue().floatValue();

            for (int i = 0; i < binCenters.length; i++) {

                float vc = binCenters[i];

                float bisectorBelow = ((i - 1) > -1) ?
                    ((binCenters[i - 1] + vc) / 2) : minValue;

                float bisectorAbove = ((i + 1) > (binCenters.length - 1)) ?
                    maxValue : ((binCenters[i + 1] + vc) / 2);

                if ((theta >= bisectorBelow) && (theta <= bisectorAbove)) {

                    //TODO: check this
                    int mappedValue = 255 - nBins + i;

                    output.setValue(p.getX(), p.getY(), mappedValue);

                    break;
                }
            }

            /*
            // if binCenters is ordered, use binary search for faster results
            int idx = Arrays.binarySearch(startBins, theta);

            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }
            int mappedValue = 255 - startBins.length + idx;

            output.setValue(p.getX(), p.getY(), mappedValue);
            */
        }
    }

    /**
     * for the given thetaValues, create a histogram of kColors bins and then
     * apply the bins to points with those values in the output image.
     * The values in the output image are populated from high value of
     * roughly 254 down to kColors in lower pixel value.
     * The range is approximate because some histograms have gaps for a color
     * bin, so those are not mapped to the final image.
     *
     * runtime complexity is O(N * lg_2(N)).
     * @param output
     * @param pixThetaMap
     * @param thetaValues
     * @param kColors
     */
    private void createAndApplyHistMapping(GreyscaleImage output,
        Map<PairInt, Float> pixThetaMap, float[] thetaValues,
        final int kColors) {

        float minValue = MiscMath.findMin(thetaValues);
        float maxValue = MiscMath.findMax(thetaValues);

        log.fine("minTheta=" + (minValue * 180./Math.PI) +
            " maxTheta=" + (maxValue * 180./Math.PI));

        int nReserved = 254 - kColors;

        HistogramHolder hist = Histogram.createSimpleHistogram(minValue,
            maxValue, (256 - nReserved - 1), thetaValues,
            Errors.populateYErrorsBySqrt(thetaValues));

        /*
        try {
            hist.plotHistogram("cie XY theta histogram", "cieXY_hist_"
                + MiscDebug.getCurrentTimeFormatted());
        } catch (Exception e) {}
        */

        int nonZeroCount = 0;
        for (int i = 0; i < hist.getXHist().length; i++) {
            int c = hist.getYHist()[i];
            if (c > 0) {
                nonZeroCount++;
            }
        }

        float[] startBins = new float[nonZeroCount];

        float halfBinWidth = (hist.getXHist()[1] - hist.getXHist()[0])/2.f;

        nonZeroCount = 0;
        for (int i = 0; i < hist.getXHist().length; i++) {
            int c = hist.getYHist()[i];
            if (c > 0) {
                startBins[nonZeroCount] = hist.getXHist()[i] - halfBinWidth;
                nonZeroCount++;
            }
        }

        Iterator<Map.Entry<PairInt, Float> > iter = pixThetaMap.entrySet().iterator();

        // O(N * lg_2(N))
        while (iter.hasNext()) {

            Map.Entry<PairInt, Float> entry = iter.next();

            PairInt p = entry.getKey();

            float theta = entry.getValue().floatValue();

            int idx = Arrays.binarySearch(startBins, theta);

            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }

            int mappedValue = 255 - startBins.length + idx;

            output.setValue(p.getX(), p.getY(), mappedValue);
        }
    }

    /**
     * NOT READY FOR USE.  STILL EXPERIMENTING.
     *
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to CIEXY Lab color space, then creates a map of the
     * CIEX and CIEY points and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters of points in CIE X, CIEY space,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * cie xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param useBlur
     * @return
     */
    public List<Set<PairInt>> calculateUsingCIEXYAndClustering(ImageExt input,
        boolean useBlur) {

        //TODO: improve the clustering results in two ways:
        // (1) for smaller ciexy clusters, merge with adjacent clusters if
        //     similar color
        // (2) any pixel with 7 neighbors of same color should be that color too

        if (useBlur) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        //TODO: consider making a segmentation method using CIEXY theta
        // for x and the frequency for y, both scaled to numerically
        // resolvable range < max of 5000.
        // this would be good to compare to the method here which
        // uses CIE XY Theta followed by histograms or KMeans++.
        // No need to specify the number of bins before use for suggested
        // version.

        //NOTE: the method needs to have gaps in the data given to it
        //    that is a lack of points for some region between the
        //    min and max of x and y data in integer space

        // max = 6250 unless reduce space complexity
        float factor = 2000;// learn this from numerical resolution

        // then subtract the minima in both cieX and cieY

        int minCIEX = Integer.MAX_VALUE;
        int minCIEY = Integer.MAX_VALUE;
        int maxCIEX = Integer.MIN_VALUE;
        int maxCIEY = Integer.MIN_VALUE;

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        // -------- debug -------
        int nGrey = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            nGrey += entry.getValue().size();
        }
        assert((points0.size() + blackPixels.size() + nGrey +
            whitePixels.size()) == input.getNPixels());
        // -------- end debug -------

        List<Set<PairInt>> greyPixelGroups = groupByPeaks(greyPixelMap);

        // ------- debug -------
        int nGrey2 = 0;
        for (Set<PairInt> set : greyPixelGroups) {
            nGrey2 += set.size();
        }
        assert(nGrey == nGrey2);
        // ------- end debug =====

        Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap0 =
            new HashMap<PairIntWithIndex, List<PairIntWithIndex>>();

        for (PairInt p : points0) {
            int idx = input.getInternalIndex(p.getX(), p.getY());
            float cx = input.getCIEX(idx);
            float cy = input.getCIEY(idx);

            int cieXInt = Math.round(factor * cx);
            int cieYInt = Math.round(factor * cy);

            PairIntWithIndex p0 = new PairIntWithIndex(cieXInt, cieYInt, idx);
            List<PairIntWithIndex> list = pointsMap0.get(p0);
            if (list == null) {
                list = new ArrayList<PairIntWithIndex>();
                pointsMap0.put(p0, list);
            }
            list.add(p0);

            if (cieXInt < minCIEX) {
                minCIEX = cieXInt;
            }
            if (cieYInt < minCIEY) {
                minCIEY = cieYInt;
            }
            if (cieXInt > maxCIEX) {
                maxCIEX = cieXInt;
            }
            if (cieYInt > maxCIEY) {
                maxCIEY = cieYInt;
            }
        }

        Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap =
            new HashMap<PairIntWithIndex, List<PairIntWithIndex>>();

        // subtract minima from the points
        for (PairIntWithIndex p : pointsMap0.keySet()) {

            int x = p.getX() - minCIEX;
            int y = p.getY() - minCIEY;

            PairIntWithIndex p2 = new PairIntWithIndex(x, y, p.getPixIndex());
            List<PairIntWithIndex> list2 = pointsMap.get(p2);
            if (list2 == null) {
                list2 = new ArrayList<PairIntWithIndex>();
                pointsMap.put(p2, list2);
            }
            // because this is a list, this will eventually be present twice:
            //list2.add(p2);

            for (PairIntWithIndex p0 : pointsMap0.get(p)) {
                PairIntWithIndex p3 = new PairIntWithIndex(
                    p0.getX() - minCIEX, p0.getY() - minCIEY, p0.getPixIndex());
                list2.add(p3);
            }
        }
        maxCIEX -= minCIEX;
        maxCIEY -= minCIEY;

        // frequency of colors:
        Map<PairIntWithIndex, Integer> freqMap = new
            HashMap<PairIntWithIndex, Integer>();
        for (Map.Entry<PairIntWithIndex, List<PairIntWithIndex>> entry :
            pointsMap.entrySet()) {
            int c = entry.getValue().size();
            freqMap.put(entry.getKey(), Integer.valueOf(c));
        }

        // ----- debug ---
        int nGreyBW = nGrey + blackPixels.size() + whitePixels.size();
        int nTot = 0;
        //Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap
        for (Entry<PairIntWithIndex, List<PairIntWithIndex>> entry : pointsMap.entrySet()) {
            nTot += entry.getValue().size();
        }
        nTot += nGreyBW;
        log.info("nTot=" + nTot + " nPixels=" + input.getNPixels());
        assert(nTot == input.getNPixels());

        // plot the points as an image to see the data first
        GreyscaleImage img = new GreyscaleImage(maxCIEX + 1, maxCIEY + 1);
        for (com.climbwithyourfeet.clustering.util.PairInt p : pointsMap.keySet()) {
            img.setValue(p.getX(), p.getY(), 255);
        }
        try {
            ImageIOHelper.writeOutputImage(
                ResourceFinder.findDirectory("bin") + "/dt_input.png", img);
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug

        //Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap

        DTClusterFinder<PairIntWithIndex> clusterFinder
            = new DTClusterFinder<PairIntWithIndex>(pointsMap.keySet(),
            maxCIEX + 1, maxCIEY + 1);

        clusterFinder.setToDebug();

        // to recover every point, set limit to 1
        clusterFinder.setMinimumNumberInCluster(1);

        clusterFinder.calculateCriticalDensity();

        clusterFinder.findClusters();

        log.info("clustering critical density=" + clusterFinder.getCriticalDensity());

        int nGroups = clusterFinder.getNumberOfClusters();

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();

        for (int k = 0; k < nGroups; ++k) {

            Set<PairIntWithIndex> group = clusterFinder.getCluster(k);

            Set<PairInt> coordPoints = new HashSet<PairInt>();

            for (PairIntWithIndex p : group) {

                int idx = p.getPixIndex();
                int xCoord = input.getCol(idx);
                int yCoord = input.getRow(idx);

                PairInt pCoord = new PairInt(xCoord, yCoord);
                coordPoints.add(pCoord);

                // include the other points of/ same color
                List<PairIntWithIndex> list = pointsMap.get(p);
                assert(list != null);
                for (PairIntWithIndex p3 : list) {
                    int idx3 = p3.getPixIndex();
                    int xCoord3 = input.getCol(idx3);
                    int yCoord3 = input.getRow(idx3);
                    pCoord = new PairInt(xCoord3, yCoord3);
                    coordPoints.add(pCoord);
                }
            }

            groupList.add(coordPoints);
        }

        // ------ debug ---------
        nTot = 0;
        for (Set<PairInt> set : groupList) {
            nTot += set.size();
        }
        nTot += nGreyBW;
        log.info("nTot=" + nTot + " nPixels=" + input.getNPixels());
        assert(nTot == input.getNPixels());
        // ------ end debug -----

        mergeOrAppendGreyWithOthers(input, greyPixelGroups, groupList,
            blackPixels, whitePixels);

        // add back in blackPixels and whitePixels
        groupList.add(blackPixels);
        groupList.add(whitePixels);

        int nTot2 = 0;
        for (Set<PairInt> groups : groupList) {
            nTot2 += groups.size();
        }

        log.info("img nPix=" + input.getNPixels() + " nTot2=" + nTot2);
        assert(nTot2 == input.getNPixels());

        return groupList;
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency)
     * and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters in that space (polar CIEXY vs frequency),
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     * The changeable parameters are the scaling of the range of polar angles
     * and the range of frequency in order to place the data between values of 0
     * and a number.
     * <code>The scale factors are double xFactor = 2000. and int yFactor = 2000
     * for example.  Smaller scale factors result in faster runtimes,
     * but must be balanced by a large enough factor to have cluster resulution.
     * </code>
     * Note that the polar angle vs frequency maps are actually partitioned into
     * 4 maps to do density based cluster finding separately.
     * partitionFreqFracs = new float[]{0.03f, 0.15f, 0.25f} is the fraction of
     * the maximum frequency defining a partition.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndClustering(ImageExt input,
        boolean useBlur) {

        return calculateUsingPolarCIEXYAndClustering(input, 2000., 2000, useBlur);
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency)
     * and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters in that space (polar CIEXY vs frequency),
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     * The changeable parameters are the scaling of the range of polar angles
     * and the range of frequency in order to place the data between values of 0
     * and a number.
     * <code>The scale factors are double xFactor = 2000. and int yFactor = 2000
     * for example.  Smaller scale factors result in faster runtimes,
     * but must be balanced by a large enough factor to have cluster resulution.
     * </code>
     * Note that the polar angle vs frequency maps are actually partitioned into
     * 4 maps to do density based cluster finding separately.
     * partitionFreqFracs = new float[]{0.03f, 0.15f, 0.25f} is the fraction of
     * the maximum frequency defining a partition.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param thetaRange range to scale the values of cie xy polar theta values
     * to (the range affects the speed because dynamic programming is used).
     * @param thetaFrequencyRange range to scale the values of frequencies of polar
     * theta values to (the range affects the speed because dynamic programming is used).
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndClustering(ImageExt input,
        double thetaRange, int thetaFrequencyRange, boolean useBlur) {

        if (useBlur) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        //making a segmentation method using CIEXY theta
        // for x and the frequency for y, both scaled to numerically
        // resolvable range < max of 5000.
        // this would be good to compare to the method here which
        // uses CIE XY Theta followed by histograms or KMeans++.
        // No need to specify the number of bins before use for suggested
        // version.

        int w = input.getWidth();
        int h = input.getHeight();

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        double[] minMaxTheta0 = findMinMaxTheta(input, points0);

        log.info("for all non-white and non-black, minTheta=" + minMaxTheta0[0]
            + " maxTheta=" + minMaxTheta0[1]);

        // -------- debug -------
        int nGrey = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            nGrey += entry.getValue().size();
        }
        int nGreyBW = + blackPixels.size() + nGrey + whitePixels.size();
        int nTot = (points0.size() + nGreyBW);
        log.info("img nPix=" + input.getNPixels() + " nTot=" + nTot);
        assert(nTot == input.getNPixels());
        // -------- end debug -------

        List<Set<PairInt>> greyPixelGroups = groupByPeaks(greyPixelMap);

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>(greyPixelGroups.size());

        if ((minMaxTheta0[1] - minMaxTheta0[0]) == 0) {
             if (points0.isEmpty()) {
                 groupList.add(points0);
             }
             if (!blackPixels.isEmpty()) {
                 groupList.add(blackPixels);
             }
             for (Set<PairInt> set : greyPixelGroups) {
                groupList.add(set);
             }
             if (!whitePixels.isEmpty()) {
                 groupList.add(whitePixels);
             }
             return groupList;
        }

        double[] minMaxTheta = new double[2];
        int[] minMaxFreq = new int[2];
        double thetaFactor0 = thetaRange/(minMaxTheta0[1] - minMaxTheta0[0]);
        Map<Integer, List<PairInt>> thetaPointMap = createThetaCIEXYMap(points0,
            input, minMaxTheta0[0], thetaFactor0, minMaxTheta, minMaxFreq);

        // ---- debug ------
        nTot = nGreyBW;
        for (Map.Entry<Integer, List<PairInt>> entry : thetaPointMap.entrySet()) {
            nTot += entry.getValue().size();
        }
        log.info("img nPix=" + input.getNPixels() + " nTot=" + nTot);
        assert(nTot == input.getNPixels());
        // ----- end debug ------

        /* ---- create frequency maps partitioned by given fractions ----
        starting w/ partitions at 3 percent (maybe discard below),
            15 percent, and 25 percent resulting in 4 maps

        For each map:
            key is pairint w/ x=theta, y=freq,
            value is all pixels having same key
        */

        final float[] partitionFreqFracs = new float[]{0.03f, 0.15f, 0.25f};

        List<Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>>> thetaFreqMaps =
            partitionIntoFrequencyMaps(input, thetaPointMap,
                partitionFreqFracs, minMaxFreq[1]);

        //---- debug, assert number of pixels ----
        nTot = nGreyBW;
        for (Map<com.climbwithyourfeet.clustering.util.PairInt,
        List<PairIntWithIndex>> map : thetaFreqMaps) {
            for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> entry : map.entrySet()) {
                nTot += entry.getValue().size();
            }
        }
        log.info("img nPix=" + input.getNPixels() + " nTot=" + nTot);
        assert(nTot == input.getNPixels());

        // TODO: handle wrap around values!
        //    if there are points at 0 and 360, and a gap elsewhere, can
        //    shift the values so the gap is at 360 instead.

        // ------ TODO: rescale each map by frequencies to span ~1000 pixels -----
        rescaleKeys(thetaFreqMaps, thetaFrequencyRange);

        int nTot2 = 0;

        for (int i = 1; i < thetaFreqMaps.size(); ++i) {

            Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> thetaFreqMapI = thetaFreqMaps.get(i);

            int[] maxXY = findMaxXY(thetaFreqMapI.keySet());

            if (maxXY[0] < 0 || maxXY[1] <= 0) {
                continue;
            }

            // ----- debug ---
            // plot the points as an image to see the data first
            GreyscaleImage img = new GreyscaleImage(maxXY[0] + 1, maxXY[1] + 1);
            for (com.climbwithyourfeet.clustering.util.PairInt p : thetaFreqMapI.keySet()) {
                img.setValue(p.getX(), p.getY(), 255);
            }
            try {
                ImageIOHelper.writeOutputImage(
                    ResourceFinder.findDirectory("bin") + "/dt2_input_" + i
                        + "_.png", img);
            } catch (IOException ex) {
                Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                    null, ex);
            }
            // --- end debug

            // Map<com.climbwithyourfeet.clustering.util.PairInt,
            //     List<PairIntWithIndex>> thetaFreqMapI

            // map w/ key=(theta, freq) value=collection of coords
            DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt>
                clusterFinder
                = new DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt>(
                    thetaFreqMapI.keySet(), maxXY[0] + 1, maxXY[1] + 1);

            clusterFinder.setToDebug();

            // to recover every point, set limit to 1
            clusterFinder.setMinimumNumberInCluster(1);

            clusterFinder.calculateCriticalDensity();

            clusterFinder.findClusters();

            int nGroups = clusterFinder.getNumberOfClusters();

            for (int k = 0; k < nGroups; ++k) {

                Set<com.climbwithyourfeet.clustering.util.PairInt> group
                    = clusterFinder.getCluster(k);

                Set<PairInt> coordPoints = new HashSet<PairInt>();

                for (com.climbwithyourfeet.clustering.util.PairInt pThetaFreq : group) {

                    // include the other points of same ciexy theta, freq
                    List<PairIntWithIndex> list = thetaFreqMapI.get(pThetaFreq);
                    assert (list != null);
                    for (PairIntWithIndex p3 : list) {
                        int idx3 = p3.getPixIndex();
                        int xCoord3 = input.getCol(idx3);
                        int yCoord3 = input.getRow(idx3);
                        PairInt pCoord = new PairInt(xCoord3, yCoord3);
                        coordPoints.add(pCoord);
                    }
                }

                nTot2 += coordPoints.size();

                groupList.add(coordPoints);
            }
        }

        mergeOrAppendGreyWithOthers(input, greyPixelGroups, groupList,
            blackPixels, whitePixels);

        // add back in blackPixels and whitePixels
        groupList.add(blackPixels);
        groupList.add(whitePixels);

        //TODO: assert npixels

        return groupList;
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction =0.03 of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.  The list is sorted by size and
     * set to values from 255 to 0.  If there are more than 255 clusters, the
     * remaining (smaller) clusters are pixels with value 0.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     * @param input
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public GreyscaleImage applyUsingPolarCIEXYAndFrequency(ImageExt input,
        boolean useBlur) {

        float fracFreqLimit = 0.03f;

        return applyUsingPolarCIEXYAndFrequency(input, fracFreqLimit, useBlur);
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction =0.03 of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.  The list is sorted by size and
     * set to values from 255 to 0.  If there are more than 255 clusters, the
     * remaining (smaller) clusters are pixels with value 0.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     * @param input
     * @param fracFreqLimit
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public GreyscaleImage applyUsingPolarCIEXYAndFrequency(ImageExt input,
        final float fracFreqLimit, boolean useBlur) {

        int w = input.getWidth();
        int h = input.getHeight();

        List<Set<PairInt>> clusters = calculateUsingPolarCIEXYAndFrequency(
            input, fracFreqLimit, useBlur);

        int n = clusters.size();
        //assert(n < 256);

        // sort indexes by set size
        int maxSize = Integer.MIN_VALUE;
        int[] indexes = new int[n];
        int[] sizes = new int[n];
        for (int i = 0; i < n; ++i) {
            indexes[i] = i;
            sizes[i] = clusters.get(i).size();
            if (sizes[i] > maxSize) {
                maxSize = sizes[i];
            }
        }
        CountingSort.sort(sizes, indexes, maxSize);

        int delta = 256/clusters.size();
        if (delta == 0) {
            delta = 1;
        }

        GreyscaleImage img2 = new GreyscaleImage(w, h);
        for (int k = 0; k < n; ++k) {

            int idx = indexes[n - k - 1];

            Set<PairInt> set = clusters.get(idx);

            int v = 255 - delta*k;
            if (v < 1) {
                continue;
            }

            for (PairInt p : set) {
                img2.setValue(p.getX(), p.getY(), v);
            }
        }

        return img2;
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction =0.03 of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     * @param input
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndFrequency(ImageExt input,
        boolean useBlur) {

        float fracFreqLimit = 0.03f;

        return calculateUsingPolarCIEXYAndFrequency(input, fracFreqLimit, useBlur);
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     * The changeable parameter is the fracFreqLimit.  Larger values exclude
     * smaller frequency peaks.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input image to find color clusters within
     * @param fracFreqLimit fraction of the maximum above which peaks will be found
     * @param useBlur if true, performs a gaussian blur of sigma=1 before finding
     * clusters.
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndFrequency(ImageExt input,
        float fracFreqLimit, boolean useBlur) {

        if (useBlur) {
            input = input.copyToImageExt();
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        //making a segmentation method using CIEXY polar theta
        // and the number of points with those colors.
        // choosing the peaks to be the cluster centers, then
        // gathering the pixels by proximity to the theta peaks
        // and when equidistant, chooses the largest peak.

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        double[] minMaxTheta0 = findMinMaxTheta(input, points0);

        log.info("for all non-white and non-black, minTheta=" + minMaxTheta0[0]
            + " maxTheta=" + minMaxTheta0[1]);

        // -------- debug -------
        int nGrey = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            nGrey += entry.getValue().size();
        }
        assert((points0.size() + blackPixels.size() + nGrey +
            whitePixels.size()) == input.getNPixels());
        // -------- end debug -------

        List<Set<PairInt>> greyPixelGroups = groupByPeaks(greyPixelMap);

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>(greyPixelGroups.size());

        if ((minMaxTheta0[1] - minMaxTheta0[0]) == 0) {
             if (points0.isEmpty()) {
                 groupList.add(points0);
             }
             if (!blackPixels.isEmpty()) {
                 groupList.add(blackPixels);
             }
             for (Set<PairInt> set : greyPixelGroups) {
                groupList.add(set);
             }
             if (!whitePixels.isEmpty()) {
                 groupList.add(whitePixels);
             }
             return groupList;
        }

        /* ----- create a map of theta and frequency ----
        need to find the peaks in frequency for frequencies larger than about
        3 percent of max frequency
        but don't want to use a spline3 to smooth, so will average every
        few pixels.
        */

        int binWidth = 3;
        Map<Integer, Collection<PairInt>> thetaPointMap = createThetaCIEXYMap(
            points0, input, binWidth);

        int n = (360/binWidth) + 1;

        int[] orderedThetaKeys = new int[n];
        for (int i = 0; i < n; ++i) {
            orderedThetaKeys[i] = i;
        }
        int maxFreq = Integer.MIN_VALUE;
        int nTot = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : thetaPointMap.entrySet()) {
            int count = entry.getValue().size();
            if (count > maxFreq) {
                maxFreq = count;
            }
            nTot += entry.getValue().size();
        }
        nTot += (blackPixels.size() + nGrey + whitePixels.size());
        assert(nTot == input.getNPixels());

        /*
        TODO: this is where the DTClusterFinder would be good to use to find
        the peaks.
        */

        PairIntArray peaks = findPeaksInThetaPointMap(orderedThetaKeys,
            thetaPointMap,
            Math.round(fracFreqLimit * maxFreq));

        /*
        // ----- debug ---
        // plot the points as an image to see the data first
        int[] minMaxXY = MiscMath.findMinMaxXY(peaks);
        int nPoints = 0;
        int maxX = Integer.MIN_VALUE;
        int maxY = Integer.MIN_VALUE;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                continue;
            }
            int y = list.size();
            nPoints++;
            if (key.intValue() > maxX) {
                maxX = key.intValue();
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        float[] xPoints = new float[nPoints];
        float[] yPoints = new float[nPoints];
        int count = 0;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                continue;
            }
            int y = list.size();
            xPoints[count] = key.intValue();
            yPoints[count] = y;
            count++;
        }
        try {
            //maxY=2000;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(0, maxX, 0, maxY, xPoints, yPoints, xPoints, yPoints, "cieXY theta vs freq");
            plotter.writeFile("_segmentation3_");
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug
        */

        if (peaks.getN() == 0) {
            for (Map.Entry<Integer, Collection<PairInt>> entry : thetaPointMap.entrySet()) {
                groupList.add(new HashSet<PairInt>(entry.getValue()));
            }
            groupList.add(blackPixels);
            groupList.add(whitePixels);
            return groupList;
        }
        for (int i = 0; i < peaks.getN(); ++i) {
            groupList.add(new HashSet<PairInt>());
        }

        /* traverse in ordered manner thetaPointMap to compare to current theta position
           w.r.t. peaks
           then place it in the groupsList.
           points before the first peak are compared with last peak too for wrap around.
        */
        int currentPeakIdx = -1;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                continue;
            }
            int idx = -1;
            if ((currentPeakIdx == -1) || (currentPeakIdx == (peaks.getN() - 1))) {
                int diffL, diffF;
                if (currentPeakIdx == -1) {
                    diffL = key.intValue() + 360 - peaks.getX(peaks.getN() - 1);
                    diffF = peaks.getX(0) - key.intValue();
                    if (diffF == 0) {
                        currentPeakIdx = 0;
                    }
                } else {
                    diffL = key.intValue() - peaks.getX(currentPeakIdx);
                    diffF = peaks.getX(0) + 360 - key.intValue();
                }
                if (diffL < diffF) {
                    idx = peaks.getN() - 1;
                } else if (diffL == diffF) {
                    int freqL = peaks.getY(peaks.getN() - 1);
                    int freqF = peaks.getY(0);
                    if (freqL < freqF) {
                        idx = peaks.getN() - 1;
                    } else {
                        idx = 0;
                    }
                } else {
                    idx = 0;
                }
            } else {
                // this has to update currentPeakIdx
                int diffP = key.intValue() - peaks.getX(currentPeakIdx);
                int diffN = peaks.getX(currentPeakIdx + 1) - key.intValue();
                if (diffN == 0) {
                    currentPeakIdx++;
                    idx = currentPeakIdx;
                } else {
                    if (diffP < diffN) {
                        idx = currentPeakIdx;
                    } else if (diffP == diffN) {
                        int freqP = peaks.getY(currentPeakIdx);
                        int freqN = peaks.getY(currentPeakIdx + 1);
                        if (freqP < freqN) {
                            idx = currentPeakIdx;
                        } else {
                            idx = currentPeakIdx + 1;
                        }
                    } else {
                        idx = currentPeakIdx + 1;
                    }
                }
            }
            assert(idx != -1);
            groupList.get(idx).addAll(list);
        }

        mergeOrAppendGreyWithOthers(input, greyPixelGroups, groupList,
            blackPixels, whitePixels);

        // add back in blackPixels and whitePixels
        /*if (!blackPixels.isEmpty()) {
            groupList.add(blackPixels);
        }
        if (!whitePixels.isEmpty()) {
            groupList.add(whitePixels);
        }

        int nTot2 = 0;
        for (Set<PairInt> groups : groupList) {
            nTot2 += groups.size();
        }

        log.info("img nPix=" + input.getNPixels() + " nTot2=" + nTot2);
        assert(nTot2 == input.getNPixels());
        */
        return groupList;
    }

    private Map<Integer, List<PairInt>> createThetaCIEXYMap(Set<PairInt>
        points0, ImageExt input, double minTheta, double thetaFactor,
        double[] outputMinMaxTheta, int[] outputMinMaxFreq) {

        CIEChromaticity cieC = new CIEChromaticity();

        // key = theta, value = pixels having that key
        Map<Integer, List<PairInt>> thetaPointMap = new HashMap<Integer, List<PairInt>>();

        double minTheta0 = minTheta;
        outputMinMaxTheta[0] = Double.MIN_VALUE;
        outputMinMaxTheta[1] = Double.MAX_VALUE;

        for (PairInt p : points0) {

            int idx = input.getInternalIndex(p.getX(), p.getY());

            float cx = input.getCIEX(idx);
            float cy = input.getCIEY(idx);

            double theta = thetaFactor * (
                (cieC.calculateXYTheta(cx, cy)*180./Math.PI) - minTheta0);

            Integer thetaCIEXY = Integer.valueOf((int)Math.round(theta));

            List<PairInt> list = thetaPointMap.get(thetaCIEXY);
            if (list == null) {
                list = new ArrayList<PairInt>();
                thetaPointMap.put(thetaCIEXY, list);
            }
            list.add(p);

            if (theta < outputMinMaxTheta[0]) {
                outputMinMaxTheta[0] = theta;
            }
            if (theta > outputMinMaxTheta[1]) {
                outputMinMaxTheta[1] = theta;
            }
        }

        outputMinMaxFreq[0] = Integer.MAX_VALUE;
        outputMinMaxFreq[1] = Integer.MIN_VALUE;
        for (Map.Entry<Integer, List<PairInt>> entry : thetaPointMap.entrySet()) {
            int count = entry.getValue().size();
            if (count < outputMinMaxFreq[0]) {
                outputMinMaxFreq[0] = count;
            }
            if (count > outputMinMaxFreq[1]) {
                outputMinMaxFreq[1] = count;
            }
        }

        return thetaPointMap;
    }

    private List<Map<com.climbwithyourfeet.clustering.util.PairInt,
    List<PairIntWithIndex>>> partitionIntoFrequencyMaps(
    ImageExt input, Map<Integer, List<PairInt>> thetaPointMap,
    float[] partitionFreqFracs, int maxFreq) {

        List<Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>>>
            mapsList =
                new ArrayList<Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>>>();

        int nMaps = partitionFreqFracs.length + 1;

        for (int i = 0; i < nMaps; ++i) {
            mapsList.add(
                new HashMap<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>>());
        }

        final int[] partitionFreqs = new int[partitionFreqFracs.length];
        for (int i = 0; i < partitionFreqs.length; ++i) {
            partitionFreqs[i] = Math.round(partitionFreqFracs[i]*maxFreq);
        }

        int[] maxXY = new int[2];
        Arrays.fill(maxXY, Integer.MIN_VALUE);

        for (Map.Entry<Integer, List<PairInt>> entry : thetaPointMap.entrySet()) {

            Integer theta = entry.getKey();

            int count = entry.getValue().size();

            List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
            for (PairInt p : entry.getValue()) {
                int pixIdx = input.getInternalIndex(p.getX(), p.getY());
                PairIntWithIndex p2 = new PairIntWithIndex(theta.intValue(),
                    count, pixIdx);
                list.add(p2);
            }

            int n = partitionFreqs.length;
            int idx = 0;
            for (int i = 0; i < n; ++i) {
                if (i == 0) {
                    if (count < partitionFreqs[0]) {
                        idx = 0;
                        break;
                    }
                } else if ((i == (n - 1)) && count >= partitionFreqs[i]) {
                    idx = n;// one past partitions is last list bin
                    break;
                } else {
                    if ((count >= partitionFreqs[i - 1]) && (count < partitionFreqs[i])) {
                        idx = i;
                        break;
                    }
                }
            }

            // this is unique to all maps, not stomping on existing key
            mapsList.get(idx).put(
                new com.climbwithyourfeet.clustering.util.PairInt(
                    theta.intValue(), count), list);

            if (theta.intValue() > maxXY[0]) {
                maxXY[0] = theta.intValue();
            }
            if (count > maxXY[1]) {
                maxXY[1] = count;
            }
        }

        // ----- temporary print of all pixels -------
        // ----- debug ---
        GreyscaleImage img = new GreyscaleImage(maxXY[0] + 1, maxXY[1] + 1);
        for (int i = 0; i < mapsList.size(); ++i) {
            for (com.climbwithyourfeet.clustering.util.PairInt p : mapsList.get(i).keySet()) {
                img.setValue(p.getX(), p.getY(), 255);
            }
        }
        try {
            ImageIOHelper.writeOutputImage(
                ResourceFinder.findDirectory("bin") + "/dt2_input.png", img);
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug

        return mapsList;
    }

    private void rescaleKeys(
        List<Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>>> thetaFreqMaps, int scaleTo) {

        // --- can remove the count after debugging ----
        int nTotBefore = 0;
        for (int i = 0; i < thetaFreqMaps.size(); ++i) {
            for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> entry : thetaFreqMaps.get(i).entrySet()) {
                nTotBefore += entry.getValue().size();
            }
        }

        for (int i = 0; i < thetaFreqMaps.size(); ++i) {

            Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> thetaFreqMap = thetaFreqMaps.get(i);

            int[] minMax = findMinMaxOfKeyYs(thetaFreqMap.keySet());

            Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> thetaFreqMap2 = rescaleKeyYs(
                    thetaFreqMap, scaleTo, minMax);

            thetaFreqMaps.set(i, thetaFreqMap2);
        }

        // --- can remove the count after debugging ----
        int nTotAfter = 0;
        for (int i = 0; i < thetaFreqMaps.size(); ++i) {
            for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> entry : thetaFreqMaps.get(i).entrySet()) {
                nTotAfter += entry.getValue().size();
            }
        }
        assert(nTotBefore == nTotAfter);
    }

    private int[] findMinMaxOfKeyYs(
        Set<com.climbwithyourfeet.clustering.util.PairInt> keySet) {

        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;

        for (com.climbwithyourfeet.clustering.util.PairInt p : keySet) {
            int y = p.getY();
            if (y < min) {
                min = y;
            }
            if (y > max) {
                max = y;
            }
        }
        return new int[]{min, max};
    }

    private Map<com.climbwithyourfeet.clustering.util.PairInt,
    List<PairIntWithIndex>> rescaleKeyYs(
    Map<com.climbwithyourfeet.clustering.util.PairInt,
    List<PairIntWithIndex>> thetaFreqMap, final int scaleTo, final int[] minMaxY) {

        Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>> scaledMap
            = new HashMap<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>>();

        if ((minMaxY[1] - minMaxY[0]) == 0) {
            return thetaFreqMap;
        }

        float factor = scaleTo/(minMaxY[1] - minMaxY[0]);

        for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>> entry : thetaFreqMap.entrySet()) {

            com.climbwithyourfeet.clustering.util.PairInt p = entry.getKey();

            int y = Math.round(factor * (p.getY() - minMaxY[0]));

            com.climbwithyourfeet.clustering.util.PairInt p2 = new
                com.climbwithyourfeet.clustering.util.PairInt(p.getX(), y);

            scaledMap.put(p2, entry.getValue());
        }

        return scaledMap;
    }

    private int[] findMaxXY(Set<com.climbwithyourfeet.clustering.util.PairInt>
        keySet) {

        int maxX = Integer.MIN_VALUE;
        int maxY = Integer.MIN_VALUE;

        for (com.climbwithyourfeet.clustering.util.PairInt p : keySet) {
            int x = p.getX();
            int y = p.getY();
            if (x > maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        return new int[]{maxX, maxY};
    }

    private Map<Integer, Collection<PairInt>> createThetaCIEXYMap(Set<PairInt>
        points, ImageExt input, int binWidth) {

        CIEChromaticity cieC = new CIEChromaticity();

        // key = theta, value = pixels having that key
        Map<Integer, Collection<PairInt>> thetaPointMap = new HashMap<Integer, Collection<PairInt>>();

        for (PairInt p : points) {

            int idx = input.getInternalIndex(p.getX(), p.getY());

            float cx = input.getCIEX(idx);
            float cy = input.getCIEY(idx);

            double thetaRadians = cieC.calculateXYTheta(cx, cy);
            double theta = thetaRadians * 180./Math.PI;

            int thetaCIEXY = (int)Math.round(theta);

            Integer binKey = Integer.valueOf(thetaCIEXY/binWidth);

            Collection<PairInt> list = thetaPointMap.get(binKey);
            if (list == null) {
                list = new ArrayList<PairInt>();
                thetaPointMap.put(binKey, list);
            }
            list.add(p);
        }

        return thetaPointMap;
    }

    /**
     * find peaks in the theta point map above lower limit.
     * @param orderedThetaKeys
     * @param thetaPointMap
     * @param limit
     * @return
     */
    protected PairIntArray findPeaksInThetaPointMap(final int[] orderedThetaKeys,
        final Map<Integer, Collection<PairInt>> thetaPointMap, final int limit) {

        int lastKey = -1;
        int lastValue = -1;
        boolean isIncr = false;
        PairIntArray peaks = new PairIntArray();
        int nInMap = 0;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                if ((nInMap > 0) && isIncr && (lastValue > limit)) {
                    peaks.add(lastKey, lastValue);
                }
                lastKey = key.intValue();
                lastValue = 0;
                isIncr = false;
                continue;
            }
            int count = list.size();
            if (nInMap == 1) {
                if (count > lastValue) {
                    isIncr = true;
                } else {
                    if (lastValue > limit) {
                        peaks.add(lastKey, lastValue);
                    }
                    isIncr = false;
                }
            } else if (nInMap != 0) {
                if (isIncr) {
                    if (count < lastValue) {
                        if (lastValue > limit) {
                            peaks.add(lastKey, lastValue);
                        }
                        isIncr = false;
                    }
                } else {
                    if (count > lastValue) {
                        isIncr = true;
                    }
                }
            }
            lastValue = count;
            lastKey = key.intValue();
            nInMap++;
        }
        if ((nInMap > 0) && isIncr && (lastValue > limit)) {
            //checking value at theta=0 to make sure this is a peak
            Integer key = orderedThetaKeys[0];
            if (key.intValue() == 0) {
                Collection<PairInt> list = thetaPointMap.get(key);
                if (list == null) {
                    peaks.add(lastKey, lastValue);
                } else if (list.size() < lastValue) {
                    peaks.add(lastKey, lastValue);
                }
            } else {
                peaks.add(lastKey, lastValue);
            }
        }

        return peaks;
    }

    private List<Set<PairInt>> groupByPeaks(
        Map<Integer, Collection<PairInt>> greyPixelMap) {

        int nTot = 0;
        int minKey = Integer.MAX_VALUE;
        int maxKey = Integer.MIN_VALUE;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            int key = entry.getKey().intValue();
            if (key < minKey) {
                minKey = key;
            }
            if (key > maxKey) {
                maxKey = key;
            }
            nTot += entry.getValue().size();
        }

        int binWidth = 8;
        greyPixelMap = binByKeys(greyPixelMap, minKey, maxKey, binWidth);

        int nTot2 = 0;
        minKey = Integer.MAX_VALUE;
        maxKey = Integer.MIN_VALUE;
        int maxFreq = Integer.MIN_VALUE;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            int key = entry.getKey().intValue();
            if (key < minKey) {
                minKey = key;
            }
            if (key > maxKey) {
                maxKey = key;
            }
            int y = entry.getValue().size();
            if (y > maxFreq) {
                maxFreq = y;
            }

            nTot2 += y;
        }
        assert(nTot == nTot2);

        int count = 0;

        /*
        // --- debug
        float[] xPoints = new float[greyPixelMap.size()];
        float[] yPoints = new float[greyPixelMap.size()];
        for (int i = minKey; i <= maxKey; ++i) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> set = greyPixelMap.get(key);
            if (set == null) {
                continue;
            }
            int y = set.size();
            xPoints[count] = key.intValue();
            yPoints[count] = y;
            count++;
        }
        try {
            //maxY=2000;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(0, maxKey, 0, maxFreq, xPoints, yPoints, xPoints, yPoints, "grey avgRGB vs freq");
            plotter.writeFile("_segmentation3_grey_");
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug
        */

        final int[] orderedKeys = new int[maxKey - minKey + 1];
        count = 0;
        for (int i = minKey; i <= maxKey; ++i) {
            orderedKeys[count] = i;
            count++;
        }
        int limit = (int)(0.03 * maxFreq);
        PairIntArray peaks = findPeaksInThetaPointMap(orderedKeys, greyPixelMap, limit);

        // if there are several peaks within small range of keys, that's noise,
        // so removing them
        filterPeaksIfNoisey(peaks);

        // ---- gather points in greyPixelMap into groups around the peaks ----
        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>(orderedKeys.length + 1);
        for (int i = 0; i < peaks.getN(); ++i) {
            groupList.add(new HashSet<PairInt>());
        }

        if (peaks.getN() == 0) {
            return groupList;
        } else if (peaks.getN() == 1) {
            for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
                Collection<PairInt> set = entry.getValue();
                groupList.get(0).addAll(set);
            }
            return groupList;
        }

        int currentPeakIdx = -1;
        for (int i : orderedKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> set = greyPixelMap.get(key);
            if (set == null) {
                continue;
            }
            int idx = -1;
            if (currentPeakIdx == -1) {
                currentPeakIdx = 0;
                idx = 0;
            } else if (currentPeakIdx == (peaks.getN() - 1)) {
                idx = currentPeakIdx;
            } else {
                // this has to update currentPeakIdx
                int diffP = key.intValue() - peaks.getX(currentPeakIdx);
                int diffN = peaks.getX(currentPeakIdx + 1) - key.intValue();
                if (diffN == 0) {
                    currentPeakIdx++;
                    idx = currentPeakIdx;
                } else {
                    if (diffP < diffN) {
                        idx = currentPeakIdx;
                    } else if (diffP == diffN) {
                        int freqP = peaks.getY(currentPeakIdx);
                        int freqN = peaks.getY(currentPeakIdx + 1);
                        if (freqP < freqN) {
                            idx = currentPeakIdx;
                        } else {
                            idx = currentPeakIdx + 1;
                        }
                    } else {
                        idx = currentPeakIdx + 1;
                    }
                }
            }
            assert(idx != -1);
            groupList.get(idx).addAll(set);
        }

        // ----- debug ----
        int nTot3 = 0;
        for (Collection<PairInt> set : groupList) {
            nTot3 += set.size();
        }
        assert(nTot == nTot3);
        // ----- end debug -----

        return groupList;
    }

    private Map<Integer, Collection<PairInt>> binByKeys(
        Map<Integer, Collection<PairInt>> greyPixelMap,
        int minKey, int maxKey, int binWidth) {

        Map<Integer, Collection<PairInt>> map2
            = new HashMap<Integer, Collection<PairInt>>();

        for (int i = minKey; i <= maxKey; ++i) {

            Integer key = Integer.valueOf(i);

            Collection<PairInt> c = greyPixelMap.get(key);

            if (c == null) {
                continue;
            }

            Integer binKey = Integer.valueOf(i/binWidth);

            Collection<PairInt> c2 = map2.get(binKey);
            if (c2 == null) {
                c2 = new HashSet<PairInt>();
                map2.put(binKey, c2);
            }
            c2.addAll(c);
        }

        return map2;
    }

    private void filterPeaksIfNoisey(PairIntArray peaks) {

        if (peaks.getN() == 0) {
            return;
        }

        int sumDeltaX = 0;
        for (int i = (peaks.getN() - 1); i > 0; --i) {
            sumDeltaX += (peaks.getX(i) - peaks.getX(i - 1));
        }
        //TODO: this may need to be revised:
        // if there are more than 1 peaks per delta x of 5 or so, re-bin by 4
        float deltaX = (float)sumDeltaX/((float)peaks.getN() - 1);
        if (deltaX < 5) {
            // re-bin by 4
            PairIntArray peaks2 = new PairIntArray();
            for (int i = 0; i < peaks.getN(); i += 4) {
                int sumX = 0;
                int sumY = 0;
                int count = 0;
                for (int j = i; j < (i + 5); ++j) {
                    sumX += peaks.getX(i);
                    sumY += peaks.getY(i);
                    count++;
                }
                sumX = Math.round((float)sumX/(float)count);
                sumY = Math.round((float)sumY/(float)count);
                peaks2.add(sumX, sumY);
            }

            peaks.removeRange(0, peaks.getN() - 1);
            peaks.addAll(peaks2);
        }

    }

    private void mergeOrAppendGreyWithOthers(ImageExt input,
        List<Set<PairInt>> greyPixelGroups, List<Set<PairInt>> groupList,
        Set<PairInt> blackPixels, Set<PairInt> whitePixels) {

        Map<PairInt, Integer> colorPixGroupMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < groupList.size(); ++i) {
            Set<PairInt> set = groupList.get(i);
            Integer groupKey = Integer.valueOf(i);
            for (PairInt p : set) {
                colorPixGroupMap.put(p, groupKey);
            }
        }

        // similarity limit for a grey pixel to join adjacent color pixel's cluster
        int limit = 40;

        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        for (Set<PairInt> greyGroup : greyPixelGroups) {

            Set<PairInt> remove = new HashSet<PairInt>();

            for (PairInt greyP : greyGroup) {
                int x = greyP.getX();
                int y = greyP.getY();

                int idx = input.getInternalIndex(x, y);
                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                // ---- check for color similarity ------
                int minDiffRGB = Integer.MAX_VALUE;
                int colorClusterIdx = -1;
                int minDiffBlack = Integer.MAX_VALUE;
                int minDiffWhite = Integer.MAX_VALUE;

                for (int i = 0; i < dxs.length; ++i) {
                    int x2 = x + dxs[i];
                    int y2 = y + dys[i];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);

                    boolean adjIsBlack = blackPixels.contains(p2);
                    boolean adjIsWhite = whitePixels.contains(p2);

                    Integer colorClusterIndex = colorPixGroupMap.get(p2);
                    if ((colorClusterIndex == null) && !adjIsBlack && !adjIsWhite) {
                        continue;
                    }

                    int idx2 = input.getInternalIndex(x2, y2);
                    int r2 = input.getR(idx2);
                    int g2 = input.getG(idx2);
                    int b2 = input.getB(idx2);

                    int diffR = Math.abs(r2 - r);
                    int diffG = Math.abs(g2 - g);
                    int diffB = Math.abs(b2 - b);

                    int diffRGB = diffR + diffG + diffB;

                    if (adjIsBlack) {
                        if (diffR == diffG && diffR == diffB) {
                            minDiffBlack = diffR;
                        }
                    } else if (adjIsWhite) {
                        if (diffR == diffG && diffR == diffB) {
                            minDiffWhite = diffR;
                        }
                    } else {
                        if ((diffRGB < minDiffRGB) && (diffRGB < limit)) {
                            minDiffRGB = diffRGB;
                            colorClusterIdx = (colorClusterIndex == null) ? -1 :
                                colorClusterIndex.intValue();
                        }
                    }
                }
                if (minDiffBlack < 75) {
                    blackPixels.add(greyP);
                    remove.add(greyP);
                } else if (minDiffWhite < 75) {
                    whitePixels.add(greyP);
                    remove.add(greyP);
                } else {
                    if (colorClusterIdx != -1) {
                        //add to color cluster and remove from grey list
                        groupList.get(colorClusterIdx).add(greyP);
                        remove.add(greyP);
                    }
                }
            }
            for (PairInt rm : remove) {
                greyGroup.remove(rm);
            }
        }

        // any remaining points in the grey list should be added as sets to
        // the color pixels list now
        for (int i = 0; i < greyPixelGroups.size(); ++i) {
            Set<PairInt> greyGroup = greyPixelGroups.get(i);
            groupList.add(greyGroup);
        }
        greyPixelGroups.clear();
    }

    protected double[] findMinMaxTheta(ImageExt input, Set<PairInt> points0) {

        CIEChromaticity cieC = new CIEChromaticity();

        double[] minMaxTheta = new double[2];

        for (PairInt p : points0) {
            int x = p.getX();
            int y = p.getY();

            double thetaRadians = cieC.calculateXYTheta(x, y);
            double theta = thetaRadians * 180. / Math.PI;

            if (theta < minMaxTheta[0]) {
                minMaxTheta[0] = theta;
            }
            if (theta > minMaxTheta[1]) {
                minMaxTheta[1] = theta;
            }
        }

        return minMaxTheta;
    }

    /**
     * Find the 4 categories of point color as black, white, grey, and color.
     * CIE XY color space is used to place a pixel in the category, with the
     * additional determination of black, grey, and white intensity limits
     * using histograms within expected intensity limits.
     *
     * @param input
     * @param points0
     * @param blackPixels
     * @param whitePixels
     * @param greyPixelMap
     */
    private void populatePixelLists(ImageExt input, Set<PairInt> points0,
        Set<PairInt> blackPixels, Set<PairInt> whitePixels,
        Map<Integer, Collection<PairInt>> greyPixelMap) {

        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        // looking for limits in peaks of (r,g,b) <= (45,45,45) and > (191,191,191)
        int[] whiteBlackLimits = findByHistogramLimitsForBlackAndWhite(input);
        // overriding:
        int whiteLimit = 245;

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int idx = input.getInternalIndex(i, j);

                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                int avg = (r + g + b)/3;
                if (avg <= whiteBlackLimits[0]) {
                    blackPixels.add(new PairInt(i, j));
                    continue;
                }

                float rgbTot = r + g + b;
                float bDivTot = (float)b/rgbTot;
                float rDivTot = (float)r/rgbTot;
                float gDivTot = (float)g/rgbTot;

                float cx = input.getCIEX(idx);
                float cy = input.getCIEY(idx);

                if (cieC.isWhite2(cx, cy) &&
                    (Math.abs(0.333f - bDivTot) < 0.02f) &&
                    (Math.abs(0.333f - rDivTot) < 0.02f) &&
                    (Math.abs(0.333f - gDivTot) < 0.02f)) {

                    //if (avg > whiteBlackLimits[1]) {
                    if ((r > whiteLimit) && (g > whiteLimit) && (b > whiteLimit)) {
                        whitePixels.add(new PairInt(i, j));
                    } else {
                        Integer avgRGB = Integer.valueOf(avg);
                        Collection<PairInt> set = greyPixelMap.get(avgRGB);
                        if (set == null) {
                            set = new HashSet<PairInt>();
                            greyPixelMap.put(avgRGB, set);
                        }
                        set.add(new PairInt(i, j));
                    }
                } else {
                    points0.add(new PairInt(i, j));
                }
            }
        }
    }

    private void mergeNoise(ImageExt input, List<Set<PairInt>> groupList) {

        Map<PairInt, Integer> pixelToGroupMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < groupList.size(); ++i) {
            Set<PairInt> set = groupList.get(i);
            for (PairInt p : set) {
                pixelToGroupMap.put(p, Integer.valueOf(i));
            }
        }

        final int w = input.getWidth();
        final int h = input.getHeight();

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        float diffLimit = 0.01f;

        for (int i = 0; i < groupList.size(); ++i) {
            Set<PairInt> set = groupList.get(i);
            Set<PairInt> remove = new HashSet<PairInt>();
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                int idx = input.getInternalIndex(x, y);

                float cieX = input.getCIEX(idx);
                float cieY = input.getCIEY(idx);

                Integer groupIndex = pixelToGroupMap.get(p);

                // key=groupIndex, value=number of pixels similar
                Map<Integer, Integer> groupSimilarCount = new HashMap<Integer, Integer>();

                // use cieXY, polar theta of cieXY, or rgb?
                for (int nIdx = 0; nIdx < dxs8.length; ++nIdx) {
                    int x2 = x + dxs8[nIdx];
                    int y2 = y + dys8[nIdx];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    int idx2 = input.getInternalIndex(x2, y2);

                    Integer groupIndex2 = pixelToGroupMap.get(new PairInt(x2, y2));

                    if (groupIndex.intValue() == groupIndex2.intValue()) {
                        continue;
                    }

                    float cieX2 = input.getCIEX(idx2);
                    float cieY2 = input.getCIEY(idx2);

                    float diffCIEX = Math.abs(cieX2 - cieX);
                    float diffCIEY = Math.abs(cieY2 - cieY);

                    if ((diffCIEX > diffLimit) || (diffCIEY > diffLimit)) {
                        continue;
                    }

                    Integer count = groupSimilarCount.get(groupIndex2);
                    if (count == null) {
                        groupSimilarCount.put(groupIndex2, Integer.valueOf(1));
                    } else {
                        groupSimilarCount.put(groupIndex2, Integer.valueOf(count.intValue() + 11));
                    }

                }
                if (groupSimilarCount.isEmpty()) {
                    continue;
                }
                for (Entry<Integer, Integer> entry : groupSimilarCount.entrySet()) {
                    if (entry.getValue() >= 6) {
                        // assign this group to pixel p
                        pixelToGroupMap.put(p, entry.getKey());
                        remove.add(p);
                        groupList.get(entry.getKey().intValue()).add(p);
                        break;
                    }
                }
            }
            for (PairInt rm : remove) {
                set.remove(rm);
            }
        }
    }

    private int[] findByHistogramLimitsForBlackAndWhite(ImageExt input) {

        //looking for limits of (r,g,b) <= (45,45,45) and > (180,180,180)
        int l0 = 45;
        int l0B = 70;
        int h0 = 245;

        List<Integer> avgL = new ArrayList<Integer>();

        List<Integer> avgH = new ArrayList<Integer>();

        for (int i = 0; i < input.getNPixels(); ++i) {
            int r = input.getR(i);
            int g = input.getG(i);
            int b = input.getB(i);
            if ((r <= l0) && (g <= l0) && (b <= l0B)) {
                int avg = (r + g + b)/3;
                avgL.add(Integer.valueOf(avg));
            } else if ((r > h0) && (g > h0) && (b > h0)) {
                int avg = (r + g + b)/3;
                avgH.add(Integer.valueOf(avg));
            }
        }

        int[] limits = new int[2];

        if (avgL.size() > 30) {
            HistogramHolder hist = Histogram.createSimpleHistogram(avgL);
            if (hist == null) {
                limits[0] = l0 - 1;
            } else {
                /*
                try {
                    hist.plotHistogram("black pixels", "black_" + MiscDebug.getCurrentTimeFormatted());
                } catch (IOException ex) {
                    Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
                }
                */
                int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
                int lastZeroIdx = MiscMath.findLastZeroIndex(hist);
                if ((lastZeroIdx != -1) && (lastZeroIdx < hist.getXHist().length)) {
                    limits[0] = Math.round(hist.getXHist()[lastZeroIdx]);
                } else if (yMaxIdx == -1) {
                    limits[0] = l0 - 1;
                } else {
                    //limits[0] = Math.round(hist.getXHist()[indexes.get(0).intValue()]);
                    limits[0] = Math.round(hist.getXHist()[yMaxIdx]);
                }
            }
        } else {
            limits[0] = l0 - 1;
        }

        if (avgH.size() > 30) {

            int[] q = ImageStatisticsHelper.getQuartiles(avgH);
            /*
            HistogramHolder hist = Histogram.createSimpleHistogram(avgH);
            if (hist == null) {
                limits[1] = h0;
            } else {
                List<Integer> indexes = MiscMath.findStrongPeakIndexes(hist, 0.1f);
                try {
                    hist.plotHistogram("hite pixels", "white_" + MiscDebug.getCurrentTimeFormatted());
                } catch (IOException ex) {
                    Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
                }
                if (indexes == null || indexes.isEmpty()) {
                    limits[1] = h0;
                } else {
                    limits[1] = Math.round(hist.getXHist()[indexes.get(0).intValue()]);
                }
            }
            */
            limits[1] = q[3];
        } else {
            limits[1] = h0;
        }

        return limits;
    }

    public Map<PairInt, Integer> calculatePolarCIEXY(ImageExt input, Set<PairInt> points) {

        Map<PairInt, Integer> map = new HashMap<PairInt, Integer>();

        CIEChromaticity cieC = new CIEChromaticity();

        for (PairInt p : points) {

            float cieX = input.getCIEX(p.getX(), p.getY());
            float cieY = input.getCIEY(p.getX(), p.getY());

            double thetaRadians = cieC.calculateXYTheta(cieX, cieY);

            int thetaDegrees = (int)Math.round(thetaRadians * 180./Math.PI);

            map.put(p, Integer.valueOf(thetaDegrees));
        }

        return map;
    }

    public float[] getValues(Map<PairInt, Integer> map) {
        float[] values = new float[map.size()];
        int count = 0;
        for (Entry<PairInt, Integer> entry : map.entrySet()) {
            values[count] = entry.getValue().floatValue();
            count++;
        }
        return values;
    }

    /**
     * NOT READY FOR USE.  STILL EXPERIMENTING.
     *
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to CIEXY Lab color space, then creates a map of the
     * CIEX and CIEY points and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters of points in CIE X, CIEY space,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * cie xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param useBlur
     */
    public void applyPolarCIEXY(ImageExt input, boolean useBlur) {

        //TODO: improve the clustering results in two ways:
        // (1) for smaller ciexy clusters, merge with adjacent clusters if
        //     similar color
        // (2) any pixel with 7 neighbors of same color should be that color too

        if (useBlur) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        Map<PairInt, Integer> clrPolarCIEXYMap = calculatePolarCIEXY(input, points0);
        float[] clrPolarCIEXY = getValues(clrPolarCIEXYMap);
        float binWidth = 1;
        HistogramHolder hist = Histogram.createSimpleHistogram(//binWidth,
            clrPolarCIEXY, Errors.populateYErrorsBySqrt(clrPolarCIEXY));
        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist,
            0.1f);
        int[] binCenters = createBinCenters360(hist, indexes);

        List<Set<PairInt>> colorPixelGroups = assignToNearestPolarCIECluster(
            clrPolarCIEXYMap, binCenters);

        List<Set<PairInt>> greyPixelGroups = groupByPeaksForGrey(input, greyPixelMap);

        ImageExt imgExt = input.copyToImageExt();

        for (PairInt p : blackPixels) {
           imgExt.setRGB(p.getX(), p.getY(), 0, 0, 0);
        }

        for (PairInt p : whitePixels) {
           imgExt.setRGB(p.getX(), p.getY(), 255, 255, 255);
        }

        int gClr = 255;
        int s = 127/colorPixelGroups.size();
        for (Set<PairInt> set : colorPixelGroups) {
            for (PairInt p : set) {
               imgExt.setRGB(p.getX(), p.getY(), 0, gClr, 0);
            }
            gClr -= s;
        }

        /*
        int rClr = 255;
        s = 127/greyPixelMap.size();
        for (Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            for (PairInt p : entry.getValue()) {
                imgExt.setRGB(p.getX(), p.getY(), rClr, 0, 0);
            }
            rClr -= s;
        }
        */
        int rClr = 255;
        s = 127/greyPixelGroups.size();
        for (Set<PairInt> set : greyPixelGroups) {
            for (PairInt p : set) {
               imgExt.setRGB(p.getX(), p.getY(), rClr, 0, 0);
            }
            rClr -= s;
        }

        input.resetTo(imgExt);
    }

    private int[] createBinCenters360(HistogramHolder hist, List<Integer> indexes) {

        if (indexes.isEmpty()) {
            return new int[0];
        }

        int n = indexes.size();

        if (n == 1) {

            int[] binCenters = new int[n + 1];

            /*
            examples for n=1
                --  90
            180
                -- 270
            ----------------
                --  90-176=  360-(176-90) = 274
            4
                -- 270-176=  94
            ----------------
                --  90+(330-180)=  240
            330
                -- 270+(330-180)=  270 + 150 - 360 = 60
            */

            int vc = Math.round(hist.getXHist()[indexes.get(0).intValue()]);

            if (vc == 180) {
                binCenters[0] = 90;
                binCenters[1] = 270;
            } else if (vc < 180) {
                int delta = 180 - vc;
                binCenters[0] = 360 - (delta - 90);
                binCenters[1] = 270 - delta;
            } else {
                int delta = vc - 180;
                binCenters[0] = 90 + delta;
                binCenters[1] = 270 + delta - 360;
            }

            return binCenters;
        }

        int[] binCenters = new int[n];

        for (int i = 0; i < n; ++i) {

            int idx = indexes.get(i);

            binCenters[i] = Math.round(hist.getXHist()[idx]);
        }

        return binCenters;
    }

    private int[] createBinCenters255(HistogramHolder hist, List<Integer> indexes) {

        if (indexes.isEmpty()) {
            return new int[0];
        }

        int n = indexes.size();

        if (n == 1) {

            int[] binCenters = new int[n + 1];

            int vc = Math.round(hist.getXHist()[indexes.get(0).intValue()]);

            binCenters[0] = vc/2;
            binCenters[1] = (255 + vc)/2;

            return binCenters;
        }

        int[] binCenters = new int[n];

        for (int i = 0; i < n; ++i) {

            int idx = indexes.get(i);

            binCenters[i] = Math.round(hist.getXHist()[idx]);
        }

        return binCenters;
    }

    private List<Set<PairInt>> groupByPeaksForGrey(ImageExt input,
        Map<Integer, Collection<PairInt>> greyPixelMap) {

        Set<PairInt> points = new HashSet<PairInt>();
        for (Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            points.addAll(entry.getValue());
        }

        Map<PairInt, Integer> polarCIEXYMap = calculatePolarCIEXY(input, points);
        float[] polarCIEXY = getValues(polarCIEXYMap);
        float binWidth = 1;
        HistogramHolder hist = Histogram.createSimpleHistogram(binWidth,
            polarCIEXY, Errors.populateYErrorsBySqrt(polarCIEXY));
        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist,
            0.1f);
        int[] binCenters = createBinCenters360(hist, indexes);

        List<Set<PairInt>> pixelGroups = assignToNearestPolarCIECluster(
            polarCIEXYMap, binCenters);

        return pixelGroups;

    }

    /**
     * a greyscale segmentation algorithm similar to the KMPP, but does not use
     * a random number generator to find seeds of intensity bins.
     * @param input
     */
    public void applyGreyscaleHistogram(GreyscaleImage input) {

        float[] values = new float[input.getNPixels()];
        for (int i = 0; i < input.getNPixels(); ++i) {
            int v = input.getValue(i);
            values[i] = v;
        }

        int[] binCenters = determineGreyscaleBinCenters(values);

        assignToNearestCluster(input, binCenters);
    }

    /**
     * a greyscale segmentation algorithm similar to the KMPP, but does not use
     * a random number generator to find seeds of intensity bins.
     *
     * @param values
     * @return
     */
    public int[] determineGreyscaleBinCenters(float[] values) {

        float binWidth = 3;
        HistogramHolder hist = Histogram.createSimpleHistogram(0, 255, binWidth,
            values, Errors.populateYErrorsBySqrt(values));
        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist, 0.05f);

        if (indexes.size() <= 3) {
            binWidth = 10;
            hist = Histogram.createSimpleHistogram(0, 255, binWidth,
                values, Errors.populateYErrorsBySqrt(values));
            indexes = MiscMath.findStrongPeakIndexesDescSort(hist, 0.05f);

            int k = 3;
            int idx = (indexes.size() < k) ? indexes.size() : k;
            indexes = indexes.subList(0, idx);
        }

        /*
        try {
            hist.plotHistogram("greyscale", "1");
        } catch (IOException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        }*/

        int[] binCenters = createBinCenters255(hist, indexes);

        return binCenters;
    }

    /**
     * segmentation algorithm using cieXY and rgb to make pixel categories,
     * and grow lists, then histograms to rescale.
     *
     * @param input
     * @return
     */
    public GreyscaleImage createGreyscale3(ImageExt input) {

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<PairInt, Float> colorPixelMap = new HashMap<PairInt, Float>();

        Set<PairInt> unassignedPixels = new HashSet<PairInt>();

        populatePixelLists2(input, blackPixels, colorPixelMap, unassignedPixels);

        int nTotC = blackPixels.size() + colorPixelMap.size() + unassignedPixels.size();
        assert(nTotC == input.getNPixels());

        // grow black pixels from unassigned pixels if within a tolerance of rgb
        growPixelsWithinRGBTolerance(input, blackPixels, unassignedPixels, 10);

        // then grow colorPixelMap from unassignedPixels if within a tolerance of rgb
        Set<PairInt> addToColor = findPixelsWithinRGBTolerance(input,
            colorPixelMap.keySet(), unassignedPixels, 10);//5

        // reassign colorPixelMap to averaged color for adjacent thetas within tolerance
        // (1) try just colorPixelMap for this step
        // (2) try adding addToColor to colorPixelMap for this step
        //    Map<PairInt, Float> colorPixelMap2 = createPolarCIEXYMap(input, addToColor);
        //    colorPixelMap.putAll(colorPixelMap2);

        int toleranceInValue = determineToleranceForIllumCorr(colorPixelMap, nTotC);

        if (toleranceInValue > -1 && !colorPixelMap.isEmpty()) {
            // use higher > 1 when bp/cp2 >> 1
            if ((blackPixels.size()/addToColor.size()) > 10) {
                correctForIllumination(input, 2, colorPixelMap);
            } else {
                correctForIllumination(input, toleranceInValue, colorPixelMap);
            }
        }
/*
ImageExt tmpInput = input.copyToImageExt();
for (PairInt p : unassignedPixels) {
tmpInput.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImage(tmpInput, "_after_illumc0_unassigned_pix" + MiscDebug.getCurrentTimeFormatted());
MiscDebug.writeImage(input, "_after_illumc0_" + MiscDebug.getCurrentTimeFormatted());
*/
        Map<PairInt, Float> colorPixelMap2 = createPolarCIEXYMap(input, addToColor);
        if (toleranceInValue > -1 && !colorPixelMap2.isEmpty()) {
            correctForIllumination(input, toleranceInValue, colorPixelMap2);
        }
/*
tmpInput = input.copyToImageExt();
for (PairInt p : colorPixelMap2.keySet()) {
tmpInput.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImage(tmpInput, "_after2_illumc0_pix" + MiscDebug.getCurrentTimeFormatted());
MiscDebug.writeImage(input, "_after2_illumc0_" + MiscDebug.getCurrentTimeFormatted());
*/
        Map<PairInt, Float> colorPixelMap3 = createPolarCIEXYMap(input, unassignedPixels);
        toleranceInValue = determineToleranceForIllumCorr(colorPixelMap3, nTotC);
        if (toleranceInValue > -1 && !colorPixelMap3.isEmpty()) {
            if ((toleranceInValue == 1) &&
                ((blackPixels.size()/colorPixelMap2.size()) > 10)) {
                if ((blackPixels.size()/colorPixelMap.size()) > 10) {
                    toleranceInValue = 10;
                } else {
                    toleranceInValue = 2;
                }
            }
            correctForIllumination(input, toleranceInValue, colorPixelMap3);
        }
/*tmpInput = input.copyToImageExt();
for (PairInt p : colorPixelMap3.keySet()) {
tmpInput.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImage(tmpInput, "_after3_illumc0_pix" + MiscDebug.getCurrentTimeFormatted());
MiscDebug.writeImage(input, "_after3_illumc0_" + MiscDebug.getCurrentTimeFormatted());
*/
        if (nTotC != blackPixels.size()) {
            //TODO: not sure will use this... might consider same method w/ rgb
            toleranceInValue = 2;
            correctForIllumination(input, blackPixels, toleranceInValue);
        }

        GreyscaleImage img = input.copyToGreyscale();

        MiscDebug.writeImage(img, "_before_greyscale_bins_" + MiscDebug.getCurrentTimeFormatted());

        float[] cValues = new float[input.getNPixels()];
        for (int i = 0; i < input.getNPixels(); ++i) {
            int v = img.getValue(i);
            cValues[i] = v;
        }

        int[] binCenters = determineGreyscaleBinCenters(cValues);

        assignToNearestCluster(img, binCenters);
        MiscDebug.writeImage(img, "_after_greyscale_bins_" + MiscDebug.getCurrentTimeFormatted());

        for (PairInt p : blackPixels) {
            img.setValue(p.getX(), p.getY(), 0);
        }

MiscDebug.writeImage(img, "_end_seg_" + MiscDebug.getCurrentTimeFormatted());

        return img;
    }

    public GreyscaleImage createCombinedWaveletBased(Image img) {
        return createCombinedWaveletBased(img.copyRedToGreyscale(),
            img.copyGreenToGreyscale(), img.copyBlueToGreyscale());
    }

    public GreyscaleImage createCombinedWaveletBased2(Image img) {
        return createCombinedWaveletBased2(img.copyRedToGreyscale(),
            img.copyGreenToGreyscale(), img.copyBlueToGreyscale());
    }

    public GreyscaleImage createCombinedWaveletBased(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg) {

        boolean use1D = true;
        GreyscaleImage rSegImg = createGreyscale5(rImg, use1D);
        GreyscaleImage gSegImg = createGreyscale5(gImg, use1D);
        GreyscaleImage bSegImg = createGreyscale5(bImg, use1D);

        GreyscaleImage combined = rSegImg.copyImage();
        for (int i = 0; i < rSegImg.getWidth(); ++i) {
            for (int j = 0; j < rSegImg.getHeight(); ++j) {
                int g = gSegImg.getValue(i, j);
                int b = bSegImg.getValue(i, j);
                if (g > 0) {
                    combined.setValue(i, j, g);
                } else if (b > 0) {
                    combined.setValue(i, j, b);
                }
            }
        }
        return combined;
    }

    public GreyscaleImage createCombinedWaveletBased2(
        GreyscaleImage rImg, GreyscaleImage gImg, GreyscaleImage bImg) {

        ATrousWaveletTransform wt = new ATrousWaveletTransform();

        GreyscaleImage coarsestCoeffR = null;
        GreyscaleImage coarsestCoeffG = null;
        GreyscaleImage coarsestCoeffB = null;

        for (int i = 0; i < 3; ++i) {
            GreyscaleImage input;
            if (i == 0) {
                input = rImg;
            } else if (i == 1) {
                input = gImg;
            } else {
                input = bImg;
            }
            List<GreyscaleImage> transformed = new ArrayList<GreyscaleImage>();
            List<GreyscaleImage> coeffs = new ArrayList<GreyscaleImage>();
            wt.calculateWithB3SplineScalingFunction(input, transformed, coeffs);
            if (i == 0) {
                coarsestCoeffR = coeffs.get(coeffs.size() - 1);
            } else if (i == 1) {
                coarsestCoeffG = coeffs.get(coeffs.size() - 1);
            } else {
                coarsestCoeffB = coeffs.get(coeffs.size() - 1);
            }
        }

        //TODO: determine top limit by frequency distr?
        int limit = 3;

        Stack<Integer> stack = new Stack<Integer>();

        GreyscaleImage coarsestCombined = coarsestCoeffB.createWithDimensions();

        for (int i = 0; i < coarsestCoeffR.getNPixels(); ++i) {
            int r = coarsestCoeffR.getValue(i);
            int g = coarsestCoeffG.getValue(i);
            int b = coarsestCoeffB.getValue(i);
            if (r > limit || g > limit || b > limit) {
                coarsestCombined.setValue(i, 250);
                stack.add(Integer.valueOf(i));
            }
        }

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        int w = coarsestCombined.getWidth();
        int h = coarsestCombined.getHeight();

        int lowerLimit = 0;
        while (limit > lowerLimit) {
            // use the canny edge 2-layer approach to pick up neighboring pixels
            // at or above a lower limit

            Set<Integer> visited = new HashSet<Integer>();
            limit--;
            while (!stack.isEmpty()) {
                Integer pixIndex = stack.pop();
                if (visited.contains(pixIndex)) {
                    continue;
                }
                int x = coarsestCombined.getCol(pixIndex.intValue());
                int y = coarsestCombined.getRow(pixIndex.intValue());
                for (int i = 0; i < dxs.length; ++i) {
                    int x2 = x + dxs[i];
                    int y2 = y + dys[i];
                    if (x2 < 0 || y2 < 0 || x2 > (w - 1) || y2 > (h - 1)) {
                        continue;
                    }
                    int r = coarsestCoeffR.getValue(x2, y2);
                    int g = coarsestCoeffG.getValue(x2, y2);
                    int b = coarsestCoeffB.getValue(x2, y2);
                    if (r > limit || g > limit || b > limit) {
                        coarsestCombined.setValue(i, 250);
                        stack.add(Integer.valueOf(i));
                    }
                }
                visited.add(pixIndex);
            }
        }

        return coarsestCombined;
    }

    /**
     * segmentation algorithm using an a trous wavelet transform.
     *
     * @param input
     * @return the segmented image holding values 0 or 250.
     */
    public GreyscaleImage createGreyscale5(GreyscaleImage input) {

        boolean use1D = false;
        return createGreyscale5(input, use1D);
    }

    /**
     * segmentation algorithm using an a trous wavelet transform.
     *
     * @param input
     * @return the segmented image holding values 0 or 250.
     */
    public GreyscaleImage createGreyscale5(GreyscaleImage input, boolean use1D) {

        ATrousWaveletTransform wt = new ATrousWaveletTransform();

        List<GreyscaleImage> transformed = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs = new ArrayList<GreyscaleImage>();

        if (use1D) {
            wt.calculateWithB3SplineScalingFunction(input, transformed, coeffs);
        } else {
            wt.calculateWithB3SplineScalingFunction2(input, transformed, coeffs);
        }

        /*
            for (int i = 0; i < coeffs.size(); ++i) {
                GreyscaleImage img = coeffs.get(i);
                String str = "coeff_" + Integer.toString(i) + "_" +
                    MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(img, str);
            }
        */

        GreyscaleImage coarsestCoeff = coeffs.get(coeffs.size() - 1);

        int limit = use1D ? 2 : 1;

        for (int i = 0; i < coarsestCoeff.getNPixels(); ++i) {
            if (coarsestCoeff.getValue(i) > limit) {
                coarsestCoeff.setValue(i, 250);
            } else {
                coarsestCoeff.setValue(i, 0);
            }
        }

        return coarsestCoeff;
    }

    /**
     * segmentation algorithm using Canny Edges.
     *
     * @param input
     * @return the segmented image holding values 0 or 250.
     */
    public GreyscaleImage createGreyscale6(GreyscaleImage input) {

        CannyEdgeFilterAdaptive filter = new CannyEdgeFilterAdaptive();
        filter.applyFilter(input.copyImage());

        OtsuThresholding ot = new OtsuThresholding();
        float otsuScaleFactor = 0.65f;//0.4f;
        double tHigh = otsuScaleFactor * ot.calculateBinaryThreshold256(
            filter.getFilterProducts().getGradientXY());
        
        double thresh = tHigh;

        GreyscaleImage img = filter.getFilterProducts().getGradientXY();

        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > thresh) {
                img.setValue(i, 250);
            } else {
                img.setValue(i, 0);
            }
        }

        return img;
    }

    /**
     * segmentation algorithm using cieXY and rgb to make segmentation for
     * the colors with CIE X or Y outside of the center region
     *
     * @param input
     * @return
     */
    public GreyscaleImage createGreyscale7(ImageExt input) {

        Map<PairInt, Integer> pixelThetaDegreesMap = populatePixelLists3(input);

        // -- scale the valus to between 0 and 255 w/ wrap around --

        int count = 0;
        float[] clrPolarCIEXY =new float[pixelThetaDegreesMap.size()];
        for (Entry<PairInt, Integer> entry : pixelThetaDegreesMap.entrySet()) {
            clrPolarCIEXY[count] = entry.getValue();
            count++;
        }

        float binWidth = 20;
        HistogramHolder hist = Histogram.createSimpleHistogram(binWidth,
            clrPolarCIEXY, Errors.populateYErrorsBySqrt(clrPolarCIEXY));
        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist, 0.1f);
        int[] binCenters = createBinCenters360(hist, indexes);

        List<Set<PairInt>> colorPixelGroups = assignToNearestPolarCIECluster(
            pixelThetaDegreesMap, binCenters);

        GreyscaleImage img = new GreyscaleImage(input.getWidth(),
            input.getHeight());

        int gClr = 255;
        int s = 127/colorPixelGroups.size();
        for (Set<PairInt> set : colorPixelGroups) {
            for (PairInt p : set) {
               img.setValue(p.getX(), p.getY(), gClr);
            }
            gClr -= s;
        }

MiscDebug.writeImage(img, "_seg_gs7_" + MiscDebug.getCurrentTimeFormatted());

        return img;
    }

    public GreyscaleImage createAWatershed(ImageExt input, String debugTag,
        int originalImageWidth, int originalImageHeight) {

        int w = input.getWidth();
        int h = input.getHeight();
        GreyscaleImage aImg = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);

        for (int i = 0; i < input.getNPixels(); ++i) {

            float[] lab = input.getCIELAB(i);

            aImg.setValue(i, Math.round(lab[1]));

        }

        ImageProcessor imageProcessor = new ImageProcessor();

        HistogramEqualization hEq = new HistogramEqualization(aImg);
        hEq.applyFilter();

        int minDimension = Math.min(originalImageWidth, originalImageHeight);
        int lowerLimitSize;
        if (minDimension > 900) {
            lowerLimitSize = 300;
        } else if (minDimension < 200) {
            lowerLimitSize = 100;
        } else {
            lowerLimitSize = 200;
        }

        imageProcessor.applyAdaptiveMeanThresholding(aImg, 1);

        GreyscaleImage ws = imageProcessor.makeWatershedFromAdaptiveMedian(aImg);

        return ws;
    }
    
    /**
     * given the list of edges, populate the output arrays with color informaion
     * from the edge points and their 8 neighbor regions.
     * note that any points in more than one output list originally because of
     * being a junction, are corrected and placed in the most similar list.
     * 
     * @param img
     * @param edges
     * @param junctions
     * @param outputPoints
     * @param outputDescripors access as [edgeListIndex][(h, s, v, nPix, cenX, cenY)]
     * @param clrSpace color space to fill the descriptors with: 0 is lab, 1 is hsv 
     */
    private void populateEdgeLists(ImageExt img,
        List<PairIntArray> edges, 
        List<Set<PairInt>> outputPoints, float[][] outputDescripors,
        int clrSpace) {
               
        int w = img.getWidth();
        int h = img.getHeight();
                
        // ----- gather edge points and their 8 neighbors into edge point sets ----
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        for (int i = 0; i < edges.size(); ++i) {
            PairIntArray edge = edges.get(i);
            Set<PairInt> set = new HashSet<PairInt>();
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                set.add(new PairInt(x, y));
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    set.add(p2);
                }
            }
            outputPoints.add(set);
        }
        
        float n = img.getNPixels();
        
        // ----- calculate descriptors of the color and location of the edge points -----
        
        populateDescriptors(img, outputPoints, outputDescripors, clrSpace);
        
        /*
        the descriptors are
             C_i = {h, s, v, nPix, cenX, cenY}  or labL, labA, labB for colorSpace = 0
        */
        
        // ======= correct for any points in more than one list ======
        
        // --- map the list indexes that a point is in --------
        Map<PairInt, Set<Integer>> pointIndexes = new HashMap<PairInt, Set<Integer>>();
        
        for (int i = 0; i < outputPoints.size(); ++i) {
            Integer key = Integer.valueOf(i);
            Set<PairInt> edgePoints = outputPoints.get(i);
            for (PairInt p : edgePoints) {
                Set<Integer> indexes = pointIndexes.get(p);
                if (indexes == null) {
                    indexes = new HashSet<Integer>();
                    pointIndexes.put(p, indexes);
                }
                indexes.add(key);
            }
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        // ---- when a point is in more than one list, choose to keep it in the
        //      list with smallest difference from it in color and remove it from others.
        
        for (Entry<PairInt, Set<Integer>> entry : pointIndexes.entrySet()) {
            Set<Integer> indexes = entry.getValue();
            if (indexes.size() == 1) {
                continue;
            }
            PairInt p = entry.getKey();
            int x = p.getX();
            int y = p.getY();
            float c1, c2, c3;
            if (clrSpace == 0) {
                float[] lab = img.getCIELAB(x, y);
                c1 = lab[0];
                c2 = lab[1];
                c3 = lab[2];
            } else {
                // hsv
                c1 = img.getHue(x, y);
                c2 = img.getSaturation(x, y);
                c3 = img.getBrightness(x, y);
            }
            double minColorDiff = Double.MAX_VALUE;
            Integer minColorDiffIndex = null;
            for (Integer index : indexes) {
                //C_i = {h, s, v, nPix, cenX, cenY}
                float[] desc = outputDescripors[index.intValue()];
                double diff;
                if (clrSpace == 0) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(c1, c2, c3, 
                        desc[0], desc[1], desc[2]));
                } else {
                    double diff1 = c1 - desc[0];
                    double diff2 = c2 - desc[1];
                    double diff3 = c3 - desc[2];
                    diff = Math.sqrt(diff1 * diff1 + diff2*diff2 + diff3*diff3);
                }
                if (diff < minColorDiff) {
                    minColorDiff = diff;
                    minColorDiffIndex = index;
                }
            }
            assert(minColorDiffIndex != null);
            
            // update the lists to remove point
            for (Integer index : indexes) {
                if (index.equals(minColorDiffIndex)) {
                    continue;
                }
                
                Set<PairInt> set = outputPoints.get(index.intValue());
                float nBefore = set.size();
                set.remove(p);
                float nAfter = set.size();
                
                ////C_i = {h, s, v, nPix, cenX, cenY}
                float[] desc = outputDescripors[index.intValue()];
                desc[0] = ((desc[0] * nBefore) - c1)/nAfter;
                desc[1] = ((desc[1] * nBefore) - c2)/nAfter;
                desc[2] = ((desc[2] * nBefore) - c3)/nAfter;
                desc[3] = nAfter;
                desc[4] = ((desc[4] * nBefore) - x)/nAfter;
                desc[5] = ((desc[5] * nBefore) - y)/nAfter;
            }
        }
    }

    /**
     * create segmented image by creating edges with phase congruence,
     * then using the edge color properties to form clusters and seeds
     * of regions to grow, then using color histograms to further merge
     * regions.
     * The algorithm follows the general outline given by
     * Jie and Peng-fei 2003, "Natural Color Image Segmentation",
       http://www-labs.iro.umontreal.ca/~mignotte/IFT6150/Articles/TRASH/ARTICLES_2010/cr1231.pdf
     * @param input
     * @return 
     */
    public List<Set<PairInt>> createColorEdgeSegmentation(ImageExt input) {

        GreyscaleImage gsImg = input.copyBlueToGreyscale();

        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 5;//10;//5;//2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.0001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = false;//true;
        
        final int w = input.getWidth();
        final int h = input.getHeight();
        final int nPix = input.getNPixels();
           
        PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseDetector.phaseCongMono(gsImg, nScale, minWavelength, mult, 
            sigmaOnf, k, increaseKIfNeeded,
            cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
        
        int[][] thinned = products.getThinned();
        {
            GreyscaleImage out2 = gsImg.createWithDimensions();
            for (int i = 0; i < thinned.length; ++i) {
                for (int j = 0; j < thinned[i].length; ++j) {
                    if (thinned[i][j] > 0) {
                        out2.setValue(j, i, 255);
                    }
                }
            }
            MiscDebug.writeImage(out2, "_EDGES_grey_");
        }
        //TODO: copy and used data from phaseDetector to release memory
        phaseDetector = null;
                
        EdgeExtractorSimple extractor = new EdgeExtractorSimple(thinned);
        extractor.extractEdges();
        List<PairIntArray> edges = extractor.getEdges();
        Set<PairInt> junctions = extractor.getJunctions();
        // put in framework of images
        for (int i = 0; i < edges.size(); ++i) {
            PairIntArray edge = edges.get(i);
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                edge.set(j, y, x);
            }
        }
        {
            Set<PairInt> tmp = new HashSet<PairInt>();
            for (PairInt p : junctions) {
                tmp.add(new PairInt(p.getY(), p.getX()));
            }
            junctions.clear();
            junctions.addAll(tmp);
        }
        
        //TODO: copy any used data from extractor to release memory
        extractor = null;
        
        int nEdges = edges.size();
        List<Set<PairInt>> clusterPoints = new ArrayList<Set<PairInt>>();
        float[][] clusterDescriptors = new float[nEdges][];
        
        // 0 for LAB, 1 for HSV
        int clrSpace = 0;
        
        populateEdgeLists(input, edges, clusterPoints, clusterDescriptors, 
            clrSpace);        
        
        int tLen = 20;
        List<Integer> longEdgeIndexes = new ArrayList<Integer>();
        List<Integer> shortEdgeIndexes = new ArrayList<Integer>();
        populateEdgeLengthLists(clusterDescriptors, tLen, 
            longEdgeIndexes, shortEdgeIndexes); 

        // merge edges in the heap for pairs with diff < tColot
        // NOTE that the moved sets modify the data structures :
        //    clusterPoints may contain empty items
        //    clusterDescriptors may contain null items
        //    both clusterPoints and clusterDescriptor non- null and non empty 
        //       items are updated for merges
        double tColor;
        if (clrSpace == 0) {
            // JND for deltaE is ~2.3
            tColor = 4.;//5.5;
        } else {
            // what is JND for HSV (a.k.a. HSB) ?
            tColor = Double.MAX_VALUE;
        }
        
        mergeEdges(clusterPoints, clusterDescriptors, clrSpace, tColor,
            longEdgeIndexes);
  
        {
            // DEBUG
            List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();
            for (Integer index : longEdgeIndexes) {
                Set<PairInt> set = clusterPoints.get(index.intValue());
                if (!set.isEmpty()) {
                    tmp.add(set);
                }
            }
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(tmp, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_longEdges_merged_" +  clrSpace);            
        }
        
        {
            // DEBUG
            List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();
            for (Integer index : shortEdgeIndexes) {
                Set<PairInt> set = clusterPoints.get(index.intValue());
                if (!set.isEmpty()) {
                    tmp.add(set);
                }
            }
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(tmp, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_shortedges_before_merge_" +  clrSpace);            
        }
        
        // ---- merge short edges (which are usually textures) ------
        mergeShortEdges(clusterPoints, clusterDescriptors, clrSpace, tColor,
            shortEdgeIndexes);
         
        {
            // DEBUG
            List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();
            for (Integer index : shortEdgeIndexes) {
                Set<PairInt> set = clusterPoints.get(index.intValue());
                if (!set.isEmpty()) {
                    tmp.add(set);
                }
            }
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(tmp, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_shortedges_merged_" +  clrSpace);            
        }

        // the number of clusters in long edges could be used in kmeans++ 

        assert(assertDescriptorCounts(clusterPoints, clusterDescriptors));  
        
        // ------ region growing -------
        growEdges(input, clusterPoints, clusterDescriptors, clrSpace, tColor,
            shortEdgeIndexes, longEdgeIndexes, tLen);
        
        // cluster descriptors not updated for unassigned pixel additions at this point
        
        {
            // DEBUG
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(clusterPoints, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_after_rgo_" +  clrSpace);            
        }
        
        // -- release some of the datastructures and condense the clusters ---
        clusterDescriptors = null;
        longEdgeIndexes = null;
        shortEdgeIndexes = null;
        
        for (int i = (clusterPoints.size() - 1); i > -1; --i) {
            Set<PairInt> set = clusterPoints.get(i);
            if (set.isEmpty()) {
                clusterPoints.remove(i);
            }
        }
        
        // ------ merge by color histograms ------
        
        clrSpace = 1;
        
        double tR;
        if (clrSpace == 0) {
            tR = 0.85*3.0;
        } else if (clrSpace == 1) {
            tR = 0.6*3.0;
        } else {
            // not tested yet
            tR = 0.6*3.0;
        }
             
        mergeByColorHistograms(input, clusterPoints, clrSpace, tR);
        
        {
            // DEBUG
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(clusterPoints, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_FINEAL_" +  clrSpace);            
        }
        
        /*
        TODO:
        ** for the clusters which are smaller than the limit tNumber,
            either merge them ith the closest cluster.        
        */
        
        int tNumber = Math.round(0.01f * nPix);
        
        throw new UnsupportedOperationException("not yet implemented");
        
    }

    /**
    if pair distance < tColor
                merge the pair and update the affected data structures.
      NOTE that the moved sets modify the dsta structures :
      clusterPoints may contain empty items, 
      clusterDescriptors may contain null items,
      both clusterPoints and clusterDescriptor non-null and non-empty 
      items are updated for merges.
      * This method uses a queue and bfs traversal pattern instead of a 
      * min fibonacci heap, so is faster than mergeEdges(...) but has the
      * possibility of a cluster's average properties wandering from
      * center of approach in mergeEdges(...).
      * 
     * @param clusterPoints
     * @param clusterDescriptors
     * @param clrSpace
     * @param tColor
     * @param longEdgesHeap
     * @param heapKeyFactor
     * @param pairEdgePindexNodes 
     */
    private void mergeEdges2(List<Set<PairInt>> clusterPoints, 
        float[][] clusterDescriptors, int clrSpace, double tColor, 
        Heap longEdgesHeap, long heapKeyFactor, 
        Map<PairInt, HeapNode> pairEdgePindexNodes) {
        
        // ---- make a map to finde and update merged data structures ------
        Map<Integer, Set<Integer>> indexToIndexMap = new HashMap<Integer, Set<Integer>>();
        for (PairInt p : pairEdgePindexNodes.keySet()) {
            
            Integer index1 = Integer.valueOf(p.getX());
            Integer index2 = Integer.valueOf(p.getY());
            
            Set<Integer> indexes = indexToIndexMap.get(index1);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                indexToIndexMap.put(index1, indexes);
            }
            indexes.add(index2);
            
            indexes = indexToIndexMap.get(index2);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                indexToIndexMap.put(index2, indexes);
            }
            indexes.add(index1);
        }
        
        
        /*
        TODO: make an alternative method which uses a queue to store the
        pair difference indexes.
           the queue is initialized with the pairs ordered by increasing 
           color difference, but thereafter is populated using a bfs
           style traversal (indexes that are merged into are added to the
           queue and visited after the smallest within threshold in the 
           initial queue are visited).
        The queue should be faster because each update would be O(1) 
            rather than O(lg2(N_pairs)) and removals are O(1) instead of
            O(lg2(N_pairs)).
        Note that the min heap pattern in this method should lead to more stable
        average properties, but not sure that the results would not be
        similar for most images.
        */
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        while (!longEdgesHeap.isEmpty()) {
            
            HeapNode node = longEdgesHeap.extractMin();
            double diff = ((double)node.getKey())/((double)heapKeyFactor);
            
            if (diff > tColor) {
                break;
            }
            
            PairInt p12 = (PairInt)node.getData();
            
            int idx1 = p12.getX();
            int idx2 = p12.getY();
            
            Set<PairInt> set1 = clusterPoints.get(idx1);
            Set<PairInt> set2 = clusterPoints.get(idx2);
            
            if (set1.isEmpty() || set2.isEmpty()) {
                continue;
            }
            
            if (set2.size() > set1.size()) {
                idx1 = p12.getY();
                idx2 = p12.getX();
                set1 = set2;
                set2 = clusterPoints.get(idx2);
            }
            
            // set1 is largest
            
            float[] desc1 = clusterDescriptors[idx1];
            float[] desc2 = clusterDescriptors[idx2];
            float n1 = set1.size();
            float n2 = set2.size();
            float nTot = n1 + n2;
            
            //{h, s, v, nPix, cenX, cenY}
            // update desc1 contents for contents in desc2
            for (int k = 0; k < desc1.length; ++k) {
                if (k == 3) {
                    desc1[k] = nTot;
                } else {
                    desc1[k] = ((desc1[k] * n1) + (desc2[k] * n2)) / nTot;
                }
            }
            clusterDescriptors[idx2] = null;
            set1.addAll(set2);
            set2.clear();
            
            n1 = set1.size();
            
            // update all pairs having set 1
            Set<Integer> indexes2 = indexToIndexMap.get(Integer.valueOf(idx1));
            for (Integer index2 : indexes2) {
                int idx3 = index2.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);            
                if (set3.isEmpty()) {
                    continue;
                }
                float[] desc3 = clusterDescriptors[idx3];
                
                //keys in pairEdgePindexNodes have smaller index in x
                PairInt p13;
                if (idx1 < idx3) {
                    p13 = new PairInt(idx1, idx3);
                } else {
                    p13 = new PairInt(idx3, idx1);
                }
                HeapNode node3 = pairEdgePindexNodes.get(p13);
                assert(node3 != null);
                
                longEdgesHeap.remove(node3);
                
                if (clrSpace == 0) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(
                        desc1[0], desc1[1], desc1[2], 
                        desc3[0], desc3[1], desc3[2]));
                } else {
                    double diff1 = desc1[0] - desc3[0];
                    double diff2 = desc1[1] - desc3[1];
                    double diff3 = desc1[2] - desc3[2];
                    diff = Math.sqrt(diff1 * diff1 + diff2*diff2 + diff3*diff3);
                }
                
                long heapKey = (long)(diff * (double)heapKeyFactor);
                node3 = new HeapNode(heapKey);
                node3.setData(p13);
                longEdgesHeap.insert(node3);
                pairEdgePindexNodes.put(p13, node3);                
            }
            
            // pairs having set2 will be skipped because of the empty set at beginning of while loop
        }        
    }
    
    private void populateEdgeLengthLists(float[][] clusterDescriptors, 
        int tLen, List<Integer> longEdgeIndexes, 
        List<Integer> shortEdgeIndexes) {
        
        //C_i = {h, s, v, nPix, cenX, cenY}
        for (int i = 0; i < clusterDescriptors.length; ++i) {
            Integer key = Integer.valueOf(i);
            float nPix = clusterDescriptors[i][3];
            if (nPix < tLen) {
                shortEdgeIndexes.add(key);
            } else {
                longEdgeIndexes.add(key);
            }
        }
    }

    private void populateColorDiffHeap(
        float[][] clusterDescriptors, int clrSpace, 
        List<Integer> longEdgeIndexes, Heap longEdgesHeap, 
        long heapKeyFactor, Map<PairInt, HeapNode> pairEdgePindexNodes) {
        
        // for heap nodes:
        //     key is the difference in color times a factor to use long instead of double
        //     data is the PairInt holding the indexes compared
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < longEdgeIndexes.size(); ++i) {
            
            int idx1 = longEdgeIndexes.get(i).intValue();
            float[] desc1 = clusterDescriptors[idx1];
            
            for (int j = (i + 1); j < longEdgeIndexes.size(); ++j) {
                
                int idx2 = longEdgeIndexes.get(j).intValue();
                float[] desc2 = clusterDescriptors[idx2];
                                
                double diff;
                if (clrSpace == 0) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(
                        desc1[0], desc1[1], desc1[2], 
                        desc2[0], desc2[1], desc2[2]));
                } else {
                    double diff1 = desc1[0] - desc2[0];
                    double diff2 = desc1[1] - desc2[1];
                    double diff3 = desc1[2] - desc2[2];
                    diff = Math.sqrt(diff1 * diff1 + diff2*diff2 + diff3*diff3);
                }
                
                // note that idx1 is always smaller than idx2
                PairInt p12 = new PairInt(idx1, idx2);
                
                long heapKey = (long)((double)heapKeyFactor * diff);
                HeapNode node = new HeapNode(heapKey);
                node.setData(p12);
                        
                longEdgesHeap.insert(node);
                
                pairEdgePindexNodes.put(p12, node);
            }
        }
    }
    
    /**
     *
     * @param input
     * @param debugLabel if null, no debug output is made, else output uses debugLabel
     * as suffix in file names
     */
    public void extractObjectEdges(ImageExt input, String debugLabel,
        int originalImageWidth, int originalImageHeight) {

        int w = input.getWidth();
        int h = input.getHeight();

        GreyscaleImage o1 = new GreyscaleImage(w, h,
            GreyscaleImage.Type.Bits32FullRangeInt);
        GreyscaleImage o2 = o1.createFullRangeIntWithDimensions();
        GreyscaleImage o3 = o1.createFullRangeIntWithDimensions();
        GreyscaleImage aImg = o1.createFullRangeIntWithDimensions();
        GreyscaleImage bImg = o1.createFullRangeIntWithDimensions();
        GreyscaleImage hueAngleImg = o1.createFullRangeIntWithDimensions();
        GreyscaleImage cieXYAngleImg = o1.createFullRangeIntWithDimensions();

        for (int i = 0; i < input.getNPixels(); ++i) {
            int r = input.getR(i);
            int g = input.getG(i);
            int b = input.getB(i);
            float[] lab = input.getCIELAB(i);
            o1.setValue(i, (int)Math.round((double)(r - g)/Math.sqrt(2)));
            o2.setValue(i, (int)Math.round((double)(r + g - 2*b)/Math.sqrt(6)));
            o3.setValue(i, (int)Math.round((double)(r + g + b)/Math.sqrt(2)));
            aImg.setValue(i, Math.round(lab[1]));
            bImg.setValue(i, Math.round(lab[2]));

            float ha;
            if (lab[1] == 0) {
                ha = 0;
            } else {
                ha = (float)(Math.atan2(lab[2], lab[1]) * 180. / Math.PI);
                if (ha < 0) {
                    ha += 360.;
                }
            }

            hueAngleImg.setValue(i, Math.round(ha));

            //TODO: replace w/ cached method
            float[] cieXY = input.getCIEXY_(i);

            float cieXYAngle = (float)(Math.atan2(cieXY[1], cieXY[0]) * 180. / Math.PI);
            if (cieXYAngle < 0) {
                cieXYAngle += 360.;
            }
            cieXYAngleImg.setValue(i, Math.round(cieXYAngle));
        }

        HistogramEqualization hEq = new HistogramEqualization(aImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(bImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(hueAngleImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(cieXYAngleImg);
        hEq.applyFilter();
        hEq = new HistogramEqualization(o1);
        hEq.applyFilter();
        hEq = new HistogramEqualization(o1);
        hEq.applyFilter();
        hEq = new HistogramEqualization(o3);
        hEq.applyFilter();
        if (debugLabel != null && !debugLabel.equals("")) {
            MiscDebug.writeImage(o1, "_o1_" + debugLabel);
            MiscDebug.writeImage(o2, "_o2_" + debugLabel);
            MiscDebug.writeImage(o3, "_o3_" + debugLabel);
            MiscDebug.writeImage(aImg, "_a_" + debugLabel);
            MiscDebug.writeImage(bImg, "_b_" + debugLabel);
            MiscDebug.writeImage(hueAngleImg, "_ha_" + debugLabel);
            MiscDebug.writeImage(cieXYAngleImg, "_ciexya_" + debugLabel);
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        imageProcessor.applyAdaptiveMeanThresholding(o1, 1);
        imageProcessor.applyAdaptiveMeanThresholding(o2, 1);
        imageProcessor.applyAdaptiveMeanThresholding(o3, 1);

        if (debugLabel != null && !debugLabel.equals("")) {
            MiscDebug.writeImage(o1, "_o1_adaptive_median_" + debugLabel);
            MiscDebug.writeImage(o2, "_o2_adaptive_median_" + debugLabel);
            MiscDebug.writeImage(o3, "_o3_adaptive_median_" + debugLabel);
        }

        /*
        imageProcessor.applyAdaptiveMeanThresholding(aImg, 1);
        if (debugLabel != null && !debugLabel.equals("")) {
            MiscDebug.writeImage(aImg, "_a_adaptive_median_" + debugLabel);
        }*/

        //TODO: revise for minimum size of contiguous pixels.
        // it should be dependent upon image resolution, that is PSF and the focal distance of objects
        // (the number of beams, that is psf diameters, across the object),
        // but the number of pixels as image size is all the information available.
        int minDimension = Math.min(originalImageWidth, originalImageHeight);
        int lowerLimitSize;
        if (minDimension > 900) {
            lowerLimitSize = 300;
        } else if (minDimension < 200) {
            lowerLimitSize = 100;
        } else {
            lowerLimitSize = 200;
        }

        //List<Set<PairInt>> maskList = imageProcessor.extractConnectedComponents(
        //    labelled, lowerLimitSize);

        imageProcessor.applyAdaptiveMeanThresholding(aImg, 1);

        GreyscaleImage aWSImg = imageProcessor.makeWatershedFromAdaptiveMedian(
            aImg);

        CannyEdgeFilterAdaptive filter = new CannyEdgeFilterAdaptive();
        filter.applyFilter(input.copyToGreyscale());
        EdgeFilterProducts ep = filter.getFilterProducts();
        MiscDebug.writeImage(ep.getGradientXY(), "_edges_" + debugLabel);
    }

    protected Map<PairInt, Float> createPolarCIEXYMap(ImageExt input,
        Set<PairInt> points) {

        Map<PairInt, Float> thetaMap = new HashMap<PairInt, Float>();

        //TODO: remove imgCp when finished debugging
        //ImageExt imgCp = input.copyToImageExt();
        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        double thetaFactor = 255./360.;

        // to set the non member colors, need to traverse all pixels.
        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; row++) {
                PairInt p = new PairInt(col, row);
                if (points.contains(p)) {
                    float cieX = input.getCIEX(col, row);
                    float cieY = input.getCIEY(col, row);
                    double thetaDegrees = cieC.calculateXYTheta(cieX, cieY) * 180. / Math.PI;
                    double theta = thetaFactor * thetaDegrees;
          //          int thetaCIEXY = (int)Math.round(theta);
                    //          imgCp.setRGB(col, row, thetaCIEXY, thetaCIEXY, thetaCIEXY);

                    thetaMap.put(p, Float.valueOf((float) theta));
                } else {
                    //         imgCp.setRGB(col, row, 0, 255, 0);
                }
            }
        }

        //MiscDebug.writeImage(imgCp, "polarciexy_" + MiscDebug.getCurrentTimeFormatted());

        return thetaMap;
    }

    private void correctForIllumination(ImageExt input, Set<PairInt> points,
        int toleranceInValue) {

        if (points.isEmpty() || (toleranceInValue == 0)) {
            return;
        }

        Map<PairInt, Float> thetaMap = createPolarCIEXYMap(input, points);

        correctForIllumination(input, toleranceInValue, thetaMap);
    }

    private void correctForIllumination(ImageExt input,
        int toleranceInValue, Map<PairInt, Float> thetaMap) {

        /*
        in CIE XY color space, reassign colors:

        for points:
        -- convert to polar theta in CIE XY space and plot in image
           as greyscale without changes.

        -- DFS traversal through the theta values of pixels in points
           to make contigous groups of pixels that are within a tolerance of
           theta value of one another.

        -- reassign the average rgb color to those pixels in a group
        */

        int w = input.getWidth();
        int h = input.getHeight();

        // find pixels similar in color that are contiguous

        DFSConnectedGroupsFinder2 groupFinder = new DFSConnectedGroupsFinder2();
        groupFinder.findConnectedPointGroups(thetaMap, 360, toleranceInValue,
            w, h, false);

        int nGroups = groupFinder.getNumberOfGroups();

        // calc avg rgb of each contiguous group and reset rgb of all in group
        // to the avg

        for (int i = 0; i < nGroups; ++i) {

            Set<PairInt> group = groupFinder.getXY(i);

            if (group.size() == 1) {
                continue;
            }

            float sumR = 0;
            float sumG = 0;
            float sumB = 0;
            for (PairInt p : group) {
                int pixIdx = input.getInternalIndex(p.getX(), p.getY());
                sumR += input.getR(pixIdx);
                sumG += input.getG(pixIdx);
                sumB += input.getB(pixIdx);
            }
            int avgR = Math.round(sumR/(float)group.size());
            int avgG = Math.round(sumG/(float)group.size());
            int avgB = Math.round(sumB/(float)group.size());

            for (PairInt p : group) {
                input.setRGB(p.getX(), p.getY(), avgR, avgG, avgB);
            }
        }
    }

    private int determineToleranceForIllumCorr(Map<PairInt, Float> thetaMap,
        int nTotCategoryMembers) {

        if (thetaMap.isEmpty()) {
            return -1;
        }

        // using a histogram of color in
        // determining whether toleranceInValue should be 1 or 4.
        // 1 is better for low contrast regions like the brown & lowe 2003 mtns
        // but 4 is better for many colors, large contrast, and large dynamic
        // range such as seen in the merton college test images

        float[] values = new float[thetaMap.size()];
        int count = 0;
        for (Entry<PairInt, Float> entry : thetaMap.entrySet()) {
            values[count] = entry.getValue().floatValue();
            count++;
        }

        int tolerance = 1;

        HistogramHolder hist = Histogram.createSimpleHistogram(values,
            Errors.populateYErrorsBySqrt(values));
        try {
            hist.plotHistogram("cieTheta", "_cietheta_");
        } catch (IOException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (hist != null) {
            List<Integer> indexes = MiscMath.findStrongPeakIndexes(hist, 0.05f);

            //TODO: for rare case that image has only one color, and the pixels
            // are all in thetaMap, we want the invoker not to remove illumination
            // because that may be the only contrast in the image.
            // an example of that happening wouldbe a blue sky and blue mountains
            // very similar in color.
            if ((indexes.size() == 1) && (nTotCategoryMembers == thetaMap.size())) {
                return -1;
            }

            // from lowest to highest peak, are there continuous values >= 1/4 lowest peak?
            if (indexes.size() > 1) {

                Collections.sort(indexes);

                int minPeak = Integer.MAX_VALUE;
                for (Integer index : indexes) {
                    int v = hist.getYHist()[index.intValue()];
                    if (v < minPeak) {
                        minPeak = v;
                    }
                }

                int nBelow = 0;
                int n = 0;
                int limit = minPeak/4;
                for (int i = indexes.get(0).intValue();
                    i < indexes.get(indexes.size() - 1).intValue(); ++i) {
                    int v = hist.getYHist()[i];
                    if (v < limit) {
                        nBelow++;
                    }
                    n++;
                }

                if (((float)nBelow/(float)n) < 0.5) {
                    tolerance = 3;
                }
            }
        }

        return tolerance;
    }

    private int correctPointsBoundaries(ImageExt input,
        Set<PairInt> points, Set<PairInt> blackPixels,
        Set<PairInt> whitePixels, Set<PairInt> greyPixels) {

         // boundary corrections:
        // iterate over points0 to find boundary points and their neighbors.
        // a boundary point has at least one non-points0 point neighbor.
        // for the boundary points and each of their neigbhors,
        // determine whether should move the point to another group
        // (or add it to its own).
        // the local rgb or ciexy can be computed and kept in a hash to avoid
        // re-calculating.
        // the local color will be the 9 point region members only regions
        //

        int imageWidth = input.getWidth();
        int imageHeight = input.getHeight();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        Set<PairInt> addP = new HashSet<PairInt>();
        Set<PairInt> rmP = new HashSet<PairInt>();
        Set<PairInt> moved = new HashSet<PairInt>();

        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();

            int nNonNeighbors = 0;
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                    (y2 > (imageHeight - 1))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (!points.contains(p2)) {
                   nNonNeighbors++;
                }
            }

            if (nNonNeighbors == 0) {
                continue;
            }

            Set<PairInt> neighbors = new HashSet<PairInt>();
            Set<PairInt> nonNB = new HashSet<PairInt>();
            Set<PairInt> nonNW = new HashSet<PairInt>();
            Set<PairInt> nonNG = new HashSet<PairInt>();

            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                    (y2 > (imageHeight - 1))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (!moved.contains(p2)) {
                    if (points.contains(p2)) {
                        neighbors.add(p2);
                    } else if (blackPixels.contains(p2)) {
                        nonNB.add(p2);
                    } else if (whitePixels.contains(p2)) {
                        nonNW.add(p2);
                    } else if (greyPixels.contains(p2)) {
                        nonNG.add(p2);
                    }
                }
            }

            // calc the avg cieX, cieY of neighbors
            // calc the avg cieX, cieY of non-neighbors
            // if p's cieXY is closer to neighbors, do not move it, but do
            // look at each of the nonNeighbors to see whether to move it...
            //float[] avgCIEXYNeighbors = ImageStatisticsHelper.calculateAvgCIEXY(
            //    input, neighbors);
            float[] avgRGBNeighbors = ImageStatisticsHelper.calculateAvgRGB(
                input, neighbors);

            //float pCIEX = input.getCIEX(x, y);
            //float pCIEY = input.getCIEY(x, y);
            int pR = input.getR(x, y);
            int pG = input.getG(x, y);
            int pB = input.getB(x, y);

            double distSq0 = Math.pow((avgRGBNeighbors[0] - pR), 2) +
                Math.pow((avgRGBNeighbors[1] - pG), 2) +
                Math.pow((avgRGBNeighbors[2] - pB), 2);

            double minDistSq = distSq0;

            // 0 = neighbors, 1=blackPixels, 2=whitePixels, 3=greyPixels
            int minDistGroup = 0;

            if (!nonNB.isEmpty()) {
                float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(
                    input, nonNB);
                double distSq = Math.pow((avgRGB[0] - pR), 2) +
                    Math.pow((avgRGB[1] - pG), 2) +
                    Math.pow((avgRGB[2] - pB), 2);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistGroup = 1;
                }
            }
            if (!nonNW.isEmpty()) {
                float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(
                    input, nonNW);
                double distSq = Math.pow((avgRGB[0] - pR), 2) +
                    Math.pow((avgRGB[1] - pG), 2) +
                    Math.pow((avgRGB[2] - pB), 2);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistGroup = 2;
                }
            }
            if (!nonNG.isEmpty()) {
                float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(
                    input, nonNG);
                double distSq = Math.pow((avgRGB[0] - pR), 2) +
                    Math.pow((avgRGB[1] - pG), 2) +
                    Math.pow((avgRGB[2] - pB), 2);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistGroup = 3;
                }
            }

            // 0 = neighbors, 1=blackPixels, 2=whitePixels, 3=greyPixels
            if (minDistGroup > 0) {
                if (minDistGroup == 1) {
                    blackPixels.add(p);
                } else if (minDistGroup == 2) {
                    whitePixels.add(p);
                } else if (minDistGroup == 3) {
                    greyPixels.add(p);
                }
                rmP.add(p);
                moved.add(p);
            }

            // determine if the non neighbors should be moved into point group
            if (!nonNB.isEmpty()) {
                for (PairInt p2 : nonNB) {
                    if (moved.contains(p2)) {
                        continue;
                    }
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    int p2R = input.getR(x2, y2);
                    int p2G = input.getG(x2, y2);
                    int p2B = input.getB(x2, y2);

                    Set<PairInt> nbN = curveHelper.findNeighbors(x2, y2, blackPixels);
                    float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(input, nbN);

                    double distSqToNeighbors = Math.pow((avgRGB[0] - p2R), 2) +
                        Math.pow((avgRGB[1] - p2G), 2) + Math.pow((avgRGB[2] - p2B), 2);

                    double distSqToP = Math.pow((pR - p2R), 2) +
                        Math.pow((pG - p2G), 2) + Math.pow((pB - p2B), 2);

                    if (distSqToP < distSqToNeighbors) {
                        blackPixels.remove(p2);
                        addP.add(p2);
                        moved.add(p2);
                    }
                }
            }
            if (!nonNW.isEmpty()) {
                for (PairInt p2 : nonNW) {
                    if (moved.contains(p2)) {
                        continue;
                    }
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    int p2R = input.getR(x2, y2);
                    int p2G = input.getG(x2, y2);
                    int p2B = input.getB(x2, y2);

                    Set<PairInt> nbN = curveHelper.findNeighbors(x2, y2, whitePixels);
                    float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(input, nonNW);

                    double distSqToNeighbors = Math.pow((avgRGB[0] - p2R), 2) +
                        Math.pow((avgRGB[1] - p2G), 2) + Math.pow((avgRGB[2] - p2B), 2);

                    double distSqToP = Math.pow((pR - p2R), 2) +
                        Math.pow((pG - p2G), 2) + Math.pow((pB - p2B), 2);

                    if (distSqToP < distSqToNeighbors) {
                        whitePixels.remove(p2);
                        addP.add(p2);
                        moved.add(p2);
                    }
                }
            }
            if (!nonNG.isEmpty()) {
                for (PairInt p2 : nonNG) {
                    if (moved.contains(p2)) {
                        continue;
                    }
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    int p2R = input.getR(x2, y2);
                    int p2G = input.getG(x2, y2);
                    int p2B = input.getB(x2, y2);

                    Set<PairInt> nbN = curveHelper.findNeighbors(x2, y2, greyPixels);
                    float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(input, nonNG);

                    double distSqToNeighbors = Math.pow((avgRGB[0] - p2R), 2) +
                        Math.pow((avgRGB[1] - p2G), 2) + Math.pow((avgRGB[2] - p2B), 2);

                    double distSqToP = Math.pow((pR - p2R), 2) +
                        Math.pow((pG - p2G), 2) + Math.pow((pB - p2B), 2);

                    if (distSqToP < distSqToNeighbors) {
                        greyPixels.remove(p2);
                        addP.add(p2);
                        moved.add(p2);
                    }
                }
            }
        }

        for (PairInt p : addP) {
            points.add(p);
        }
        for (PairInt p : rmP) {
            points.remove(p);
        }

        return moved.size();
    }

    /**
     * an algorithm to traverse the image and store black pixels in a list
     * then use cieXY to place pixels in colorPixelMap if they are not in the
     * large white central space, else place them in unassignedPixels.
     * @param input
     * @param blackPixels
     * @param colorPixelMap
     * @param unassignedPixels
     */
    private void populatePixelLists2(ImageExt input, Set<PairInt> blackPixels,
        Map<PairInt, Float> colorPixelMap, Set<PairInt> unassignedPixels) {

        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        // looking for limits in peaks of (r,g,b) <= (45,45,45) and > (191,191,191)
        int[] whiteBlackLimits = findByHistogramLimitsForBlackAndWhite(input);
        // overriding:

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int idx = input.getInternalIndex(i, j);

                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                float c = r + g + b;
                int avg = (int)(c/3);

                /*
                float rd = (float)r/c;
                float gd = (float)g/c;
                float bd = (float)b/c;
                // r and b 20 percent or more below and g 20 percent or more above
                float diff = 0.1f * 0.33333f;
                boolean veryGreen = (rd < (0.333f - diff)) && (bd < (0.333f - diff)) && (gd > (0.333f + diff));
                */
                //log.info(String.format(
                //    "rgb=(%3d, %3d, %3d)  rgbdiv=(%.2f, %.2f, %.2f) cie=(%.3f, %.3f) xy=(%3d, %3d)",
                //    r, g, b, rd, gd, bd, cieX, cieY, i, j));

                if ((avg <= whiteBlackLimits[0]) /*&& !veryGreen*/) {
                    blackPixels.add(new PairInt(i, j));
                    continue;
                }

                float cieX = input.getCIEX(idx);
                float cieY = input.getCIEY(idx);

                if (cieC.isInLargeWhiteCenter(cieX, cieY) /*&& !veryGreen*/) {
                    unassignedPixels.add(new PairInt(i, j));
                    continue;
                }

                double thetaRadians = cieC.calculateXYTheta(cieX, cieY);

                colorPixelMap.put(new PairInt(i, j),
                    Float.valueOf((float)thetaRadians));
            }
        }
    }

    private Map<PairInt, Integer> populatePixelLists3(ImageExt input) {

        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        Map<PairInt, Integer> pixelCIETheta = new HashMap<PairInt, Integer>();

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int idx = input.getInternalIndex(i, j);

                float cieX = input.getCIEX(idx);
                float cieY = input.getCIEY(idx);

                if (!cieC.isInLargeWhiteCenter(cieX, cieY) /*&& !veryGreen*/) {

                    double thetaRadians = cieC.calculateXYTheta(cieX, cieY);

                    double thetaDegrees = thetaRadians * 180./Math.PI;

                    int thetaDegreesInt = (int)Math.round(thetaDegrees);

                    pixelCIETheta.put(new PairInt(i, j), thetaDegreesInt);
                }
            }
        }

        return pixelCIETheta;
    }

    private void growPixelsWithinRGBTolerance(Image img,
        Set<PairInt> assignedPixels, Set<PairInt> unassignedPixels, int tolRGB) {

        Set<PairInt> addTo = findPixelsWithinRGBTolerance(img,
            assignedPixels, unassignedPixels, tolRGB);

        assignedPixels.addAll(addTo);
    }

    /**
     * traverses assignedPixels and finds adjacent unassignedPixels that are
     * within tolRGB of r, g, b colors and removes those points from
     * unassignedPixels and returns them as a set.
     * @param input
     * @param assignedPixels
     * @param unassignedPixels
     * @param tolRGB
     * @return
     */
    private Set<PairInt> findPixelsWithinRGBTolerance(Image img,
        Set<PairInt> assignedPixels, Set<PairInt> unassignedPixels, int tolRGB) {

        Set<PairInt> addToAssigned = new HashSet<PairInt>();

        int w = img.getWidth();
        int h = img.getHeight();

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        Set<PairInt> visited = new HashSet<PairInt>();

        Stack<PairInt> stack = new Stack<PairInt>();
        stack.addAll(assignedPixels);

        while (!stack.isEmpty()) {
            PairInt uPoint = stack.pop();
            if (visited.contains(uPoint)) {
                continue;
            }
            int x = uPoint.getX();
            int y = uPoint.getY();

            int r = img.getR(x, y);
            int g = img.getG(x, y);
            int b = img.getB(x, y);

            for (int i = 0; i < dxs8.length; ++i) {
                int x2 = x + dxs8[i];
                int y2 = y + dys8[i];
                if ((x2 < 0) || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                    continue;
                }
                PairInt vPoint = new PairInt(x2, y2);
                if (!unassignedPixels.contains(vPoint)) {
                    continue;
                }
                int r2 = img.getR(x2, y2);
                int g2 = img.getG(x2, y2);
                int b2 = img.getB(x2, y2);
                int diffR = Math.abs(r2 - r);
                int diffG = Math.abs(g2 - g);
                int diffB = Math.abs(b2 - b);

                if ((diffR < tolRGB) && (diffG < tolRGB) && (diffB < tolRGB)) {
                    unassignedPixels.remove(vPoint);
                    addToAssigned.add(vPoint);
                    stack.add(vPoint);
                }
            }

            visited.add(uPoint);
        }

        return addToAssigned;
    }

    public void createContrastImages(ImageExt input) {

        int n = input.getNPixels();
        double[] luma = new double[n];
        for (int i = 0; i < n; ++i) {
            int r = input.getR(i);
            int g = input.getG(i);
            int b = input.getB(i);
            double lumaI = (0.256*r) - (-0.148*g) + (0.439*b);
            luma[i] = lumaI;
        }

        // create contrast as (avgLuma - luma[i])/luma[i]
        GreyscaleImage lumaAvg = new GreyscaleImage(input.getWidth(), input.getHeight());
        for (int i = 0; i < n; ++i) {
            lumaAvg.setValue(i, (int)Math.round(luma[i]));
        }

        // range of luma is 0 to 215.  adding 1 to avoid divide by zero
        double[] contrast = new double[n];
        MedianSmooth medSmooth = new MedianSmooth();
        lumaAvg = medSmooth.calculate(lumaAvg, 2, 2);
        //ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.applyAdaptiveMeanThresholding(lumaAvg, 1);
        //medSmooth.calculate(lumaAvg, 1, 1);
        //medSmooth.calculate(lumaAvg, 1, 1);
        for (int i = 0; i < n; ++i) {
            int vAvg = lumaAvg.getValue(i);
            double vI = luma[i] + 1;
            contrast[i] = ((double)vAvg - vI)/vI;
        }
        int[] contrastInt = MiscMath.rescale(contrast, 0, 255);

        double[] blueDivContrast = new double[n];
        double[] redDivContrast = new double[n];

        for (int i = 0; i < n; ++i) {
            int r = input.getR(i);
            //int g = input.getG(i);
            int b = input.getB(i);
            double c = contrastInt[i] + 1; // plus 1 to void divide by zero
            blueDivContrast[i] = (double)b/c;
            redDivContrast[i] = (double)r/c;
        }

        // rescale to be between 0 and 255

        int[] blueDivContrastInt = MiscMath.rescale(blueDivContrast, 0, 255);
        int[] redDivContrastInt = MiscMath.rescale(redDivContrast, 0, 255);

        GreyscaleImage contrastImg = new GreyscaleImage(input.getWidth(), input.getHeight());
        GreyscaleImage blueDivContrastImg = new GreyscaleImage(input.getWidth(), input.getHeight());
        GreyscaleImage redDivContrastImg = new GreyscaleImage(input.getWidth(), input.getHeight());

        for (int i = 0; i < n; ++i) {
            contrastImg.setValue(i, contrastInt[i]);
            blueDivContrastImg.setValue(i, blueDivContrastInt[i]);
            redDivContrastImg.setValue(i, redDivContrastInt[i]);
        }

        long ts = MiscDebug.getCurrentTimeFormatted();

        MiscDebug.writeImage(contrastImg, "_contrast_" + ts);
        MiscDebug.writeImage(blueDivContrastImg, "_blue_div_contrast_" + ts);
        MiscDebug.writeImage(redDivContrastImg, "_red_div_contrast_" + ts);
    }

    /**
     * makes a lower resolution image and uses segmentation based upon
     * cie xy polar theta to find contiguous points of same value, then
     * finds the perimeters of those points and orders them counter-clockwise
     * into the resulting curves.  Note, the resulting curves should all be
     * closed and that can be tested with (edge instance of PairIntWithColor).
     * Note also that straight line segments in the curves have been removed
     * where simple to do so, so that the resulting polynomial can be better
     * used with "point in polygon" tests.
     *
     * runtime complexity is
     * ~ O(N)
     * + O(N_small) + O(N_small * lg_2(N_small)) where N_small is 75X75 pix^2
     * + 16 * O(N_small)
     * + (N_blobs * O(N_perimeter_pts * lg_2(N_perimeter_pts) which is much smaller than O(N_small * lg_2(N_small)))
     * so total is O(N) + constant factor * O(N_small) + O(N_small * lg_2(N_small)).
     * @param img
     * @return
     */
    public BoundingRegions extractBlobBoundsFromLowRes(ImageExt img,
        boolean debug, String debugTag, IntensityClrFeatures features) {

        int w = img.getWidth();
        int h = img.getHeight();
        float maxDimension = 100;//75;

        int binFactor = (int) Math.ceil(Math.max((float)w/maxDimension,
            (float)h/maxDimension));

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageExt img2 = imageProcessor.binImage(img, binFactor);

        ImageExt img2Cp = img2.copyToImageExt();

        // runtime complexity: ~ O(N_small) + O(N_small * lg_2(N_small)) where N_small is 75X75
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        GreyscaleImage segImg = imageSegmentation
            .applyUsingCIEXYPolarThetaThenHistEq(img2, 16, false);

        if (debug) {
            MiscDebug.writeImage(segImg, "seg_cluster_" + MiscDebug.getCurrentTimeFormatted());
        }

        int smallestGroupLimit = 15;
        int largestGroupLimit = Integer.MAX_VALUE;
        boolean filterOutImageBoundaryBlobs = false;
        boolean filterOutZeroPixels = false;
        // runtime complexity is N_freq * O(N) where N_freq is at most 16 and
        // the O(N) term may be as high as O(N*8) if highly connected.
        List<Set<PairInt>> blobs =  BlobsAndPerimeters.extractBlobsFromSegmentedImage(
            segImg, smallestGroupLimit, largestGroupLimit,
            filterOutImageBoundaryBlobs, filterOutZeroPixels, debugTag);

        //--------- begin section to log colors to look at selecting matchable bounds by color ------
        CIEChromaticity cieC = new CIEChromaticity();

        List<Double> labLAvg = new ArrayList<Double>();
        List<Double> labAAvg = new ArrayList<Double>();
        List<Double> labBAvg = new ArrayList<Double>();
        List<Double> rAvg = new ArrayList<Double>();
        List<Double> gAvg = new ArrayList<Double>();
        List<Double> bAvg = new ArrayList<Double>();

        for (int i = 0; i < blobs.size(); ++i) {
            double redSum = 0;
            double greenSum = 0;
            double blueSum = 0;
            for (PairInt p : blobs.get(i)) {
                int x = p.getX();
                int y = p.getY();
                int red = img.getR(x, y);
                int green = img.getG(x, y);
                int blue = img.getB(x, y);
                redSum += red;
                greenSum += green;
                blueSum += blue;
            }
            double n = (double)blobs.get(i).size();
            redSum /= n;
            greenSum /= n;
            blueSum /= n;
            float[] avgLAB = cieC.rgbToCIELAB((int)Math.round(redSum),
                (int)Math.round(greenSum), (int)Math.round(blueSum));

            labLAvg.add(Double.valueOf(avgLAB[0]));
            labAAvg.add(Double.valueOf(avgLAB[1]));
            labBAvg.add(Double.valueOf(avgLAB[2]));

            rAvg.add(Double.valueOf(redSum));
            gAvg.add(Double.valueOf(greenSum));
            bAvg.add(Double.valueOf(blueSum));
        }

        // less than O(N)
        List<Set<PairInt>> borderPixelSets = BlobsAndPerimeters.extractBlobPerimeterAsPoints(
            blobs, segImg.getWidth(), segImg.getHeight());

        assert(blobs.size() == borderPixelSets.size());

        List<PairIntArray> perimetersList = new ArrayList<PairIntArray>();

        float srchRadius = (float)Math.sqrt(2) * (float)binFactor;

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        // scale blobs and borderPixelSets back to current frame
        for (int i = 0; i < borderPixelSets.size(); ++i) {
            Set<PairInt> blob = blobs.get(i);
            Set<PairInt> tmp = new HashSet<PairInt>();
            for (PairInt p : blob) {
                int x = p.getX() * binFactor;
                int y = p.getY() * binFactor;
                tmp.add(new PairInt(x, y));
            }
            blob.clear();
            blob.addAll(tmp);

            Set<PairInt> borderPixels = borderPixelSets.get(i);
            tmp = new HashSet<PairInt>();
            for (PairInt p : borderPixels) {
                int x = p.getX() * binFactor;
                int y = p.getY() * binFactor;
                tmp.add(new PairInt(x, y));
            }
            borderPixels.clear();
            borderPixels.addAll(tmp);
        }

        // create a map for reverse look-ups later
        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < blobs.size(); ++i) {
            Integer key = Integer.valueOf(i);
            Set<PairInt> blob = blobs.get(i);
            for (PairInt p : blob) {
                pointIndexMap.put(p, key);
            }
        }

        BlobMedialAxes bma = new BlobMedialAxes(blobs, labLAvg, labAAvg, labBAvg,
            rAvg, gAvg, bAvg);

        for (int i = 0; i < borderPixelSets.size(); ++i) {

            Set<PairInt> blob = blobs.get(i);
            Set<PairInt> borderPixels = borderPixelSets.get(i);

            // approx O(N_perimeter), but has factors during searches that could be improved
            PairIntArray orderedPerimeter = perimeterFinder.orderThePerimeter(
                borderPixels, blob, srchRadius,
                bma, i);

            /*Image imgCp = img.copyImage();
            ImageIOHelper.addCurveToImage(orderedPerimeter, imgCp, 2, 255, 0, 0);
            MiscDebug.writeImage(imgCp, "_" + i + "_" + MiscDebug.getCurrentTimeFormatted());
            */

            // runtime complexity is O(N_perimeter_pts * lg_2(N_perimeter_pts)
            // remove straight line segments except their endpoints to make simpler
            // polynomial for "point in polygon" tests
            makeStraightLinesHollow(orderedPerimeter, img.getWidth(),
                img.getHeight(), srchRadius);

            /*
            Image imgCp = img.copyImage();
            ImageIOHelper.addCurveToImage(orderedPerimeter, imgCp, 2, 255, 0, 0);
            MiscDebug.writeImage(imgCp, "_" + i + "_" + MiscDebug.getCurrentTimeFormatted());
            */

            perimetersList.add(orderedPerimeter);
        }

        /*
        TODO: consider expanding the bounds to the nearest smaller segmentation
        as refinement of their rougher locations.
        */

        BoundingRegions br = new BoundingRegions(perimetersList, bma, pointIndexMap);

        //plotOrientation(features, br, img, debugTag);

        return br;
    }

    private double countFreq(PairIntArray sortedFreqL) {

        long n = 0;
        for (int i = 0; i < sortedFreqL.getN(); ++i) {
            n += sortedFreqL.getY(i);
        }
        return n;
    }

    private List<Float> countFraction(GreyscaleImage img, int value,
        int colCellSize, int rowCellSize) {

        int w = img.getWidth();
        int h = img.getHeight();

        List<Float> fractionZeros = new ArrayList<Float>();

        int col = 0;
        while (col < w) {
            int row = 0;
            while (row < h) {
                //look at next 50x50 region
                int nZ = 0;
                int nNonZ = 0;
                for (int i = col; i < (col + colCellSize); ++i) {
                    if (i > (w - 1)) {
                        break;
                    }
                    for (int j = row; j < (row + rowCellSize); ++j) {
                        if (j > (h - 1)) {
                            break;
                        }
                        int v = img.getValue(i, j);
                        if (v == value) {
                            nZ++;
                        } else {
                            nNonZ++;
                        }
                    }
                }
                fractionZeros.add(Float.valueOf((float)nZ/(float)(nZ + nNonZ)));
                row += 50;
            }
            col += 50;
        }

        return fractionZeros;
    }

    private List<Float> countFractionNonZeros(GreyscaleImage img,
        int colCellSize, int rowCellSize) {

        int w = img.getWidth();
        int h = img.getHeight();

        List<Float> fraction = new ArrayList<Float>();

        int col = 0;
        while (col < w) {
            int row = 0;
            while (row < h) {
                //look at next 50x50 region
                int n1 = 0;
                int n2 = 0;
                for (int i = col; i < (col + colCellSize); ++i) {
                    if (i > (w - 1)) {
                        break;
                    }
                    for (int j = row; j < (row + rowCellSize); ++j) {
                        if (j > (h - 1)) {
                            break;
                        }
                        int v = img.getValue(i, j);
                        if (v == 0) {
                            n2++;
                        } else {
                            n1++;
                        }
                    }
                }
                fraction.add(Float.valueOf((float)n1/(float)(n1 + n2)));
                row += 50;
            }
            col += 50;
        }

        return fraction;
    }

    private GreyscaleImage keepContigAboveLimit(GreyscaleImage img, int value,
        int sizeLimit) {

        DFSContiguousValueFinder cf = new DFSContiguousValueFinder(img);
        cf.findGroups(value);
        GreyscaleImage tmpImg2 = img.createWithDimensions();
        for (int i = 0; i < cf.getNumberOfGroups(); ++i) {
            PairIntArray zeros = cf.getXY(i);
            if (zeros.getN() > sizeLimit) {
                for (int j = 0; j < zeros.getN(); ++j) {
                    tmpImg2.setValue(zeros.getX(j), zeros.getY(j), 255);
                }
            }
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        imageProcessor.applyAdaptiveMeanThresholding(tmpImg2, 1);

        return tmpImg2;
    }

    private GreyscaleImage expandBy1(GreyscaleImage img, int value) {

        int w = img.getWidth();
        int h = img.getHeight();

        GreyscaleImage tmpImg2 = img.copyImage();
        int[] dxs0 = Misc.dx8;
        int[] dys0 = Misc.dy8;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = img.getValue(i, j);
                if (v != value) {
                    continue;
                }
                for (int k = 0; k < dxs0.length; ++k) {
                    int x1 = i + dxs0[k];
                    int y1 = j + dys0[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    tmpImg2.setValue(x1, y1, v);
                }
            }
        }

        return tmpImg2;
    }

    public GreyscaleImage fillInGapsOf1(GreyscaleImage img,
        Set<PairInt> outputAddedGaps, int value) {

        int w = img.getWidth();
        int h = img.getHeight();

        /*
        0  1  2
        7     3
        6  5  4
        fill in !value if these pairs are filled in:
            0:3, 0:4, 0:5
            1:4, 1:5, 1:6
            2:5, 2:6, 2:7
            3:6, 3:7, 3:0
            4:7
        so a +1 and -1 in x or y and a +1 or -1 in y or x
        */
        int[] dxs0 = new int[]{-1, -1, -1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1};
        int[] dys0 = new int[]{+1, +1, +1,  1,  1,  1,  1,  1,  1,  0,  0,  0, -1};
        int[] dxs1 = new int[]{1,  +1,  0,  1,  0, -1,  0, -1, -1, -1, -1, -1, -1};
        int[] dys1 = new int[]{0,  -1, -1, -1, -1, -1, -1, -1,  0, -1,  0,  1,  0};

        GreyscaleImage tmpImg2 = img.copyImage();

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int v = img.getValue(i, j);

                if (v == value) {
                    continue;
                }

                for (int k = 0; k < dxs0.length; ++k) {
                    int x1 = i + dxs0[k];
                    int y1 = j + dys0[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    int v1 = img.getValue(x1, y1);
                    if (v1 != value) {
                        continue;
                    }
                    int x2 = i + dxs1[k];
                    int y2 = j + dys1[k];
                    if (x2 < 0 || (x2 > (w - 1)) || y2 < 0 || (y2 > (h - 1))) {
                        continue;
                    }
                    int v2 = img.getValue(x2, y2);
                    if (v2 != value) {
                        continue;
                    }
                    tmpImg2.setValue(i, j, value);
                    outputAddedGaps.add(new PairInt(i, j));
                    break;
                }
            }
        }

        return tmpImg2;
    }

    private GreyscaleImage fillInCompleteGapsOf1(GreyscaleImage img,
        Set<PairInt> outputAddedGaps,int value) {

        int w = img.getWidth();
        int h = img.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        GreyscaleImage tmpImg2 = img.copyImage();

        int nIter = 0;
        int nChanged = 0;
        int nMaxIter = 5;
        while ((nIter == 0) || ((nChanged > 0) && (nIter < nMaxIter))) {
            nChanged = 0;
            for (int i = 0; i < w; ++i) {
                for (int j = 0; j < h; ++j) {
                    int v = img.getValue(i, j);
                    if (v == value) {
                        continue;
                    }
                    int count = 0;
                    int neighborCount = 0;
                    for (int k = 0; k < dxs.length; ++k) {
                        int x2 = i + dxs[k];
                        int y2 = j + dys[k];
                        if (x2 < 0 || (x2 > (w - 1)) || y2 < 0 || (y2 > (h - 1))) {
                            continue;
                        }
                        count++;
                        int v1 = img.getValue(x2, y2);
                        if (v1 != value) {
                            continue;
                        }
                        neighborCount++;
                    }
                    if (count == neighborCount) {
                        tmpImg2.setValue(i, j, value);
                        outputAddedGaps.add(new PairInt(i, j));
                        nChanged++;
                    }
                }
            }
            nIter++;
        }

        return tmpImg2;
    }

    GreyscaleImage shrinkBy1(GreyscaleImage img, int edgeValue, int nonEdgeValue) {

        int w = img.getWidth();
        int h = img.getHeight();

        /*
         any pixel with neighbors that are not edgeValue can be removed
         */
        PairInt[][] neighborCoordOffsets
            = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                0, 0);

        GreyscaleImage tmpImg2 = img.copyImage();
        int[] dxs0 = Misc.dx8;
        int[] dys0 = Misc.dy8;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int v = img.getValue(i, j);
                if (v != edgeValue) {
                    continue;
                }
                int nEmptyNeigbhors = 0;
                for (int k = 0; k < dxs0.length; ++k) {
                    int x1 = i + dxs0[k];
                    int y1 = j + dys0[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    int v2 = img.getValue(x1, y1);
                    if (v2 != edgeValue) {
                        nEmptyNeigbhors++;
                    }
                }

                if (nEmptyNeigbhors > 0
                    && !doesDisconnect(tmpImg2, neighborCoordOffsets, i, j, edgeValue)) {
                    tmpImg2.setValue(i, j, nonEdgeValue);
                }
            }
        }

        return tmpImg2;
    }

    private Set<PairInt> createZerosSet(List<Set<PairInt>> segmentedCellList,
        int w, int h, Set<PairInt> mask) {

        Set<PairInt> points = new HashSet<PairInt>();

        for (Set<PairInt> set : segmentedCellList) {
            points.addAll(set);
        }

        Set<PairInt> nonPoints = new HashSet<PairInt>();

        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                if (!points.contains(p) && !mask.contains(p)) {
                    nonPoints.add(p);
                }
            }
        }

        return nonPoints;
    }

    /**
     * continues the segmentation by placing the unassigned pixels into adjacent
     * cells if color is similar to the average color of the cell.
     * @param input
     * @param segmentedCellList
     * @param unassigned
     * @param useAvgCellColor uses the average color of the cell if true, else
     * uses the color of the adjacent pixel.
     */
    private void placeUnassignedByGrowingCells(ImageExt input,
        List<Set<PairInt>> segmentedCellList, Set<PairInt> unassigned,
        Map<PairInt, Integer> pointIndexMap, double deltaELimit,
        boolean useDeltaE2000) {

        boolean useAvgCellColor = true;

        placeUnassignedByGrowingCells(input, segmentedCellList, unassigned,
            pointIndexMap, deltaELimit, useDeltaE2000, useAvgCellColor);
    }

    /**
     * continues the segmentation by placing the unassigned pixels into adjacent
     * cells if color is similar to the average color of the cell.
     * @param input
     * @param segmentedCellList
     * @param unassigned
     * @param useAvgCellColor uses the average color of the cell if true, else
     * uses the color of the adjacent pixel.
     */
    private void placeUnassignedByGrowingCells(ImageExt input,
        List<Set<PairInt>> segmentedCellList, Set<PairInt> unassigned,
        Map<PairInt, Integer> pointIndexMap, double deltaELimit,
        boolean useDeltaE2000, boolean useAvgCellColor) {

        long t0 = System.currentTimeMillis();

        CIEChromaticity cieC = new CIEChromaticity();

        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();

        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            queue.addAll(set);
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        Map<Integer, Colors> segmentedCellAvgLabColors = useAvgCellColor ?
            new HashMap<Integer, Colors>() : null;

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        Set<PairInt> visited = new HashSet<PairInt>();

        while (!queue.isEmpty()) {

            PairInt p0 = queue.pop();

            if (visited.contains(p0)) {
                continue;
            }

            int x = p0.getX();
            int y = p0.getY();

            Integer listIndex = pointIndexMap.get(p0);

            float[] lab0 = null;

            for (int i = 0; i < dxs.length; ++i) {

                int x2 = x + dxs[i];
                int y2 = y + dys[i];

                PairInt p2 = new PairInt(x2, y2);
                if (!unassigned.contains(p2)) {
                    continue;
                }

                if (lab0 == null) {
                    if (useAvgCellColor) {
                        Colors colors0 = segmentedCellAvgLabColors.get(listIndex);
                        if (colors0 == null) {
                            colors0 = imageProcessor.calculateAverageLAB(input,
                                segmentedCellList.get(listIndex.intValue()));
                            segmentedCellAvgLabColors.put(listIndex, colors0);
                        }
                        lab0 = colors0.getColors();
                    } else {
                        lab0 = input.getCIELAB(x, y);
                    }
                }

                assert(!pointIndexMap.containsKey(p2));

                float[] lab2 = input.getCIELAB(x2, y2);

                double deltaE;
                if (useDeltaE2000) {
                    deltaE = Math.abs(cieC.calcDeltaECIE2000(lab0, lab2));
                } else {
                    deltaE = Math.abs(cieC.calcDeltaECIE94(lab0, lab2));
                }

                // jnd ~ 2.3
                // 4 is good and can be continued w/ labelling
                // 5 is fine if goal is finding objects
                if (deltaE > deltaELimit) {
                    continue;
                }

                segmentedCellList.get(listIndex.intValue()).add(p2);

                pointIndexMap.put(p2, listIndex);

                unassigned.remove(p2);

                queue.add(p2);
            }

            visited.add(p0);
        }

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec to place " + unassigned.size()
            + " points by growing cells");
    }

    /**
     *
     * @param input
     * @param segmentedCellList
     * @param zeros
     * @param deltaELimit
     * @param radius distance to search around a point for the nearest points
     * within sets in segmentedCellList. The search runtime is roughly
     * O(N * radius * lg2(N)) where N is the total number of points in
     * segmentedCellList, so keep radius as small as possible (and less than 10).
     */
    private void placeUnassignedUsingNearest(ImageExt input,
        List<Set<PairInt>> segmentedCellList, Set<PairInt> zeros,
        double deltaELimit, float radius, boolean useDeltaE2000) {

        /*
        for the points which are 0's in the greyGradient,
           finding the nearest points from segmented cells within a radius
           of 28 or so, limited by max for deltaE.
           for each returned,
               calculate dist = xyDist + abs(deltaE)
               and determine minDist
           assign the 0 point to the list of the minimum distance
           and update the np
        */

        Set<PairInt> assignedZeros = new HashSet<PairInt>();

        CIEChromaticity cieC = new CIEChromaticity();

        long t0 = System.currentTimeMillis();
        NearestPointsInLists np = new NearestPointsInLists(segmentedCellList);

        boolean useDistAndColor = true;

        // since some of the results will not be contiguous, need to research
        // the edited lists with dfs when finished to split set if there are gaps
        Set<Integer> addedTo = new HashSet<Integer>();
        for (PairInt p : zeros) {

            Map<Integer, PairInt> listPointMap = np.findNeighbors(p.getX(),
                p.getY(), radius);

            if (listPointMap.isEmpty()) {
                continue;
            }

            float[] lab0 = input.getCIELAB(p.getX(), p.getY());

            double minDist = Double.MAX_VALUE;
            Integer minDistIdx = null;

            for (Entry<Integer, PairInt> entry : listPointMap.entrySet()) {

                PairInt pClosest = entry.getValue();

                float[] lab1 = input.getCIELAB(pClosest.getX(), pClosest.getY());

                double deltaE;
                if (useDeltaE2000) {
                    deltaE = Math.abs(cieC.calcDeltaECIE2000(lab0, lab1));
                } else {
                    deltaE = Math.abs(cieC.calcDeltaECIE94(lab0, lab1));
                }

                //the Just Noticeable Difference is ~2.3

                if (deltaE > deltaELimit) {
                    continue;
                }

                double dist = deltaE;

                if (useDistAndColor) {
                    int diffX = p.getX() - pClosest.getX();
                    int diffY = p.getY() - pClosest.getY();
                    dist += Math.sqrt(diffX*diffX + diffY*diffY);
                }

                if (dist < minDist) {
                    minDist = dist;
                    minDistIdx = entry.getKey();
                }
            }

            if (minDistIdx != null) {
                addedTo.add(minDistIdx);
                segmentedCellList.get(minDistIdx.intValue()).add(p);
                np.addPoint(p, minDistIdx);
                assignedZeros.add(p);
            }
        }

        zeros.removeAll(assignedZeros);

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec to place " + zeros.size()
            + " points in nearest cells by color");
        t0 = System.currentTimeMillis();

        List<Set<PairInt>> segmentedCellList2 = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Integer index = Integer.valueOf(i);
            Set<PairInt> points = segmentedCellList.get(i);
            if (!addedTo.contains(index)) {
                segmentedCellList2.add(points);
                continue;
            }
            DFSConnectedGroupsFinder finder2 = new DFSConnectedGroupsFinder();
            finder2.setMinimumNumberInCluster(1);
            finder2.findConnectedPointGroups(points);
            for (int ii = 0; ii < finder2.getNumberOfGroups(); ++ii) {
                Set<PairInt> group = finder2.getXY(ii);
                segmentedCellList2.add(group);
            }
        }

        t1 = System.currentTimeMillis();
        t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec to refine clustered point sets");

        segmentedCellList.clear();
        segmentedCellList.addAll(segmentedCellList2);
    }

    private float calculateAverageHueAngle(ImageExt input, Set<PairInt> points) {

        double labA = 0;
        double labB = 0;

        for (PairInt p : points) {
            float[] lab = input.getCIELAB(p.getX(), p.getY());
            labA += lab[1];
            labB += lab[2];
        }
        labA /= (double)points.size();
        labB /= (double)points.size();

        float ha;
        if (labA == 0) {
            ha = 0;
        } else {
            ha = (float) (Math.atan2(labB, labA) * 180. / Math.PI);
            if (ha < 0) {
                ha += 360.;
            }
        }

        return ha;
    }

    private boolean areAdjacent(Set<PairInt> setA, Set<PairInt> setB) {

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        for (PairInt p : setA) {
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = p.getX() + dxs[i];
                int y2 = p.getY() + dys[i];
                PairInt p2 = new PairInt(x2, y2);
                if (setB.contains(p2)) {
                    return true;
                }
            }
        }
        return false;
    }

    private void invertImage(GreyscaleImage img) {
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v = img.getValue(i);
            img.setValue(i, 255 - v);
        }
    }
    private void setAllNonZeroTo255(GreyscaleImage img) {
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v = img.getValue(i);
            if (v > 0) {
                img.setValue(i, 255);
            }
        }
    }
    private void setAllNon255To0(GreyscaleImage img) {
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v = img.getValue(i);
            if (v < 255) {
                img.setValue(i, 0);
            }
        }
    }

    private void mergeEmbeddedIfSimilar(ImageExt input, List<Set<PairInt>>
        segmentedCellList, Map<PairInt, Integer> pointIndexMap,
        double deltaELimit, boolean useDeltaE2000) {

        Map<Integer, Colors> segmentedCellAvgLabColors = new HashMap<Integer, Colors>();

        ImageProcessor imageProcessor = new ImageProcessor();

        CIEChromaticity cieC = new CIEChromaticity();

        long t0 = System.currentTimeMillis();

        /*
        for each set:
           for all points, if a neighbor is in an adjacent cell, store the
           cell numbers and if there is only one, then the
           set is embedded (or adjacent).
           merge if they are similar.
        */

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int nIter = 0;
        int nChanged = 0;

        while ((nIter == 0) || (nChanged > 0)) {

            nChanged = 0;

            for (int i = 0; i < segmentedCellList.size(); ++i) {

                Set<PairInt> set = segmentedCellList.get(i);

                if (set.isEmpty()) {
                    continue;
                }

                Colors colors1 = segmentedCellAvgLabColors.get(Integer.valueOf(i));
                if (colors1 == null) {
                    colors1 = imageProcessor.calculateAverageLAB(input, set);
                    segmentedCellAvgLabColors.put(Integer.valueOf(i), colors1);
                }

                Integer adjacentListIndex = null;

                boolean cannotMerge = false;

                for (PairInt p : set) {
                    int x = p.getX();
                    int y = p.getY();

                    for (int k = 0; k < dxs.length; ++k) {
                        int x2 = x + dxs[k];
                        int y2 = y + dys[k];
                        PairInt p2 = new PairInt(x2, y2);
                        if (set.contains(p2)) {
                            continue;
                        }
                        Integer listIndex2 = pointIndexMap.get(p2);
                        if (listIndex2 == null) {
                            continue;
                        }
                        if ((adjacentListIndex != null)
                            && !adjacentListIndex.equals(listIndex2)) {
                            cannotMerge = true;
                            break;
                        } else if (adjacentListIndex == null) {
                            adjacentListIndex = listIndex2;
                        }
                    }
                    if (cannotMerge) {
                        break;
                    }
                }
                if (cannotMerge || (adjacentListIndex == null)) {
                    continue;
                }

                Set<PairInt> set2 = segmentedCellList.get(adjacentListIndex.intValue());

                Colors colors2 = segmentedCellAvgLabColors.get(adjacentListIndex);
                if (colors2 == null) {
                    colors2 = imageProcessor.calculateAverageLAB(input, set2);
                    segmentedCellAvgLabColors.put(adjacentListIndex, colors2);
                }

                double deltaE;
                if (useDeltaE2000) {
                    deltaE = Math.abs(cieC.calcDeltaECIE2000(colors1.getColors(),
                        colors2.getColors()));
                } else {
                    deltaE = Math.abs(cieC.calcDeltaECIE94(colors1.getColors(),
                        colors2.getColors()));
                }

                if (deltaE > deltaELimit) {
                    continue;
                }

                nChanged++;

                /*choose to merge such that the set that keeps same index is
                   the set with the largest number of points so that the cached
                   color remains closer to representing the set's points without
                   recalculating
                */

                int n1 = set.size();
                int n2 = set2.size();

                boolean merge2Into1 = true;

                if (n1 == n2) {
                    if (i > adjacentListIndex.intValue()) {
                        merge2Into1 = false;
                    }
                } else if (n1 < n2) {
                    merge2Into1 = false;
                }

                if (merge2Into1) {
                    set.addAll(set2);
                    for (PairInt p2 : set2) {
                        pointIndexMap.put(p2, Integer.valueOf(i));
                    }
                    set2.clear();
                } else {
                    set2.addAll(set);
                    for (PairInt p1 : set) {
                        pointIndexMap.put(p1, adjacentListIndex);
                    }
                    set.clear();
                    break;
                }
            }

            nIter++;
        }

        int count = 0;
        for (int i = (segmentedCellList.size() - 1); i > -1; --i) {
            Set<PairInt> set0 = segmentedCellList.get(i);
            if (set0.size() == 0) {
                segmentedCellList.remove(i);
                count++;
            }
        }

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, merged " + count + " embedded similar cells "
            + " nIter=" + nIter);
    }

    private void mergeAdjacentIfSimilar(ImageExt input, List<Set<PairInt>>
        segmentedCellList, double deltaELimit, boolean useDeltaE2000,
        String debugTag) {

        SegmentedCellMerger scm = new SegmentedCellMerger(input,
            segmentedCellList, useDeltaE2000, (float)deltaELimit, debugTag);

        scm.merge();

        segmentedCellList.clear();
        segmentedCellList.addAll(scm.getSegmentedCellList());
    }

    /**
     * a merge algorithm that looks at the colors of the adjacent pixels individually,
     * then performs stats on all adjacent for two sets to determine if the
     * edge is similar enough to merge.  This differs from the other method
     * which uses the average color of all points in a cell and the difference
     * between those to decide on merging.
     * @param input
     * @param segmentedCellList
     * @param pointIndexMap
     * @param deltaELimit
     * @param useDeltaE2000
     * @param debugTag
     */
    private void mergeAdjacentIfSimilar2(ImageExt input, List<Set<PairInt>>
        segmentedCellList, Map<PairInt, Integer> pointIndexMap,
        double deltaELimit, boolean useDeltaE2000,
        String debugTag) {

        long t0 = System.currentTimeMillis();

        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int count = 0;

        CIEChromaticity cieC = new CIEChromaticity();

        Stack<Integer> stack = new Stack<Integer>();
        for (int i = (segmentedCellList.size() - 1); i > -1; --i) {
            stack.add(Integer.valueOf(i));
        }

        while (!stack.isEmpty()) {

            Integer index = stack.pop();

            Set<PairInt> set = segmentedCellList.get(index.intValue());

            if (set.size() == 0) {
                continue;
            }

            // only comparing bordering points
            // storing all then averaging and comparing to deltaELimit
            Map<Integer, List<Double>> listIndexDeltaEMap = new HashMap<Integer, List<Double>>();

            for (PairInt p : set) {

                int x = p.getX();
                int y = p.getY();
                Integer listIndex = pointIndexMap.get(p);
                assert(listIndex.intValue() == index.intValue());

                float[] lab = input.getCIELAB(x, y);

                for (int k = 0; k < dxs.length; ++k) {

                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];

                    PairInt p2 = new PairInt(x2, y2);
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    Integer listIndex2 = pointIndexMap.get(p2);
                    if (listIndex2 == null || listIndex.equals(listIndex2)) {
                        continue;
                    }

                    Set<PairInt> set2 = segmentedCellList.get(listIndex2.intValue());

                    int n2 = set2.size();

                    if (n2 == 0) {
                        continue;
                    }

                    float[] lab2 = input.getCIELAB(x2, y2);

                    double deltaE;
                    if (useDeltaE2000) {
                        deltaE = Math.abs(cieC.calcDeltaECIE2000(lab, lab2));
                    } else {
                        deltaE = Math.abs(cieC.calcDeltaECIE94(lab, lab2));
                    }

                    List<Double> dEs = listIndexDeltaEMap.get(listIndex2);
                    if (dEs == null) {
                        dEs = new ArrayList<Double>();
                        listIndexDeltaEMap.put(listIndex2, dEs);
                    }
                    dEs.add(Double.valueOf(deltaE));
                }
            }

            boolean changed = false;

            for (Entry<Integer, List<Double>> entry : listIndexDeltaEMap.entrySet()) {

                Integer listIndex2 = entry.getKey();

                List<Double> deltaEs = entry.getValue();

                double binSize = 0.25;
                if ((deltaELimit/2.) < binSize) {
                    binSize = deltaELimit/2.;
                }

                // {mean, stdev of mean, mode}
                double[] mnStDevMode = MiscStats.calculateMeanStDevAndMode(
                    deltaEs, 0.25);

                /*
                TODO: this inequality in doMerge suggests that this part of
                a larger segmentation algorithm could be solved by an
                integer programming layer of neural network,
                that is, converting min cost to integer programming.
                The other layer(s) would primarily be strict comparison to
                a deltaE limit to be optimized for merging or growing.
                Need to re-do the labelling needed for such an algorithm...
                */

                boolean doMerge = (mnStDevMode[0] <= deltaELimit) ||
                    (
                    (Math.abs(mnStDevMode[1] - deltaELimit) < (0.4*deltaELimit))
                    && (mnStDevMode[0] > mnStDevMode[2]) && (mnStDevMode[0] > mnStDevMode[1])
                    && (Math.abs(mnStDevMode[0] - deltaELimit) < (0.25*mnStDevMode[1]))
                    && (Math.abs(mnStDevMode[2] - deltaELimit) < (0.25*mnStDevMode[1]))
                    );

                if (!doMerge) {
                    continue;
                }

                /*
                if (mnStDevMode[0] > deltaELimit) {
                    continue;
                }*/

                ++count;

                Set<PairInt> set2 = segmentedCellList.get(listIndex2.intValue());
                set.addAll(set2);
                for (PairInt p3 : set2) {
                    pointIndexMap.put(p3, index);
                }
                set2.clear();

                changed = true;
            }

            if (changed) {
                stack.add(index);
            }
        }

        count = 0;
        for (int i = (segmentedCellList.size() - 1); i > -1; --i) {
            Set<PairInt> set0 = segmentedCellList.get(i);
            if (set0.size() == 0) {
                segmentedCellList.remove(i);
                count++;
            }
        }

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec to merge  " + count + " cells");
    }

    private GreyscaleImage createUncrossableEdges(ImageExt img,
        GreyscaleImage greyImg, String debugTag) {

        GreyscaleImage greyImgC = greyImg.copyImage();

        ImageProcessor imageProcessor = new ImageProcessor();

        imageProcessor.highPassIntensityFilter(greyImg, 0.3);//0.1
        invertImage(greyImg);
        setAllNon255To0(greyImg);

        MiscDebug.writeImage(greyImg, "_tmp_grey_cuts_0_" + debugTag);
        MedianSmooth smooth = new MedianSmooth();
        greyImg = smooth.calculate(greyImg, 4, 4);
        greyImg = smooth.calculate(greyImg, 2, 2);
        MiscDebug.writeImage(greyImg, "_tmp_grey_cuts_1_" + debugTag);
        greyImg = fillInGapsOf1(greyImg, new HashSet<PairInt>(), 0);
        MiscDebug.writeImage(greyImg, "_tmp_grey_cuts_2_" + debugTag);
        imageProcessor.applyAdaptiveMeanThresholding(greyImg, 1);

        MiscDebug.writeImage(greyImg, "_tmp_grey_cuts_3_" + debugTag);

        /*
        bRImg = s.calculate(bRImg, 3, 3);
        invertImage(bRImg);
//NOTE: because the 0's are placed in algorithm below, consider not using adaptive mean here
        imageProcessor.applyAdaptiveMeanThresholding(bRImg, 1);
        */

        return greyImg;
    }

    private List<Set<PairInt>> findContiguousCells(int value, GreyscaleImage img,
        Set<PairInt> mask) {

        long t0 = System.currentTimeMillis();

        // uses the 4-neighbor region in search
        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img, mask);
        finder.setMinimumNumberInCluster(2);
        finder.findGroups(value);

        List<Set<PairInt>> segmentedCellList = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            PairIntArray group = finder.getXY(i);
            Set<PairInt> set = new HashSet<PairInt>();
            for (int j = 0; j < group.getN(); ++j) {
                set.add(new PairInt(group.getX(j), group.getY(j)));
            }
            segmentedCellList.add(set);
        }

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec to make sets of contiguous pixels");
        t0 = System.currentTimeMillis();

        return segmentedCellList;
    }

    private void placeUnassignedAndMergeEmbedded(ImageExt input,
        List<Set<PairInt>> segmentedCellList, double deltaELimit,
        Set<PairInt> mask, boolean useDeltaE2000, float radius) {

        int w = input.getWidth();
        int h = input.getHeight();

        Set<PairInt> unassigned = createZerosSet(segmentedCellList, w, h, mask);

        // for placing points that are near more than one boundary:
        placeUnassignedUsingNearest(input, segmentedCellList, unassigned,
            deltaELimit, radius, useDeltaE2000);

        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }

        placeUnassignedByGrowingCells(input, segmentedCellList, unassigned,
            pointIndexMap, deltaELimit, useDeltaE2000);

        mergeEmbeddedIfSimilar(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000);
    }

    private void extractSetsThatShouldNotBeMerged(ImageExt input,
        List<Set<PairInt>> segmentedCellList,
        List<Set<PairInt>> outputMaskList) {

        ImageProcessor imageProcessor = new ImageProcessor();

        for (int i = 0; i < segmentedCellList.size(); ++i) {

            Set<PairInt> set = segmentedCellList.get(i);

            Colors colors = imageProcessor.calculateAverageRGB(input, set);

            if (isBrightBlueOrGrey(colors.getColors())) {
                /*
                might need to widen this, but hue angle near 51 and near 275
                seem to be one of the areas that deltaE2000 improves deltaE94,
                grey and blueish...
                */
                //double hueAngle = calculateAverageHueAngle(input, set);
                //if ((Math.abs(hueAngle - 51) < 25)
                //    || (Math.abs(hueAngle - 275) < 25)) {
                    outputMaskList.add(new HashSet<PairInt>(set));
                    set.clear();
                //}
            }
        }

        int count = 0;
        for (int i = (segmentedCellList.size() - 1); i > -1; --i) {
            Set<PairInt> set0 = segmentedCellList.get(i);
            if (set0.size() == 0) {
                segmentedCellList.remove(i);
                count++;
            }
        }
    }

    private boolean isBrightBlueOrGrey(float[] rgb) {

        if ((rgb[1] < 85) && (rgb[2] < 85)) {
            return false;
        }

        if (((rgb[0] < rgb[1]) || ((rgb[0] - rgb[1]) <= 5)) && (rgb[1] < rgb[2])) {
            return true;
        }

        int limit = 10;
        if ((Math.abs(rgb[0] - rgb[1]) < limit) &&
            (Math.abs(rgb[0] - rgb[2]) < limit) &&
            (Math.abs(rgb[1] - rgb[2]) < limit)) {
            return true;
        }

        return false;
    }

    private void findDarkPixels(ImageExt img, Set<PairInt> mask) {

        for (int i = 0; i < img.getNPixels(); ++i) {

            int r = img.getR(i);
            int g = img.getG(i);
            int b = img.getB(i);

            if ((r < 85) && (g < 85) && (b < 85)) {
                mask.add(new PairInt(img.getCol(i), img.getRow(i)));
            }
        }
    }

    private void reassignSmallestGroups(ImageExt input, List<Set<PairInt>>
        segmentedCellList, Map<PairInt, Integer> pointIndexMap,
        boolean useDeltaE2000, String debugTag) {

        long t0 = System.currentTimeMillis();

        int n = segmentedCellList.size();

        Map<Integer, Colors> segmentedCellAvgLabColors = new HashMap<Integer, Colors>();

        ImageProcessor imageProcessor = new ImageProcessor();

        CIEChromaticity cieC = new CIEChromaticity();

        int[] sizes = new int[n];
        int[] indexes = new int[n];
        for (int i = 0; i < n; ++i) {
            Set<PairInt> group = segmentedCellList.get(i);
            sizes[i] = group.size();
            indexes[i] = i;
        }
        MultiArrayMergeSort.sortByDecr(sizes, indexes);

        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int limit = 12;
        for (int i = (n - 1); i > -1; --i) {

            int idx = indexes[i];

            Set<PairInt> group = segmentedCellList.get(idx);
            if (group.size() > limit) {
                break;
            }

            Integer listIndex = Integer.valueOf(idx);

            Colors colors = segmentedCellAvgLabColors.get(listIndex);
            if (colors == null) {
                colors = imageProcessor.calculateAverageLAB(input, group);
                segmentedCellAvgLabColors.put(listIndex, colors);
            }

            double minDeltaE = Double.MAX_VALUE;
            int minDeltaEIdx = -1;

            for (PairInt p : group) {

                int x = p.getX();
                int y = p.getY();

                for (int k = 0; k < dxs.length; ++k) {

                    int x1 = x + dxs[k];
                    int y1 = y + dys[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    Integer listIndex1 = pointIndexMap.get(new PairInt(x1, y1));
                    if ((listIndex1 == null) || listIndex1.equals(listIndex)) {
                        continue;
                    }

                    Set<PairInt> set1 = segmentedCellList.get(listIndex1.intValue());

                    Colors colors1 = segmentedCellAvgLabColors.get(listIndex1);
                    if (colors1 == null) {
                        colors1 = imageProcessor.calculateAverageLAB(input, set1);
                        segmentedCellAvgLabColors.put(listIndex1, colors1);
                    }

                    double deltaE;
                    if (useDeltaE2000) {
                        deltaE = Math.abs(cieC.calcDeltaECIE2000(colors.getColors(),
                            colors1.getColors()));
                    } else {
                        deltaE = Math.abs(cieC.calcDeltaECIE94(colors.getColors(),
                            colors1.getColors()));
                    }

                    if (deltaE < minDeltaE) {
                        minDeltaE = deltaE;
                        minDeltaEIdx = listIndex1.intValue();
                    }
                }
            }
            if (minDeltaEIdx > -1) {
                Set<PairInt> set1 = segmentedCellList.get(minDeltaEIdx);
                group.addAll(set1);
                for (PairInt p1 : set1) {
                    pointIndexMap.put(p1, listIndex);
                }
                set1.clear();
            }
        }
        int count = 0;
        int nAssigned = 0;
        for (int i = (segmentedCellList.size() - 1); i > -1; --i) {
            Set<PairInt> set0 = segmentedCellList.get(i);
            if (set0.size() == 0) {
                segmentedCellList.remove(i);
                count++;
            } else {
                nAssigned += set0.size();
            }
        }

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, reassigned " + count + " small cells ");
        log.info(" nAssigned=" + nAssigned + " nUnassigned="
            + (input.getNPixels() - nAssigned));

    }

    private void assignedRemainingUnassigned(ImageExt input,
        List<Set<PairInt>> segmentedCellList,
        Map<PairInt, Integer> pointIndexMap, boolean useDeltaE2000,
        String debugTag) {

        long t0 = System.currentTimeMillis();

        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int count = 0;

        double deltaELimit = 8;

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                PairInt p = new PairInt(i, j);
                if (pointIndexMap.containsKey(p)) {
                    continue;
                }

                float[] lab = input.getCIELAB(i, j);

                double minDeltaE = Double.MAX_VALUE;
                int minDeltaEIdx = -1;

                for (int k = 0; k < dxs.length; ++k) {

                    int x1 = i + dxs[k];
                    int y1 = j + dys[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    Integer listIndex1 = pointIndexMap.get(new PairInt(x1, y1));
                    if (listIndex1 == null) {
                        continue;
                    }

                    float[] lab1 = input.getCIELAB(x1, y1);

                    double deltaE;
                    if (useDeltaE2000) {
                        // max value is ~19
                        deltaE = Math.abs(cieC.calcDeltaECIE2000(lab, lab1));
                    } else {
                        // max value is ~29
                        deltaE = Math.abs(cieC.calcDeltaECIE94(lab, lab1));
                    }

                    //double deltaL = Math.pow(Math.abs(lab[0] - lab1[0]), 2);
                    //deltaE += deltaL;

                    if ((deltaE < minDeltaE) && (deltaE < deltaELimit)) {
                        minDeltaE = deltaE;
                        minDeltaEIdx = listIndex1.intValue();
                    }
                }

                if (minDeltaEIdx > -1) {
                    segmentedCellList.get(minDeltaEIdx).add(p);
                    pointIndexMap.put(p, Integer.valueOf(minDeltaEIdx));
                    count++;
                }
            }
        }

        Set<PairInt> unassigned = createZerosSet(segmentedCellList, w, h,
            new HashSet<PairInt>());

        int count2 = 0;
        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(unassigned);
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            Set<PairInt> group = finder.getXY(i);
            count2 += group.size();
            Integer index = Integer.valueOf(segmentedCellList.size());
            for (PairInt p : group) {
                pointIndexMap.put(p, index);
            }
            segmentedCellList.add(group);
        }

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, assigned " + count + " unassigned points" +
            " and added new groups for " + count2 + " points");
    }

    private GreyscaleImage exploreCombinedColorDifferences(
        GreyscaleImage o1Img, GreyscaleImage bGImg, GreyscaleImage bRImg,
        String debugTag) {

        GreyscaleImage o1ImgCp = o1Img.copyImage();
        GreyscaleImage bGImgCp = bGImg.copyImage();
        GreyscaleImage bRImgCp = bRImg.copyImage();

        CannyEdgeFilterLite filter = new CannyEdgeFilterLite();
        filter.setToUseSobel();
        filter.applyFilter(o1ImgCp);
        filter = new CannyEdgeFilterLite();
        filter.setToUseSobel();
        filter.applyFilter(bGImgCp);
        filter = new CannyEdgeFilterLite();
        filter.setToUseSobel();
        filter.applyFilter(bRImgCp);

        GreyscaleImage combinedImg = o1Img.createWithDimensions();
        for (int i = 0; i < o1Img.getNPixels(); ++i) {
            int v = o1ImgCp.getValue(i) + bGImgCp.getValue(i) + bRImgCp.getValue(i);
            combinedImg.setValue(i, v);
        }

        HistogramEqualization histEq = new HistogramEqualization(combinedImg);
        histEq.applyFilter();

        //MiscDebug.writeImage(combinedImg, "_tmp_pre_thresh_" + debugTag);

        // --- making a combined image from thresholded at 127 and at 200 -----

        GreyscaleImage combinedImg200 = combinedImg.copyImage();
        for (int i = 0; i < combinedImg200.getNPixels(); ++i) {
            int v = combinedImg200.getValue(i);
            if (v <= 200) {
                combinedImg200.setValue(i, 0);
            } else {
                combinedImg200.setValue(i, 255);
            }
        }
        combinedImg200 = fillInGapsOf1(combinedImg200, new HashSet<PairInt>(), 255);
        removeIsolatedPixels(combinedImg200, 255, 0);

        return combinedImg200;
    }

    private GreyscaleImage edgesForLowContrastSmoothSkyColorDifferences(
        ImageExt img, GreyscaleImage o1Img, GreyscaleImage bGImg,
        GreyscaleImage bRImg, String debugTag) {

        long t0 = System.currentTimeMillis();

        GreyscaleImage greyImg = img.copyToGreyscale();

        GreyscaleImage combinedImg = greyImg.copyImage();

        CannyEdgeFilterLite filter = new CannyEdgeFilterLite();
        filter.setToUseSobel();
        filter.overrideHighThreshold(0.5f);
        filter.applyFilter(combinedImg);
        for (int i = 0; i < combinedImg.getNPixels(); ++i) {
            int v = combinedImg.getValue(i);
            if (v > 0) {
                combinedImg.setValue(i, 255);
            }
        }
        invertImage(combinedImg);
        removeEdgesSmallerThanLimit(combinedImg, 0, 255, 2);

        //MiscDebug.writeImage(combinedImg, "_tmp_color_low_contrast_2" + debugTag);

        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, first combined color gradient image");

        Set<PairInt> edgePoints = new HashSet<PairInt>();
        for (int col = 0; col < combinedImg.getWidth(); ++col) {
            for (int row = 0; row < combinedImg.getHeight(); ++row) {
                if (combinedImg.getValue(col, row) == 0) {
                    edgePoints.add(new PairInt(col, row));
                }
            }
        }

        t0 = System.currentTimeMillis();

        DistanceTransform dt = new DistanceTransform();
        int[][] dtTr = dt.applyMeijsterEtAl(edgePoints, combinedImg.getWidth(),
            combinedImg.getHeight());

        t1 = System.currentTimeMillis();
        t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, distance transform");

        t0 = System.currentTimeMillis();

        float[] values = new float[combinedImg.getNPixels()];
        int count = 0;
        for (int col = 0; col < combinedImg.getWidth(); ++col) {
            for (int row = 0; row < combinedImg.getHeight(); ++row) {
                values[count] = dtTr[col][row];
                count++;
            }
        }

        // if the histogram only has points at small x, then the spacings are
        // small and the gradient image is probably saturated, so need to create
        // the gradient with a threshold of 200 instead

        float max = 500.f;

        HistogramHolder hist = Histogram.createSimpleHistogram(0.f, max, 5.f,
            values, Errors.populateYErrorsBySqrt(values));
        /*
        try {
            for (int i = 0; i < hist.getYHistFloat().length; ++i) {
                float v = (float)Math.log(hist.getYHistFloat()[i]);
                hist.getYHistFloat()[i] = v;
                hist.getYHist()[i] = Math.round(v);
            }
            hist.plotHistogram("edge dist transform", debugTag);
        } catch (IOException e) {
        }
        */

        t1 = System.currentTimeMillis();
        t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, histogram created");

        int lastXIdx = MiscMath.findLastNonZeroIndex(hist);
        float lastX = (lastXIdx > -1) ? hist.getXHist()[lastXIdx] : Float.MAX_VALUE;

        // ------- analyze the distance transform, and if too noisey,
        //  create gradient image thresholded by value 200
        int avgDim = (o1Img.getWidth() + o1Img.getHeight())/2;
        double maxPossibleY = Math.log(o1Img.getNPixels());
        boolean createThresh200 = false;

        StringBuilder sb = new StringBuilder();
        sb.append(debugTag + " avgDim=" + avgDim + " maxPossibleY=" + maxPossibleY);
        sb.append(" lastX=" + lastX);
        if ((lastX < (avgDim/3)) || (avgDim <= 50)) {
            // the comparison should be proportional to image size. wanting to
            // see if points are clustered at small distances
            createThresh200 = true;
        } else {
            // for a 300x300 image, looking at slope between x=50 and x=250
            int y50 = MiscMath.findYforX(hist, 50);
            if (y50 < 1) {
                createThresh200 = true;
            } else {
                double deltaLogY = maxPossibleY - Math.log(y50);
                double f = deltaLogY/maxPossibleY;
                if (f > 0.55) {
                    createThresh200 = true;
                } else {
                    // look at the slope between 50 and 150.  if steep drop,
                    // then create thresh 200
                    int y150 = MiscMath.findYforX(hist, 150);
                    if (y150 < 1) {
                        createThresh200 = true;
                    } else {
                        double deltaLogY2 = Math.log(y50) - Math.log(y150);
                        double f2 = deltaLogY2/maxPossibleY;
                        if (f2 > 0.35) {
                            createThresh200 = true;
                        }
                        sb.append(" deltaLogY2=" + deltaLogY2 + " 0.35*max="
                            + 0.3*maxPossibleY + " f2=" + f2);
                    }
                }
                sb.append(" deltaLogY=" + deltaLogY + " 0.55*max="
                    + 0.55*maxPossibleY + " f=" + f);
            }
        }
        sb.append(" createThresh200=" + createThresh200);
        log.info(sb.toString());

        if (createThresh200) {

            t0 = System.currentTimeMillis();

            combinedImg = exploreCombinedColorDifferences(o1Img, bGImg, bRImg,
                debugTag);

            invertImage(combinedImg);

            // edges are now 0's

            //MiscDebug.writeImage(combinedImg, "_tmp_color_low_contrast_3" + debugTag);

            t1 = System.currentTimeMillis();
            t1Sec = (t1 - t0)/1000;
            log.info(t1Sec + " sec, thresh 200");
        }

        t0 = System.currentTimeMillis();

        filter = new CannyEdgeFilterLite();
        filter.setToUseSobel();
        filter.overrideHighThreshold(3.5f);
        filter.overrideLowThreshold(0.5f); //0.5
        filter.applyFilter(greyImg);
        for (int i = 0; i < greyImg.getNPixels(); ++i) {
            if (greyImg.getValue(i) > 0) {
                greyImg.setValue(i, 255);
            }
        }
        greyImg = fillInGapsOf1(greyImg, new HashSet<PairInt>(), 255);
        invertImage(greyImg);

        //MiscDebug.writeImage(greyImg, "_tmp_greyImg_" + debugTag);

        t1 = System.currentTimeMillis();
        t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, grey_low_contrast");

        //NOTE: this method needs to be improved.  if greyImg has many small contig regions, this takes along time to finish
        t0 = System.currentTimeMillis();
        //List<Set<PairInt>> perimeterLists = findPerimeters(greyImg, 255);
        List<Set<PairInt>> perimeterLists = new ArrayList<Set<PairInt>>();
        t1 = System.currentTimeMillis();
        t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, to find perimeters in greyImg");

        t0 = System.currentTimeMillis();

        addAdjacentEdges0(combinedImg, new GreyscaleImage[]{greyImg}, 0);
        //addAdjacentEdges(combinedImg, new GreyscaleImage[]{greyImg}, 0);

        t1 = System.currentTimeMillis();
        t1Sec = (t1 - t0)/1000;
        log.info(t1Sec + " sec, added adjacent edges from other images");

        return combinedImg;
    }

    private int[] countFractionZeros(GreyscaleImage img) {
        return countFraction(img, 0);
    }

    private int[] countFraction(GreyscaleImage img, int value) {

        List<Float> fraction = countFraction(img, value, 50, 50);
        // edges are 0's, so when many fractions are near 0.5 or higher,
        // should not use this image
        int nHigh2 = 0;
        int nHigh3 = 0;
        int nHigh4 = 0;
        for (Float fracZ : fraction) {
            if (fracZ.floatValue() >= 0.2f) {
                nHigh2++;
            }
            if (fracZ.floatValue() >= 0.3f) {
                nHigh3++;
            }
            if (fracZ.floatValue() >= 0.4f) {
                nHigh4++;
            }
        }

        return new int[]{nHigh2, nHigh3, nHigh4, fraction.size()};
    }

    private int[] countFractionNonZeros(GreyscaleImage img) {

        List<Float> fraction = countFractionNonZeros(img, 50, 50);
        // edges are 0's, so when many fractions are near 0.5 or higher,
        // should not use this image
        int nHigh2 = 0;
        int nHigh3 = 0;
        int nHigh4 = 0;
        for (Float fracZ : fraction) {
            if (fracZ.floatValue() >= 0.2f) {
                nHigh2++;
            }
            if (fracZ.floatValue() >= 0.3f) {
                nHigh3++;
            }
            if (fracZ.floatValue() >= 0.4f) {
                nHigh4++;
            }
        }

        return new int[]{nHigh2, nHigh3, nHigh4, fraction.size()};
    }

    private void removeIsolatedPixels(GreyscaleImage img, int pixValue,
        int pixNullValue) {

        removeIsolatedPixels(img, pixValue, pixNullValue, true);
    }

    private void removeIsolatedPixels(GreyscaleImage img, int pixValue,
        int pixNullValue, boolean use8Neighbors) {

        int[] dxs, dys;
        if (use8Neighbors) {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        } else {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int v = img.getValue(x, y);
                if (v != pixValue) {
                    continue;
                }
                int count = 0;
                int nSame = 0;
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) ||
                        (y2 > (h - 1))) {
                        continue;
                    }
                    count++;
                    int v2 = img.getValue(x2, y2);
                    if (v2 == pixValue) {
                        nSame++;
                        break;
                    }
                }
                if (nSame == 0) {
                    img.setValue(x, y, pixNullValue);
                }
            }
        }

    }

    /**
     * NOT READY FOR USE
     * @param input
     * @param greyGradient
     * @param mt
     * @param debugTag
     * @return 
     */
    public List<Set<PairInt>> performSegmentationWithColorEdges(ImageExt input,
        GreyscaleImage greyGradient, SegmentationMergeThreshold mt,
        String debugTag) {

        boolean fineDebug = true;

        int w = input.getWidth();
        int h = input.getHeight();
        
        Set<PairInt> mask = new HashSet<PairInt>();

        List<Set<PairInt>> segmentedCellList = findContiguousCells(0,
            greyGradient, mask);
 
        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_before_" + debugTag);
        }

        Map<PairInt, Integer> pointIndexMap;
        Set<PairInt> unassigned;
        boolean useDeltaE2000;
        double deltaELimit;
        boolean useAvgCellColor;

        unassigned = createZerosSet(segmentedCellList, w, h, mask);
        
        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }

// TODO: edit performSegmentationWithColorEdges
        
        useAvgCellColor = true;
        useDeltaE2000 = true;
        deltaELimit = 4;
        placeUnassignedByGrowingCells(input, segmentedCellList, unassigned,
            pointIndexMap, deltaELimit, useDeltaE2000, useAvgCellColor);
        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_0_" + debugTag);
        }

        deltaELimit = 10.0;//6.0
        mergeEmbeddedIfSimilar(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000);

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_1_" + debugTag);
        }

        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }

        assignedRemainingUnassigned(input, segmentedCellList, pointIndexMap,
            useDeltaE2000, debugTag);

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_2_" + debugTag);
        }

        int n2 = segmentedCellList.size();

        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
        
        deltaELimit = 0.25;
        mergeAdjacentIfSimilar2(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000, debugTag);

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_3_" + debugTag);
        }
        
        if (mt.equals(SegmentationMergeThreshold.EXTREMELY_LOW_CONTRAST)) {
            return segmentedCellList;
        }

        int n3 = segmentedCellList.size();
        float f23 = (float)n2/(float)n3;
        //log.info(debugTag + " n2-n3=" + (n2 - n3) + " n2/n3=" + f23);
        
        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
        int[] n2Andn8 = getEdgeProperties(input, segmentedCellList, 
            pointIndexMap, debugTag);
        float div = (n2Andn8[0] == 0) ? Float.MAX_VALUE : (float)n2Andn8[1]/(float)n2Andn8[0];
        log.info(debugTag + " n2DeltaE < 2 = " + n2Andn8[0] 
            + " n2DeltaE > 2 and < 8 = " + n2Andn8[1] + " div=" + div + 
            " n2-n3=" + (n2 - n3) + " n2/n3=" + f23);
        
        deltaELimit = 0.25;
        if ((n2Andn8[1] > 0) && ((n2Andn8[0]/n2Andn8[1]) > 3)) {
            deltaELimit = 1;
        }
        mergeAdjacentIfSimilar2(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000, debugTag);
        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_4_" + debugTag);
        }
        
        // -------------
        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
        n2Andn8 = getEdgeProperties(input, segmentedCellList, 
            pointIndexMap, debugTag);
        div = (n2Andn8[0] == 0) ? Float.MAX_VALUE : (float)n2Andn8[1]/(float)n2Andn8[0];
        log.info(debugTag + " n4DeltaE < 2 = " + n2Andn8[0] 
            + " n4DeltaE > 2 and < 8 = " + n2Andn8[1] + " div=" +div);
        
        //TODO: may need a setting for moderate contrast and then set deltaELimit=0.5 here
        deltaELimit = 0.25;//1.;
        mergeAdjacentIfSimilar(input, segmentedCellList, deltaELimit,
            useDeltaE2000, debugTag);

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_5_" + debugTag);
        }
        
        // -----------
        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
        deltaELimit = 5.0;//6.0
        mergeEmbeddedIfSimilar(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000);
        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_6_" + debugTag);
        }
        
        if (mt.equals(SegmentationMergeThreshold.DEFAULT)) {
            return segmentedCellList;
        }        

        int n5 = segmentedCellList.size();

        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
                    
        n2Andn8 = getEdgeProperties(input, segmentedCellList, 
            pointIndexMap, debugTag);
        div = (n2Andn8[0] == 0) ? Float.MAX_VALUE : (float)n2Andn8[1]/(float)n2Andn8[0];
        log.info(debugTag + " n7DeltaE < 2 = " + n2Andn8[0] 
            + " n7DeltaE > 2 and < 8 = " + n2Andn8[1] + " div=" + div);
        deltaELimit = 0.25;
        boolean lower = true;
        if ((n2Andn8[0] > n2Andn8[1]) && (div <= 0.5)) {
            lower = false;
            deltaELimit = 1;
        }

        mergeAdjacentIfSimilar(input, segmentedCellList, deltaELimit,
            useDeltaE2000, debugTag);

        int n6 = segmentedCellList.size();

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_7_" + debugTag);
        }

        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
        deltaELimit = 0.25;
        /*
        if ((n2Andn8[1] > 0) && ((n2Andn8[0]/n2Andn8[1]) > 3)) {
            deltaELimit = 1.;
        }
        */
        /*if (!lower) {
            
            n2Andn8 = getEdgeProperties(input, segmentedCellList, 
                pointIndexMap, debugTag);
            div = (n2Andn8[0] == 0) ? Float.MAX_VALUE : (float)n2Andn8[1]/(float)n2Andn8[0];
            log.info(debugTag + " n8DeltaE < 2 = " + n2Andn8[0] 
                + " n8DeltaE > 2 and < 8 = " + n2Andn8[1] + " div=" + div);
            
            if (n2Andn8[0] > n2Andn8[1]) {
                deltaELimit = 1.;
            } else if (div < 3) {
                deltaELimit = 1.; //0.75
                lower = false;
            }
        }*/
        
        mergeAdjacentIfSimilar2(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000, debugTag);

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_8_" + debugTag);
        }
        
        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
        deltaELimit = 5.0;//6.0
        /*if (!lower) {
            n2Andn8 = getEdgeProperties(input, segmentedCellList, 
                pointIndexMap, debugTag);
        
            log.info(debugTag + " n_DeltaE < 2 = " + n2Andn8[0] 
                + " n_DeltaE > 2 and < 8 = " + n2Andn8[1] + " div=" +
                (n2Andn8[1]/n2Andn8[0]));
            
            if (n2Andn8[0] > n2Andn8[1]) {
                deltaELimit = 6.0;
            } else if ((n2Andn8[1]/n2Andn8[0]) < 3) {
                deltaELimit = 6.0;
            }
        }*/      
        mergeEmbeddedIfSimilar(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000);
        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_9_" + debugTag);
        }

        pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            Integer key = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, key);
            }
        }
                    
        n2Andn8 = getEdgeProperties(input, segmentedCellList, 
            pointIndexMap, debugTag);
        div = (n2Andn8[0] == 0) ? Float.MAX_VALUE : (float)n2Andn8[1]/(float)n2Andn8[0];
        log.info(debugTag + " n9DeltaE < 2 = " + n2Andn8[0] 
            + " n9DeltaE > 2 and < 8 = " + n2Andn8[1] + " div=" + div);
        
        deltaELimit = 2.0;//2.5
        mergeAdjacentIfSimilar2(input, segmentedCellList, pointIndexMap,
            deltaELimit, useDeltaE2000, debugTag);

        if (fineDebug && debugTag != null && !debugTag.equals("")) {
            MiscDebug.writeAlternatingColor(input.copyImage(),
                segmentedCellList, "_tmp_10_" + debugTag);
        }

        return segmentedCellList;
    }

    private int findLargestContiguous(GreyscaleImage img, int value) {

        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);
        finder.setMinimumNumberInCluster(2);
        finder.findGroups(value);

        int maxN = Integer.MIN_VALUE;

        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            PairIntArray group = finder.getXY(i);
            if (group.getN() > maxN) {
                maxN = group.getN();
            }
        }

        return maxN;
    }

    /**
     * add to greyGradient, edges from addImages if the edges are adjacent
     * to an edge in greyGradient.
     * @param img
     * @param addImages
     * @param edgeValue
     */
    private void addAdjacentEdges(GreyscaleImage img, int imgEdgeValue,
        GreyscaleImage[] addImages, int edgeValue) {

        GreyscaleImage imgCp = img.copyImage();

        int nImages = addImages.length;

        DFSContiguousValueFinder[] finders = new DFSContiguousValueFinder[nImages];
        List<Map<PairInt, Integer>> pointListIndexes = new ArrayList<Map<PairInt, Integer>>();

        for (int i = 0; i < nImages; ++i) {

            DFSContiguousValueFinder finder = new DFSContiguousValueFinder(addImages[i]);
            finder.setMinimumNumberInCluster(1);
         finder.setToUse8Neighbors();
            finder.findGroups(edgeValue);
            finders[i] = finder;

            Map<PairInt, Integer> pointIndexes = new HashMap<PairInt, Integer>();
            for (int j = 0; j < finder.getNumberOfGroups(); ++j) {
                PairIntArray group = finder.getXY(j);
                Integer key = Integer.valueOf(j);
                for (int ii = 0; ii < group.getN(); ++ii) {
                    PairInt p = new PairInt(group.getX(ii), group.getY(ii));
                    pointIndexes.put(p, key);
                }
            }
            pointListIndexes.add(pointIndexes);
        }

        int w = imgCp.getWidth();
        int h = imgCp.getHeight();

        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                
                int v = imgCp.getValue(x, y);
                if (v != imgEdgeValue) {
                    continue;
                }

                for (int dx = -1; dx <= +1; ++dx) {
                    int x2 = x + dx;
                    if ((x2 < 0) || (x2 > (w - 1))) {
                        continue;
                    }
                    for (int dy = -1; dy <= +1; ++dy) {
                        int y2 = y + dy;
                        if ((y2 < 0) || (y2 > (h - 1))) {
                            continue;
                        }
                        PairInt p2 = new PairInt(x2, y2);

                        // if p2 is in any of the point index lists,
                        // add that edge to greyGradient and remove it from
                        // the maps to make next traversal faster

                        for (int ii = 0; ii < nImages; ++ii) {
                            Map<PairInt, Integer> pointIndexes = pointListIndexes.get(ii);
                            Integer listIndex = pointIndexes.get(p2);
                            if (listIndex != null) {
                                PairIntArray edge2 = finders[ii].getXY(listIndex.intValue());
                                for (int jj = 0; jj < edge2.getN(); ++jj) {
                                    img.setValue(edge2.getX(jj),
                                        edge2.getY(jj), imgEdgeValue);
                                }
                                pointIndexes.remove(p2);
                            }
                        }
                    }
                }
            }
        }
    }

    private void addAdjacentEdges0(GreyscaleImage img,
        GreyscaleImage[] addImages, int edgeValue) {

        int nImages = addImages.length;

        int w = img.getWidth();
        int h = img.getHeight();

        for (int ii = 0; ii < nImages; ++ii) {
            GreyscaleImage img2 = addImages[ii];
            for (int x = 0; x < w; ++x) {
                for (int y = 0; y < h; ++y) {
                    if (img2.getValue(x, y) == edgeValue) {
                        img.setValue(x, y, edgeValue);
                    }
                }
            }
        }
    }

    private List<Set<PairInt>> findPerimeters(GreyscaleImage greyImg,
        int edgeValue) {

        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        int imageMaxColumn = greyImg.getWidth() - 1;
        int imageMaxRow = greyImg.getHeight() - 1;
        int[] rowMinMax = new int[2];

        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(greyImg);
        finder.setMinimumNumberInCluster(1);
        finder.setToUse8Neighbors();
        finder.findGroups(edgeValue);

        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {

            PairIntArray edge = finder.getXY(i);
            Set<PairInt> points = Misc.convert(edge);

            Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();

            Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
                points, rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);

            if (!outputEmbeddedGapPoints.isEmpty()) {
                // update the perimeter for "filling in" embedded points
                perimeterFinder.updateRowColRangesForAddedPoints(rowColRanges,
                    rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
            }

            Set<PairInt> perimeterPoints = perimeterFinder.getBorderPixels0(
                rowColRanges, rowMinMax, imageMaxColumn, imageMaxRow);

            output.add(perimeterPoints);
        }

        return output;
    }

    private void plotHistograms(ImageExt input,
        List<Set<PairInt>> segmentedCellList, String debugTag) {

        int n = segmentedCellList.size();

        ImageProcessor imageProcessor = new ImageProcessor();

        float[] values = new float[n];

        for (int i = 0; i < n; ++i) {

            float ha = imageProcessor.calculateAverageHueAngle(input,
                segmentedCellList.get(i));

            values[i] = ha;
        }

        HistogramHolder hist = Histogram.createSimpleHistogram(0.f, 360.f, 361,
            values, Errors.populateYErrorsBySqrt(values));

        try {
            hist.plotHistogram(debugTag, debugTag);
        } catch (IOException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private int[] getEdgeProperties(ImageExt input, 
        List<Set<PairInt>> segmentedCellList, 
        Map<PairInt, Integer> pointIndexMap, String debugTag) {
        
        /*
        for adjacent segmented cells, determining these properties either on
        a point by point basis of the adjacent points, or using the average of the
        entire cell for adjacent cells:
           -- deltaE
           -- hueAngle1 vs hueAngle2
           -- (360 - hueAngle1) vs hueAngle2
           -- hueAngle1 vs (360 - hueAngle2)
        */
        List<Double> deltaE = new ArrayList<Double>();
        List<Integer> hueAngle1 = new ArrayList<Integer>();
        List<Integer> hueAngle2 = new ArrayList<Integer>();
        
        boolean useCellAverages = false;//true;
        
        if (useCellAverages) {
            populateAdjacentCellAverages(input, segmentedCellList, pointIndexMap,
                deltaE, hueAngle1, hueAngle2);
        } else {
            populateAdjacentCellPoints(input, segmentedCellList, pointIndexMap,
                deltaE, hueAngle1, hueAngle2);
        }
        
        int nDELT2 = 0;
        int nDELT28 = 0;

        //---------
        int n = deltaE.size();
        for (int i = 0; i < n; ++i) {
            float x = deltaE.get(i).floatValue();
            if (x <= 2.) {
                nDELT2++;
            } else if (x <= 8.) {
                nDELT28++;
            }
        }
        nDELT28 /= 3;
        return new int[]{nDELT2, nDELT28};
        /*
        try { 
            float[] xPolygon = null;
            float[] yPolygon = null;
            int n = deltaE.size();
            float[] x = new float[n];
            float[] y = new float[n];
            for (int i = 0; i < n; ++i) {
                x[i] = deltaE.get(i).floatValue();
                y[i] = hueAngle1.get(i).floatValue();
                if (x[i] <= 2.) {
                    nDELT2++;
                } else if (x[i] <= 8.) {
                    nDELT28++;
                }
            }
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            float minX = MiscMath.findMin(x);
            float maxX = MiscMath.findMax(x);
            plotter.addPlot(minX, maxX, 0.f, 361.f, x, y,xPolygon, yPolygon, 
                "dE vs ha1 ");
            //------
            plotter.writeFile(debugTag);
            
            log.info(debugTag + " nDeltaE < 2 = " + nDELT2 
                + " nDeltaE > 2 and < 8 = " + nDELT28);
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        */
    }

    private void populateAdjacentCellPoints(ImageExt input, 
        List<Set<PairInt>> segmentedCellList, 
        Map<PairInt, Integer> pointIndexMap, 
        List<Double> deltaE, List<Integer> hueAngle1, List<Integer> hueAngle2) {
        
        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        CIEChromaticity cieC = new CIEChromaticity();

        for (int i = 0; i < segmentedCellList.size(); ++i) {

            Integer index = Integer.valueOf(i);

            Set<PairInt> set = segmentedCellList.get(index.intValue());

            // collecting bordering points
            // storing all then averaging and comparing to deltaELimit
            Map<Integer, List<Double>> listIndexDeltaEMap = new HashMap<Integer, List<Double>>();
            Map<Integer, List<Integer>> listIndexHA1Map = new HashMap<Integer, List<Integer>>();
            Map<Integer, List<Integer>> listIndexHA2Map = new HashMap<Integer, List<Integer>>();

            for (PairInt p : set) {

                int x = p.getX();
                int y = p.getY();
                Integer listIndex = pointIndexMap.get(p);
                assert(listIndex.intValue() == index.intValue());

                float[] lab = input.getCIELAB(x, y);

                for (int k = 0; k < dxs.length; ++k) {

                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];

                    PairInt p2 = new PairInt(x2, y2);
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    Integer listIndex2 = pointIndexMap.get(p2);
                    if (listIndex2 == null || listIndex.equals(listIndex2)) {
                        continue;
                    }

                    float[] lab2 = input.getCIELAB(x2, y2);

                    double dE = Math.abs(cieC.calcDeltaECIE2000(lab, lab2));
                    
                    List<Double> dEs = listIndexDeltaEMap.get(listIndex2);
                    if (dEs == null) {
                        dEs = new ArrayList<Double>();
                        listIndexDeltaEMap.put(listIndex2, dEs);
                    }
                    dEs.add(Double.valueOf(dE));
                    
                    float ha1;
                    if (lab[1] == 0) {
                        ha1 = 0;
                    } else {
                        ha1 = (float) (Math.atan2(lab[2], lab[1]) * 180. / Math.PI);
                        if (ha1 < 0) {
                            ha1 += 360.;
                        }
                    }
                    
                    List<Integer> ha1s = listIndexHA1Map.get(listIndex2);
                    if (ha1s == null) {
                        ha1s = new ArrayList<Integer>();
                        listIndexHA1Map.put(listIndex2, ha1s);
                    }
                    ha1s.add(Integer.valueOf(Math.round(ha1)));
                    
                    float ha2;
                    if (lab2[1] == 0) {
                        ha2 = 0;
                    } else {
                        ha2 = (float) (Math.atan2(lab2[2], lab2[1]) * 180. / Math.PI);
                        if (ha2 < 0) {
                            ha2 += 360.;
                        }
                    }
                    
                    List<Integer> ha2s = listIndexHA2Map.get(listIndex2);
                    if (ha2s == null) {
                        ha2s = new ArrayList<Integer>();
                        listIndexHA2Map.put(listIndex2, ha2s);
                    }
                    ha2s.add(Integer.valueOf(Math.round(ha2)));
                }
            }

            for (Entry<Integer, List<Double>> entry : listIndexDeltaEMap.entrySet()) {

                Integer listIndex2 = entry.getKey();

                List<Double> deltaEs = entry.getValue();

                double sumDeltaE = 0;
                for (int ii = 0; ii < deltaEs.size(); ++ii) {
                    sumDeltaE += deltaEs.get(ii).doubleValue();
                }
                sumDeltaE /= (double)deltaEs.size();
                
                List<Integer> ha1s = listIndexHA1Map.get(listIndex2.intValue());
                int[] hueAngles1 = new int[ha1s.size()];
                List<Integer> ha2s = listIndexHA2Map.get(listIndex2.intValue());
                int[] hueAngles2 = new int[ha2s.size()];
                for (int ii = 0; ii < ha1s.size(); ++ii) {
                    hueAngles1[ii] = ha1s.get(ii).intValue();
                    hueAngles2[ii] = ha2s.get(ii).intValue();
                }
                
                float avgHA1 = AngleUtil.calculateAverageWithQuadrantCorrections(
                    hueAngles1, hueAngles1.length - 1);
                
                float avgHA2 = AngleUtil.calculateAverageWithQuadrantCorrections(
                    hueAngles2, hueAngles2.length - 1);
                
                if (avgHA1 < 0) {
                    avgHA1 += 360;
                } else if (avgHA1 > 359) {
                    avgHA1 -= 360;
                }
                if (avgHA2 < 0) {
                    avgHA2 += 360;
                } else if (avgHA2 > 359) {
                    avgHA2 -= 360;
                }
                
                deltaE.add(Double.valueOf(sumDeltaE));
                hueAngle1.add(Integer.valueOf(Math.round(avgHA1)));
                hueAngle2.add(Integer.valueOf(Math.round(avgHA2)));
            }
        }

    }
    
    private void populateAdjacentCellAverages(ImageExt input, 
        List<Set<PairInt>> segmentedCellList, 
        Map<PairInt, Integer> pointIndexMap, 
        List<Double> deltaE, List<Integer> hueAngle1, List<Integer> hueAngle2) {
        
        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Set<PairInt> addedPairIndexes = new HashSet<PairInt>();

        Map<Integer, Colors> segmentedCellAvgLabColors 
            = new HashMap<Integer, Colors>();
        
        ImageProcessor imageProcessor = new ImageProcessor();
         
        CIEChromaticity cieC = new CIEChromaticity();

        for (int i = 0; i < segmentedCellList.size(); ++i) {

            Integer index = Integer.valueOf(i);

            Set<PairInt> set = segmentedCellList.get(index.intValue());
            
            Set<Integer> list2Indexes = new HashSet<Integer>();

            for (PairInt p : set) {

                int x = p.getX();
                int y = p.getY();
                Integer listIndex = pointIndexMap.get(p);
                assert(listIndex.intValue() == index.intValue());

                for (int k = 0; k < dxs.length; ++k) {

                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];

                    PairInt p2 = new PairInt(x2, y2);
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }

                    Integer listIndex2 = pointIndexMap.get(p2);
                    if (listIndex2 == null || listIndex.equals(listIndex2)) {
                        continue;
                    }

                    list2Indexes.add(listIndex2);
                }
            }
            
            for (Integer index2 : list2Indexes) {
                int idx1, idx2;
                if (index.intValue() < index2.intValue()) {
                    idx1 = index.intValue();
                    idx2 = index2.intValue();
                } else {
                    idx2 = index.intValue();
                    idx1 = index2.intValue();
                }
                PairInt p12 = new PairInt(idx1, idx2);
                if (addedPairIndexes.contains(p12)) {
                    continue;
                }
                addedPairIndexes.add(p12);
                
                Colors colors1 = segmentedCellAvgLabColors.get(Integer.valueOf(idx1));
                if (colors1 == null) {
                    colors1 = imageProcessor.calculateAverageLAB(input,
                        segmentedCellList.get(idx1));
                    segmentedCellAvgLabColors.put(Integer.valueOf(idx1), colors1);
                }
                float[] lab1 = colors1.getColors();
                
                Colors colors2 = segmentedCellAvgLabColors.get(Integer.valueOf(idx2));
                if (colors2 == null) {
                    colors2 = imageProcessor.calculateAverageLAB(input,
                        segmentedCellList.get(idx2));
                    segmentedCellAvgLabColors.put(Integer.valueOf(idx2), colors2);
                }
                float[] lab2 = colors2.getColors();
                
                double dE = Math.abs(cieC.calcDeltaECIE2000(lab1, lab2));
                    
                float ha1;
                if (lab1[1] == 0) {
                    ha1 = 0;
                } else {
                    ha1 = (float) (Math.atan2(lab1[2], lab1[1]) * 180. / Math.PI);
                    if (ha1 < 0) {
                        ha1 += 360.;
                    }
                }
                
                float ha2;
                if (lab2[1] == 0) {
                    ha2 = 0;
                } else {
                    ha2 = (float) (Math.atan2(lab2[2], lab2[1]) * 180. / Math.PI);
                    if (ha2 < 0) {
                        ha2 += 360.;
                    }
                }
                
                deltaE.add(Double.valueOf(dE));
                hueAngle1.add(Integer.valueOf(Math.round(ha1)));
                hueAngle2.add(Integer.valueOf(Math.round(ha2)));
            }
        }
    }

    private GreyscaleImage exploreCombiningImages(GreyscaleImage o1Img, 
        GreyscaleImage labAImg, GreyscaleImage labBImg, 
        GreyscaleImage greyGradient, String debugTag) {
        
        //NOTE:  not keeping this as is... loses disconnected edges not in
        //   the base of intesection2 or intersection3
        
        o1Img = o1Img.copyImage();
        labAImg = labAImg.copyImage();
        labBImg = labBImg.copyImage();
        greyGradient = greyGradient.copyImage();
        
        // labA and labB have already been scaled.  greyGradient too, differently
        HistogramEqualization hEq = new HistogramEqualization(o1Img);
        hEq.applyFilter();
      
        CannyEdgeFilterLite cannyFilter = new CannyEdgeFilterLite();
        cannyFilter.setToUseSobel();
        cannyFilter.applyFilter(o1Img);
        
        cannyFilter = new CannyEdgeFilterLite();
        cannyFilter.setToUseSobel();
        cannyFilter.applyFilter(labAImg);
        
        cannyFilter = new CannyEdgeFilterLite();
        cannyFilter.setToUseSobel();
        cannyFilter.applyFilter(labBImg);
        
        
        //ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.blur(o1Img, SIGMA.ONE);
        //imageProcessor.blur(labAImg, SIGMA.ONE);
        //imageProcessor.blur(labBImg, SIGMA.ONE);
        
        //MiscDebug.writeImage(o1Img, "_canny_o1_" + debugTag);
        //MiscDebug.writeImage(labAImg, "_canny_labA_" + debugTag);
        //MiscDebug.writeImage(labBImg, "_canny_labB_" + debugTag);
        
        int n = greyGradient.getNPixels();
        int w = greyGradient.getWidth();
        int h = greyGradient.getHeight();
        
        int n2 = 0;
        int n4 = 0;
        GreyscaleImage intersection4 = new GreyscaleImage(w, h);
        for (int i = 0; i < n; ++i) {
            if ((o1Img.getValue(i) == 0) || (labAImg.getValue(i) == 0) ||
                (labBImg.getValue(i) == 0) || (greyGradient.getValue(i) == 255)) {
                continue;
            }
            intersection4.setValue(i, 255);
            n4++;
        }
        GreyscaleImage intersection2 = new GreyscaleImage(w, h);
        for (int i = 0; i < n; ++i) {
            if ((labAImg.getValue(i)== 0) || (greyGradient.getValue(i) == 255)) {
                continue;
            }
            intersection2.setValue(i, 255);
            n2++;
        }
        GreyscaleImage intersection3 = new GreyscaleImage(w, h);
        for (int i = 0; i < n; ++i) {
            if ((labAImg.getValue(i)== 0) || (labBImg.getValue(i)== 0)) {
                continue;
            }
            intersection3.setValue(i, 255);
            n2++;
        }
        int nI24 = 0;
        for (int i = 0; i < n; ++i) {
            if ((intersection2.getValue(i)== 0) || (intersection4.getValue(i) == 0)) {
                continue;
            }
            nI24++;
        }
        
        MiscDebug.writeImage(intersection4, "_intersection4_" + debugTag);
        MiscDebug.writeImage(intersection3, "_intersection3_" + debugTag);
        MiscDebug.writeImage(intersection2, "_intersection2_" + debugTag);
        
        intersection3 = fillInGapsOf1(intersection3, new HashSet<PairInt>(), 255);
        MiscDebug.writeImage(intersection3, "_intersection3_ext" + debugTag);
        
        float n2F = (float)n2/(float)n;
        float n4F = (float)n4/(float)n;
        float n4F2 = (float)n4/(float)n2;
        float nI24FI4 = (float)nI24/(float)n4;
        float nI24FI2 = (float)nI24/(float)n2;
        log.info(debugTag + " n2F=" + n2F + " n4F=" + n4F + " nFdivN2=" 
            + n4F2 +" n=" + n);
        
        boolean choseInter4 = false;
        GreyscaleImage baseImg;
        if (n4F < 0.0075) {
            baseImg = intersection2;
        } else {
            baseImg = intersection4;
            choseInter4 = true;
        }
        MiscDebug.writeImage(baseImg, "_combined_pre_" + debugTag);

        addAdjacent(baseImg, 255, greyGradient, 0);
        
        MiscDebug.writeImage(baseImg, "_combined_pre_2_" + debugTag);
        
        if (choseInter4) {
            //setAllNonZeroTo255(labAImg);
            //addAdjacent(intersection2, 255, labAImg, 255);
            //MiscDebug.writeImage(intersection2, "_canny_labA_ext_" + debugTag);
            
            addAdjacent(baseImg, 255, intersection2, 255);
            MiscDebug.writeImage(baseImg, "_combined_pre_3_" + debugTag);
            addAdjacent(baseImg, 255, intersection3, 255);
            MiscDebug.writeImage(baseImg, "_combined_pre_4_" + debugTag);
            
            //setAllNonZeroTo255(o1Img);
            //addAdjacent(baseImg, 255, o1Img, 255);
            
            //MiscDebug.writeImage(baseImg, "_combined_pre_5_" + debugTag);
            //baseImg = fillInGapsOf1(baseImg, new HashSet<PairInt>(), 255);
        }
                
        /*
        //ImageProcessor imageProcessor = new ImageProcessor();
        
        // build a stack from baseImg and then add connected from all 4 images
        Stack<Integer> stack = new Stack<Integer>();
        for (int i = 0; i < n; ++i) {
            if (baseImg.getValue(i) == 255) {
                stack.add(Integer.valueOf(i));
            }
        }
        
        Set<Integer> visited = new HashSet<Integer>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        while (!stack.isEmpty()) {
            Integer pixIndex = stack.pop();
            if (visited.contains(pixIndex)) {
                continue;
            }
            int pixIdx = pixIndex.intValue();
            int x = baseImg.getCol(pixIdx);
            int y = baseImg.getRow(pixIdx);

            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                    continue;
                }
                
                int pixIdx2 = baseImg.getIndex(x2, y2);
                if (baseImg.getValue(pixIdx2) == 255) {
                    continue;
                }
                if (
                    (o1Img.getValue(pixIdx2) > 0) || 
                    (labAImg.getValue(pixIdx2) > 0) ||
                    (labBImg.getValue(pixIdx2) > 0)) {
                    stack.add(Integer.valueOf(pixIdx2));
                    baseImg.setValue(pixIdx2, 255);
                }
            }
            
            visited.add(pixIndex);
        }
        */
        //imageProcessor.applyAdaptiveMeanThresholding(baseImg, 1);
        
        //baseImg = fillInGapsOf1(baseImg, new HashSet<PairInt>(), 255);
                
        MiscDebug.writeImage(baseImg, "_combined_" + debugTag);
        
        return baseImg;
                
    }

    /**
     * add perimeters of largest contiguous non-edge values in addToImg to
     * img.
     * @param img
     * @param edgeValue
     * @param addToImage
     * @param nonEdgeValue 
     */
    private void addLargestContiguous(GreyscaleImage img, int edgeValue, 
        GreyscaleImage addToImg, int nonEdgeValue) {
        
        int nPix = addToImg.getNPixels();
        int minSize = (int)Math.round(0.02 * nPix);
        
        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(addToImg);
        finder.findGroups(nonEdgeValue);
        
        Set<PairInt> contig = new HashSet<PairInt>();
        
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            PairIntArray group = finder.getXY(i);
            if (group.getN() > minSize) {
                for (int j = 0; j < group.getN(); ++j) {
                    contig.add(new PairInt(group.getX(j), group.getY(j)));
                }
            }
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        int w = addToImg.getWidth();
        int h = addToImg.getHeight();
        
        for (PairInt p : contig) {
            int x = p.getX();
            int y = p.getY();
            // any neighbor that is not a non-edge value can be added to img
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                    continue;
                }
                if (addToImg.getValue(x2, y2) != nonEdgeValue) {
                    img.setValue(x2, y2, edgeValue);
                }
            }
        }
    }

    private void addAdjacent(GreyscaleImage img, int edgeValue, 
        GreyscaleImage addToImg, int addToImgEdgeValue) {
        
        int n = img.getNPixels();
        int w = img.getWidth();
        int h = img.getHeight();
        
        // build a stack from baseImg and then add connected from all 4 images
        Stack<Integer> stack = new Stack<Integer>();
        for (int i = 0; i < n; ++i) {
            if (img.getValue(i) == edgeValue) {
                stack.add(Integer.valueOf(i));
            }
        }
        
        Set<Integer> visited = new HashSet<Integer>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        while (!stack.isEmpty()) {
            Integer pixIndex = stack.pop();
            if (visited.contains(pixIndex)) {
                continue;
            }
            int pixIdx = pixIndex.intValue();
            int x = img.getCol(pixIdx);
            int y = img.getRow(pixIdx);

            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                    continue;
                }
                
                int pixIdx2 = img.getIndex(x2, y2);
                if (img.getValue(pixIdx2) == edgeValue) {
                    continue;
                }
                if (addToImg.getValue(pixIdx2) == addToImgEdgeValue) {
                    stack.add(Integer.valueOf(pixIdx2));
                    img.setValue(pixIdx2, edgeValue);
                }
            }
            
            visited.add(pixIndex);
        }
    }

    private boolean neighborIsInSet(Set<PairInt> set, PairInt p) {
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        for (int i = 0; i < dxs.length; ++i) {
            int x2 = p.getX() + dxs[i];
            int y2 = p.getY() + dys[i];
            PairInt p2 = new PairInt(x2, y2);
            if (set.contains(p2)) {
                return true;
            }
        }
        return false;
    }

    /**
     * check whether a neighbor belongs to another edge (NOTE that the bounds
     * have to have been checked before this)
     * @param x
     * @param y
     * @param edgeIndexMap
     * @param hN
     * @return 
     */
    private boolean foundAdjacentEdge(int x, int y, Map<PairInt, Integer> 
        edgeIndexMap, Integer index, int hN) {
        
        for (int dy = -hN; dy <= hN; ++dy) {
            int y2 = y + dy;
            for (int dx = -hN; dx <= hN; ++dx) {
                int x2 = x + dx;
                PairInt p2 = new PairInt(x2, y2);
                Integer index2 = edgeIndexMap.get(p2);
                if ((index2 != null) && !index2.equals(index)) {
                    return true;
                }
            }
        }
        return false;
    }
    
    private boolean aMemberIsOutOfBounds(int x, int y, int hN, int w, int h) {
        
        for (int dy = -hN; dy <= hN; ++dy) {
            int y2 = y + dy;
            if ((y2 < 0) || (y2 > (h - 1))) {
                return true;
            }
            for (int dx = -hN; dx <= hN; ++dx) {
                int x2 = x + dx;
                if ((x2 < 0) || (x2 > (w - 1))) {
                    return true;
                }
            }
        }
        
        return false;
    }

    /**
     * 
     * @param img
     * @param clrIdx expected to be 0 or 1 where 0 indicates to calculate for O1
     * and 1 indicates to calculate for labA.
     * @param x
     * @param y
     * @param sparseValueMap
     * @param orientationMap
     * @param gradientXMap
     * @param gradientYMap
     * @param gradientMap
     * @param hN larger neighbor region half width which covers at least the
     * half length of the kernel.
     * @param hN2 smaller neighbor region half width which may cover less than
     * hN, but still needs to cover the half kernel.
     * @param kernel 
     */
    private void calculateGradientAndOrientation(
        ImageExt img, int clrIdx,
        int x, int y, 
        Map<PairInt, Float> sparseValueMap, 
        Map<PairInt, Integer> orientationMap, 
        Map<PairInt, Float> gradientXMap, 
        Map<PairInt, Float> gradientYMap, 
        Map<PairInt, Float> gradientMap, int hN, int hN2, float[] kernel) {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // calculate the labA and o1 associated values for the surrounding
        // hN x hN region
        for (int dy = -hN; dy <= hN; ++dy) {
            int y2 = y + dy;
            for (int dx = -hN; dx <= hN; ++dx) {
                int x2 = x + dx;
                PairInt p2 = new PairInt(x2, y2);
                if (sparseValueMap.containsKey(p2)) {
                    continue;
                }
                float value;
                if (clrIdx == 0) {
                    value = img.getR(x2, y2) - img.getG(x2, y2);
                } else {
                    value = img.getCIELAB(x2, y2)[1];
                }
                sparseValueMap.put(p2, Float.valueOf(value));
            }
        }

        boolean calcForX = true;
        for (int dy = -hN2; dy <= hN2; ++dy) {
            int y2 = y + dy;
            for (int dx = -hN2; dx <= hN2; ++dx) {
                int x2 = x + dx;
                PairInt p2 = new PairInt(x2, y2);
                if (gradientXMap.containsKey(p2)) {
                    continue;
                }
                float convolvedLabA = imageProcessor.convolve1D(sparseValueMap,
                    x2, y2, kernel, calcForX);
                gradientXMap.put(p2, Float.valueOf(convolvedLabA));
            }
        }
        calcForX = false;
        for (int dy = -hN2; dy <= hN2; ++dy) {
            int y2 = y + dy;
            //sb.append(String.format("         %d : ", y2));
            for (int dx = -hN2; dx <= hN2; ++dx) {
                int x2 = x + dx;
                PairInt p2 = new PairInt(x2, y2);
                if (orientationMap.containsKey(p2)) {
                    continue;
                }
                float convolvedY = imageProcessor.convolve1D(
                    sparseValueMap, x2, y2, kernel, calcForX);
                gradientYMap.put(p2, Float.valueOf(convolvedY));
                float convolvedX = gradientXMap.get(p2).intValue();
                double gradient = Math.sqrt(convolvedX * convolvedX
                    + convolvedY * convolvedY);
                gradientMap.put(p2, Float.valueOf((float) gradient));

                double orientation = Math.atan2(convolvedY, convolvedX);
                if (orientation < 0) {
                    orientation += Math.PI;
                }
                orientation *= (180. / Math.PI);

                int or = Integer.valueOf((int) Math.round(orientation));
                orientationMap.put(p2, or);
            }
        }
    }

    private int countNeighbors(Map<PairInt, Integer> pointMap, int x, int y) {
        
        int count = 0;
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
                
        for (int jj = 0; jj < dxs.length; ++jj) {
            PairInt p3 = new PairInt(x + dxs[jj], y + dys[jj]);
            if (pointMap.containsKey(p3)) {
                count++;
            }
        }
        
        return count;
    }

    private void mergeEdges(List<Set<PairInt>> clusterPoints, 
        float[][] clusterDescriptors, int clrSpace, double tColor,
        List<Integer> edgeIndexes) {
        
        System.out.println(edgeIndexes.size() + " edges");
        
        final long heapKeyFactor = 1000000l;
        Heap heap = new Heap();  
        Map<PairInt, HeapNode> pairEdgePindexNodes = new HashMap<PairInt, HeapNode>();
        populateColorDiffHeap(clusterDescriptors, clrSpace,
            edgeIndexes, heap, heapKeyFactor, pairEdgePindexNodes);

        // ---- make a map to find and update merged data structures ------
        Map<Integer, Set<Integer>> indexToIndexMap = new HashMap<Integer, Set<Integer>>();
        for (PairInt p : pairEdgePindexNodes.keySet()) {
            
            Integer index1 = Integer.valueOf(p.getX());
            Integer index2 = Integer.valueOf(p.getY());
            
            Set<Integer> indexes = indexToIndexMap.get(index1);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                indexToIndexMap.put(index1, indexes);
            }
            indexes.add(index2);
            
            indexes = indexToIndexMap.get(index2);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                indexToIndexMap.put(index2, indexes);
            }
            indexes.add(index1);
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        while (!heap.isEmpty()) {
            
            HeapNode node = heap.extractMin();
            double diff = ((double)node.getKey())/((double)heapKeyFactor);
            
            if (diff > tColor) {
                break;
            }
            
            PairInt p12 = (PairInt)node.getData();
            
            int idx1 = p12.getX();
            int idx2 = p12.getY();
            
            Set<PairInt> set1 = clusterPoints.get(idx1);
            Set<PairInt> set2 = clusterPoints.get(idx2);
            
            if (set1.isEmpty() || set2.isEmpty()) {
                continue;
            }
            
            if (set2.size() > set1.size()) {
                idx1 = p12.getY();
                idx2 = p12.getX();
                set1 = set2;
                set2 = clusterPoints.get(idx2);
            }
            
            // set1 is largest
            
            float[] desc1 = clusterDescriptors[idx1];
            float[] desc2 = clusterDescriptors[idx2];
            float n1 = set1.size();
            float n2 = set2.size();
            float nTot = n1 + n2;
           
            assert(Math.abs(n1 - desc1[3]) < 0.1);
            assert(Math.abs(n2 - desc2[3]) < 0.1);
            
            //{h, s, v, nPix, cenX, cenY}
            // update desc1 contents for contents in desc2
            for (int k = 0; k < desc1.length; ++k) {
                if (k == 3) {
                    desc1[k] = nTot;
                } else {
                    desc1[k] = ((desc1[k] * n1) + (desc2[k] * n2)) / nTot;
                }
            }
            clusterDescriptors[idx2] = null;
            set1.addAll(set2);
            set2.clear();
            
            n1 = set1.size();
                        
            // update all pairs having set 1
            Set<Integer> indexes2 = indexToIndexMap.get(Integer.valueOf(idx1));
            for (Integer index2 : indexes2) {
                int idx3 = index2.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);            
                if (set3.isEmpty()) {
                    continue;
                }
                //keys in pairEdgePindexNodes have smaller index in x
                PairInt p13;
                if (idx1 < idx3) {
                    p13 = new PairInt(idx1, idx3);
                } else {
                    p13 = new PairInt(idx3, idx1);
                }
                
                float[] desc3 = clusterDescriptors[idx3];
               
                HeapNode node3 = pairEdgePindexNodes.get(p13);
                assert(node3 != null);
                
                heap.remove(node3);
            
                double diffUpdated;
                if (clrSpace == 0) {
                    diffUpdated = Math.abs(cieC.calcDeltaECIE2000(
                        desc1[0], desc1[1], desc1[2], 
                        desc3[0], desc3[1], desc3[2]));
                } else {
                    double diff1 = desc1[0] - desc3[0];
                    double diff2 = desc1[1] - desc3[1];
                    double diff3 = desc1[2] - desc3[2];
                    diffUpdated = Math.sqrt(diff1 * diff1 + diff2*diff2 + diff3*diff3);
                }
                
                long heapKey = (long)((double)heapKeyFactor * diffUpdated);
                node3 = new HeapNode(heapKey);
                node3.setData(p13);
                heap.insert(node3);
                pairEdgePindexNodes.put(p13, node3);                
            }
            
            // pairs having set2 will be skipped because of the empty set at beginning of while loop
        }
        
        {
            //DEBUG
            int nEdges = 0;
            for (Integer index : edgeIndexes) {
                Set<PairInt> set = clusterPoints.get(index.intValue());
                if (set.size() > 0) {
                    nEdges++;
                }
            }
            System.out.println(nEdges + " edges after merge edges");
        }
    }

    /**
     * populate doubly linked circular list with nodes holding the differences
     * between pairs of edge descriptor colors, ordered by increasing difference.
     * The nodes in the list have a data field which contains the indexes
     * of the pair as x an y and the key holds the value of the difference 
     * multiplied by keyFactor.
     * @param clusterDescriptors
     * @param clrSpace
     * @param longEdgeIndexes
     * @param keyFactor a large number that is multiplied by the floating point
     * difference of colors to result in a long key for the node.
     * @return 
     */
    private DoubleLinkedCircularList populateQueueByAscendingPairDifferences(
        float[][] clusterDescriptors, int clrSpace, 
        List<Integer> longEdgeIndexes, final long keyFactor) {

        int n = longEdgeIndexes.size();
        int nComb = (n * (n-1))/2;// n!/(2!*(n-2)! = (n * (n-1))/2
        
        PairInt[] indexes = new PairInt[nComb];
        float[] diffs = new float[nComb];
        
        CIEChromaticity cieC = new CIEChromaticity();

        int count = 0;
        for (int i = 0; i < n; ++i) {
            
            int idx1 = longEdgeIndexes.get(i).intValue();
            float[] desc1 = clusterDescriptors[idx1];
            
            for (int j = (i + 1); j < n; ++j) {
                
                int idx2 = longEdgeIndexes.get(j).intValue();
                float[] desc2 = clusterDescriptors[idx2];
                
                double diff;
                if (clrSpace == 0) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(
                        desc1[0], desc1[1], desc1[2],
                        desc2[0], desc2[1], desc2[2]));
                } else {
                    double diff1 = desc1[0] - desc2[0];
                    double diff2 = desc1[1] - desc2[1];
                    double diff3 = desc1[2] - desc2[2];
                    diff = Math.sqrt(diff1 * diff1 + diff2 * diff2 + diff3 * diff3);
                }
                
                indexes[count] = new PairInt(idx1, idx2);
                diffs[count] = (float)diff;
                count++;
            }
        }
        
        assert(nComb == count);
        
        QuickSort.sortBy1stArg(diffs, indexes);
        
        DoubleLinkedCircularList queue = new DoubleLinkedCircularList();
        
        for (int i = 0; i < indexes.length; ++i) {
            PairInt p12 = indexes[i];
            float diff = diffs[i];
            long key = (long)(diff * (double)keyFactor);
            HeapNode node = new HeapNode(key);
            node.setData(p12);
            queue.insert(node);
        }

        return queue;
    }

    private void growEdges(ImageExt img,
        List<Set<PairInt>> clusterPoints, 
        float[][] clusterDescriptors, int clrSpace, double tColor, 
        List<Integer> shortEdgeIndexes, List<Integer> longEdgeIndexes,
        int tLen) {
        
        // for each edge, use bfs traversal of neighbors to add those with diff < tColor
        //    if an adjacent pixel is part of a short edge cluster,
        //    then all of that short edge is added to the cluster
        
        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < clusterPoints.size(); ++i) {
            Set<PairInt> set = clusterPoints.get(i);            
            Integer index = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, index);
            }
        }
        Set<PairInt> shortEdgePoints = new HashSet<PairInt>();
        for (Integer index : shortEdgeIndexes) {
            int idx = index.intValue();
            Set<PairInt> set = clusterPoints.get(idx);
            shortEdgePoints.addAll(set);
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int width = img.getWidth();
        int height = img.getHeight();
        
        for (int idx1 = 0; idx1 < clusterPoints.size(); ++idx1) {
            
            Set<PairInt> set1 = clusterPoints.get(idx1);
            
            if (set1.isEmpty()) {
                continue;
            }
            
            float[] desc1 = clusterDescriptors[idx1];
            
            Integer index1 = Integer.valueOf(idx1);
            
            // BFS traversal for each edge until no more pixels are added
            
            // NOTE: the descriptors are updated w/ each addition
            //   but may want to compare results to not updating the
            //   descriptors
            
            ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();
            for (PairInt p : set1) {
                queue.add(p);
            }
                                    
            while (!queue.isEmpty()) {
                
                PairInt p = queue.poll();
                                
                int x = p.getX();
                int y = p.getY();
                
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if (x2 < 0 || y2 < 0 || (x2 > (width - 1)) || 
                        (y2 > (height - 1))) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    Integer index2 = pointIndexMap.get(p2);
                    if (index2 != null) {
                        // if this is a short edge and not the same as index1, do not skip
                        if (!(shortEdgePoints.contains(p2) && !index2.equals(index1))) {
                            continue;
                        }
                    }
                    
                    float[] clrs1 = getColors(img, x2, y2, clrSpace);
                                        
                    double diff;
                    if (clrSpace == 0) {
                        diff = Math.abs(cieC.calcDeltaECIE2000(
                            clrs1[0], clrs1[1], clrs1[2], desc1[0], desc1[1], desc1[2]));
                    } else {
                        double diff1 = clrs1[0] - desc1[0];
                        double diff2 = clrs1[1] - desc1[1];
                        double diff3 = clrs1[2] - desc1[2];
                        diff = Math.sqrt(diff1 * diff1 + diff2 * diff2 + diff3 * diff3);
                    }
                    
                    if (diff > tColor) {
                        continue;
                    }
                    
                    if (shortEdgePoints.contains(p2)) {
                        int idx2 = index2.intValue();
                        Set<PairInt> set2 = clusterPoints.get(idx2);
                        float[] desc2 = clusterDescriptors[idx2];
                        float n1 = set1.size();
                        float n2 = set2.size();
                        float nTot = n1 + n2;
                        
                        //{h, s, v, nPix, cenX, cenY}
                        // update desc1 contents for contents in desc2
                        for (int ii = 0; ii < desc1.length; ++ii) {
                            if (ii == 3) {
                                desc1[ii] = nTot;
                            } else {
                                desc1[ii] = ((desc1[ii] * n1) + (desc2[ii] * n2)) / nTot;
                            }
                        }
                        
                        for (PairInt p3 : set2) {
                            shortEdgePoints.remove(p3);
                            pointIndexMap.put(p3, index1);
                            set1.add(p3);
                            queue.offer(p3);
                        }
                        clusterDescriptors[idx2] = null;
                        set2.clear();
                            
                    } else {
                        float n1 = set1.size();
                        float nTot = n1 + 1;
                        //{h, s, v, nPix, cenX, cenY}
                        desc1[0] = ((desc1[0] * n1) + clrs1[0])/nTot;
                        desc1[1] = ((desc1[1] * n1) + clrs1[1])/nTot;
                        desc1[2] = ((desc1[2] * n1) + clrs1[2])/nTot;
                        desc1[3]++;
                        desc1[4] = ((desc1[4] * n1) + x2)/nTot;
                        desc1[5] = ((desc1[5] * n1) + y2)/nTot;
                            
                        pointIndexMap.put(p2, index1);
                        set1.add(p2);
                        queue.offer(p2);
                    }
                }
            }
        }
        
        // for pixels that are in one region and adjacent to another,
        //    evaluate it's min distance color to all it is adjacent too
        //    and relocate it if necessary.
        
        // --- initialize with the first points on boundaries ----
        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();
        for (Entry<PairInt, Integer> entry : pointIndexMap.entrySet()) {
            PairInt p = entry.getKey();
            Integer index1 = entry.getValue();
            int x = p.getX();
            int y = p.getY();
            Set<Integer> neighboringIndexes = new HashSet<Integer>();
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                Integer index2 = pointIndexMap.get(p2);
                if (index2 == null) {
                    continue;
                }
                if (!index1.equals(index2)) {
                    neighboringIndexes.add(index2);
                }
            }
            if (neighboringIndexes.size() > 0) {
                queue.add(p);
            }
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        while (!queue.isEmpty()) {
            
            PairInt p = queue.poll();
            
            if (visited.contains(p)) {
                continue;
            }

            visited.add(p);
            
            Integer index1 = pointIndexMap.get(p);
            
            int x = p.getX();
            int y = p.getY();

            Set<Integer> neighboringIndexes = new HashSet<Integer>();
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                Integer index2 = pointIndexMap.get(p2);
                if (index2 == null) {
                    continue;
                }
    
                if (!index1.equals(index2)) {
                    neighboringIndexes.add(index2);
                }
            }
            if (neighboringIndexes.isEmpty()) {
                continue;
            }
            
            neighboringIndexes.add(index1);
                        
            // find min diference between this point and the indexes and index1
            float[] clrs1 = getColors(img, x, y, clrSpace);
            int minDiffIdx = -1;
            double minDiff = Double.MAX_VALUE;
            for (Integer index2 : neighboringIndexes) {
                double diff;
                float[] desc2 = clusterDescriptors[index2.intValue()];
                if (clrSpace == 0) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(
                        clrs1[0], clrs1[1], clrs1[2], desc2[0], desc2[1], desc2[2]));
                } else {
                    double diff1 = clrs1[0] - desc2[0];
                    double diff2 = clrs1[1] - desc2[1];
                    double diff3 = clrs1[2] - desc2[2];
                    diff = Math.sqrt(diff1 * diff1 + diff2 * diff2 + diff3 * diff3);
                }
                if (diff > tColor) {
                    continue;
                }
                if (diff < minDiff) {
                    minDiff = diff;
                    minDiffIdx = index2.intValue();
                }
            }
            
            if (minDiffIdx == -1 || (minDiffIdx == index1.intValue())) {
                continue;
            }
            
            // change point p to minDiffIdx cluster
            float[] desc = clusterDescriptors[minDiffIdx];
            Set<PairInt> set = clusterPoints.get(minDiffIdx);
            float n1 = set.size();
            float nTot = n1 + 1;
            //{h, s, v, nPix, cenX, cenY}
            desc[0] = ((desc[0] * n1) + clrs1[0]) / nTot;
            desc[1] = ((desc[1] * n1) + clrs1[1]) / nTot;
            desc[2] = ((desc[2] * n1) + clrs1[2]) / nTot;
            desc[3]++;
            desc[4] = ((desc[4] * n1) + x) / nTot;
            desc[5] = ((desc[5] * n1) + y) / nTot;
            pointIndexMap.put(p, Integer.valueOf(minDiffIdx));
            set.add(p);
            
            desc = clusterDescriptors[index1.intValue()];
            set = clusterPoints.get(index1.intValue());
            n1 = set.size();
            nTot = n1 - 1;
            //{h, s, v, nPix, cenX, cenY}
            desc[0] = ((desc[0] * n1) - clrs1[0]) / nTot;
            desc[1] = ((desc[1] * n1) - clrs1[1]) / nTot;
            desc[2] = ((desc[2] * n1) - clrs1[2]) / nTot;
            desc[3]--;
            desc[4] = ((desc[4] * n1) - x) / nTot;
            desc[5] = ((desc[5] * n1) - y) / nTot;
            set.remove(p);
            
            //add to queue, neighbors of p which are not in minDiffIdx
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                Integer index2 = pointIndexMap.get(p2);
                if (index2 == null || (index2.intValue() == minDiffIdx)) {
                    continue;
                }
                queue.offer(p2);
            }
        }
        
        // --------- handle the unassigned pixels ------
        List<Set<PairInt>> unassignedClusters = new ArrayList<Set<PairInt>>();
        
        float[][] unassignedDescriptors = 
            findAndMergeUnassignedPixels(img, pointIndexMap, 
                unassignedClusters, clrSpace, tLen, tColor);
        
        {   // DEBUG
            int nExtraForDot = 1;
            Image img2 = img.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(unassignedClusters, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_unassigned_clustered_" +  clrSpace);            
        }
                
        // merge by proximity, but if there is more than one adjacent cluster
        // from pointIndexMap, choose then closest in color.
        for (int idx = 0; idx < unassignedClusters.size(); ++idx) {
            Set<PairInt> set = unassignedClusters.get(idx);
            
            float[] desc1 = unassignedDescriptors[idx];
            
            double minDiff = Double.MAX_VALUE;
            Integer minDiffIndex = null;
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    Integer index2 = pointIndexMap.get(p2);
                    if (index2 == null) {
                        continue;
                    }
                    float[] desc2 = clusterDescriptors[index2.intValue()];
                    assert(desc2 != null);
                    double  diff;
                    if (clrSpace == 0) {
                        diff = Math.abs(cieC.calcDeltaECIE2000(
                            desc1[0], desc1[1], desc1[2], desc2[0], desc2[1], desc2[2]));
                    } else {
                        double diff1 = desc1[0] - desc2[0];
                        double diff2 = desc1[1] - desc2[1];
                        double diff3 = desc1[2] - desc2[2];
                        diff = Math.sqrt(diff1 * diff1 + diff2 * diff2 + diff3 * diff3);
                    }
                    if (diff < minDiff) {
                        minDiff = diff;
                        minDiffIndex = index2;
                    }
                }
            }
            if (minDiffIndex == null) {
                clusterPoints.add(set);
 
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
System.out.println("unassigned embedded? " 
+ Arrays.toString(curveHelper.calculateXYCentroids(set)));

            } else {
                for (PairInt p3 : set) {
                    pointIndexMap.put(p3, minDiffIndex);
                }
                // add set to cluster minDiffIndex, but do not update the descriptors
                clusterPoints.get(minDiffIndex.intValue()).addAll(set);
            }
        }
        
        {   // DEBUG
            int nExtraForDot = 1;
            Image img2 = img.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(clusterPoints, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_all_merged_" +  clrSpace);   
            
            System.out.println("pointIndexMap.size() = " + pointIndexMap.size()
                + " nPix=" + img.getNPixels());
        }
    }

    private void mergeShortEdges(List<Set<PairInt>> clusterPoints, 
        float[][] clusterDescriptors, int clrSpace, double tColor, 
        List<Integer> shortEdgeIndexes) {
        
        /*
        the paper suggests:
            "The pair of lines where the distance of centroids between them is 
            nearest is always merged into one if their color difference is not 
            exceeding the predefined threshold Tc."
        
        this is only merging a single pair at most for every short edge
        */
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < shortEdgeIndexes.size(); ++i) {
            
            Integer index1 = shortEdgeIndexes.get(i);
            int idx1 = index1.intValue();
            Set<PairInt> set1 = clusterPoints.get(idx1);
            if (set1.isEmpty()) {
                continue;
            }
            
            float[] desc1 = clusterDescriptors[idx1];
            
            assert(Math.abs(set1.size() - desc1[3]) < 0.1);
            
            double minDistSq = Double.MAX_VALUE;
            int minDistIdx = -1;
            
            for (int j = (i + 1); j < shortEdgeIndexes.size(); ++j) {
                
                Integer index2 = shortEdgeIndexes.get(j);
                int idx2 = index2.intValue();
                Set<PairInt> set2 = clusterPoints.get(idx2);
                if (set2.isEmpty()) {
                    continue;
                }
                
                //{h, s, v, nPix, cenX, cenY}
                float[] desc2 = clusterDescriptors[idx2];
                
                assert(Math.abs(set2.size() - desc2[3]) < 0.1);
                
                double diff;
                if (clrSpace == 0) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(
                        desc1[0], desc1[1], desc1[2], 
                        desc2[0], desc2[1], desc2[2]));
                } else {
                    double diff1 = desc1[0] - desc2[0];
                    double diff2 = desc1[1] - desc2[1];
                    double diff3 = desc1[2] - desc2[2];
                    diff = Math.sqrt(diff1 * diff1 + diff2*diff2 + diff3*diff3);
                }
                
                if (diff >= tColor) {
                    continue;
                }
                
                float diffX = desc1[4] - desc2[4];
                float diffY = desc1[5] - desc2[5];
                double distSq = diffX * diffX + diffY * diffY;
                
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistIdx = idx2;
                }
            }
            
            if (minDistIdx == -1) {
                continue;
            }
            
            // merge set2 with set1 and update associated data structures
            Set<PairInt> set2 = clusterPoints.get(minDistIdx);
            float[] desc2 = clusterDescriptors[minDistIdx];
            
            float n1 = set1.size();
            float n2 = set2.size();
            float nTot = n1 + n2;

            //{h, s, v, nPix, cenX, cenY}
            // update desc1 contents for contents in desc2
            for (int k = 0; k < desc1.length; ++k) {
                if (k == 3) {
                    desc1[k] = nTot;
                } else {
                    desc1[k] = ((desc1[k] * n1) + (desc2[k] * n2)) / nTot;
                }
            }
            clusterDescriptors[minDistIdx] = null;
            set1.addAll(set2);
            set2.clear();
        }        
    }

    /**
     * get cie lab or hsv colors from img for coordinates (x, y)
     * @param img
     * @param x
     * @param y
     * @param clrSpace
     * @return 
     */
    private float[] getColors(ImageExt img, int x, int y, int clrSpace) {

        if (clrSpace == 0) {
            float[] lab2 = img.getCIELAB(x, y);
            return lab2;
        } else {
            float[] hsv = new float[3];
            hsv[0] = img.getHue(x, y);
            hsv[1] = img.getSaturation(x, y);
            hsv[2] = img.getBrightness(x, y);
            return hsv;
        }
    }

    private boolean assertDescriptorCounts(List<Set<PairInt>> clusterPoints, 
        float[][] clusterDescriptors) {
        
        for (int i = 0; i < clusterPoints.size(); ++i) {
            int n = clusterPoints.get(i).size();
            if (n == 0) { continue;}
            float diff = Math.abs(n - clusterDescriptors[i][3]);
            assert(diff < 0.1);
        }
        
        return true;
    }

    private void mergeByColorHistograms(ImageExt input, 
        List<Set<PairInt>> clusterPoints, int clrSpace, double tR) {
                
        int[][][] colorHistograms = calculateColorHistograms(input, 
            clusterPoints, clrSpace);
/*        
int t0 = -1;
int t1 = -1;
for (int i = 0; i < clusterPoints.size(); ++i ) {
    Set<PairInt> set = clusterPoints.get(i);
    if (set.contains(new PairInt(172, 158))) {
        t0 = i;
    }
    if (set.contains(new PairInt(156, 170))) {
        t1 = i;
    }
}
System.out.println("---> " + t0 + " ----> " + t1);       
*/      
        Map<Integer, Set<Integer>> adjacencyMap = createAdjacencyMap(
            clusterPoints);
        
        System.out.println(adjacencyMap.size() + " adjacenct clusters");
        
        // key is index1, index2 where index1 < index2
        Map<PairInt, HeapNode> nodesMap = new HashMap<PairInt, HeapNode>();
        
        Heap heap = new Heap();
        
        ColorHistogram ch = new ColorHistogram();
        
        // the histogram intersection range of values is 0 : nColors * 1
        // so for 3 colors, expect that max similarity is 3.0.
        // need to merge by higher similarity, so need to invert the keys.
        // 3 - similairty bcomes the new key.
        // a tR of 0.7*3.0 = 2.1 becomes 0.9 and any values larger than
        //    that are less similar.
        double tRInv = 3.0 - tR;
        
        long heapKeyFactor = input.getNPixels();

        for (Entry<Integer, Set<Integer>> entry : adjacencyMap.entrySet()) {
            
            Integer index1 = entry.getKey();
            int idx1 = index1.intValue();
            
            int[][] hist1 = colorHistograms[idx1];
            assert(hist1 != null);
            
            Set<Integer> indexes2 = entry.getValue();
            
            for (Integer index2 : indexes2) {
                
                int idx2 = index2.intValue();
                
                assert(idx1 != idx2);
                PairInt p12;
                if (idx1 < idx2) {
                    p12 = new PairInt(idx1, idx2);
                } else {
                    p12 = new PairInt(idx2, idx1);
                }

                if (nodesMap.containsKey(p12)) {
                    continue;
                }
                
                int[][] hist2 = colorHistograms[index2.intValue()];
                assert(hist2 != null);
                
                float similarity = 3.0f - ch.intersection(hist1, hist2);
/*
//  14: (85.8, 36.1)  9: (17, 166)   diff=2.14
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
double[] xyCen1 = curveHelper.calculateXYCentroids(clusterPoints.get(idx1));
double[] xyCen2 = curveHelper.calculateXYCentroids(clusterPoints.get(idx2));
System.out.println( String.format("%d, %d  (%.1f, %.1f) (%.1f, %.1f)", idx1, idx2,
(float)xyCen1[0], (float)xyCen1[1],  (float)xyCen2[0], (float)xyCen2[1])
+ " diff=" + similarity + " tRInv=" + tRInv);               
*/      
                long key = (long)(similarity * (double)heapKeyFactor);
                HeapNode node = new HeapNode(key);
                node.setData(p12);
                heap.insert(node);
                nodesMap.put(p12, node);

                assert(heap.getNumberOfNodes() == nodesMap.size());
            }
        }
        
        int nMerged = 0;
        
        while(!heap.isEmpty()) {
            
            HeapNode node = heap.extractMin();
            
            nodesMap.remove((PairInt)node.getData());
            
            // this is 3.0 - similarity
            double diff = ((double)node.getKey())/((double)heapKeyFactor);

            if (diff > tRInv) {
                break;
            }
            
            PairInt p12 = (PairInt)node.getData();
            
            int idx1 = p12.getX();
            int idx2 = p12.getY();
            
            Set<PairInt> set1 = clusterPoints.get(idx1);
            Set<PairInt> set2 = clusterPoints.get(idx2);

            if (set1.isEmpty() || set2.isEmpty()) {
                continue;
            }

            if (set2.size() > set1.size()) {
                idx1 = p12.getY();
                idx2 = p12.getX();
                set1 = set2;
                set2 = clusterPoints.get(idx2);
            }

            // set1 is largest
            
            int[][] hist1 = colorHistograms[idx1];
            ch.add2To1(hist1, colorHistograms[idx2]);
            colorHistograms[idx2] = null;
            
            float n1 = set1.size();
            float n2 = set2.size();
            float nTot = n1 + n2;
            
            set1.addAll(set2);
            set2.clear();
            n1 = set1.size();
            
            Integer index1 = Integer.valueOf(idx1);
            Integer index2 = Integer.valueOf(idx2);
            
            // update all pairs having set 1
            Set<Integer> indexes3 = adjacencyMap.get(index1);
            for (Integer index3 : indexes3) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                if (set3.isEmpty() || idx1 == idx3) {
                    continue;
                }
                int[][] hist3 = colorHistograms[idx3];
                assert(hist3 != null);
                
                PairInt p13;
                if (idx1 < idx3) {
                    p13 = new PairInt(idx1, idx3);
                } else {
                    p13 = new PairInt(idx3, idx1);
                }
                
                HeapNode node3 = nodesMap.get(p13);
                assert(node3 != null);
                
                heap.remove(node3);

                float similarity3 = 3.0f - ch.intersection(hist1, hist3);
                
                long key3 = (long)(similarity3 * (double)heapKeyFactor);
                node3 = new HeapNode(key3);
                node3.setData(p13);
                
                heap.insert(node3);
                nodesMap.put(p13, node3);
                
                assert(heap.getNumberOfNodes() == nodesMap.size());
             
            }

            //update adjacency map
            Set<Integer> adj2 = new HashSet<Integer>(adjacencyMap.get(index2));
            for (Integer index3 : adj2) {
                if (!index3.equals(index1)) {
                    adjacencyMap.get(index3).remove(index2);
                    adjacencyMap.get(index3).add(index1);
                }
            }
            indexes3.addAll(adj2);
            indexes3.remove(index2);
            adjacencyMap.remove(index2);
            
            nMerged++;
        }
        
        for (int i = (clusterPoints.size() - 1); i > -1; --i) {
            if (clusterPoints.get(i).isEmpty()) {
                clusterPoints.remove(i);
            }
        }
       
        System.out.println("color histogram nMerged=" + nMerged);
    }

    private int[][][] calculateColorHistograms(ImageExt input, 
        List<Set<PairInt>> clusterPoints, int clrSpace) {
        
        //0 == cie lab,  1 = hsv, 2 = rgb
        
        int n = clusterPoints.size();
        
        int[][][] hist = new int[n][][];
        
        ColorHistogram ch = new ColorHistogram();
        
        for (int i = 0; i < n; ++i) {
            
            Set<PairInt> set = clusterPoints.get(i);
            
            if (set.isEmpty()) {
                continue;
            }
            
            if (clrSpace == 0) {
                hist[i] = ch.histogramCIELAB(input, set);
            } else if (clrSpace == 1) {
                hist[i] = ch.histogramHSV(input, set);
            } if (clrSpace == 2) {
                hist[i] = ch.histogramRGB(input, set);
            }
        }
        
        return hist;
    }

    private Map<Integer, Set<Integer>> createAdjacencyMap(List<Set<PairInt>> 
        clusterPoints) {
        
        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        
        for (int i = 0; i < clusterPoints.size(); ++i) {
            Set<PairInt> set = clusterPoints.get(i);
            Integer index = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, index);
            }
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Map<Integer, Set<Integer>> adjMap = new HashMap<Integer, Set<Integer>>();
        
        for (int i = 0; i < clusterPoints.size(); ++i) {
            
            Set<PairInt> set = clusterPoints.get(i);
            
            Set<Integer> indexes2 = new HashSet<Integer>();
            
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    Integer index2 = pointIndexMap.get(p2);
                    if (index2 == null || index2.intValue() == i) {
                        continue;
                    }                    
                    indexes2.add(index2);
                }
            }
            
            if (indexes2.isEmpty()) {
                continue;
            }
            
            // add these to all point sets in adjacency map
            
            indexes2.add(Integer.valueOf(i));
            
            for (Integer key : indexes2) {
                Set<Integer> v = new HashSet<Integer>(indexes2);
                v.remove(key);
                
                Set<Integer> mapV = adjMap.get(key);
                if (mapV == null) {
                    adjMap.put(key, v);
                } else {
                    mapV.addAll(v);
                }
            }
        }
        
        return adjMap;
    }

    private Set<PairInt> findUnassigned(ImageExt img, Map<PairInt, Integer>
        pointIndexMap) {
        
        Set<PairInt> unassigned = new HashSet<PairInt>();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                PairInt p = new PairInt(i, j);
                if (!pointIndexMap.containsKey(p)) {
                    unassigned.add(p);
                }
            }
        }
        
        return unassigned;
    }

    private List<Set<PairInt>> findContiguous(Set<PairInt> points) {
        
        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(points);
        
        List<Set<PairInt>>  sets = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            Set<PairInt> set = finder.getXY(i);
            sets.add(set);
        }
        
        return sets;
    }

    private float[][] findAndMergeUnassignedPixels(ImageExt img, 
        Map<PairInt, Integer> pointIndexMap, 
        List<Set<PairInt>> outputUnassignedClusters, int clrSpace, int tLen,
        double tColor) {
        
        // assign unassigned pixels to adjacent clusters.
        Set<PairInt> unassigned = findUnassigned(img, pointIndexMap);
        
        outputUnassignedClusters.addAll(findContiguous(unassigned));
        
        {   // DEBUG
            int nExtraForDot = 1;
            Image img2 = img.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(outputUnassignedClusters, 0, 0, 
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_unassigned_clustered_0_" +  clrSpace);            
        }
        
        // treating the unassigned as did the long and short edges
        //the descriptors are
        //     C_i = {h, s, v, nPix, cenX, cenY}  or labL, labA, labB for colorSpace = 0
        
        float[][] descriptors = new float[outputUnassignedClusters.size()][];
        populateDescriptors(img, outputUnassignedClusters, descriptors, clrSpace);
        
        List<Integer> largeSetIndexes = new ArrayList<Integer>();
        List<Integer> smallSetIndexes = new ArrayList<Integer>();
        populateEdgeLengthLists(descriptors, tLen, 
            largeSetIndexes, smallSetIndexes); 

        mergeEdges(outputUnassignedClusters, descriptors, clrSpace, tColor,
            largeSetIndexes);
  
        mergeShortEdges(outputUnassignedClusters, descriptors, clrSpace, tColor,
            smallSetIndexes);
        
        //seperate the clusters into contiguous groups
        List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < outputUnassignedClusters.size(); ++i) {
            Set<PairInt> set = outputUnassignedClusters.get(i);
            if (set.isEmpty()) {
                continue;
            }
            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(set);
            for (int j = 0; j < finder.getNumberOfGroups(); ++j) {
                Set<PairInt> set2 = finder.getXY(j);
                tmp.add(set2);
            }
        }
        outputUnassignedClusters.clear();
        outputUnassignedClusters.addAll(tmp);
        
        descriptors = new float[outputUnassignedClusters.size()][];
        populateDescriptors(img, outputUnassignedClusters, descriptors, clrSpace);
        
        return descriptors;
    }

    private void populateDescriptors(ImageExt img, 
        List<Set<PairInt>> outputPoints, 
        float[][] outputDescripors, int clrSpace) {
        
        /*
        the descriptors are
             C_i = {h, s, v, nPix, cenX, cenY}  or labL, labA, labB for colorSpace = 0
        */
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < outputPoints.size(); ++i) {
            
            Set<PairInt> edgePoints = outputPoints.get(i);
            
            outputDescripors[i] = new float[6];
                        
            double[] xyCen = curveHelper.calculateXYCentroids(edgePoints);
            
            double c1Sum = 0;
            double c2Sum = 0;
            double c3Sum = 0;
            for (PairInt p : edgePoints) {
                
                int x = p.getX();
                int y = p.getY();
                
                if (clrSpace == 0) {
                    float[] lab = img.getCIELAB(x, y);
                    c1Sum += lab[0];
                    c2Sum += lab[1];
                    c3Sum += lab[2];
                } else {
                    // hsv
                    c1Sum += img.getHue(x, y);
                    c2Sum += img.getSaturation(x, y);
                    c3Sum += img.getBrightness(x, y);
                }
            }
            c1Sum /= (float)edgePoints.size();
            c2Sum /= (float)edgePoints.size();
            c3Sum /= (float)edgePoints.size();
            
            //C_i = {h, s, v, nPix, cenX, cenY}
            outputDescripors[i][0] = (float)c1Sum;
            outputDescripors[i][1] = (float)c2Sum;
            outputDescripors[i][2] = (float)c3Sum;
            outputDescripors[i][3] = edgePoints.size();
            outputDescripors[i][4] = (float)xyCen[0];
            outputDescripors[i][5] = (float)xyCen[1];
        }        
    }

    public static class BoundingRegions {
        private final List<PairIntArray> perimeterList;
        private final BlobMedialAxes bma;
        private final Map<PairInt, Integer> pointIndexMap;
        public BoundingRegions(List<PairIntArray> perimeters, BlobMedialAxes
            skeletons, Map<PairInt, Integer> pointIndexMap) {
            this.perimeterList = perimeters;
            this.bma = skeletons;
            this.pointIndexMap = pointIndexMap;
        }
        public List<PairIntArray> getPerimeterList() {
            return perimeterList;
        }
        public BlobMedialAxes getBlobMedialAxes() {
            return bma;
        }
        public Map<PairInt, Integer> getPointIndexMap() {
            return pointIndexMap;
        }

        /**
         * update the internal datasets
         * @param removeIndexes an ascending list of unique indexes to remove
         */
        public void removeIndexes(final List<Integer> removeIndexes) {

            /*
            updating: Map<PairInt, Integer> pointIndexMap
            convert to List<Set<PairInt>>
            perform removal of indexes,
            then re-populate pointIndexMap
            */
            List<Set<PairInt>> indexPoints = new ArrayList<Set<PairInt>>();
            for (int i = 0; i < perimeterList.size(); ++i) {
                indexPoints.add(new HashSet<PairInt>());
            }
            for (Entry<PairInt, Integer> entry : pointIndexMap.entrySet()) {
                PairInt p = entry.getKey();
                int idx = entry.getValue().intValue();
                indexPoints.get(idx).add(p);
            }
            for (int i = (removeIndexes.size() - 1); i > -1; --i) {
                int rmIdx = removeIndexes.get(i);
                indexPoints.remove(rmIdx);
            }
            pointIndexMap.clear();
            for (int i = 0; i < indexPoints.size(); ++i) {
                Integer key = Integer.valueOf(i);
                Set<PairInt> points = indexPoints.get(i);
                for (PairInt p : points) {
                    pointIndexMap.put(p, key);
                }
            }

            for (int i = (removeIndexes.size() - 1); i > -1; --i) {
                int rmIdx = removeIndexes.get(i);
                perimeterList.remove(rmIdx);
            }

            bma.removeIndexes(removeIndexes);
        }
    }

    /**
     * runtime complexity is approx O(N_perimeter_pts * lg_2(N_perimeter_pts))
     *
     * @param orderedPerimeter
     * @param width
     * @param height
     * @param srchRadius
     */
    public void makeStraightLinesHollow(PairIntArray orderedPerimeter,
        int width, int height, float srchRadius) {

        HoughTransform ht = new HoughTransform();

        //O(N_edge_pts), but includes transcendental operations
        Map<PairInt, Set<PairInt>> trPointsMap = ht.calculateLineGivenEdge(
            orderedPerimeter, width, height);

        //O(N_edge_pts * lg_2(N_edge_pts))
        List<PairInt> outSortedKeys = ht.sortByVotes(trPointsMap);

        int thetaTol = 1;
        int radiusTol = (int) Math.ceil(srchRadius);

        //runtime complexity is approx O(N_pix * lg_2(N_pix)).
        // === find indiv lines within the edge ====
        HoughTransform.HoughTransformLines htl
            = ht.createPixTRMapsFromSorted(outSortedKeys, trPointsMap, thetaTol,
                radiusTol);

        Map<PairInt, PairInt> pixToTRMap = htl.getPixelToPolarCoordMap();

        Set<PairInt> remove = new HashSet<PairInt>();
        for (Set<PairInt> line : htl.getSortedLineGroups()) {
            if (line.size() < 3) {
                continue;
            }
            PairInt tr = pixToTRMap.get(line.iterator().next());

            if ((Math.abs(tr.getX() - 180) < 20) || (Math.abs(tr.getX() - 0) < 20)
                || (Math.abs(tr.getX() - 360) < 20)) {
                // these are nearly vertical lines, so endpoints are min and max Y
                int minY = Integer.MAX_VALUE;
                int maxY = Integer.MIN_VALUE;
                for (PairInt p : line) {
                    int y = p.getY();
                    if (y < minY) {
                        minY = y;
                    }
                    if (y > maxY) {
                        maxY = y;
                    }
                }
                for (PairInt p : line) {
                    int y = p.getY();
                    if (y > minY && y < maxY) {
                        remove.add(p);
                    }
                }

            } else {
                // endpoints are min and max X
                int minX = Integer.MAX_VALUE;
                int maxX = Integer.MIN_VALUE;
                for (PairInt p : line) {
                    int x = p.getX();
                    if (x < minX) {
                        minX = x;
                    }
                    if (x > maxX) {
                        maxX = x;
                    }
                }
                for (PairInt p : line) {
                    int x = p.getX();
                    if (x > minX && x < maxX) {
                        remove.add(p);
                    }
                }
            }
        }
        PairIntArray output = new PairIntArray();
        for (int i = 0; i < orderedPerimeter.getN(); ++i) {
            int x = orderedPerimeter.getX(i);
            int y = orderedPerimeter.getY(i);
            PairInt p = new PairInt(x, y);
            if (!remove.contains(p)) {
                output.add(x, y);
            }
        }
        if (output.getN() > 0) {
            orderedPerimeter.swapContents(output);
        }
    }

    /**
     * a quick look at the feature orientations and a knowledge of "inward"
     * w.r.t. the encapsulating segmentation region.
     * @param features
     * @param blobs
     * @param img
     * @param lbl
     */
    private void plotOrientation(IntensityClrFeatures features,
        BoundingRegions boundingRegion, ImageExt img, String lbl) {

        /*
        -- pass into method the features so can determine orientation,
        -- then for each point on the bounds, plot the orientation first
           point in red/white and subsequent short line in red, then
           add a black/red/white dot for the segment endpoint pointing inward.
        */

        List<PairIntArray> perimetersList = boundingRegion.getPerimeterList();
        BlobMedialAxes bma = boundingRegion.getBlobMedialAxes();

        ImageExt imgCp = img.copyToImageExt();

        // make skeleton for each blob for detailed perimeter "inward" directions
        for (int i = 0; i < perimetersList.size(); ++i) {

            PairIntArray perimeter = perimetersList.get(i);

            for (int ii = 0; ii < perimeter.getN(); ++ii) {

                int x = perimeter.getX(ii);
                int y = perimeter.getY(ii);

                int rotD;
                try {
                    rotD = features.calculateOrientation(x, y);
                } catch (CornerRegion.CornerRegionDegneracyException e) {
                    continue;
                }

                // find the closest skeleton point then find the
                // neighbor closest to it in counter clockwise direction.
                PairInt xySkel = bma.findClosestPoint(i, x, y);

                PairInt xyCen = bma.getOriginalBlobXYCentroid(i);

                /*
                          90
                   135    |    45
                          |
                180 ---------------  0
                          |
                   225    |    315
                         270
                */
                // direction away from skeleton or centroid
                int thetaOut;
                if ((x != xySkel.getX()) || (y != xySkel.getY())) {
                    double theta = Math.atan2(y - xySkel.getY(), x - xySkel.getX());
                    // transform to 0 to 2*pi radians
                    if (theta < 0) {
                        theta += 2. * Math.PI;
                    }
                    thetaOut = (int)Math.round(theta * 180./Math.PI);
                } else {
                    double theta = Math.atan2(y - xyCen.getY(), x - xyCen.getX());
                    // transform to 0 to 2*pi radians
                    if (theta < 0) {
                        theta += 2. * Math.PI;
                    }
                    thetaOut = (int)Math.round(theta * 180./Math.PI);
                }

                //log.info(String.format("rotation=%d  x,y=(%d,%d) skel=(%d,%d) centroid=(%d,%d) theta outward=%d",
                //    rotD, x, y, xySkel.getX(), xySkel.getY(), (int)Math.round(xyCen[0]),
                //    (int)Math.round(xyCen[1]), thetaOut));

                /*
                if orientation and thetaOut are closer than the 180 opposite
                configuration:
                    draw a large white box with black center, then segment of 3 pixels
                else
                    draw 3 pixel segmen
                */
                boolean same = Math.abs(AngleUtil.getAngleDifference(rotD, thetaOut)) < 90;
                if ((Math.abs(rotD - 90) < 20) || (Math.abs(rotD - 270) < 20)) {
                    // use dy steps
                    for (int dy = -1; dy <= 1; ++dy) {
                        int y2 = y + dy;
                        if (y2 < 0 || y2 > (imgCp.getHeight() - 1)) {
                            continue;
                        }
                        for (int dx = -1; dx <= 1; ++dx) {
                            int x2 = x + dx;
                            if (x2 < 0 || x2 > (imgCp.getWidth() - 1)) {
                                continue;
                            }
                            imgCp.setRGB(x2, y2, 255, 255, 255);
                        }
                    }
                    if (same){
                        imgCp.setRGB(x, y, 0, 0, 0);
                    } else {
                        imgCp.setRGB(x, y, 255, 0, 0);
                    }
                    double slope = Math.tan(rotD);
                    // write 2 more pixels in red
                    //y1 - y0 = slope*(x1 - x0);  ---> (1./slope)*(y1-y0) + x0 = x1
                    int dy;
                    if (rotD <= 180) {
                        dy = 1;
                    } else {
                        dy = -1;
                    }
                    int y1 = y + dy;
                    double x1 = (1./slope)*((double)dy) + (double) x;
                    if (x1 < 0 || x1 > (imgCp.getWidth() - 1) || y1 < 0 ||
                        y1 > (imgCp.getHeight() - 1)) {
                        continue;
                    }
                    imgCp.setRGB((int)Math.round(x1), y1, 255, 0, 0);
                    y1 = y + 2*dy;
                    x1 = (1./slope)*(2.*dy) + (double) x;
                    if (x1 < 0 || x1 > (imgCp.getWidth() - 1) || y1 < 0 ||
                        y1 > (imgCp.getHeight() - 1)) {
                        continue;
                    }
                    imgCp.setRGB((int)Math.round(x1), y1, 255, 0, 0);
                } else {
                    // use dx steps
                    for (int dy = -1; dy <= 1; ++dy) {
                        int y2 = y + dy;
                        if (y2 < 0 || y2 > (imgCp.getHeight() - 1)) {
                            continue;
                        }
                        for (int dx = -1; dx <= 1; ++dx) {
                            int x2 = x + dx;
                            if (x2 < 0 || x2 > (imgCp.getWidth() - 1)) {
                                continue;
                            }
                            imgCp.setRGB(x2, y2, 255, 255, 255);
                        }
                    }
                    if (same) {
                        imgCp.setRGB(x, y, 0, 0, 0);
                    } else {
                        imgCp.setRGB(x, y, 255, 0, 0);
                    }
                    double slope = Math.tan(rotD);
                    // write 2 more pixels in red
                    //y1 - y0 = slope*(x1 - x0);
                    int dx;
                    if (rotD >= 90 && rotD <= 270) {
                        dx = -1;
                    } else {
                        dx = 1;
                    }
                    int x1 = x + dx;
                    double y1 = (slope * (double)dx) + (double)y;
                    if (x1 < 0 || x1 > (imgCp.getWidth() - 1) || y1 < 0 ||
                        y1 > (imgCp.getHeight() - 1)) {
                        continue;
                    }
                    imgCp.setRGB(x1, (int)Math.round(y1), 255, 0, 0);
                    x1 = x + 2 * dx;
                    y1 = (slope * 2. * dx) + (double)y;
                    if (x1 < 0 || x1 > (imgCp.getWidth() - 1) || y1 < 0 ||
                        y1 > (imgCp.getHeight() - 1)) {
                        continue;
                    }
                    imgCp.setRGB(x1, (int)Math.round(y1), 255, 0, 0);
                }
                /*
                          90
                   135    |    45
                          |
                180 ---------------  0
                          |
                   225    |    315
                         270
                */
            }
            MiscDebug.writeImage(imgCp, "_orientation_" + lbl);

            // when !same, would compare half descriptor from the
            // "bottom" half of the oriented descriptor.
            int z = 1;
        }

    }

    /**
     * for the full 8 neighbor region, determine whether nulling the pixel
     * at (col, row) would disconnect the remaining line.  Note that the
     * boolean logic is embedded in the comments.  One should be able to
     * combine the rules for multiple pixel tests to reduce the redundant
     * comparisons for the regions in common.
     *
     * Note, that row and col are expected to be at least 1 pixel distant
     * from the image borders.
     *
     * @param input
     * @param col
     * @param row
     * @return
     */
    public static boolean doesDisconnect(final GreyscaleImage input,
        PairInt[][] neighborCoords, int col, int row, int edgeValue) {

        int w = input.getWidth();
        int h = input.getHeight();

        if (((col - 1) < 0) || ((row - 1) < 0) || ((col + 1) > (w - 1)) ||
            ((row + 1) > (h - 1))) {
            // general rule so that invoker doesn't disconnect a line that is
            // connected to image boundaries
            return true;
        }

        /*
        coordinates of the 8 neighbors as already created PairInts without
        bound checks.
        indexes are found as +1 of the difference relative to center,
        for example, a point for (col-1, row-1) is found as neighborCoords[0][0]
        */

         /*
            6  7  8      +1  2      transformed by 90 rot:     15  11  6
           11 *C* 12     0   1                                 16  C*  7
           15  16 17     -1  0                                 17  12  8

           -1  0   1
            0  1   2

        disconnects:
           -- if (6) && (8) && !(7) && (!(11) || !(16) || !(12))
           -- if (6) && (12) && !(7) && (!(11) || !(16))
           -- if (6) && (15) && !(11) && (!(16) || !(12) || !(7))
           -- if (6) && (16) && !(7) && !(11)
           -- if (6) && (17) && ( (!(7) || !(12)) && (!(11) || !(16)) )
           -- if (7) && (15) && !(11) && (!(12) || !(16))
           -- if (7) && (17) && !(12) && (!(11) || !(16))
           -- if (7) && (16) && !(11) && !(12)
           -- if (8) && (11) && !(7) && (!(12) || !(16))
           -- if (8) && (17) && !(12) && (!(16) || !(11) || !(7))
           -- if (8) && (16) && !(7) && !(12)
           -- if (8) && (15) && ( (!(7) || !(11)) && (!(12) || !(16)) )
           -- if (11) && (12) && !(7) && !(16)
           -- if (11) && (17) && !(16) && (!(7) || !(12))
           -- if (12) && (15) && !(16) && (!(7) || !(11))
           -- if (15) && (17) && !(16) && (!(11) || !(7) || !(12))

        does not disconnect
           -- if (6 || 7 || 8) && !(15) && !(16) && !(17) && !(11) && !(12)

        then rotate 90 and test, then rotate 90 and test, then rotate 90 and test
        */

        boolean t6 = (input.getValue(neighborCoords[0][2].getX() + col,
            neighborCoords[0][2].getY() + row) == edgeValue);
        boolean t7 = (input.getValue(neighborCoords[1][2].getX() + col,
            neighborCoords[1][2].getY() + row) == edgeValue);
        boolean t8 = (input.getValue(neighborCoords[2][2].getX() + col,
            neighborCoords[2][2].getY() + row) == edgeValue);
        boolean t11 = (input.getValue(neighborCoords[0][1].getX() + col,
            neighborCoords[0][1].getY() + row) == edgeValue);
        boolean t12 = (input.getValue(neighborCoords[2][1].getX() + col,
            neighborCoords[2][1].getY() + row) == edgeValue);
        boolean t15 = (input.getValue(neighborCoords[0][0].getX() + col,
            neighborCoords[0][0].getY() + row) == edgeValue);
        boolean t16 = (input.getValue(neighborCoords[1][0].getX() + col,
            neighborCoords[1][0].getY() + row) == edgeValue);
        boolean t17 = (input.getValue(neighborCoords[2][0].getX() + col,
            neighborCoords[2][0].getY() + row) == edgeValue);

       if ((t6) && (t8) && !(t7) && (!(t11) || !(t16) || !(t12))) {
            return true;
        } else if ((t6) && (t12) && !(t7) && (!(t11) || !(t16))) {
            return true;
        } else if ((t6) && (t15) && !(t11) && (!(t16) || !(t12) || !(t7))) {
            return true;
        } else if ((t6) && (t16) && !(t7) && !(t11)) {
            return true;
        } else if ((t6) && (t17) && ( (!(t7) || !(t12)) && (!(t11) || !(t16)) )) {
            return true;
        } else if ((t7) && (t15) && !(t11) && (!(t12) || !(t16))) {
            return true;
        } else if ((t7) && (t17) && !(t12) && (!(t11) || !(t16))) {
            return true;
        } else if ((t7) && (t16) && !(t11) && !(t12)) {
            return true;
        } else if ((t8) && (t11) && !(t7) && (!(t12) || !(t16))) {
            return true;
        } else if ((t8) && (t17) && !(t12) && (!(t16) || !(t11) || !(t7))) {
            return true;
        } else if ((t8) && (t16) && !(t7) && !(t12)) {
            return true;
        } else if ((t8) && (t15) && ( (!(t7) || !(t11)) && (!(t12) || !(t16)) )) {
            return true;
        } else if ((t11) && (t12) && !(t7) && !(t16)) {
            return true;
        } else if ((t11) && (t17) && !(t16) && (!(t7) || !(t12))) {
            return true;
        } else if ((t12) && (t15) && !(t16) && (!(t7) || !(t11))) {
            return true;
        } else if ((t15) && (t17) && !(t16) && (!(t11) || !(t7) || !(t12))) {
            return true;
        }

        return false;
    }
    
    /**
     * for input with zeros for non-neighbor pixels else any value,
     * look within the neighborhood of point (col, row) to see if there are 
     * edges points to either side of the point that would be disconnected
     * if this one were removed.   A non-edge point is defined as having value 0.
     * 
     * @param input
     * @param neighborCoords
     * @param col
     * @param row
     * @return 
     */
    public static boolean doesDisconnect(final GreyscaleImage input,
        PairInt[][] neighborCoords, int col, int row) {

        int w = input.getWidth();
        int h = input.getHeight();

        if (((col - 1) < 0) || ((row - 1) < 0) || ((col + 1) > (w - 1)) ||
            ((row + 1) > (h - 1))) {
            // general rule so that invoker doesn't disconnect a line that is
            // connected to image boundaries
            return true;
        }

        /*
        coordinates of the 8 neighbors as already created PairInts without
        bound checks.
        indexes are found as +1 of the difference relative to center,
        for example, a point for (col-1, row-1) is found as neighborCoords[0][0]
        */

         /*
            6  7  8      +1  2      transformed by 90 rot:     15  11  6
           11 *C* 12     0   1                                 16  C*  7
           15  16 17     -1  0                                 17  12  8

           -1  0   1
            0  1   2

        disconnects:
           -- if (6) && (8) && !(7) && (!(11) || !(16) || !(12))
           -- if (6) && (12) && !(7) && (!(11) || !(16))
           -- if (6) && (15) && !(11) && (!(16) || !(12) || !(7))
           -- if (6) && (16) && !(7) && !(11)
           -- if (6) && (17) && ( (!(7) || !(12)) && (!(11) || !(16)) )
           -- if (7) && (15) && !(11) && (!(12) || !(16))
           -- if (7) && (17) && !(12) && (!(11) || !(16))
           -- if (7) && (16) && !(11) && !(12)
           -- if (8) && (11) && !(7) && (!(12) || !(16))
           -- if (8) && (17) && !(12) && (!(16) || !(11) || !(7))
           -- if (8) && (16) && !(7) && !(12)
           -- if (8) && (15) && ( (!(7) || !(11)) && (!(12) || !(16)) )
           -- if (11) && (12) && !(7) && !(16)
           -- if (11) && (17) && !(16) && (!(7) || !(12))
           -- if (12) && (15) && !(16) && (!(7) || !(11))
           -- if (15) && (17) && !(16) && (!(11) || !(7) || !(12))

        does not disconnect
           -- if (6 || 7 || 8) && !(15) && !(16) && !(17) && !(11) && !(12)

        then rotate 90 and test, then rotate 90 and test, then rotate 90 and test
        */

        boolean t6 = (input.getValue(neighborCoords[0][2].getX() + col,
            neighborCoords[0][2].getY() + row) > 0);
        boolean t7 = (input.getValue(neighborCoords[1][2].getX() + col,
            neighborCoords[1][2].getY() + row) > 0);
        boolean t8 = (input.getValue(neighborCoords[2][2].getX() + col,
            neighborCoords[2][2].getY() + row) > 0);
        boolean t11 = (input.getValue(neighborCoords[0][1].getX() + col,
            neighborCoords[0][1].getY() + row) > 0);
        boolean t12 = (input.getValue(neighborCoords[2][1].getX() + col,
            neighborCoords[2][1].getY() + row) > 0);
        boolean t15 = (input.getValue(neighborCoords[0][0].getX() + col,
            neighborCoords[0][0].getY() + row) > 0);
        boolean t16 = (input.getValue(neighborCoords[1][0].getX() + col,
            neighborCoords[1][0].getY() + row) > 0);
        boolean t17 = (input.getValue(neighborCoords[2][0].getX() + col,
            neighborCoords[2][0].getY() + row) > 0);

       if ((t6) && (t8) && !(t7) && (!(t11) || !(t16) || !(t12))) {
            return true;
        } else if ((t6) && (t12) && !(t7) && (!(t11) || !(t16))) {
            return true;
        } else if ((t6) && (t15) && !(t11) && (!(t16) || !(t12) || !(t7))) {
            return true;
        } else if ((t6) && (t16) && !(t7) && !(t11)) {
            return true;
        } else if ((t6) && (t17) && ( (!(t7) || !(t12)) && (!(t11) || !(t16)) )) {
            return true;
        } else if ((t7) && (t15) && !(t11) && (!(t12) || !(t16))) {
            return true;
        } else if ((t7) && (t17) && !(t12) && (!(t11) || !(t16))) {
            return true;
        } else if ((t7) && (t16) && !(t11) && !(t12)) {
            return true;
        } else if ((t8) && (t11) && !(t7) && (!(t12) || !(t16))) {
            return true;
        } else if ((t8) && (t17) && !(t12) && (!(t16) || !(t11) || !(t7))) {
            return true;
        } else if ((t8) && (t16) && !(t7) && !(t12)) {
            return true;
        } else if ((t8) && (t15) && ( (!(t7) || !(t11)) && (!(t12) || !(t16)) )) {
            return true;
        } else if ((t11) && (t12) && !(t7) && !(t16)) {
            return true;
        } else if ((t11) && (t17) && !(t16) && (!(t7) || !(t12))) {
            return true;
        } else if ((t12) && (t15) && !(t16) && (!(t7) || !(t11))) {
            return true;
        } else if ((t15) && (t17) && !(t16) && (!(t11) || !(t7) || !(t12))) {
            return true;
        }

        return false;
    }
    
    /**
     * for input with zeros for non-neighbor pixels else any value,
     * look within the neighborhood of point (col, row) to see if there are 
     * edges points to either side of the point that would be disconnected
     * if this one were removed.   A non-edge point is defined as having value 0.
     * 
     * @param input
     * @param neighborCoords
     * @param col
     * @param row
     * @param w width of image
     * @param h height of image
     * @return 
     */
    public static boolean doesDisconnect(final Set<PairInt> input,
        PairInt[][] neighborCoords, int col, int row, int w, int h) {

        if (((col - 1) < 0) || ((row - 1) < 0) || ((col + 1) > (w - 1)) ||
            ((row + 1) > (h - 1))) {
            // general rule so that invoker doesn't disconnect a line that is
            // connected to image boundaries
            return true;
        }

        /*
        coordinates of the 8 neighbors as already created PairInts without
        bound checks.
        indexes are found as +1 of the difference relative to center,
        for example, a point for (col-1, row-1) is found as neighborCoords[0][0]
        */

         /*
            6  7  8      +1  2      transformed by 90 rot:     15  11  6
           11 *C* 12     0   1                                 16  C*  7
           15  16 17     -1  0                                 17  12  8

           -1  0   1
            0  1   2

        disconnects:
           -- if (6) && (8) && !(7) && (!(11) || !(16) || !(12))
           -- if (6) && (12) && !(7) && (!(11) || !(16))
           -- if (6) && (15) && !(11) && (!(16) || !(12) || !(7))
           -- if (6) && (16) && !(7) && !(11)
           -- if (6) && (17) && ( (!(7) || !(12)) && (!(11) || !(16)) )
           -- if (7) && (15) && !(11) && (!(12) || !(16))
           -- if (7) && (17) && !(12) && (!(11) || !(16))
           -- if (7) && (16) && !(11) && !(12)
           -- if (8) && (11) && !(7) && (!(12) || !(16))
           -- if (8) && (17) && !(12) && (!(16) || !(11) || !(7))
           -- if (8) && (16) && !(7) && !(12)
           -- if (8) && (15) && ( (!(7) || !(11)) && (!(12) || !(16)) )
           -- if (11) && (12) && !(7) && !(16)
           -- if (11) && (17) && !(16) && (!(7) || !(12))
           -- if (12) && (15) && !(16) && (!(7) || !(11))
           -- if (15) && (17) && !(16) && (!(11) || !(7) || !(12))

        does not disconnect
           -- if (6 || 7 || 8) && !(15) && !(16) && !(17) && !(11) && !(12)

        then rotate 90 and test, then rotate 90 and test, then rotate 90 and test
        */

        boolean t6 = input.contains(
            new PairInt(neighborCoords[0][2].getX() + col,
            neighborCoords[0][2].getY() + row));
        boolean t7 = input.contains(
            new PairInt(neighborCoords[1][2].getX() + col,
            neighborCoords[1][2].getY() + row));
        boolean t8 = input.contains(
            new PairInt(neighborCoords[2][2].getX() + col,
            neighborCoords[2][2].getY() + row));
        boolean t11 = input.contains(
            new PairInt(neighborCoords[0][1].getX() + col,
            neighborCoords[0][1].getY() + row));
        boolean t12 = input.contains(
            new PairInt(neighborCoords[2][1].getX() + col,
            neighborCoords[2][1].getY() + row));
        boolean t15 = input.contains(
            new PairInt(neighborCoords[0][0].getX() + col,
            neighborCoords[0][0].getY() + row));
        boolean t16 = input.contains(
            new PairInt(neighborCoords[1][0].getX() + col,
            neighborCoords[1][0].getY() + row));
        boolean t17 = input.contains(
            new PairInt(neighborCoords[2][0].getX() + col,
            neighborCoords[2][0].getY() + row));

       if ((t6) && (t8) && !(t7) && (!(t11) || !(t16) || !(t12))) {
            return true;
        } else if ((t6) && (t12) && !(t7) && (!(t11) || !(t16))) {
            return true;
        } else if ((t6) && (t15) && !(t11) && (!(t16) || !(t12) || !(t7))) {
            return true;
        } else if ((t6) && (t16) && !(t7) && !(t11)) {
            return true;
        } else if ((t6) && (t17) && ( (!(t7) || !(t12)) && (!(t11) || !(t16)) )) {
            return true;
        } else if ((t7) && (t15) && !(t11) && (!(t12) || !(t16))) {
            return true;
        } else if ((t7) && (t17) && !(t12) && (!(t11) || !(t16))) {
            return true;
        } else if ((t7) && (t16) && !(t11) && !(t12)) {
            return true;
        } else if ((t8) && (t11) && !(t7) && (!(t12) || !(t16))) {
            return true;
        } else if ((t8) && (t17) && !(t12) && (!(t16) || !(t11) || !(t7))) {
            return true;
        } else if ((t8) && (t16) && !(t7) && !(t12)) {
            return true;
        } else if ((t8) && (t15) && ( (!(t7) || !(t11)) && (!(t12) || !(t16)) )) {
            return true;
        } else if ((t11) && (t12) && !(t7) && !(t16)) {
            return true;
        } else if ((t11) && (t17) && !(t16) && (!(t7) || !(t12))) {
            return true;
        } else if ((t12) && (t15) && !(t16) && (!(t7) || !(t11))) {
            return true;
        } else if ((t15) && (t17) && !(t16) && (!(t11) || !(t7) || !(t12))) {
            return true;
        }

        return false;
    }
    
    /**
     * for input with zeros for non-neighbor pixels else any value,
     * look within the neighborhood of point (col, row) to see if there are 
     * edges points to either side of the point that would be disconnected
     * if this one were removed.   A non-edge point is defined as having value 0.
     * 
     * @param input
     * @param neighborCoords
     * @param col
     * @param row
     * @return 
     */
    public static boolean doesDisconnect(final double[][] input,
        PairInt[][] neighborCoords, int col, int row) {

        int w = input.length;
        int h = input[0].length;

        if (((col - 1) < 0) || ((row - 1) < 0) || ((col + 1) > (w - 1)) ||
            ((row + 1) > (h - 1))) {
            // general rule so that invoker doesn't disconnect a line that is
            // connected to image boundaries
            return true;
        }

        /*
        coordinates of the 8 neighbors as already created PairInts without
        bound checks.
        indexes are found as +1 of the difference relative to center,
        for example, a point for (col-1, row-1) is found as neighborCoords[0][0]
        */

         /*
            6  7  8      +1  2      transformed by 90 rot:     15  11  6
           11 *C* 12     0   1                                 16  C*  7
           15  16 17     -1  0                                 17  12  8

           -1  0   1
            0  1   2

        disconnects:
           -- if (6) && (8) && !(7) && (!(11) || !(16) || !(12))
           -- if (6) && (12) && !(7) && (!(11) || !(16))
           -- if (6) && (15) && !(11) && (!(16) || !(12) || !(7))
           -- if (6) && (16) && !(7) && !(11)
           -- if (6) && (17) && ( (!(7) || !(12)) && (!(11) || !(16)) )
           -- if (7) && (15) && !(11) && (!(12) || !(16))
           -- if (7) && (17) && !(12) && (!(11) || !(16))
           -- if (7) && (16) && !(11) && !(12)
           -- if (8) && (11) && !(7) && (!(12) || !(16))
           -- if (8) && (17) && !(12) && (!(16) || !(11) || !(7))
           -- if (8) && (16) && !(7) && !(12)
           -- if (8) && (15) && ( (!(7) || !(11)) && (!(12) || !(16)) )
           -- if (11) && (12) && !(7) && !(16)
           -- if (11) && (17) && !(16) && (!(7) || !(12))
           -- if (12) && (15) && !(16) && (!(7) || !(11))
           -- if (15) && (17) && !(16) && (!(11) || !(7) || !(12))

        does not disconnect
           -- if (6 || 7 || 8) && !(15) && !(16) && !(17) && !(11) && !(12)

        then rotate 90 and test, then rotate 90 and test, then rotate 90 and test
        */

        boolean t6 = (input[neighborCoords[0][2].getX() + col][
            neighborCoords[0][2].getY() + row] > 0);
        boolean t7 = (input[neighborCoords[1][2].getX() + col][
            neighborCoords[1][2].getY() + row] > 0);
        boolean t8 = (input[neighborCoords[2][2].getX() + col][
            neighborCoords[2][2].getY() + row] > 0);
        boolean t11 = (input[neighborCoords[0][1].getX() + col][
            neighborCoords[0][1].getY() + row] > 0);
        boolean t12 = (input[neighborCoords[2][1].getX() + col][
            neighborCoords[2][1].getY() + row] > 0);
        boolean t15 = (input[neighborCoords[0][0].getX() + col][
            neighborCoords[0][0].getY() + row] > 0);
        boolean t16 = (input[neighborCoords[1][0].getX() + col][
            neighborCoords[1][0].getY() + row] > 0);
        boolean t17 = (input[neighborCoords[2][0].getX() + col][
            neighborCoords[2][0].getY() + row] > 0);

       if ((t6) && (t8) && !(t7) && (!(t11) || !(t16) || !(t12))) {
            return true;
        } else if ((t6) && (t12) && !(t7) && (!(t11) || !(t16))) {
            return true;
        } else if ((t6) && (t15) && !(t11) && (!(t16) || !(t12) || !(t7))) {
            return true;
        } else if ((t6) && (t16) && !(t7) && !(t11)) {
            return true;
        } else if ((t6) && (t17) && ( (!(t7) || !(t12)) && (!(t11) || !(t16)) )) {
            return true;
        } else if ((t7) && (t15) && !(t11) && (!(t12) || !(t16))) {
            return true;
        } else if ((t7) && (t17) && !(t12) && (!(t11) || !(t16))) {
            return true;
        } else if ((t7) && (t16) && !(t11) && !(t12)) {
            return true;
        } else if ((t8) && (t11) && !(t7) && (!(t12) || !(t16))) {
            return true;
        } else if ((t8) && (t17) && !(t12) && (!(t16) || !(t11) || !(t7))) {
            return true;
        } else if ((t8) && (t16) && !(t7) && !(t12)) {
            return true;
        } else if ((t8) && (t15) && ( (!(t7) || !(t11)) && (!(t12) || !(t16)) )) {
            return true;
        } else if ((t11) && (t12) && !(t7) && !(t16)) {
            return true;
        } else if ((t11) && (t17) && !(t16) && (!(t7) || !(t12))) {
            return true;
        } else if ((t12) && (t15) && !(t16) && (!(t7) || !(t11))) {
            return true;
        } else if ((t15) && (t17) && !(t16) && (!(t11) || !(t7) || !(t12))) {
            return true;
        }

        return false;
    }
    
    /**
     * given a greyscale image, makes edges (0's are edges and the background
     * is 255).
     * @param img
     * @param debugTag 
     */
    public void createEdges02(GreyscaleImage img, String debugTag) {
                
        GreyscaleImage greyGradient2 = img.copyImage();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(greyGradient2, SIGMA.ONE);
        
//TODO: an adaptive gradient might help here
        
        CannyEdgeFilterAdaptive fl = new CannyEdgeFilterAdaptive();
        fl.applyFilter(greyGradient2);    
        removeIsolatedPixels(greyGradient2, 0, 255, true);
        removeIsolatedPixels(greyGradient2, 255, 0, true);
        MedianSmooth s = new MedianSmooth();
        GreyscaleImage tmp2 = s.calculate(greyGradient2, 2, 2);
        greyGradient2 = tmp2;
        for (int i = 0; i < greyGradient2.getNPixels(); ++i) {
            int v = greyGradient2.getValue(i);
            if (v > 1) {
                img.setValue(i, 0);
            } else {
                img.setValue(i, 255);
            }
        }
        removeEdgesSmallerThanLimit(img, 0, 255, 2);
        //removeIsolatedPixels(img, 0, 255, true);
    }

    public void createEdges03(GreyscaleImage img, String debugTag) {
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();
        createEdges02(img, debugTag);
    }
    
    public void createEdges01(GreyscaleImage img, String debugTag) {

        ImageProcessor imageProcessor = new ImageProcessor();
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();
        CannyEdgeFilterAdaptive cannyFilter = new CannyEdgeFilterAdaptive();
        cannyFilter.applyFilter(img);
        
        //MiscDebug.writeImage(img, "_canny_" + debugTag);
        
        setAllNonZeroTo255(img);
        removeIsolatedPixels(img, 0, 255, false);
        removeIsolatedPixels(img, 255, 0, true);

        MedianSmooth s = new MedianSmooth();
        GreyscaleImage tmp2 = s.calculate(img, 3, 3);

        //MiscDebug.writeImage(tmp2, "tmp_edges01_2_" + debugTag);

        removeIsolatedPixels(tmp2, 255, 0, true);
        removeIsolatedPixels(tmp2, 0, 255, true);
        invertImage(tmp2);
        imageProcessor.applyAdaptiveMeanThresholding(tmp2, 1);
        img.resetTo(tmp2);

        //MiscDebug.writeImage(img, "tmp_edges01_3_" + debugTag);
    }

    public void removeSmallBubblesFromEdges(GreyscaleImage img, int edgeValue,
        int nonEdgeValue, String label) {

        removeSmallBubblesFromEdges(img, edgeValue, nonEdgeValue, -1, -1, label);
    }

     /**
     *
     * @param img
     * @param edgeValue
     * @param nonEdgeValue
     * @param limit edges less than or equal to this size will be removed
     */
    public void removeEdgesSmallerThanLimit(GreyscaleImage img, int edgeValue,
        int nonEdgeValue, int limit) {

        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);
        finder.setMinimumNumberInCluster(1);
        finder.setToUse8Neighbors();
        finder.findGroups(edgeValue);

        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            PairIntArray edge = finder.getXY(i);
            if (edge.getN() > limit) {
                continue;
            }
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                img.setValue(x, y, nonEdgeValue);
            }
        }
    }

    /**
     *
     * @param img
     * @param edgeValue
     * @param nonEdgeValue
     * @param overrideEdgesLimit0 if larger than -1, this is used as the size
     * limit of contiguous bounded white space or edges to remove, else a
     * conservative estimate is made using histogram of sizes.
     * @param label
     */
    public void removeSmallBubblesFromEdges(GreyscaleImage img, int edgeValue,
        int nonEdgeValue, int overrideEdgesLimit0, int overrideEdgesLimit1,
        String label) {

        DFSContiguousValueFinder cf0 = new DFSContiguousValueFinder(img);
        cf0.setMinimumNumberInCluster(1);
        cf0.setToUse8Neighbors();
        cf0.findGroups(edgeValue);
        Map<PairInt, Integer> nonNeighborEdges = new HashMap<PairInt, Integer>();
        for (int i = 0; i < cf0.getNumberOfGroups(); ++i) {
            PairIntArray group = cf0.getXY(i);
            Integer index = Integer.valueOf(i);
            for (int j = 0; j < group.getN(); ++j) {
                PairInt p2 = new PairInt(group.getX(j), group.getY(j));
                nonNeighborEdges.put(p2, index);
            }
        }

        DFSContiguousValueFinder cf = new DFSContiguousValueFinder(img);
        cf.setMinimumNumberInCluster(1);
        cf.setToUse8Neighbors();
        cf.findGroups(nonEdgeValue);
        int n = cf.getNumberOfGroups();

        float xLimit;
        if (overrideEdgesLimit0 > -1) {
            xLimit = overrideEdgesLimit0;
        } else {
            float[] sizes = new float[n];
            for (int i = 0; i < n; ++i) {
                PairIntArray group = cf.getXY(i);
                sizes[i] = group.getN();
            }
            float xMin = 0.f;
            float xMax = 500.f;
            int nBins = 20;
            HistogramHolder hist = Histogram.createSimpleHistogram(xMin, xMax,
                nBins, sizes, Errors.populateYErrorsBySqrt(sizes));
            /*
            try {
                hist.plotHistogram(label, label);
            } catch (IOException ex) {
                Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
            }
            */
            int firstMinIdx = -1;
            int firstNonZeroIdx = -1;
            for (int i = 0; i < hist.getXHist().length; ++i) {
                float y = hist.getYHist()[i];
                if (firstNonZeroIdx == -1) {
                    if (y > 0) {
                        firstNonZeroIdx = i;
                    }
                } else {
                    if (y > hist.getYHist()[i - 1]) {
                        firstMinIdx = i;
                        break;
                    }
                }
            }

            //TODO: this may need to be changed
            xLimit = (firstMinIdx > -1) ? hist.getXHist()[firstMinIdx] : 125;
            xLimit /= 5.f;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        List<Set<PairInt>> neighborLists = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < n; ++i) {
            neighborLists.add(new HashSet<PairInt>());
            PairIntArray group = cf.getXY(i);
            if (group.getN() > xLimit) {
                continue;
            }
            Set<PairInt> set = Misc.convert(group);
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    if (img.getValue(x2, y2) != edgeValue) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    if (!set.contains(p2)) {
                        neighborLists.get(i).add(p2);
                        nonNeighborEdges.remove(p2);
                    }
                }
            }
        }

        // for each set in neighborLists, if neighbors are all in neighborLists, can null the pixels
        for (int i = 0; i < neighborLists.size(); ++i) {
            Set<PairInt> set = neighborLists.get(i);
            if (set.size() == 0) {
                continue;
            }
            boolean canNull = true;
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    if (img.getValue(x2, y2) != edgeValue) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    if (nonNeighborEdges.containsKey(p2)) {
                        canNull = false;
                        break;
                    }
                }
                if (!canNull) {
                    break;
                }
            }
            if (canNull) {
                for (PairInt p : set) {
                    int x = p.getX();
                    int y = p.getY();
                    img.setValue(x, y, nonEdgeValue);
                }
            }
        }

        overrideEdgesLimit1 = (overrideEdgesLimit1 == -1) ? 12 :
            overrideEdgesLimit1;

        // null small sets in cf0
        for (int i = 0; i < cf0.getNumberOfGroups(); ++i) {
            PairIntArray group = cf0.getXY(i);
            if (group.getN() > overrideEdgesLimit1) {
                continue;
            }
            for (int j = 0; j < group.getN(); ++j) {
                int x = group.getX(j);
                int y = group.getY(j);
                img.setValue(x, y, nonEdgeValue);
            }
        }
    }
    
    /**
     * NOT READY FOR USE.
     * create edges for img using phase congruency edges and then sparse color 
     * gradients (b-g, g-b, and r-b) to complete the curves.  The results 
     * contain closed curves, but need to be followed by merging similar color
     * cells.
     * Note that the returned result contains values 0 or 255 for
     * easier display.  This returned format in the future will likely hold 
     * just binary 0 or 1 in a more compact internal structure in the 
     * GreyscaleImage.
     * @param img
     * @return edge image holding values 0 or 255
     */
    public GreyscaleImage createColorEdges(Image img) {
      
        GreyscaleImage gsImg = img.copyBlueToGreyscale();

        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 10;//5;//2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.0001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = false;//true;
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        final int w = img.getWidth();
        final int h = img.getHeight();
           
        // half width of neighbor region
        final int hN = 2;//3;

        PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseDetector.phaseCongMono(gsImg, nScale, minWavelength, mult, 
            sigmaOnf, k, increaseKIfNeeded,
            cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
        
        int[][] thinned = products.getThinned();
        {
            GreyscaleImage out2 = gsImg.createWithDimensions();
            for (int i = 0; i < thinned.length; ++i) {
                for (int j = 0; j < thinned[i].length; ++j) {
                    if (thinned[i][j] > 0) {
                        out2.setValue(j, i, 255);
                    }
                }
            }
            MiscDebug.writeImage(out2, "_EDGES_grey_");
        }
                
        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage edgeImg = null;
       
        /*
        sparse color gradients:
            0  r-g
            1  b-g
            2  r-b
        */
        for (int clrIdx = 0; clrIdx < 3; ++clrIdx) {
        
            EdgeExtractorSimple extractor = new EdgeExtractorSimple(thinned);
            extractor.extractEdges();
            List<PairIntArray> edgeList = extractor.getEdges();
            Set<PairInt> junctions = extractor.getJunctions();

            Map<PairInt, Integer> edgeIndexMap = new HashMap<PairInt, Integer>();
            for (int i = 0; i < edgeList.size(); ++i) {
                PairIntArray curve = edgeList.get(i);
                Integer index = Integer.valueOf(i);
                for (int j = 0; j < curve.getN(); ++j) {
                    int x = curve.getX(j);
                    int y = curve.getY(j);
                    curve.set(j, y, x);

                    PairInt p = new PairInt(y, x);
                    edgeIndexMap.put(p, index);
                }
            }
            {
                // junctions are using format a[row][col] so change to GreyscaleImage col, row
                Set<PairInt> tmp = new HashSet<PairInt>();
                for (PairInt p : junctions) {
                    int x = p.getX();
                    int y = p.getY();
                    tmp.add(new PairInt(y, x));
                }
                junctions.clear();
                junctions.addAll(tmp);
            }

            Set<PairInt> added = new HashSet<PairInt>();

            // to avoid processing some of the noise, there's a minimum line 
            // length for curves used to initialize the stack,
            // but because the edgeextractorsimple prefers to keep curves with
            // junctions separated by the junction, those curves adjacent to 
            // junctions, even when small in length, should be processed.

            Set<Integer> adjToJunction = new HashSet<Integer>();
            for (int i = 0; i < edgeList.size(); ++i) {
                PairIntArray curve = edgeList.get(i);
                if ((curve instanceof PairIntArrayWithColor)
                    && ((PairIntArrayWithColor) curve).isClosedCurve()) {
                    continue;
                }
                int n = curve.getN();
                boolean foundJunction = false;
                for (int j = 0; j < 2; ++j) {
                    PairInt p = (j == 0) ? 
                        new PairInt(curve.getX(0), curve.getY(0)) :
                        new PairInt(curve.getX(n - 1), curve.getY(n - 1));
                    for (int dIdx = 0; dIdx < dxs.length; ++dIdx) {
                        int x2 = p.getX() + dxs[dIdx];
                        int y2 = p.getY() + dys[dIdx];
                        PairInt p2 = new PairInt(x2, y2);                    
                        if (junctions.contains(p2)) {
                            adjToJunction.add(Integer.valueOf(i));
                            foundJunction = true;
                            break;
                        } else {
                            Integer index2 = edgeIndexMap.get(p2);
                            if (index2 != null && index2.intValue() != i) {
                                adjToJunction.add(Integer.valueOf(i));
                                foundJunction = true;
                                break;
                            }
                        }
                    }
                    if (foundJunction) {
                        break;
                    }
                }
            }
        
            Set<PairInt> visited = new HashSet<PairInt>();
            Stack<PairInt> stack = new Stack<PairInt>();
            for (int i = 0; i < edgeList.size(); ++i) {
                PairIntArray curve = edgeList.get(i);
                int n = curve.getN(); 

                if ((curve instanceof PairIntArrayWithColor)
                    && ((PairIntArrayWithColor) curve).isClosedCurve()) {
                    continue;
                }
                if ((n > 4) || adjToJunction.contains(Integer.valueOf(i))) {
                    PairInt p = new PairInt(curve.getX(0), curve.getY(0));
                    stack.add(p);
                    if (n > 1) {
                        p = new PairInt(curve.getX(n - 1), curve.getY(n - 1));
                        stack.add(p);
                    }
                }
            }

            GreyscaleImage outputGradients = gsImg.createFullRangeIntWithDimensions();
            outputGradients.fill(-1);

            int sz = 3;
            int kHL = 3;

            while (!stack.isEmpty()) {
                PairInt p = stack.pop();

                if (visited.contains(p)) {
                    continue;
                }
                visited.add(p);
                Integer index = edgeIndexMap.get(p);
                final int x = p.getX();
                final int y = p.getY();
                // quick check to see if this point is adjacent to another edge.
                // and if so, do not process further.
                // also, if any point is out of bounds, can skip processing

                if (aMemberIsOutOfBounds(x, y, hN, w, h)) {
                    continue;
                }
                if (foundAdjacentEdge(x, y, edgeIndexMap, index, 1)) {
                    continue;
                }

                if ((x < (sz + kHL)) || (y < (sz + kHL)) || 
                    (x > (img.getWidth() - 1 - (sz + kHL))) ||
                    (y > (img.getHeight() - 1 - (sz + kHL)))) {
                    // TODO: need to adjust the sparse gradient method 
                    // to handle points ner image boundaries
                    continue;
                }

                CannyEdgeFilterLite cnf = new CannyEdgeFilterLite();

                int len = (2 * (sz + kHL)) + 1;
                float[] values = new float[len * len];
                int count = 0;
                int startX = (x - sz - kHL);
                int endX = (x + sz + kHL);
                int startY = (y - sz - kHL);
                int endY = (y + sz + kHL);
                for (int xp = startX; xp <= endX; ++xp) {
                    for (int yp = startY; yp <= endY; ++yp) {
                        if (clrIdx == 0) {
                            values[count] = img.getR(xp, yp) - img.getG(xp, yp);
                        } else if (clrIdx == 1) {
                            values[count] = img.getB(xp, yp) - img.getG(xp, yp);
                        } else if (clrIdx == 2) {
                            values[count] = img.getR(xp, yp) - img.getB(xp, yp);
                        }

                        count++;
                    }
                }
                values = MiscMath.rescale(values, 0, 255);
                Map<PairInt, Integer> scaledValues = new HashMap<PairInt, Integer>();
                count = 0;
                for (int xp = startX; xp <= endX; ++xp) {
                    for (int yp = startY; yp <= endY; ++yp) {
                        int v = Math.round(values[count]);
                        scaledValues.put(new PairInt(xp, yp), v);
                        count++;
                    }
                }

                cnf.applyFilterToRegion(scaledValues, outputGradients,
                    x - sz, y - sz, x + sz, y + sz);
              
                Stack<PairInt> stack2 = new Stack<PairInt>();
                stack2.add(new PairInt(x, y));

                PairInt lastPAdded = null;

                while (!stack2.isEmpty()) {
                    PairInt pG = stack2.pop();
                    int max = Integer.MIN_VALUE;
                    PairInt maxP = null;
                    for (int dIdx = 0; dIdx < dxs.length; ++dIdx) {
                        int x2 = pG.getX() + dxs[dIdx];
                        int y2 = pG.getY() + dys[dIdx];
                        PairInt p2 = new PairInt(x2, y2);
                        if (edgeIndexMap.containsKey(p2)) {
                            continue;
                        }
                        if (countNeighbors(edgeIndexMap, x2, y2) > 2) {
                            continue;
                        }
                        int v2 = outputGradients.getValue(x2, y2);
                        if ((v2 > 0) && (v2 > max)) {
                            max = v2;
                            maxP = p2;
                        }
                    }
                    if (maxP != null) {
                        edgeIndexMap.put(maxP, index);
                        added.add(maxP);
                        lastPAdded = maxP;

                        // if a point on the boundary of gradient region is found, exit now
                        if ((Math.abs(maxP.getX() - x) == (sz - 1)) || 
                            (Math.abs(maxP.getY() - y) == (sz - 1))) {
                            break;
                        }
                        stack2.add(maxP);
                    }
                }
                if (lastPAdded != null) {
                    stack.add(lastPAdded);
                }
            }
            System.out.println("number of points added = " + added.size());

            for (PairInt p : added) {
                int x = p.getX();
                int y = p.getY();
                thinned[y][x] = 1;
            }
            
            {
                edgeImg = gsImg.createWithDimensions();
                for (int i = 0; i < thinned.length; ++i) {
                    for (int j = 0; j < thinned[i].length; ++j) {
                        if (thinned[i][j] > 0) {
                            edgeImg.setValue(j, i, 255);
                        }
                    }
                }
                MiscDebug.writeImage(edgeImg, "_EDGES_2_" + clrIdx + "_");
            }
        }
        
        return edgeImg;
    }
    
    /**
     * NOT READY FOR USE.
     * create edges for img using phase congruency edges for grey 
     * (b-g, g-b, and r-b).
     * 
     * Note that the returned result contains values 0 or 255 for
     * easier display.  This returned format in the future will likely hold 
     * just binary 0 or 1 in a more compact internal structure in the 
     * GreyscaleImage.
     * @param img
     * @param cellSize if 0 or less is full image, else cellSize is the 
     * length in x and the length in y to perform phase congruency on in a 
     * tiled, overlapping manner.
     * @return edge image holding values 0 or 255
     */
    public GreyscaleImage createColorEdges_1(Image img, int cellSize) {
        
        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;//4
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 10;//2;
        float g = 10;
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = false;
        
        int width = img.getWidth();
        int height = img.getHeight();
        int x0 = 0;
        int x1 = width;
        int y0 = 0;
        int y1 = height;

        int szX, szY, dXOff, dYOff;
        if (cellSize < 1) {
            szX = width;
            szY = height;
            dXOff = 0;
            dYOff = 0;
        } else {
            szX = cellSize;
            szY = cellSize;
            dXOff = 15;
            dYOff = 15;
        }
        
        GreyscaleImage combined = new GreyscaleImage(width, height);
        for (int xOff = 0; xOff < width; xOff += szX) {
            if (xOff > 0) {
                xOff -= dXOff;
            }
            for (int yOff = 0; yOff < height; yOff += szY) {
                if (yOff > 0) {
                    yOff -= dYOff;
                }
                for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {
                    float[] values = null;
                    int count = 0;
                    if (clrIdx == 0) {
                        values = read_R_G_B(img, x0 + xOff, x0 + xOff + szX, 
                            y0 + yOff, y0 + yOff + szY);
                    } else if (clrIdx == 1) {
                        values = read_R_minus_G(img, x0 + xOff, x0 + xOff + szX, 
                            y0 + yOff, y0 + yOff + szY);
                    } else if (clrIdx == 2) {
                        values = read_B_minus_G(img, x0 + xOff, x0 + xOff + szX, 
                            y0 + yOff, y0 + yOff + szY);
                    } else if (clrIdx == 3) {
                        values = read_R_minus_B(img, x0 + xOff, x0 + xOff + szX, 
                            y0 + yOff, y0 + yOff + szY);
                    }
                    values = MiscMath.rescale(values, 0, 255);

                    GreyscaleImage img2 = new GreyscaleImage(szX, szY);
                    count = 0;
                    int e0 = szX;
                    if ((x0 + xOff + szX) > x1) {
                        e0 = x1 - x0 - xOff;//sz - (x0 + xOff + sz - x1);
                    }
                    int e1 = szY;
                    if ((y0 + yOff + szY) > y1) {
                        e1 = szY - (y0 + yOff + szY - y1);
                    }
                    for (int i = 0; i < e0; ++i) {
                        for (int j = 0; j < e1; ++j) {
                            int v = Math.round(values[count]);
                            img2.setValue(i, j, v);
                            count++;
                        }
                    }

                    PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
                    PhaseCongruencyDetector.PhaseCongruencyProducts products
                        = phaseDetector.phaseCongMono(img2, nScale, minWavelength, mult,
                            sigmaOnf, k, increaseKIfNeeded,
                            cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
                    int[][] thinned = products.getThinned();
                    assert(thinned.length == szX);
                    assert(thinned[0].length == szY);
                    for (int j = 0; j < thinned.length; ++j) {
                        for (int i = 0; i < thinned[j].length; ++i) {
                            if (((i + xOff) > (combined.getWidth() - 1)) ||
                                ((j + yOff) > (combined.getHeight() - 1))) {
                                continue;
                            }
                            if (thinned[j][i] > 0) {
                                combined.setValue(i + xOff, j + yOff, 255);
                            }
                        }
                    }
                }
            }
        }
        
        return combined;
    }
    
    /**
     * NOT READY FOR USE.
     * create edges for img the maximum of sobel edges that are grey 
     * (b-g, g-b, and r-b).
     * 
     * Note that the returned result contains values 0 or 255 for
     * easier display.  This returned format in the future will likely hold 
     * just binary 0 or 1 in a more compact internal structure in the 
     * GreyscaleImage.
     * @param img
     * @return edge image holding values 0 or 255
     */
    public GreyscaleImage createColorEdges_2(Image img) {
        
        GreyscaleImage[] gradients = new GreyscaleImage[4];
          
        /*
        0 grey    min:    0    max: 255 
        1 r-g     min: -255    max: 255 
        2 b-g       "
        3 r-b       "
        */
        for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {
            gradients[clrIdx] = new GreyscaleImage(img.getWidth(), 
                img.getHeight(), GreyscaleImage.Type.Bits32FullRangeInt);
            if (clrIdx == 0) {
                for (int i = 0; i < img.getNPixels(); ++i) {
                    int v = (img.getR(i) + img.getG(i) + img.getB(i))/3;
                    gradients[clrIdx].setValue(i, v);
                }
            } else if (clrIdx == 1) {
                for (int i = 0; i < img.getNPixels(); ++i) {
                    int v = img.getR(i) - img.getG(i);
                    v = (v + 255)/2;
                    gradients[clrIdx].setValue(i, v);
                }
            } else if (clrIdx == 2) {
                for (int i = 0; i < img.getNPixels(); ++i) {
                    int v = img.getB(i) - img.getG(i);
                    v = (v + 255)/2;
                    gradients[clrIdx].setValue(i, v);
                }
            } else if (clrIdx == 3) {
                for (int i = 0; i < img.getNPixels(); ++i) {
                    int v = img.getR(i) - img.getB(i);
                    v = (v + 255)/2;
                    gradients[clrIdx].setValue(i, v);
                }
            }               
            CannyEdgeFilterLite filter = new CannyEdgeFilterLite();
            filter.setToUseSobel();
            filter.setToUseOtsu();
            filter.applyFilter(gradients[clrIdx]);            
        }

        GreyscaleImage combined =  new GreyscaleImage(img.getWidth(), 
            img.getHeight(), GreyscaleImage.Type.Bits32FullRangeInt);
        for (int i = 0; i < img.getNPixels(); ++i) {
            int gradientMax = Integer.MIN_VALUE;
            int maxClrIdx = -1;
            for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {
                int v = gradients[clrIdx].getValue(i);
                if (Math.abs(v) > Math.abs(gradientMax)) {
                    gradientMax = v;
                    maxClrIdx = clrIdx;
                }
            }
            combined.setValue(i, gradientMax);
        }
        
        return combined;
    }
    
    private float[] read_R_G_B(Image img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int r = img.getR(x, y);
                int g = img.getG(x, y);
                int b = img.getB(x, y);
                values[count] = r + g + b;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private float[] read_R_minus_B(Image img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int r = img.getR(x, y);
                int b = img.getB(x, y);
                values[count] = r - b;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private float[] read_B_minus_G(Image img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int g = img.getG(x, y);
                int b = img.getB(x, y);
                values[count] = b - g;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private float[] read_R_minus_G(Image img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int r = img.getR(x, y);
                int g = img.getG(x, y);
                values[count] = r - g;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private int[][] copy(int[][] a) {
        
        int[][] b = new int[a.length][];
        for (int i = 0; i < b.length; ++i) {
            b[i] = Arrays.copyOf(a[i], a[i].length);
        }
        
        return b;
    }
}
