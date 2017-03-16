package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusColor;
import algorithms.connected.ConnectedValuesFinder;
import algorithms.imageProcessing.util.GroupAverageColors;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.features.UnsupervisedTextureFinder;
import algorithms.imageProcessing.features.UnsupervisedTextureFinder.TexturePatchesAndResponse;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MedianSmooth;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor1D;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.ResourceFinder;
import algorithms.util.TwoDIntArray;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTree;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

/**
 * Many methods in here will be removed soon.
 * Meanwhile, see MSEREdges.java for segmentation.
 * 
 * class holding several different image segmentation methods.  Note that
 * some other techniques involving contrast for example, are elsewhere.
 *
 * A few of the methods use a density based clustering algorithm from
       http://nking.github.io/two-point-correlation/
       which has an MIT license
      ---- begin nking copyright ----
      The MIT License (MIT)
      Copyright (c) 2013-* Nichole King
      http://nking.github.io/two-point-correlation/

        Permission is hereby granted, free of charge, to any person obtaining 
        a copy of this software and associated documentation files 
        (the "Software"), to deal in the Software without restriction, 
        including without limitation the rights to use, copy, modify, merge, 
        publish, distribute, sublicense, and/or sell copies of the Software, 
        and to permit persons to whom the Software is furnished to do so, 
        subject to the following conditions:

        The above copyright notice and this permission notice shall be included 
        in all copies or substantial portions of the Software.
        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
        OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
        MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
        IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
        CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
        SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
     ---- end nking copyright ---- 
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

        KMeansPlusPlus kmpp = new KMeansPlusPlus();
        kmpp.computeMeans(kBands, input);

        int[] seeds = kmpp.getCenters();

        int[] imgSeedIndexAssignments = kmpp.getImgPixelSeedIndexes();

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
                    if (x2 < 0 || y2 < 0 || 
                        (x2 > (w - 1)) || (y2 > (h - 1))) {
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
     *
     * NOTE: parameters in this algorithm are sensitive to
     * the PSF.
     *
     * NOTE: this implementation doesn't reproduce their results as
     * precisely so needs some improvements.
     *
     * @param input
     * @return
     */
    public List<Set<PairInt>> createColorEdgeSegmentation(ImageExt input,
        String debugTag) {
        
        // 0 is CIE LAB, 1 is HSV
        final int clrSpace = 0;

        boolean reduceNoise = false;

        double tColor;
        int tLen;
        double tR;
        double tSmallMerge;
        if (clrSpace == 0) {
            // JND for deltaE is ~2.3, so tColor must be that or larger
            tColor = 2.8;//4.0;//5.5;
            tR = 0.8;//1.0;
            tLen = 1;
            tSmallMerge = 0.02;//0.095;
        } else {
            // what is JND for HSV (a.k.a. HSB) ?  each range of values is 0:1
            tColor =  0.125;//0.125;  between 0.1 and 0.175
            tR = 1.5;
            tLen = 5;
            tSmallMerge = 0.02;
        }

        return createColorEdgeSegmentation(input, clrSpace, tLen, tColor, tR,
            reduceNoise, tSmallMerge, debugTag);
    }

    public List<PairIntArray> extractEdges(Image img,
        boolean reduceNoise, String debugTag) {

        if (debugTag == null) {
            debugTag = "";
        }

        GreyscaleImage gsImg = img.copyBlueToGreyscale();

        PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
        if (reduceNoise) {
            phaseDetector.setK(5);
        } else {
            phaseDetector.setK(2);
        }
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseDetector.phaseCongMono(gsImg);

        // thinned is in row major format
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
            MiscDebug.writeImage(out2, "_EDGES_grey_" + debugTag);
        }

        EdgeExtractorSimple extractor = new EdgeExtractorSimple(thinned);
        extractor.extractEdges();
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        // put in framework of images
        for (int i = 0; i < extractor.getEdges().size(); ++i) {
            PairIntArray edge = extractor.getEdges().get(i).copy();
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                edge.set(j, y, x);
            }
            edges.add(edge);
        }

        return edges;
    }
    
    public List<PairIntArray> extractEdges2(Image img,
        String debugTag) {

        if (debugTag == null) {
            debugTag = "";
        }

        CannyEdgeFilterAdaptiveDeltaE2000 canny = 
            new CannyEdgeFilterAdaptiveDeltaE2000();
        canny.applyFilter(img.copyToImageExt());

        GreyscaleImage gXY = canny.getFilterProducts().getGradientXY();
        
        int w = gXY.getWidth();
        int h = gXY.getHeight();
        
        assert(img.getWidth() == w);
        assert(img.getHeight() == h);
        
        int[][] thinned = new int[w][];
        for (int i = 0; i < w; ++i) {
            thinned[i] = new int[h];
            for (int j = 0; j < h; ++j) {
                thinned[i][j] = gXY.getValue(i, j);
            }
        }

        EdgeExtractorSimple extractor = new EdgeExtractorSimple(thinned);
        extractor.extractEdges();
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        // put in framework of images
        for (int i = 0; i < extractor.getEdges().size(); ++i) {
            PairIntArray edge = extractor.getEdges().get(i).copy();
            for (int j = 0; j < edge.getN(); ++j) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                edge.set(j, x, y);
            }
            edges.add(edge);
        }

        return edges;
    }

    /**
     * create segmented image by creating edges with phase congruence,
     * then using the edge color properties to form clusters and seeds
     * of regions to grow, then using color histograms to further merge
     * regions.
     * The algorithm follows the general outline given by
     * Jie and Peng-fei 2003, "Natural Color Image Segmentation",
       http://www-labs.iro.umontreal.ca/~mignotte/IFT6150/Articles/TRASH/ARTICLES_2010/cr1231.pdf

     NOTE: this implementation doesn't reproduce their results as
     * precisely so needs some improvements.
     * 
     * @param input
     * @return
     */
    public List<Set<PairInt>> createColorEdgeSegmentation(ImageExt input,
        int clrSpace, int tLen, double tColor, double tR, boolean reduceNoise,
        double tSmallMerge, String debugTag) {

        List<PairIntArray> edges = extractEdges2(input, debugTag);

        //List<PairIntArray> edges = extractEdges(input, reduceNoise, debugTag);
        
        return createColorEdgeSegmentation(input, edges,
            clrSpace, tLen, tColor, tR, reduceNoise, tSmallMerge, debugTag);
    }

    /**
     * create segmented image by creating edges with phase congruence,
     * then using the edge color properties to form clusters and seeds
     * of regions to grow, then using color histograms to further merge
     * regions.
     * The algorithm follows the general outline given by
     * Jie and Peng-fei 2003, "Natural Color Image Segmentation",
       http://www-labs.iro.umontreal.ca/~mignotte/IFT6150/Articles/TRASH/ARTICLES_2010/cr1231.pdf

     NOTE: this implementation doesn't reproduce their results as
     * precisely so needs some improvements.
     * 
     * @param input
     * @return
     */
    public List<Set<PairInt>> createColorEdgeSegmentation(ImageExt input,
        List<PairIntArray> edges,
        int clrSpace, int tLen, double tColor, double tR, boolean reduceNoise,
        double tSmallMerge, String debugTag) {

        boolean doPlot = false;

        if (debugTag == null) {
            debugTag = "";
        }

        final int w = input.getWidth();
        final int h = input.getHeight();
        final int nPix = input.getNPixels();

        int nEdges = edges.size();
        List<Set<PairInt>> clusterPoints = new ArrayList<Set<PairInt>>();

        if (nEdges == 0) {
            // add all picels to one set
             Set<PairInt> set = new HashSet<PairInt>();
            for (int i = 0; i < w; ++i) {
                for (int j = 0; j < h; ++j) {
                    set.add(new PairInt(i, j));
                }
            }
            clusterPoints.add(set);
            return clusterPoints;
        }

        float[][] clusterDescriptors = new float[nEdges][];

        populateEdgeLists(input, edges, clusterPoints, clusterDescriptors,
            clrSpace);

        assert(clusterPoints.size() == clusterDescriptors.length);

        List<Integer> longEdgeIndexes = new ArrayList<Integer>();
        List<Integer> shortEdgeIndexes = new ArrayList<Integer>();
        populateEdgeLengthLists(clusterDescriptors, tLen, longEdgeIndexes,
            shortEdgeIndexes);

        assert(clusterPoints.size() == clusterDescriptors.length);

        // ----------- merge long edges ----------

        // NOTE that the moved sets modify the data structures :
        //    clusterPoints may contain empty items
        //    clusterDescriptors may contain null items
        //    both clusterPoints and clusterDescriptor non- null and non empty
        //       items are updated for merges

        // the authors consider this algorithm of min-heap merging within a
        // radius of color, a kmeans method as it updates the descriptors upon
        // each merge, but the minimum distance ordering is an improvement over
        // standard kmeans if one can use it as one can here
        // (pairs results in outer loop iteration of approx O(N^2),
        // specifically (N*(N-1)/2)), while kmeans ordering by index uses O(N))
        mergeEdges(clusterPoints, clusterDescriptors, clrSpace, tColor,
            longEdgeIndexes);

        //TODO: consider a number limit to use an alternate here when
        //  n edges is a large number.  determine a fixed k and use kmeans.
        //  can roughly determine a fixed k from
        //  a color histogram with bin size being color tolerance
        //  and counting the number of peaks.

        assert(clusterPoints.size() == clusterDescriptors.length);

        if (doPlot) {
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
            MiscDebug.writeImage(img2, "_longEdges_merged_" +  debugTag + "_"
                + clrSpace);
        }

        // ---- merge short edges (which are usually textures) ------

        mergeShortEdges(clusterPoints, clusterDescriptors, clrSpace, tColor,
            shortEdgeIndexes);

        assert(clusterPoints.size() == clusterDescriptors.length);

        if (doPlot) {
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
            MiscDebug.writeImage(img2, "_shortedges_merged_" +  debugTag +
                "_" + clrSpace);
        }

        assert(assertDescriptorCounts(clusterPoints, clusterDescriptors));

        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < clusterPoints.size(); ++i) {
            Set<PairInt> set = clusterPoints.get(i);
            Integer index = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, index);
            }
        }

        // ------ region growing -------
        growEdges(input, clusterPoints, clusterDescriptors, pointIndexMap,
            clrSpace, tColor, shortEdgeIndexes, longEdgeIndexes);

        assert(clusterPoints.size() == clusterDescriptors.length);

        longEdgeIndexes = null;
        shortEdgeIndexes = null;

        if (doPlot) {
            // DEBUG
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(clusterPoints, 0, 0,
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_after_rgo_" +  debugTag + "_" + clrSpace);
        }

        // ------ merge by color histograms ------

        clusterDescriptors = condenseAndUpdate(clusterPoints,
            clusterDescriptors, pointIndexMap);

        assert(clusterPoints.size() == clusterDescriptors.length);

        Map<Integer, Set<Integer>> adjacencyMap = createAdjacencyMap(
            clusterPoints);

        mergeByColorHistograms(input, clusterPoints, adjacencyMap,
            clrSpace, tR);

        assert(clusterPoints.size() == clusterDescriptors.length);

        if (doPlot) {
            // DEBUG
            int nExtraForDot = 1;
            Image img2 = input.copyImage().copyToGreyscale().copyToColorGreyscale();
            ImageIOHelper.addAlternatingColorPointSetsToImage(clusterPoints, 0, 0,
                nExtraForDot, img2);
            MiscDebug.writeImage(img2, "_FINAL_" +  debugTag + "_" + clrSpace);
        }

        // ----- merge smallest clusters into adjacent larger --------
        int tNumber = (int)Math.round(tSmallMerge * nPix);

        clusterDescriptors = condenseAndUpdate(clusterPoints,
            clusterDescriptors, pointIndexMap);

        assert(clusterPoints.size() == clusterDescriptors.length);

        mergeSmallClusters(input, clusterPoints, clusterDescriptors,
            clrSpace, tNumber, debugTag);

        assert(clusterPoints.size() == clusterDescriptors.length);

        return clusterPoints;
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
     * expecting binary image input where a pixel of value "1" is significant.
     * if any of the pixels in gapsOf1 are completely surrounded by
     * pixels of value "1" in their 8 pixel neighborhood,
     * those gap pixels are set to "0".
     * 
     * @param img binary image
     * @param gapsOf1 pixel indexes of filled in gaps (pixels set to "1")
     * @param value 
     */
    public void restoreGapsOf1WhereSurrounded(GreyscaleImage img,
        TIntSet gapsOf1, int value) {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        int w = img.getWidth();
        int h = img.getHeight();
        TIntSet reset = new TIntHashSet();
        
        TIntIterator iter = gapsOf1.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/w;
            int x = pixIdx - (y * w);
            int n1s = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || x2 >= w || y2 >= h) {
                    continue;
                }
                if (img.getValue(x2, y2) == 1) {
                    n1s++;
                }
            }
            if (n1s == dxs.length) {
                reset.add(pixIdx);
            }
        }
        iter = reset.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            img.setValue(pixIdx, 0);
        } 
        
    }

    /**
     * expecting binary image input where a pixel of value "1" is significant.
     * The algorithm is similar to dilation, in that if any pixel has
     * a gap of size 1 pixel in between itself and another, that gap is
     * set to value 0 and stored in outputAddedGaps.
     * 
     * @param img
     * @param outputAddedGaps
     * @param value
     * @return 
     */
    public GreyscaleImage fillInGapsOf1(GreyscaleImage img,
        TIntSet outputAddedGaps, int value) {

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
                    int pixIdx = (j * w) + i;
                    outputAddedGaps.add(pixIdx);
                    break;
                }
            }
        }

        return tmpImg2;
    }

    public GreyscaleImage fillInCompleteGapsOf1(GreyscaleImage img,
        TIntSet outputAddedGaps,int value) {

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
                        int pixIdx = (j * w) + i;
                        outputAddedGaps.add(pixIdx);
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

        log.fine(edgeIndexes.size() + " edges");

        if (edgeIndexes.isEmpty()) {
            return;
        }

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

            Integer index1 = Integer.valueOf(idx1);
            Integer index2 = Integer.valueOf(idx2);

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

            // remove the idx1 --> set<integer> pairs from map and heap
            // remove the idx2 --> set<integer> pairs from map and heap

            Set<Integer> indexes3 = indexToIndexMap.get(index1);
            for (Integer index3 : indexes3) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                if (set3.isEmpty() || idx1 == idx3) {
                    continue;
                }
                //keys in pairEdgePindexNodes have smaller index in x
                PairInt p13;
                if (idx1 < idx3) {
                    p13 = new PairInt(idx1, idx3);
                } else {
                    p13 = new PairInt(idx3, idx1);
                }
                HeapNode node3 = pairEdgePindexNodes.get(p13);
                assert(node3 != null);
                heap.remove(node3);
                pairEdgePindexNodes.remove(p13);
            }

            indexes3 = indexToIndexMap.get(index2);
            for (Integer index3 : indexes3) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                if (set3.isEmpty() || idx1 == idx3) {
                    continue;
                }
                PairInt p23;
                if (idx2 < idx3) {
                    p23 = new PairInt(idx2, idx3);
                } else {
                    p23 = new PairInt(idx3, idx2);
                }
                HeapNode node3 = pairEdgePindexNodes.get(p23);
                assert(node3 != null);
                heap.remove(node3);
                pairEdgePindexNodes.remove(p23);
            }

            // update the indexToIndexMap
            Set<Integer> iim = new HashSet<Integer>(indexToIndexMap.get(index2));
            for (Integer index3 : iim) {
                if (!index3.equals(index1)) {
                    Set<Integer> indexes4 = indexToIndexMap.get(index3);
                    if (indexes4 != null) {
                        indexes4.remove(index2);
                        indexes4.add(index1);
                    }
                }
            }
            indexes3.addAll(iim);
            indexes3.remove(index2);
            indexToIndexMap.remove(index2);

            // add node for updated idx1 ---> set<integer> to map and node
            indexes3 = indexToIndexMap.get(Integer.valueOf(idx1));
            for (Integer index3 : indexes3) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                if (set3.isEmpty() || idx1 == idx3) {
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
                HeapNode node3 = new HeapNode(heapKey);
                node3.setData(p13);
                heap.insert(node3);
                pairEdgePindexNodes.put(p13, node3);
            }

            // pairs having set2 will be skipped because of the empty set at beginning of while loop
        }
    }

    private void growEdges(ImageExt img,
        List<Set<PairInt>> clusterPoints,
        float[][] clusterDescriptors,
        Map<PairInt, Integer> pointIndexMap, int clrSpace, double tColor,
        List<Integer> shortEdgeIndexes, List<Integer> longEdgeIndexes) {

        if (pointIndexMap.isEmpty()) {
            return;
        }

        /*
        traverse all image points to make a map of unassigned pixels and the
            cluster indexes they are adjacent to if any.
            --> O(N)

        initialize an outer queue with the unassigned which have nIndexes > 0
           sorted by descending number of adjacenct indexes
           --> O(N_i * lg2(N_i))

        create an inner queue

        --> O(|V| + |E|)
        visit each outer queue member, adding it to adjacent cluster
           which has most similar color.
           update the cluster's descriptor
           add each of the 8 neighbors which aren't assigned to the inner
              queue.
           add the point to the visited set.
        when all outer queue members have been visited, fill the outer queue
            with the inner queue and empty the inner queue.
            (could sort again here for more precise growing)
        continue in this manner until the inner queue is empty and hence outer
            queue is empty.
            assert that visited.size == unassigned map.size

        The same pattern should be applied elsewhere too
        */

        // for each edge, add neighbors with diff < tColor
        //    if an adjacent pixel is part of a short edge cluster,
        //    then all of that short edge is added to the cluster

        // removing short edge points from pointIndexMap and creating
        // shortPointIndexMap to more easily add them as a whole
        Map<PairInt, Integer> shortPointIndexMap = new HashMap<PairInt, Integer>();
        for (Integer index : shortEdgeIndexes) {
            int idx = index.intValue();
            Set<PairInt> set = clusterPoints.get(idx);
            for (PairInt p : set) {
                shortPointIndexMap.put(p, index);
                pointIndexMap.remove(p);
            }
        }

        CIEChromaticity cieC = new CIEChromaticity();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int width = img.getWidth();
        int height = img.getHeight();

        /*
        on first iteration, the edge regions are grown to include adjacent points
        that are within tColor tolerance and short edges which are ajacent
        regardless of color difference.

        on second iteration, unassigned pixels are added to adjacent indexes
        most similar in color.
        */

        int lastInnerQ = 0;

        for (int nIter = 0; nIter < 2; ++nIter) {

            Map<PairInt, Set<Integer>> unassignedAndIndexes =
                findUnassignedPixelsAndAdjacentIndexes(img, pointIndexMap);

            if (unassignedAndIndexes.isEmpty()) {
                return;
            }

            assert(img.getNPixels() == (unassignedAndIndexes.size() +
                pointIndexMap.size()));

            int count = 0;
            for (Entry<PairInt, Set<Integer>> entry : unassignedAndIndexes.entrySet()) {
                if (entry.getValue().size() > 0) {
                    ++count;
                }
            }
            PairInt[] unassigned = new PairInt[count];
            int[] nAdjIndexes = new int[unassigned.length];
            count = 0;
            for (Entry<PairInt, Set<Integer>> entry :
                unassignedAndIndexes.entrySet()) {
                if (entry.getValue().size() > 0) {
                    unassigned[count] = entry.getKey();
                    nAdjIndexes[count] = entry.getValue().size();
                    ++count;
                }
            }
            QuickSort.sortBy1stArg(nAdjIndexes, unassigned);
            ArrayDeque<PairInt> outerQueue = new ArrayDeque<PairInt>();
            for (int i = (unassigned.length - 1); i > -1; --i) {
                outerQueue.add(unassigned[i]);
            }

            ArrayDeque<PairInt> innerQueue = new ArrayDeque<PairInt>();
            while (true) {
                while (!outerQueue.isEmpty()) {
                    PairInt p = outerQueue.poll();
                    if (pointIndexMap.containsKey(p)) {
                        continue;
                    }
                    int x = p.getX();
                    int y = p.getY();
                    boolean isAShortEdge = shortPointIndexMap.containsKey(p);
                    float[] clrs1;
                    if (isAShortEdge) {
                        clrs1 = clusterDescriptors[shortPointIndexMap.get(p).intValue()];
                    } else {
                        clrs1 = getColors(img, x, y, clrSpace);
                    }

                    double minDiff = Double.MAX_VALUE;
                    Integer minDiffIndex = null;
                    for (int k = 0; k < dxs.length; ++k) {
                        int x2 = x + dxs[k];
                        int y2 = y + dys[k];
                        if (x2 < 0 || y2 < 0 || (x2 > (width - 1))
                            || (y2 > (height - 1))) {
                            continue;
                        }
                        PairInt p2 = new PairInt(x2, y2);
                        Integer index2 = pointIndexMap.get(p2);
                        if (index2 == null) {
                            continue;
                        } else {
                            if (isAShortEdge && index2.equals(shortPointIndexMap.get(p))) {
                                continue;
                            }
                            if (clusterPoints.get(index2.intValue()).isEmpty()) {
                                continue;
                            }
                        }
                        float[] desc2 = clusterDescriptors[index2.intValue()];
                        double diff;
                        if (clrSpace == 0) {
                            diff = Math.abs(cieC.calcDeltaECIE2000(
                                clrs1[0], clrs1[1], clrs1[2], desc2[0], desc2[1], desc2[2]));
                        } else {
                            double diff1 = clrs1[0] - desc2[0];
                            double diff2 = clrs1[1] - desc2[1];
                            double diff3 = clrs1[2] - desc2[2];
                            diff = Math.sqrt(diff1 * diff1 + diff2 * diff2 + diff3 * diff3);
                        }
                        if ((nIter == 0) && !isAShortEdge && (diff > tColor)) {
                            continue;
                        }
                        if (diff < minDiff) {
                            minDiff = diff;
                            minDiffIndex = index2;
                        }
                    }
                    if (minDiffIndex == null) {
                        continue;
                    }
                    int idx2 = minDiffIndex.intValue();
                    Set<PairInt> set2 = clusterPoints.get(idx2);
                    float[] desc2 = clusterDescriptors[idx2];
                    if (isAShortEdge) {
                        Integer index1 = shortPointIndexMap.get(p);
                        int idx1 = index1.intValue();
                        Set<PairInt> set1 = clusterPoints.get(idx1);
                        float[] desc1 = clusterDescriptors[idx1];
                        float n1 = set1.size();
                        float n2 = set2.size();
                        float nTot = n1 + n2;
                        //{h, s, v, nPix, cenX, cenY}
                        for (int k = 0; k < 6; ++k) {
                            if (k == 3) {
                                desc2[k] = nTot;
                            } else {
                                desc2[k] = ((desc2[k] * n2) + desc1[k] * n1)/nTot;
                            }
                        }
                        //add unassigned perimeter of short edge to innerqueue
                        for (PairInt p3 : set1) {
                            shortPointIndexMap.remove(p3);
                            pointIndexMap.put(p3, minDiffIndex);
                            set2.add(p3);
                            for (int k = 0; k < dxs.length; ++k) {
                                int x4 = p3.getX() + dxs[k];
                                int y4 = p3.getY() + dys[k];
                                if (x4 < 0 || y4 < 0 || (x4 > (width - 1))
                                    || (y4 > (height - 1))) {
                                    continue;
                                }
                                PairInt p4 = new PairInt(x4, y4);
                                if (!shortPointIndexMap.containsKey(p4) &&
                                    !pointIndexMap.containsKey(p4)) {
                                    innerQueue.offer(p4);
                                }
                            }
                        }
                        clusterDescriptors[idx1] = null;
                        set1.clear();
                    } else {
                        float n2 = set2.size();
                        float nTot = n2 + 1;
                        //{h, s, v, nPix, cenX, cenY}
                        desc2[0] = ((desc2[0] * n2) + clrs1[0])/nTot;
                        desc2[1] = ((desc2[1] * n2) + clrs1[1])/nTot;
                        desc2[2] = ((desc2[2] * n2) + clrs1[2])/nTot;
                        desc2[3] = nTot;
                        desc2[4] = ((desc2[4] * n2) + x)/nTot;
                        desc2[5] = ((desc2[5] * n2) + y)/nTot;

                        pointIndexMap.put(p, minDiffIndex);
                        set2.add(p);
                        for (int k = 0; k < dxs.length; ++k) {
                            int x4 = x + dxs[k];
                            int y4 = y + dys[k];
                            if (x4 < 0 || y4 < 0 || (x4 > (width - 1))
                                || (y4 > (height - 1))) {
                                continue;
                            }
                            PairInt p4 = new PairInt(x4, y4);
                            if (!shortPointIndexMap.containsKey(p4) &&
                                !pointIndexMap.containsKey(p4)) {
                                innerQueue.offer(p4);
                            }
                        }
                    }
                }
                if (innerQueue.isEmpty()) {
                    if (nIter == 0) {
                        break;
                    }
                    if (!shortPointIndexMap.isEmpty()) {
                        if (shortPointIndexMap.size() == lastInnerQ) {
                            return;
                        }
                        innerQueue.addAll(shortPointIndexMap.keySet());
                        lastInnerQ = shortPointIndexMap.size();
                    } else {
                        lastInnerQ = 0;
                        break;
                    }
                }
                outerQueue.addAll(innerQueue);
                innerQueue.clear();
            }
        }

        assert(assertShortEdgesAreEmpty(shortEdgeIndexes, clusterPoints));

        assert(pointIndexMap.size() == img.getNPixels());
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

        //TODO: consider improvements of this for large shortEdgeIndexes.size()

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
        List<Set<PairInt>> clusterPoints,
        Map<Integer, Set<Integer>> adjacencyMap,
        int clrSpace, double tR) {

        int[][][] colorHistograms = calculateColorHistograms(input,
            clusterPoints, clrSpace);

        // key is index1, index2 where index1 < index2
        Map<PairInt, HeapNode> nodesMap = new HashMap<PairInt, HeapNode>();

        Heap heap = new Heap();

        ColorHistogram ch = new ColorHistogram();

        // the histogram intersection range of values
        //   is 0 : nColors * 1
        // so for 3 colors, expect that max similarity is 3.0.
        // need to merge by higher similarity, so need to invert
        //   the keys.
        // 3 - similairty bcomes the new key.
        // a tR of 0.7*3.0 = 2.1 becomes 0.9 and any values larger than
        //    that are less similar...smalled values are more similar
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

            PairInt p12 = (PairInt)node.getData();

            nodesMap.remove(p12);

            // this is 3.0 - similarity
            double diff = ((double)node.getKey())/((double)heapKeyFactor);

            if (diff > tRInv) {
                break;
            }

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

            Integer index1 = Integer.valueOf(idx1);
            Integer index2 = Integer.valueOf(idx2);

            // remove the idx1 --> set<integer> pairs from map and heap
            // remove the idx2 --> set<integer> pairs from map and heap
            // update the adjacencyMap
            //   add node for updated idx1 ---> set<integer> to map and node

            Set<Integer> indexes1 = adjacencyMap.get(index1);
            for (Integer index3 : indexes1) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                if (set3.isEmpty() || idx1 == idx3) {
                    continue;
                }
                PairInt p13;
                if (idx1 < idx3) {
                    p13 = new PairInt(idx1, idx3);
                } else {
                    p13 = new PairInt(idx3, idx1);
                }
                HeapNode node3 = nodesMap.get(p13);
                assert(node3 != null);
                heap.remove(node3);
                nodesMap.remove(p13);
                assert(heap.getNumberOfNodes() == nodesMap.size());
            }
            Set<Integer> indexes2 = adjacencyMap.get(index2);
            for (Integer index3 : indexes2) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                if (set3.isEmpty() || idx1 == idx3 || idx2 == idx3) {
                    continue;
                }
                PairInt p23;
                if (idx2 < idx3) {
                    p23 = new PairInt(idx2, idx3);
                } else {
                    p23 = new PairInt(idx3, idx2);
                }
                HeapNode node3 = nodesMap.get(p23);
                assert(node3 != null);
                heap.remove(node3);
                nodesMap.remove(p23);
                assert(heap.getNumberOfNodes() == nodesMap.size());
            }

            //update adjacency map
            for (Integer index3 : indexes2) {
                if (!index3.equals(index1)) {
                    Set<Integer> indexes4 = adjacencyMap.get(index3);
                    if (indexes4 != null) {
                        indexes4.remove(index2);
                        indexes4.add(index1);
                    }
                }
            }
            indexes1.addAll(indexes2);
            indexes1.remove(index1);
            indexes1.remove(index2);
            adjacencyMap.remove(index2);

            // add nodes back into heap and map for the updated idx1 --> set<integer>
            for (Integer index3 : indexes1) {
                int idx3 = index3.intValue();
                Set<PairInt> set3 = clusterPoints.get(idx3);
                assert(idx1 != idx3);
                if (set3.isEmpty()) {
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

                float similarity3 = 3.0f - ch.intersection(hist1, hist3);

                long key3 = (long)(similarity3 * (double)heapKeyFactor);
                HeapNode node3 = new HeapNode(key3);
                node3.setData(p13);

                heap.insert(node3);
                nodesMap.put(p13, node3);
                assert(heap.getNumberOfNodes() == nodesMap.size());
            }

            nMerged++;
        }

        log.fine("color histogram nMerged=" + nMerged);
    }

    public int[][][] calculateColorHistograms(ImageExt input,
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

    private Map<PairInt, Set<Integer>> findUnassignedPixelsAndAdjacentIndexes(
        ImageExt img, Map<PairInt, Integer> pointIndexMap) {

        Map<PairInt, Set<Integer>> unassignedAndIndexes =
            new HashMap<PairInt, Set<Integer>>();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int w = img.getWidth();
        int h = img.getHeight();

        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                if (pointIndexMap.containsKey(p)) {
                    continue;
                }
                Set<Integer> adjIndexes = new HashSet<Integer>();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    Integer index2 = pointIndexMap.get(p2);
                    if (index2 != null) {
                        adjIndexes.add(index2);
                    }
                }
                unassignedAndIndexes.put(p, adjIndexes);
            }
        }

        assert(pointIndexMap.size() + unassignedAndIndexes.size() == img.getNPixels());

        return unassignedAndIndexes;
    }

    private void populateDescriptors(ImageExt img,
        List<Set<PairInt>> pointSets,
        float[][] outputDescripors, int clrSpace) {

        /*
        the descriptors are
             C_i = {h, s, v, nPix, cenX, cenY}  or labL, labA, labB for colorSpace = 0
        */
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        for (int i = 0; i < pointSets.size(); ++i) {

            Set<PairInt> edgePoints = pointSets.get(i);

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

    private float[][] condenseAndUpdate(List<Set<PairInt>> clusterPoints,
        float[][] clusterDescriptors, Map<PairInt, Integer> pointIndexMap) {

        int nNonNull = 0;
        for (Set<PairInt> set : clusterPoints) {
            if (!set.isEmpty()) {
                nNonNull++;
            }
        }

        List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();

        float[][] outputDescriptors = new float[nNonNull][];

        for (int i = 0; i < clusterPoints.size(); ++i) {
            Set<PairInt> set = clusterPoints.get(i);
            if (!set.isEmpty()) {
                outputDescriptors[tmp.size()] = clusterDescriptors[i];
                tmp.add(set);
            }
        }
        clusterPoints.clear();
        clusterPoints.addAll(tmp);

        pointIndexMap.clear();
        for (int i = 0; i < clusterPoints.size(); ++i) {
            Set<PairInt> set = clusterPoints.get(i);
            Integer index = Integer.valueOf(i);
            for (PairInt p : set) {
                pointIndexMap.put(p, index);
            }
        }

        return outputDescriptors;
    }

    private void mergeSmallClusters(ImageExt img,
        List<Set<PairInt>> clusterPoints, float[][] clusterDescriptors,
        int clrSpace, int tNumber, String debugTag) {

        if (clusterPoints.isEmpty()) {
            return;
        }

        Map<Integer, Set<Integer>> adjacencyMap = createAdjacencyMap(
            clusterPoints);

        CIEChromaticity cieC = new CIEChromaticity();

        for (Entry<Integer, Set<Integer>> entry : adjacencyMap.entrySet()) {

            Integer index1 = entry.getKey();

            int idx1 = index1.intValue();

            Set<PairInt> set1 = clusterPoints.get(idx1);

            if ((set1.size() > tNumber) || set1.isEmpty()) {
                continue;
            }

            Set<Integer> indexes2 = entry.getValue();

            // merge with closest in color

            float[] desc1 = clusterDescriptors[idx1];

            double minDiff = Double.MAX_VALUE;
            Integer minDiffIndex = null;

            for (Integer index2 : indexes2) {
                int idx2 = index2.intValue();
                Set<PairInt> set2 = clusterPoints.get(idx2);
                if (set2.isEmpty()) {
                    continue;
                }
                float[] desc2 = clusterDescriptors[idx2];
                double diff;
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

            if (minDiffIndex == null) {
                continue;
            }
            int idx2 = minDiffIndex.intValue();
            Set<PairInt> set2 = clusterPoints.get(idx2);
            assert(!set2.isEmpty());
            float[] desc2 = clusterDescriptors[idx2];
            int n2 = set2.size();
            int n1 = set1.size();
            set2.addAll(set1);
            int nTot = set2.size();
            set1.clear();

            //{h, s, v, nPix, cenX, cenY}
            for (int ii = 0; ii < desc1.length; ++ii) {
                desc2[ii] = ((desc2[ii] * n2) - desc1[ii] * n1) / nTot;
            }
            desc2[3] = nTot;

            clusterDescriptors[idx1] = null;
        }
    }

    /**
     * order the indexes by decreasing number of neighbors
     * within pixIndexes.
     * runtime complexity is max(O(N*log_2(N)), O(N*9)).
     * @param pixIndexes
     * @param imgWidth
     * @param imgHeight
     * @return
     */
    private ArrayDeque<Integer> populateByNumberOfNeighbors(
        TIntSet pixIndexes, int imgWidth, int imgHeight) {

        int n = pixIndexes.size();

        ArrayDeque<Integer> output = new ArrayDeque<Integer>(n);
        if (n == 0) {
            return output;
        }

        int[] idxs = new int[n];
        int[] nn = new int[n];

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        int count = 0;
        TIntIterator iter = pixIndexes.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            int y = idx/imgWidth;
            int x = idx - (y * imgWidth);
            int nc = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || x2 > (imgWidth - 1)) {
                    continue;
                }
                if (y2 < 0 || y2 > (imgHeight - 1)) {
                    continue;
                }
                int pixIdx2 = (y2 * imgWidth) + x2;
                if (pixIndexes.contains(pixIdx2)) {
                    nc++;
                }
            }
            idxs[count] = idx;
            nn[count] = nc;
            count++;
        }

        MultiArrayMergeSort.sortByDecr(nn, idxs);

        for (int i = 0; i < nn.length; ++i) {
            output.add(Integer.valueOf(idxs[i]));
        }

        return output;
    }

    private List<GroupPixelRGB> calculateRGB(ImageExt input,
        List<Set<PairInt>> sets) {

        List<GroupPixelRGB> out = new ArrayList<GroupPixelRGB>();

        for (Set<PairInt> set : sets) {
            GroupPixelRGB gp = new GroupPixelRGB(set, input, 0, 0);
            out.add(gp);
        }

        return out;
    }

    private List<Set<PairInt>> copy(
        List<Set<PairInt>> contiguousSets) {

        List<Set<PairInt>> c = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < contiguousSets.size(); ++i) {
            Set<PairInt> set = contiguousSets.get(i);
            Set<PairInt> set2 = new HashSet<PairInt>(set);
            c.add(set2);
        }

        return c;
    }

    public void mergeSmallSegments(ImageExt img,
        int[] labels, int sizeLimit, ColorSpace clrSpace) {

        CIEChromaticity cieC = null;
        if (ColorSpace.CIELAB.equals(clrSpace)) {
            cieC = new CIEChromaticity();
        }

        //key = label, value = label indexes
        TIntObjectMap<TIntSet> labelToIndexMap =
            LabelToColorHelper.createLabelIndexMap(labels);

        // key = label, value = adjacent label indexes
        TIntObjectMap<TIntSet> adjacencyMap =
            LabelToColorHelper.createAdjacencyLabelMap(
                img, labels, false);

        //key = label, value = average color in clrSpace
        TIntObjectMap<Colors> labelColorMap =
            new TIntObjectHashMap<Colors>();

        TIntObjectMap<Colors> labelRGBMap = new
            TIntObjectHashMap<Colors>();

        TIntSet merged = new TIntHashSet();

        //labels
        Stack<Integer> stack = new Stack<Integer>();
        TIntSet labelSet = labelToIndexMap.keySet();
        TIntIterator iter = labelSet.iterator();
        while (iter.hasNext()) {
            stack.add(Integer.valueOf(iter.next()));
        }

        int nSmall = 0;
        while (!stack.isEmpty()) {
            int label1 = stack.pop().intValue();
            if (merged.contains(label1)) {
                continue;
            }

            TIntSet set = labelToIndexMap.get(label1);
            if (set.size() < sizeLimit) {
                nSmall++;
                // merge w/ closest adjacent above sizeLimit
                Colors clrs = labelColorMap.get(label1);
                Colors rgbClrs = labelRGBMap.get(label1);
                int minLabel = -1;
                double minDiff = Double.MAX_VALUE;
                Colors minDiffClrs = null;
                Colors minDiffRGBClrs = null;

                TIntSet adj = adjacencyMap.get(label1);
                if (adj != null) {
                    TIntIterator iter2 = adj.iterator();
                    while (iter2.hasNext()) {
                        int label3 = iter2.next();
                        if (merged.contains(label3) ||
                            labelToIndexMap.get(label3).size() < sizeLimit) {
                            continue;
                        }
                        if (rgbClrs == null) {
                            rgbClrs = calculateSetColor(set, img,
                                ColorSpace.RGB);
                            clrs = calculateSetColor(rgbClrs, clrSpace);
                        }
                        Colors rgbClrs3 = labelRGBMap.get(label3);
                        Colors clrs3 = labelColorMap.get(label3);
                        if (rgbClrs3 == null) {
                            TIntSet set3 = labelToIndexMap.get(label3);
                            if (set3.size() < sizeLimit) {
                                continue;
                            }
                            rgbClrs3 = calculateSetColor(set3, img,
                                ColorSpace.RGB);
                            labelRGBMap.put(label3, rgbClrs3);
                            clrs3 = calculateSetColor(rgbClrs3, clrSpace);
                            labelColorMap.put(label3, clrs3);
                        }

                        double diff = colorDiff(clrs, clrs3, clrSpace);
                        if (diff < 0) {
                            diff *= -1;
                        }
                        if (diff < minDiff) {
                            minDiff = diff;
                            minLabel = label3;
                            minDiffClrs = clrs3;
                            minDiffRGBClrs = rgbClrs3;
                        }
                    }
                }

                if (minLabel > -1) {

                    // merge this label's data into minLabel's data

                    //update adjacency maps
                    TIntSet adjLabels1 = adjacencyMap.get(label1);
                    TIntIterator iter1 = adjLabels1.iterator();
                    while (iter1.hasNext()) {
                        // change to adjacent to minLabel instead of label1
                        int label3 = iter1.next();
                        if (label3 == minLabel || label3 == label1) {
                            continue;
                        }
                        TIntSet adjLabels3 = adjacencyMap.get(label3);
                        adjLabels3.remove(label1);
                        adjLabels3.add(minLabel);
                    }
                    adjacencyMap.get(minLabel).addAll(adjLabels1);
                    adjacencyMap.get(minLabel).remove(label1);

                    // reassign label1 indexes to minLabel
                    TIntSet set1 = labelToIndexMap.get(label1);
                    float n1 = set.size();
                    float n2 = set1.size();
                    float nTot = n1 + n2;
                    iter1 = set1.iterator();
                    while (iter1.hasNext()) {
                        int lIdx = iter1.next();
                        labels[lIdx] = minLabel;
                    }

                    // update color maps
                    int r = Math.round((rgbClrs.getColors()[0]*n1 +
                        minDiffRGBClrs.getColors()[0]*n2)/nTot);
                    int g = Math.round((rgbClrs.getColors()[1]*n1 +
                        minDiffRGBClrs.getColors()[1]*n2)/nTot);
                    int b = Math.round((rgbClrs.getColors()[2]*n1 +
                        minDiffRGBClrs.getColors()[2]*n2)/nTot);
                    minDiffRGBClrs.getColors()[0] = r;
                    minDiffRGBClrs.getColors()[1] = g;
                    minDiffRGBClrs.getColors()[2] = b;
                    if (ColorSpace.RGB.equals(clrSpace)) {
                        minDiffClrs.getColors()[0] = r;
                        minDiffClrs.getColors()[1] = g;
                        minDiffClrs.getColors()[2] = b;
                    } else if (ColorSpace.HSV.equals(clrSpace)) {
                        Color.RGBtoHSB(r, g, b, minDiffClrs.getColors());
                    } else if (ColorSpace.CIELAB.equals(clrSpace)) {
                        float[] lab = cieC.rgbToCIELAB(r, g, b);
                        System.arraycopy(lab, 0, minDiffClrs.getColors(),
                            0, lab.length);
                    }

                    TIntSet set0 = labelToIndexMap.get(minLabel);
                    set0.addAll(set);
                    set.clear();
                    set0.remove(label1);
                    set0.remove(minLabel);

                    merged.add(label1);
                    stack.add(minLabel);
                }
            }
        }
        log.fine("number of small merges=" + merged.size()
           + " nSmall=" + nSmall + " nSets=" + labelToIndexMap.size());
    }

    private Colors calculateSetColor(Colors rgbClrs,
        ColorSpace clrSpace) {

        CIEChromaticity cieC = null;
        if (ColorSpace.CIELAB.equals(clrSpace)) {
            cieC = new CIEChromaticity();
        }

        float r = rgbClrs.getColors()[0];
        float g = rgbClrs.getColors()[1];
        float b = rgbClrs.getColors()[2];

        if (ColorSpace.RGB.equals(clrSpace)) {
            return new Colors(new float[]{r, g, b});
        } else if (ColorSpace.HSV.equals(clrSpace)) {
            float[] hsb = new float[3];
            Color.RGBtoHSB(Math.round(r), Math.round(g),
                Math.round(b), hsb);
            return new Colors(hsb);
        } else if (ColorSpace.CIELAB.equals(clrSpace)) {
            float[] lab = cieC.rgbToCIELAB(
                Math.round(r), Math.round(g), Math.round(b));
            return new Colors(lab);
        }

        throw new IllegalArgumentException("ColorSpace " +
            clrSpace + " calc is not implemented.");
    }

    private double colorDiff(Colors clrs1, Colors clrs2,
        ColorSpace clrSpace) {

        if (ColorSpace.CIELAB.equals(clrSpace)) {
            CIEChromaticity cieC = new CIEChromaticity();
            double d = cieC.calcDeltaECIE2000(clrs1.getColors(),
                clrs2.getColors());
            return d;
        } else if (ColorSpace.HSV.equals(clrSpace) ||
            ColorSpace.RGB.equals(clrSpace)) {

            double diff = 0;
            for (int i = 0; i < clrs1.getColors().length; ++i) {
                float d = clrs1.getColors()[i] -
                    clrs2.getColors()[i];
                diff = (d * d);
            }

            return Math.sqrt(diff);
        }

        throw new IllegalArgumentException("ColorSpace " +
            clrSpace + " calc is not implemented.");
    }

    private Colors calculateSetColor(
        TIntSet indexesSet, ImageExt img,
        ColorSpace clrSpace) {

        CIEChromaticity cieC = null;
        if (ColorSpace.CIELAB.equals(clrSpace)) {
            cieC = new CIEChromaticity();
        }

        GroupPixelRGB0 gpb = new GroupPixelRGB0();
        gpb.calculateColors(indexesSet, img, 0, 0);

        int r = Math.round(gpb.getAvgRed());
        int g = Math.round(gpb.getAvgGreen());
        int b = Math.round(gpb.getAvgBlue());
        if (ColorSpace.RGB.equals(clrSpace)) {
            return new Colors(new float[]{r, g, b});
        } else if (ColorSpace.HSV.equals(clrSpace)) {
            float[] hsb = new float[3];
            Color.RGBtoHSB(r, g, b, hsb);
            return new Colors(hsb);
        } else if (ColorSpace.CIELAB.equals(clrSpace)) {
            float[] lab = cieC.rgbToCIELAB(r, g, b);
            return new Colors(lab);
        }

        throw new IllegalArgumentException("ColorSpace " +
            clrSpace + " calc is not implemented.");
    }

    public float[][] createSobelGradient(GreyscaleImage img) {

        int nPix = img.getNPixels();
        int width = img.getWidth();
        int height = img.getHeight();

        float[][] gradient = new float[width][];
        for (int i = 0; i < width; ++i) {
            gradient[i] = new float[height];
        }

        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();

        int h = (kernel.length - 1) >> 1;

        for (int i = 0; i < nPix; ++i) {
            final int x1 = img.getCol(i);
            final int y1 = img.getRow(i);

            float xSum = 0;
            float ySum = 0;

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
                xSum += gg * img.getValue(x2, y2);

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
                ySum += gg * img.getValue(x2, y2);
            }

            double c = Math.sqrt(xSum * xSum + ySum * ySum);

            gradient[x1][y1] = (float) c;
        }

        return gradient;
    }

    /**
     *
     * @param img
     * @param gradientMethod
     * 0=CannyEdgeFilterAdaptiveDeltaE2000,
     * 1=CannyEdgeFilterAdaptive,
     * 2=PhaseCongruencyDetector
     * @param ts timestamp used in debugging image name
     * @return
     */
    public EdgeFilterProducts createGradient(Image img,
        int gradientMethod, long ts) {

        EdgeFilterProducts products = null;

        if (gradientMethod == 0) {

            ImageExt imgCp = img.copyToImageExt();
            
            CannyEdgeFilterAdaptiveDeltaE2000 canny =
                new CannyEdgeFilterAdaptiveDeltaE2000();
            canny.setOtsuScaleFactor(0.3f);
            canny.setToUseSingleThresholdIn2LayerFilter();
            canny.applyFilter(imgCp);

            products = canny.getFilterProducts();

        } else if (gradientMethod == 1) {
            
            CannyEdgeFilterAdaptive canny2 = new CannyEdgeFilterAdaptive();
            canny2.overrideToNotUseLineThinner();
            //canny2.setOtsuScaleFactor(0.3f);
            canny2.setToUseSingleThresholdIn2LayerFilter();
            canny2.applyFilter(img.copyToGreyscale2());

            products = canny2.getFilterProducts();

        } else if (gradientMethod == 2) {
            
            products = createPhaseCongruencyGradient(
                img.copyBlueToGreyscale());

        }

        return products;
    }

    private int[] calcStdDev(List<GroupAverageColors> listOfColors,
        List<Integer> indexes) {

        float avgR = 0;
        float avgG = 0;
        float avgB = 0;
        for (Integer index : indexes) {
            GroupAverageColors clrs = listOfColors.get(index.intValue());
            avgR += clrs.getR();
            avgG += clrs.getG();
            avgB += clrs.getB();
        }
        float length = (float)indexes.size();
        avgR /= length;
        avgG /= length;
        avgB /= length;

        float stdvR = 0;
        float stdvG = 0;
        float stdvB = 0;
        for (Integer index : indexes) {
            GroupAverageColors clrs = listOfColors.get(index.intValue());
            float diffR = clrs.getR() - avgR;
            float diffG = clrs.getG() - avgG;
            float diffB = clrs.getB() - avgB;
            stdvR += (diffR * diffR);
            stdvG += (diffG * diffG);
            stdvB += (diffB * diffB);
        }
        int stdDevR = (int)Math.round(Math.sqrt(stdvR/(length - 1.0f)));
        int stdDevG = (int)Math.round(Math.sqrt(stdvG/(length - 1.0f)));
        int stdDevB = (int)Math.round(Math.sqrt(stdvB/(length - 1.0f)));

        return new int[]{stdDevR, stdDevG, stdDevB};
    }

    public EdgeFilterProducts packageToEdgeProduct(
        PhaseCongruencyDetector.PhaseCongruencyProducts pr) {
        
        EdgeFilterProducts eProduct = new EdgeFilterProducts();

        int nCols = pr.getThinned()[0].length;//img.getWidth();
        int nRows = pr.getThinned().length;//img.getHeight();

        GreyscaleImage pcImg = new GreyscaleImage(nCols, nRows);
        double[][] pc = pr.getPhaseCongruency();
        for (int i = 0; i < pr.getThinned().length; ++i) {
            for (int j = 0; j < pr.getThinned()[i].length; ++j) {
                if (pr.getThinned()[i][j] > 0) {
                    int v = (int)Math.round(255. * pc[i][j]);
                    pcImg.setValue(j, i, v);
                }
            }
        }

        eProduct.setGradientXY(pcImg);

        GreyscaleImage paImg = new GreyscaleImage(nCols, nRows,
            GreyscaleImage.Type.Bits32FullRangeInt);
        // range -pi to pi
        double[][] pa = pr.getPhaseAngle();
        for (int i = 0; i < pa.length; ++i) {
            for (int j = 0; j < pa[i].length; ++j) {
                double v = pa[i][j];
                int d = (int)Math.round(v * 180./Math.PI);
                paImg.setValue(j, i, d);
            }
        }
        eProduct.setPhaseAngle(paImg);

        GreyscaleImage orImg = new GreyscaleImage(nCols, nRows);
        double[][] or = pr.getOrientation();
        // orientation is already in range 0 to 180
        for (int i = 0; i < or.length; ++i) {
            for (int j = 0; j < or[i].length; ++j) {
                double v = or[i][j];
                orImg.setValue(j, i, (int)Math.round(v));
            }
        }
        eProduct.setTheta(orImg);
        
        return eProduct;
    }
    
    public EdgeFilterProducts createPhaseCongruencyGradient(
        GreyscaleImage img) {
        
        PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
        phaseDetector.setK(2);
        
        PhaseCongruencyDetector.PhaseCongruencyProducts pr =
            phaseDetector.phaseCongMono(img);

        EdgeFilterProducts eProduct = packageToEdgeProduct(pr);
        
        return eProduct;
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
     * @param w width of image
     * @param h height of image
     * @return
     */
    public static boolean doesDisconnect(final TIntSet input,
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
            ((neighborCoords[0][2].getY() + row) * w) +
            neighborCoords[0][2].getX() + col);
        boolean t7 = input.contains(
            ((neighborCoords[1][2].getY() + row) * w) +
            neighborCoords[1][2].getX() + col);
        boolean t8 = input.contains(
            ((neighborCoords[2][2].getY() + row) * w) +
            neighborCoords[2][2].getX() + col);
        boolean t11 = input.contains(
            ((neighborCoords[0][1].getY() + row) * w) +
            neighborCoords[0][1].getX() + col);
        boolean t12 = input.contains(
            ((neighborCoords[2][1].getY() + row) * w) +
            neighborCoords[2][1].getX() + col);
        boolean t15 = input.contains(
            ((neighborCoords[0][0].getY() + row) * w) +
            neighborCoords[0][0].getX() + col);
        boolean t16 = input.contains(
            ((neighborCoords[1][0].getY() + row) * w) +
            neighborCoords[1][0].getX() + col);
        boolean t17 = input.contains(
            ((neighborCoords[2][0].getY() + row) * w) +
            neighborCoords[2][0].getX() + col);

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

     /**
     *
     * @param img
     * @param edgeValue
     * @param nonEdgeValue
     * @param limit edges less than or equal to this size will be removed
     */
    public void removeEdgesSmallerThanLimit(GreyscaleImage img, int edgeValue,
        int nonEdgeValue, int limit) {

        ConnectedValuesFinder finder = new ConnectedValuesFinder(img);
        finder.setMinimumNumberInCluster(1);
        finder.setToUse8Neighbors();
        finder.findGroups(edgeValue);

        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            TIntSet edge = finder.getXY(i);
            if (edge.size() > limit) {
                continue;
            }
            TIntIterator iter = edge.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                int y = pixIdx/img.getWidth();
                int x = pixIdx - (y * img.getWidth());
            
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

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        final int w = img.getWidth();
        final int h = img.getHeight();

        // half width of neighbor region
        final int hN = 2;//3;

        PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
        
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseDetector.phaseCongMono(gsImg);

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

    private int[][] copy(int[][] a) {

        int[][] b = new int[a.length][];
        for (int i = 0; i < b.length; ++i) {
            b[i] = Arrays.copyOf(a[i], a[i].length);
        }

        return b;
    }

    private boolean assertShortEdgesAreEmpty(List<Integer> indexes,
        List<Set<PairInt>> clusterSets) {

        for (Integer index : indexes) {
            Set<PairInt> set = clusterSets.get(index.intValue());
            assert(set.isEmpty());
        }

        return true;
    }

    /**
     * apply a super-pixels algorithm followed by normalized cuts to obtain
     * labels that are still slightly over-segmented, but should be helpful
     * for object identification.
     * Note that images with average width and length less than 100 and larger
     * than 500 have not been tested yet for the automatic setting of number of
     * super pixels (large images, espec. might need internal handling in
     * block sizes near 512 to keep the number of super pixels less than 1000
     * to reduce the eigenvector calculations).
     * (Note the algorithms implemented are those of
       "SLIC Superpixels Compared to State-of-the-Art Superpixel Methods"
       by Achanta, Appu Shaji,Smith,  Lucchi, Fua, and Suestrun
       and "Normalized Cuts and Image Segmentation"
       by Jianbo Shi and Jitendra Malik)
     * @param img
     */
    public void applySuperPixelsAndNormalizedCuts(ImageExt img) {

        int[] labels = calcSuperPixelsAndNormalizedCutsLabels(img);

        LabelToColorHelper.applyLabels(img, labels);
    }

    /**
     * apply a super-pixels algorithm followed by normalized cuts to obtain
     * labels that are still slightly over-segmented, but should be helpful
     * for object identification.
     * Note that images with average width and length less than 100 and larger
     * than 500 have not been tested yet for the automatic setting of number of
     * super pixels (large images, espec. might need internal handling in
     * block sizes near 512 to keep the number of super pixels less than 1000
     * to reduce the eigenvector calculations).
     * (Note the algorithms implemented are those of
       "SLIC Superpixels Compared to State-of-the-Art Superpixel Methods"
       by Achanta, Appu Shaji,Smith,  Lucchi, Fua, and Suestrun
       and "Normalized Cuts and Image Segmentation"
       by Jianbo Shi and Jitendra Malik)
     * @param img
     */
    public void applySuperPixelsAndNormalizedCuts(ImageExt img, int nIter) {

        int[] labels = calcSuperPixelsAndNormalizedCutsLabels(img, nIter);

        LabelToColorHelper.applyLabels(img, labels);
    }

    /**
     * Apply a super-pixels algorithm followed by normalized cuts to obtain
     * labels that are still slightly over-segmented, but should be helpful
     * for object identification.
     * Note that images with average width and length less than 100 and larger
     * than 500 have not been tested yet for the automatic setting of number of
     * super pixels (large images, espec. might need internal handling in
     * block sizes near 512 to keep the number of super pixels less than 1000
     * to reduce the eigenvector calculations).
       (Note the algorithms implemented are those of
       "SLIC Superpixels Compared to State-of-the-Art Superpixel Methods"
       by Achanta, Appu Shaji,Smith,  Lucchi, Fua, and Suestrun
       and "Normalized Cuts and Image Segmentation"
       by Jianbo Shi and Jitendra Malik)
     * @param img
     */
    public int[] calcSuperPixelsAndNormalizedCutsLabels(ImageExt img) {
        int kCell;

        int avgDimension = (img.getWidth() + img.getHeight()) / 2;

        if (avgDimension < 25) {
            kCell = 2 * avgDimension;
        } else if (avgDimension < 100) {
            //kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 10);
            kCell = 200;
        } else if (avgDimension < 200) {
            //kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 10);
            kCell = 200;
        } else if (avgDimension < 301) {
            //kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 100);
            kCell = 1050;
        } else if (avgDimension < 400) {
            //kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 350);
            kCell = 750;
        } else if (avgDimension < 500) {
            kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 1000);
            kCell *= 2;  // creates a more segmented defined labelling
        } else {
            // this section has not been tested well yet
            kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 1000);
            if (kCell > 2000) {
                kCell = 2000;
            }
        }

        System.out.println("kCell=" + kCell + " avgDim=" + avgDimension);

        kCell = 200;
        int clNorm = 1;

        SLICSuperPixels slic
            = new SLICSuperPixels(img, kCell, clNorm);
        slic.calculate();
        int[] labels = slic.getLabels();

        //ImageExt img2 = img.copyToImageExt();
        //ImageIOHelper.addAlternatingColorLabelsToRegion(img2, labels);
        //MiscDebug.writeImage(img2, "_slic_" + trainingData[i].imgFileName);

        NormalizedCuts normCuts = new NormalizedCuts();
        int[] labels2 = normCuts.normalizedCut(img, labels);

        return labels2;
    }

    /**
     * Apply a super-pixels algorithm followed by normalized cuts to obtain
     * labels that are still slightly over-segmented, but should be helpful
     * for object identification.
     * Note that images with average width and length less than 100 and larger
     * than 500 have not been tested yet for the automatic setting of number of
     * super pixels (large images, espec. might need internal handling in
     * block sizes near 512 to keep the number of super pixels less than 1000
     * to reduce the eigenvector calculations).
       (Note the algorithms implemented are those of
       "SLIC Superpixels Compared to State-of-the-Art Superpixel Methods"
       by Achanta, Appu Shaji,Smith,  Lucchi, Fua, and Suestrun
       and "Normalized Cuts and Image Segmentation"
       by Jianbo Shi and Jitendra Malik)
     * @param img
     */
    public int[] calcSuperPixelsAndNormalizedCutsLabels(ImageExt img, int nIter) {

        int kCell;

        int avgDimension = (img.getWidth() + img.getHeight()) / 2;

        if (avgDimension < 25) {
            kCell = 2 * avgDimension;
        } else if (avgDimension < 100) {
            kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 100);
        } else if (avgDimension < 500) {
            kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 1000);
            //kCell *= 2;  // creates a more segmented defined labelling
        } else {
            // this section has not been tested well yet
            kCell = Math.round((float) (img.getWidth() * img.getHeight()) / 1000);
            if (kCell > 2000) {
                kCell = 2000;
            }
        }

        System.out.println("kCell=" + kCell);

        SLICSuperPixels slic = new SLICSuperPixels(img, kCell);
        slic.calculate();
        int[] labels = slic.getLabels();

        //ImageExt img2 = img.copyToImageExt();
        //ImageIOHelper.addAlternatingColorLabelsToRegion(img2, labels);
        //MiscDebug.writeImage(img2, "_slic_" + trainingData[i].imgFileName);

        for (int i = 0; i < nIter; ++i) {

            NormalizedCuts normCuts = new NormalizedCuts();
            labels = normCuts.normalizedCut(img, labels);

        }

        return labels;
    }

    public static class DecimatedData {
        // decimated comparison size is the closest to 128
        // for that have labels, decimated image, perimeter.
        // for decimated sizes 256 and 512, have
        //    decimated images and labels also.
        // that allows for a range of scale up to 4 in
        // feature comparisons
        public TIntObjectMap<TIntSet> fullLabels;

        // note that a small image may have nulls for
        // dimensions larger than it's image.
        // the 128, 256, and 512 decimated images
        public ImageExt[] dImages = new ImageExt[3];
        public int[] dBinFactors = new int[3];

        // first list indexes are for 128, 256, or 512
        // then map indexes are the labels of the segments
        //  and the values of the maps are the characteristics
        //  of those segments
        public List<TIntObjectMap<TIntSet>> dLabeledIndexes
            = new ArrayList<TIntObjectMap<TIntSet>>();
        public List<TIntObjectMap<PairInt>> dLabelCentroids
            = new ArrayList<TIntObjectMap<PairInt>>();
    }

    /**
     * using default spacing of 6 pixels for a feature.
     * @param dd
     * @param decimatedImageIndex
     * @param imgWidth
     * @param imgHeight
     * @return
     */
    public TIntObjectMap<List<PairInt>> calculateKeyPoints(
        DecimatedData dd, int decimatedImageIndex,
        int imgWidth, int imgHeight) {

        TIntObjectMap<List<PairInt>> output =
            new TIntObjectHashMap<List<PairInt>>();

        /*
        Filling out keypoints across object shape to cover all
        points (including exterior points is allowed).

        Using a BFS search pattern within dSpace limits of
        current center keyPoint.
        */

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        // -dSpace to +dSpace is the feature range along a col or row
        int dSpace = 6;

        TIntObjectMap<TIntSet> labeledIndexes =
            dd.dLabeledIndexes.get(decimatedImageIndex);

        TIntObjectIterator<TIntSet> iter =
            labeledIndexes.iterator();

        for (int sIdx = 0; sIdx < labeledIndexes.size(); ++sIdx) {

            iter.advance();

            int label = iter.key();

            List<PairInt> keyPoints = new ArrayList<PairInt>();
            output.put(label, keyPoints);

            TIntSet indexes = iter.value();

            // assign keyPoints to cover all of indexes:
            // BFS search within indexes, adding neighbors
            //    within dSpace.

            // key = pixIdx, value = index of keyPoints list
            TIntIntMap indexKPMap = new TIntIntHashMap();

            // if need speed, can replace q0 with any order of indexes
            ArrayDeque<Integer> q0 = populateByNumberOfNeighbors(
                indexes, imgWidth, imgHeight);
            ArrayDeque<Integer> q1 = new ArrayDeque<Integer>();

            int nIter = 0;
            Set<Integer> visited = new HashSet<Integer>();
                while (!q0.isEmpty() || !q1.isEmpty()) {
                Integer uIndex;
                if (!q1.isEmpty()) {
                    // draw from q1 first if any because it searches
                    // current keypoint
                    uIndex = q1.poll();
                } else {
                    uIndex = q0.poll();
                }
                if (visited.contains(uIndex)) {
                    continue;
                }
                visited.add(uIndex);
                int uIdx = uIndex.intValue();
                int uY = uIdx / imgWidth;
                int uX = uIdx - (uY * imgWidth);
                int kpIdx;
                // lookup center point of keyPoint or assign a new
                int xc, yc;
                if (indexKPMap.containsKey(uIdx)) {
                    kpIdx = indexKPMap.get(uIdx);
                    PairInt kp = keyPoints.get(kpIdx);
                    xc = kp.getX();
                    yc = kp.getY();
                } else {
                    PairInt kp = new PairInt(uX, uY);
                    kpIdx = keyPoints.size();
                    indexKPMap.put(uIdx, kpIdx);
                    keyPoints.add(kp);
                    xc = uX;
                    yc = uY;
                }
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = uX + dxs[k];
                    if (x2 < 0 || x2 > (imgWidth - 1)) {
                        continue;
                    }
                    int diffX = Math.abs(x2 - xc);
                    if (diffX > dSpace) {
                        continue;
                    }
                    int y2 = uY + dys[k];
                    if (y2 < 0 || y2 > (imgHeight - 1)) {
                        continue;
                    }
                    int diffY = Math.abs(y2 - yc);
                    if (diffY > dSpace) {
                        continue;
                    }
                    int vIdx = (y2 * imgWidth) + x2;
                    if (!indexKPMap.containsKey(vIdx)) {
                        indexKPMap.put(vIdx, kpIdx);
                        q1.add(Integer.valueOf(vIdx));
                    }
                }
                nIter++;
            }
        }

        return output;
    }

    public void replaceSinglePixelLabelsCIELAB(int[] labels,
        ImageExt img) {

        // ----- replace single pixels w/ adjacent nearest in color -----
        int[] dx2 = Misc.dx4;
        int[] dy2 = Misc.dy4;
        // single pixels should join closest,,,
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                int v = labels[pixIdx];
                boolean oneIsSame = false;
                for (int z = 0; z < dx2.length; ++z) {
                    int x2 = i + dx2[z];
                    int y2 = j + dy2[z];
                    if (x2 < 0 || y2 < 0 || (x2 > (img.getWidth() - 1))
                        || (y2 > (img.getHeight() - 1))
                        ) {
                        continue;
                    }
                    int pixIdx2 = img.getInternalIndex(x2, y2);
                    int v2 = labels[pixIdx2];
                    if (v2 == v) {
                        oneIsSame = true;
                        break;
                    }
                }
                if (!oneIsSame) {
                    float[] lab = img.getCIELAB(pixIdx);
                    double minDiff = Double.MAX_VALUE;
                    int minIdx = -1;
                    for (int z = 0; z < dx2.length; ++z) {
                        int x2 = i + dx2[z];
                        int y2 = j + dy2[z];
                        if (x2 < 0 || y2 < 0 || (x2 > (img.getWidth() - 1))
                            || (y2 > (img.getHeight() - 1))) {
                            continue;
                        }
                        int pixIdx2 = img.getInternalIndex(x2, y2);
                        float[] lab2 = img.getCIELAB(pixIdx2);
                        
                        double dClrSq = 0;
                        for (int i2 = 0; i2 < 3; ++i2) {
                            float diff = lab[i2] - lab2[i2];
                            dClrSq += (diff * diff);
                        }
                        
                        if (dClrSq < minDiff) {
                            minDiff = dClrSq;
                            minIdx = pixIdx2;
                        }
                    }
                    labels[pixIdx] = labels[minIdx];
                }
            }
        }
    }

    public boolean filterByLUVDeltaE(ImageExt templateImage, 
        Set<PairInt> templateSet, ImageExt img, 
        List<Set<PairInt>> pointSets, float luvDeltaELimit) {
        
        GroupPixelCIELUV luvTemplate = new GroupPixelCIELUV(
            templateSet, templateImage);
        luvTemplate.calculateColors(templateSet, templateImage, 0, 0);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        TIntList rm = new TIntArrayList();
                
        for (int i = 0; i < pointSets.size(); ++i) {
            
            Set<PairInt> set = pointSets.get(i);
            
            GroupPixelCIELUV luv = new GroupPixelCIELUV(set, img);
            luv.calculateColors(set, img, 0, 0);
            
            double deltaE = cieC.calcDeltaECIE2000(
                luvTemplate.getAvgL(), luvTemplate.getAvgU(),
                luvTemplate.getAvgV(),
                luv.getAvgL(), luv.getAvgU(), luv.getAvgV());
             
            // TODO: this may need to be a higher limit
            if (Math.abs(deltaE) > luvDeltaELimit) {
                rm.add(i);
            }
        }
        
        for (int i = (rm.size() - 1); i > -1; --i) {
            int rmIdx = rm.get(i);
            pointSets.remove(rmIdx);
        }
                
        return !rm.isEmpty();
    }
    
    public boolean filterByCIETheta(ImageExt templateImage, 
        Set<PairInt> templateSet, ImageExt img, 
        List<Set<PairInt>> pointSets) {
     
        ImageProcessor imageProcessor = new ImageProcessor();
        GreyscaleImage theta1 = 
            imageProcessor.createCIELUVTheta(templateImage, 255);
        
        GreyscaleImage theta2 = 
            imageProcessor.createCIELUVTheta(img, 255);
        
        //TODO:  may need special handling for colors near center,
        //   that is grey-ish colors, because the polar direction
        //   error in that small radius is large.
        //  looking at bigger picture changes to use deltaE2000
        //  again...
        
        ColorHistogram ch = new ColorHistogram();
        
        int[] tHist = ch.histogram1D(theta1, templateSet, 255);
        
        TIntList rm = new TIntArrayList();
        
        float limit = 0.5f;//0.2f;
        
        for (int i = 0; i < pointSets.size(); ++i) {
            Set<PairInt> set = pointSets.get(i);
            int[] hist = ch.histogram1D(theta2, set, 255);
            float intersection = ch.intersection(tHist, hist);
            if (intersection < limit) {
                rm.add(i);
            }
        }
        
        for (int i = (rm.size() - 1); i > -1; --i) {
            int rmIdx = rm.get(i);
            pointSets.remove(rmIdx);
        }
                
        return !rm.isEmpty();
    }
    
    public boolean filterByCIECH(ImageExt templateImage, 
        Set<PairInt> templateSet, ImageExt img, 
        List<Set<PairInt>> pointSets, float lowerLimit) {
     
        ColorHistogram ch = new ColorHistogram();
        
        int[] tHist = ch.histogramCIECH64(templateImage, 
            templateSet);
        
        TIntList rm = new TIntArrayList();
                
        for (int i = 0; i < pointSets.size(); ++i) {
            Set<PairInt> set = pointSets.get(i);
            int[] hist = ch.histogramCIECH64(img, set);
            float intersection = ch.intersection(tHist, hist);
            if (intersection < lowerLimit) {
                rm.add(i);
            }
        }
        
        for (int i = (rm.size() - 1); i > -1; --i) {
            int rmIdx = rm.get(i);
            pointSets.remove(rmIdx);
        }
                
        return !rm.isEmpty();
    }
    
}
