package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.imageProcessing.features.PhaseCongruencyDetector;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.set.TIntSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

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

    private Map<PairInt, Set<Integer>> 
    findUnassignedPixelsAndAdjacentIndexes(
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

    private boolean assertShortEdgesAreEmpty(List<Integer> indexes,
        List<Set<PairInt>> clusterSets) {

        for (Integer index : indexes) {
            Set<PairInt> set = clusterSets.get(index.intValue());
            assert(set.isEmpty());
        }

        return true;
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
    
}
