package algorithms.imageProcessing.segmentation;

import algorithms.disjointSets.UnionFind;
import algorithms.imageProcessing.*;
import algorithms.misc.MiscMath0;

import java.util.*;

public class MergeLabels {

    protected static int[][] offsets = new int[][]{
            {-1, 0}, {0, 1}, {1, 0}, {-1, 0},
            {-1,1}, {-1,-1}, {1,-1}, {1, 1}
    };

    public static enum METHOD {
        MEAN, MODE, MIN_GRADIENT
    }

    /**
     * merge adjacent labeled regions if the deltaE2000 color difference between two adjacent regions is
     * less than the given thresh.
     *
     * @param img
     * @param labels input labels, they are updated in place with merged labels
     * @param thresh threshold for the deltaE2000 color difference.  The JND (just noticeable difference) threshold is
     *               about 2.3 and the maximum value is 20?  (TODO: calculate this).
     *               Caveat is that the if the method is MEAN of color histograms of regions are taken as the representative
     *               color for the region, and so a threshold between 0 and 1.0 performs better.
     *               For MODE,
     * @param method mean or mode
     * @return number of unique labels
     */
    public static int mergeUsingDeltaE2000(ImageExt img, int[] labels, double thresh, METHOD method) {

        if (method == null) {
            throw new IllegalArgumentException("method cannot be null");
        }

        if (method == METHOD.MIN_GRADIENT) {
            return mergeWithMinGrad(img, labels, thresh);
        }

        // merging adjacent labels if their deltaE2000 are < thresh
        //
        // will use the mode of each labeled region
        //
        // datastructures:
        //  hashmap w/ key = label, value = color histogram
        //  hashmap w/ key = label, value = set of adjacent labels
        //
        // will use BFS traversal of adjacency map and when deltaE2000 of parent and child are < thresh
        // will merge the two using union find.
        //
        // the final components are the parents in union find.
        // for each original label, there will be a new label and total number of unique labels is <= original number.

        CIEChromaticity cieC = new CIEChromaticity();

        int r1, c1, r2, c2, i2, label2;
        int w = img.getWidth();
        int h = img.getHeight();
        Map<Integer, Set<Integer>> adjLabelMap = new HashMap<>();
        Map<Integer, int[][]> labelHistMap = new HashMap<>();

        int[][] tmp = new int[3][51];

        for (int i1 = 0; i1 < labels.length; ++i1) {
            int label1 = labels[i1];
            adjLabelMap.putIfAbsent(label1, new HashSet<>());

            r1 = img.getRow(i1);
            c1 = img.getCol(i1);
            for (int[] offset : offsets) {
                r2 = r1 + offset[0];
                c2 = c1 + offset[1];
                if (r2 < 0 || c2 < 0 || r2 == h || c2 == w) {
                    continue;
                }
                i2 = img.getInternalIndex(c2, r2);
                label2 = labels[i2];

                if (label1 == label2) {
                    continue;
                }

                adjLabelMap.putIfAbsent(label2, new HashSet<>());

                adjLabelMap.get(label1).add(label2);
                adjLabelMap.get(label2).add(label1);
            }

            // add to histogram entry.  use bin sizes of 5. 255/5 = 51.
            labelHistMap.putIfAbsent(label1, new int[3][51]);

            int[][] lHist = labelHistMap.get(label1);
            lHist[0][img.getR(i1)/51]++;
            lHist[1][img.getG(i1)/51]++;
            lHist[2][img.getB(i1)/51]++;
        }

        int nLabels = adjLabelMap.size();
        UnionFind uf = new UnionFind(nLabels);

        // BFS traverse adj map
        Set<Integer> visited = new HashSet<>();
        Deque<Integer> q = new ArrayDeque<>();
        q.add(labels[0]);
        int u;
        int[][] uHist;
        int[][] vHist;
        while (!q.isEmpty()) {
            u = q.poll();
            visited.add(u);
            uHist = labelHistMap.get(u);
            for (int v : adjLabelMap.get(u)) {
                if (visited.contains(v)) {
                    continue;
                }
                q.add(v);
                if (uf.find(u) == uf.find(v)) {
                    continue;
                }
                vHist = labelHistMap.get(v);
                double diff = calcColorDiff(cieC, uHist, vHist, method);
                if (diff < thresh) {
                    uf.union(u, v);
                    // update uHist and vHist
                    add(uHist, vHist, tmp);
                }
            }
        }

        Map<Integer, Set<Integer>> mergedLabelMap = uf.getComponents();

        int j = 0;
        // reverse the map so given previous label, we have new label, but renumber new label from 0 to size()-1
        int[] rev = new int[nLabels];
        for (Map.Entry<Integer, Set<Integer>> entry : mergedLabelMap.entrySet()) {
            for (int oldLabel : entry.getValue()) {
                rev[oldLabel] = j;
            }
            ++j;
        }

        for (int i = 0; i < labels.length; ++i) {
            labels[i] = rev[labels[i]];
        }

        return mergedLabelMap.size();
    }

    /**
     * the color descriptor of each label is determined from the point which has the min gradient for each label region
     * as was done in the SuperPixels algorithm.
     * @param img
     * @param labels
     * @param thresh
     * @return the number of unique labels after merging
     */
    protected static int mergeWithMinGrad(ImageExt img, int[] labels, double thresh) {

        ImageSegmentation seg = new ImageSegmentation();
        EdgeFilterProducts edgeProducts = seg.createGradient(img, 0, System.currentTimeMillis());
        GreyscaleImage grad = edgeProducts.getGradientXY();

        CIEChromaticity cieC = new CIEChromaticity();

        int r1, c1, r2, c2, i2, label2;
        int w = img.getWidth();
        int h = img.getHeight();
        Map<Integer, Set<Integer>> adjLabelMap = new HashMap<>();

        // descriptor is [L, A, B, minGrad]
        Map<Integer, float[]> labelDescMap = new HashMap<>();

        for (int i1 = 0; i1 < labels.length; ++i1) {
            int label1 = labels[i1];
            adjLabelMap.putIfAbsent(label1, new HashSet<>());

            r1 = img.getRow(i1);
            c1 = img.getCol(i1);
            for (int[] offset : offsets) {
                r2 = r1 + offset[0];
                c2 = c1 + offset[1];
                if (r2 < 0 || c2 < 0 || r2 == h || c2 == w) {
                    continue;
                }
                i2 = img.getInternalIndex(c2, r2);
                label2 = labels[i2];

                if (label1 == label2) {
                    continue;
                }

                adjLabelMap.putIfAbsent(label2, new HashSet<>());

                adjLabelMap.get(label1).add(label2);
                adjLabelMap.get(label2).add(label1);
            }

            labelDescMap.putIfAbsent(label1, new float[]{0.f, 0.f, 0.f, Integer.MAX_VALUE});

            int gradient = grad.getValue(c1, r1);
            if (gradient < labelDescMap.get(label1)[3]) {
                float[] lab = img.getCIELAB(c1, r1);
                float[] tmp = labelDescMap.get(label1);
                System.arraycopy(lab, 0, tmp, 0, lab.length);
                tmp[3] = gradient;
            }
        }

        int nLabels = adjLabelMap.size();
        UnionFind uf = new UnionFind(nLabels);

        // BFS traverse adj map
        Set<Integer> visited = new HashSet<>();
        Deque<Integer> q = new ArrayDeque<>();
        q.add(labels[0]);
        int u;
        float[] uDesc;
        float[] vDesc;
        while (!q.isEmpty()) {
            u = q.poll();
            visited.add(u);
            uDesc = labelDescMap.get(u);
            for (int v : adjLabelMap.get(u)) {
                if (visited.contains(v)) {
                    continue;
                }
                q.add(v);
                if (uf.find(u) == uf.find(v)) {
                    continue;
                }
                vDesc = labelDescMap.get(v);
                double deltaE = cieC.calcDeltaECIE2000( uDesc[0], uDesc[1], uDesc[2], vDesc[0], vDesc[1], vDesc[2]);
                if (deltaE < thresh) {
                    uf.union(u, v);
                    // update uDesc and vDesc
                    if (uDesc[3] <= vDesc[3]) {
                        System.arraycopy(uDesc, 0, vDesc, 0, uDesc.length);
                    } else {
                        System.arraycopy(vDesc, 0, uDesc, 0, vDesc.length);
                    }
                }
            }
        }

        Map<Integer, Set<Integer>> mergedLabelMap = uf.getComponents();

        int j = 0;
        // reverse the map so given previous label, we have new label, but renumber new label from 0 to size()-1
        int[] rev = new int[nLabels];
        for (Map.Entry<Integer, Set<Integer>> entry : mergedLabelMap.entrySet()) {
            for (int oldLabel : entry.getValue()) {
                rev[oldLabel] = j;
            }
            ++j;
        }

        for (int i = 0; i < labels.length; ++i) {
            labels[i] = rev[labels[i]];
        }

        return mergedLabelMap.size();

    }

    private static void add(int[][] uHist, int[][] vHist, int[][] tmp) {
        for (int c = 0; c < 3; ++c) {
            System.arraycopy(uHist[c], 0, tmp[c], 0, uHist[c].length);
        }
        // u += v
        for (int c = 0; c < 3; ++c) {
            for (int i = 0; i < uHist[c].length; ++i) {
                uHist[c][i] += vHist[c][i];
            }
        }
        // v += tmp
        for (int c = 0; c < 3; ++c) {
            for (int i = 0; i < vHist[c].length; ++i) {
                vHist[c][i] += tmp[c][i];
            }
        }
    }

    private static double calcColorDiff(CIEChromaticity cieC, int[][] uHist, int[][] vHist,
                                        METHOD method) {
        // ideally, the histograms would have single modes, but cannot guarantee that, so will use mean instead.

        int[] u, v;
        if (method == METHOD.MEAN) {
            u = meanRGB(uHist);
            v = meanRGB(vHist);
        } else {
            u = modeRGB(uHist);
            v = modeRGB(vHist);
        }

        float[] lab1 = cieC.cieXYZToCIELAB(cieC._rgbToCIEXYZ(u[0], u[1], u[2]));
        float[] lab2 = cieC.cieXYZToCIELAB(cieC._rgbToCIEXYZ(v[0], v[1], v[2]));
        return cieC.calcDeltaECIE2000(lab1, lab2);
    }

    protected static int[] meanRGB(int[][] colorHist) {
        int[] mean = new int[3];
        for (int c = 0; c < colorHist.length; ++c) {
            float total = 0;
            float sum = 0;
            float x = 2.5f;
            for (int i = 0; i < colorHist[c].length; ++i, x += 5.0f) {
                sum += (x * colorHist[c][i]);
                total += colorHist[c][i];
            }
            mean[c] = Math.round(sum/total);
        }
        return mean;
    }

    protected static int[] modeRGB(int[][] colorHist) {
        int[] mode = new int[3];
        for (int i = 0; i < 3; ++i) {
            int idx = MiscMath0.findMax(colorHist[i]);
            mode[i] = Math.round((idx * 5)+2.5f);
        }
        return mode;
    }

}
