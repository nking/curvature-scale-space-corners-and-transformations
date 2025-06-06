package algorithms.imageProcessing.segmentation;

import algorithms.disjointSets.UnionFind;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.Image;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.*;

public class MSEREdgesToLabelsHelper {

    /**
     * populate outLabels given the list of edges and a color image.
     * @param img
     * @param edgeList
     * @param outLabels
     * @return the number of label components
     */
    public static int createLabels(Image img, List<TIntSet> edgeList, int[] outLabels) {

        int n = img.getNPixels();
        if (outLabels.length != n) {
            throw new IllegalArgumentException("outLabels.length must == img.getNPixels()");
        }

        TIntSet edgePoints = new TIntHashSet();
        for (TIntSet set : edgeList) {
            edgePoints.addAll(set);
        }

        if (edgePoints.size() == n) {
            Arrays.fill(outLabels, 0);
            return 1;
        }

        ArrayDeque<Integer> q = new ArrayDeque<>();
        for (int i = 0; i < n; ++i) {
            if (!edgePoints.contains(i)) {
                q.add(i);
            }
        }

        int w = img.getWidth();
        int h = img.getHeight();

        UnionFind uf = new UnionFind(n);

        // r,c offsets.  q is populated along columns then next row
        int[][] offsets0 = new int[][]{
                //{-1, 0}, {0, 1}, {1, 0}, {-1, 0}
                {0, 1}, {-1,0},{-1,1}
        };

        // join neighbors
        int r, c, r2, c2, i2;
        while (!q.isEmpty()) {
            int i = q.pollFirst();
            r = img.getRow(i);
            c = img.getCol(i);
            for (int[] offset : offsets0) {
                r2 = r + offset[0];
                c2 = c + offset[1];
                if (r2 < 0 || c2 < 0 || r2 == h || c2 == w) {
                    continue;
                }
                i2 = img.getInternalIndex(c2, r2);
                if (!edgePoints.contains(i2)) {
                    if (uf.find(i) != uf.find(i2)) {
                        uf.union(i, i2);
                    }
                }
            }
        }

        Map<Integer, Set<Integer>> components = uf.getComponents();

        // remove edgePoints.  they are sets of size 1
        Set<Integer> rm = new HashSet<>();
        for (Map.Entry<Integer, Set<Integer>> entry : components.entrySet()) {
            if (entry.getValue().size() == 1 && edgePoints.contains(entry.getValue().iterator().next())) {
                rm.add(entry.getKey());
            }
        }
        for (int label : rm) {
            components.remove(label);
        }

        final int nMax = Math.max(w, h);
        // array of {index, nIter}
        ArrayDeque<int[]> q2= new ArrayDeque<>();
        TIntIterator iter = edgePoints.iterator();
        // assign parent of edgepoints to -1
        while (iter.hasNext()) {
            int i = iter.next();
            q2.add(new int[]{i, 0});
            uf.getParent()[i] = -1;
        }

        int[][] offsets = new int[][]{
                {-1, 0}, {0, 1}, {1, 0}, {-1, 0},
                {-1,1}, {-1,-1}, {1,-1}, {1, 1}
        };

        // assign edge points to closest component in terms of color
        CIEChromaticity cieC = new CIEChromaticity();
        float[] f, f2;
        while (!q2.isEmpty()) {
            int[] a = q2.pollFirst();
            int i = a[0];
            r = img.getRow(i);
            c = img.getCol(i);

            f = cieC._rgbToCIEXYZ(img.getR(i), img.getG(i), img.getB(i));
            f = cieC.cieXYZToCIELAB(f);

            double closestDiff = Double.POSITIVE_INFINITY;
            int closestI = -1;

            for (int[] offset : offsets) {
                r2 = r + offset[0];
                c2 = c + offset[1];
                if (r2 < 0 || c2 < 0 || r2 == h || c2 == w) {
                    continue;
                }
                i2 = img.getInternalIndex(c2, r2);
                // if it is an assigned edgePoint, proceed, else skip
                if (edgePoints.contains(i2) && uf.getParent()[i2] == -1) {
                    continue;
                }

                // using delta CIE 2000 as color difference
                f2 = cieC._rgbToCIEXYZ(img.getR(i2), img.getG(i2), img.getB(i2));
                f2 = cieC.cieXYZToCIELAB(f2);

                double diff = cieC.calcDeltaECIE2000(f, f2);

                if (diff < closestDiff) {
                    closestI = i2;
                    closestDiff = diff;
                }
            }

            if (closestI != -1) {
                int label = uf.find(closestI);
                components.get(label).add(i);
                uf.getParent()[i] = label;
            } else if (a[1] < nMax) {
                q2.add(new int[]{i, a[1] + 1});
            }
        }

        Arrays.fill(outLabels, -1);

        // extract groups and assign labels

        int label = 0;
        for (Map.Entry<Integer, Set<Integer>> entry : components.entrySet()) {
            for (int i : entry.getValue()) {
                outLabels[i] = label;
            }
            ++label;
        }

        return components.size();
    }

    public static boolean allAreNonNegative(int[] labels) {
        for (int label : labels) {
            if (label < 0) {
                return false;
            }
        }
        return true;
    }
}
