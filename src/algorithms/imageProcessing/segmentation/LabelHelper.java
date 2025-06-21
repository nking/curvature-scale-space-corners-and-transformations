package algorithms.imageProcessing.segmentation;

import algorithms.compGeometry.PerimeterFinder3;
import algorithms.disjointSets.UnionFind;
import algorithms.misc.MiscMath0;

import java.util.*;

/**
 * TODO: consider moving all the miscelaneous label helper methods into this class
 */
public class LabelHelper {

    protected static int[][] offsets4 = new int[][]{
            {0, -1}, {-1, 0}, {0, 1}, {1, 0}
    };
    protected static int[][] offsets8 = new int[][]{
            {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1,1}, {1, 0}, {1, -1},
    };

    /**
     * given the labels, further divide each lable by connectedness.  e.g. if
     * 2 regions share the same label but have no adjacent pixels, they will be
     * relabeled to have 2 separate labels.
     * @param labels indexes are pixel indexes following the convention idx = row * width + col
     * @param imgWidth width of image
     * @param imgHeight height of image
     * @param use4Neighors if true, considers a neighbor to be in 4 neighbor neighborhood, else uses 8 neighbor region
     @return the number of unique labels in labels now
     */
    public static int resolveByConnectedness(int[] labels, int imgWidth, int imgHeight, boolean use4Neighors) {

        Map<Integer, Set<Integer>> pointsMap = PerimeterFinder3.getLabeledPointsMap(labels);

        int maxLabel = MiscMath0.findMax(labels);

        int[][] offsets;
        if (use4Neighors) {
            offsets = offsets4;
        } else {
            offsets = offsets8;
        }

        for (Map.Entry<Integer, Set<Integer>> entry : pointsMap.entrySet()) {
            final int label = entry.getKey();
            Set<Integer> points = entry.getValue();
            UnionFind uf = new UnionFind(labels.length);
            Deque<Integer> q = new ArrayDeque<>(points);
            Set<Integer> visited = new HashSet<>();
            int r1, r2, c1, c2, idx1, idx2;
            while (!q.isEmpty()) {
                idx1 = q.pollFirst();
                if (visited.contains(idx1)) {
                    continue;
                }
                visited.add(idx1);
                r1 = idx1 / imgWidth;
                c1 = idx1 % imgWidth;
                for (int[] offset : offsets) {
                    r2 = offset[0] + r1;
                    c2 = offset[1] + c1;
                    if (r2 < 0 || c2 < 0 || r2 == imgHeight || c2 == imgWidth) {
                        continue;
                    }
                    idx2 = r2 * imgWidth + c2;
                    if (!points.contains(idx2)) {
                        continue;
                    }
                    if (uf.find(idx1) != uf.find(idx2)) {
                        uf.union(idx1, idx2);
                    }
                }
            }
            // extract the connected regions within this label
            Map<Integer, Set<Integer>> labelMap = new HashMap<>();
            for (int idx : points) {
                int p = uf.find(uf.getParent()[idx]);
                labelMap.putIfAbsent(p, new HashSet<>());
                labelMap.get(p).add(idx);
            }
            //relabel them if there is more than 1 region
            if (labelMap.size() > 1) {
                int nI = -1;
                for (Map.Entry<Integer, Set<Integer>> entry2 : labelMap.entrySet()) {
                    ++nI;
                    if (nI == 0) {
                        continue; // keeps default label
                    }
                    ++maxLabel;
                    for (int idx : entry2.getValue()) {
                        labels[idx] = maxLabel;
                    }
                }
            }
        }
        return maxLabel;
    }
}
