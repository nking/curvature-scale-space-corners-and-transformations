package algorithms.compGeometry;

import algorithms.disjointSets.UnionFind;
import algorithms.disjointSets.UnionFindExt;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;

import java.util.*;

//TODO: consider adding a minimum size argument to some of the methods or overload to not solve for borders if number
// of points is fewer thn a minimum size.
public class PerimeterFinder3 {

    // using convention {col, row}
    protected static int[][] offsets8 = new int[][]{
            {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1,1}, {1, 0}, {1, -1},
    };
    protected static int[][] offsets16 = new int[][]{
            {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1,1}, {1, 0}, {1, -1},
            {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1,1}, {1, 0}, {1, -1},
    };
    protected static Map<PairInt, Integer> offsetIndexMap = new HashMap<>();
    static {
        offsetIndexMap.put(new PairInt(0, -1), 0);
        offsetIndexMap.put(new PairInt(-1, -1), 1);
        offsetIndexMap.put(new PairInt(-1, 0), 2);
        offsetIndexMap.put(new PairInt(-1,1), 3);
        offsetIndexMap.put(new PairInt(0,1), 4);
        offsetIndexMap.put(new PairInt(1,1), 5);
        offsetIndexMap.put(new PairInt(1,0), 6);
        offsetIndexMap.put(new PairInt(1,-1), 7);
    }

    /**
     * given labels for each pixel image where idx = row * imageWidth + col, extract each region's
     * border pixels in a connected clockwise manner to form a closed curve, missing the last point which equals
     * the first.
     * @param labels labels of each image pixel
     * @param imgWidth width of image
     * @param imgHeight height of image
     * @return each region's bounding pixels in a connected clockwise manner to forma closed curve, missing the last point which equals
     * the first.
     */
    public static Map<Integer, PairIntArray> extractOrderedBorders(int[] labels, int imgWidth, int imgHeight) {

        Map<Integer, Set<Integer>> filledRegions = regionFill(labels, imgWidth, imgHeight);

        Map<Integer, PairIntArray> out = new HashMap<>();

        for (Map.Entry<Integer, Set<Integer>> entry : filledRegions.entrySet()) {

            int label = entry.getKey();

            Set<Integer> points = entry.getValue();

            // using Moore boundary tracing algorithm with Jacob's stopping criteria
            // https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/moore.html
            PairIntArray bounds = mooreTracingWithJacob(points, imgWidth, imgHeight);

            out.put(label, bounds);
        }

        return out;
    }

    /**
     * given labels for each pixel image where idx = row * imageWidth + col,
     * for each region, fill in any interior "holes" where a "hole" is any other labeled region.
     * @param labels labels of each image pixel
     * @param imgWidth width of image
     * @param imgHeight height of image
     * @return a map of labeled regions where each region has been in-filled so that all
     * interior points are present as a mask where the map holds key=label, value = in-filled region points.
     */
    public static Map<Integer, Set<Integer>> regionFill(int[] labels, int imgWidth, int imgHeight) {

        /*
        (1) create hashmap w/ key = label, value = indexes of points having that label
        (2) determine which labels touch any image boundaries by scanning the image boundary pixels
        (3) for each labeled region, a region-fill algorithm fills in all interior
            labeled regions, even for complex cases like a region interior to another region.
            (3.a) determine all labeled regions that touch an image boundary and keep sets of
                  those that do and those that do not.
            (3.b) make an adjacency map for the labels.
            (3.c) for each labeled region:
                 - create a UnionFind initialized with all boundary labels having same parent.
                 - take the non-boundary set of labels and subtract current label from it.
                 - traverse the non-boundary set and union all adjacent pairs excepting current label.
                 - extract components from union find.  those that are not in the boundary touching
                   group are internal to the current label and can be merged for the infill.

         */

        Map<Integer, Set<Integer>> labelPointsMap = getLabeledPointsMap(labels);

        int nLabels = labelPointsMap.size();

        // determine which labels touch any image boundaries by scanning bounds
        Set<Integer> boundaryTouchingLabels = new HashSet<>();
        int r, c, idx;
        for (r = 0; r < imgHeight; ++r) {
            c = 0;
            idx = r * imgWidth + c;
            boundaryTouchingLabels.add(labels[idx]);
            c = imgWidth - 1;
            idx = r * imgWidth + c;
            boundaryTouchingLabels.add(labels[idx]);
        }
        for (c = 0; c < imgWidth; ++c) {
            r = 0;
            idx = r * imgWidth + c;
            boundaryTouchingLabels.add(labels[idx]);
            r = imgHeight - 1;
            idx = r * imgWidth + c;
            boundaryTouchingLabels.add(labels[idx]);
        }

        // make a set of all boundary touching labels and remove a set as needed
        Set<Integer> notBoundaryConnectedLabels = new HashSet<>();
        for (int label : labelPointsMap.keySet()) {
            if (!boundaryTouchingLabels.contains(label)) {
                notBoundaryConnectedLabels.add(label);
            }
        }

        // union find initialized for all boundary connected labels
        UnionFind uf = new UnionFind(nLabels);
        int label0 = -1;
        for (int label : boundaryTouchingLabels) {
            if (label0 == -1) {
                label0 = label;
                continue;
            }
            uf.union(label0, label);
        }
        assert(label0 != -1);

        Map<Integer, Set<Integer>> adjLabelMap = getLabelAdjMap(labels, imgWidth, imgHeight);

        Map<Integer, Set<Integer>> filledRegions = new HashMap<>();

        for (Map.Entry<Integer, Set<Integer>> entry : labelPointsMap.entrySet()) {
            int label = entry.getKey();

            // region filling algorithm:
            // UF init w/ boundary touching labels:
            UnionFind uf2 = new UnionFindExt(uf);
            // non-boundary touching labels minus current label
            Set<Integer> notBoundingLabels = new HashSet<>(notBoundaryConnectedLabels);
            notBoundingLabels.remove(label);
            // traverse notBoundingLabels and use adj map to merge adjacent labels,
            //  excluding current label.
            Set<Integer> visited = new HashSet<>();
            Deque<Integer> q = new ArrayDeque<>();
            q.addAll(notBoundingLabels);
            int u;
            while (!q.isEmpty()) {
                u = q.pollFirst();
                visited.add(u);
                if (!adjLabelMap.containsKey(u)) {
                    continue;
                }
                for (int v : adjLabelMap.get(u)) {
                    if (visited.contains(v) || v == label) {
                        continue;
                    }
                    if (uf2.find(u) != uf2.find(v)) {
                        uf2.union(u, v);
                    }
                }
            }
            // and labels connected to the boundary touching labels in uf2 are not internal
            // to current label region, and all else are
            Map<Integer, Set<Integer>> labelComponents = uf2.getComponents();

            int boundaryTouchingParent = uf2.find(label0);

            Set<Integer> points = new HashSet<>(entry.getValue());
            for (Map.Entry<Integer, Set<Integer>> entry2 : labelComponents.entrySet()) {
                if (entry2.getKey() == boundaryTouchingParent) {
                    continue;
                }
                for (int label2 : entry2.getValue()) {
                    points.addAll(labelPointsMap.get(label2));
                }
            }

            filledRegions.put(label, points);
        }

        return filledRegions;
    }

    /**
     * given a set of points containing no holes in it, find the clockwise ordered bounding closed curve.
     * Note that the last point has been removed to be consistent with PerimiterFinder2 in not repeating the point.
     * @param points
     * @param imgWidth
     * @param imgHeight
     * @return
     */
    public static PairIntArray mooreTracingWithJacob(Set<Integer> points, int imgWidth, int imgHeight) {

        /*Moore-Neighbor boundary tracing algorithm
        -- identifying a starting pixel on the border as the bottommost leftmost by scanning up
        every row for colomn 0, then every row for column 1, etc until find the first
        pixel in the labeled region, store that position as the first point and save the position of the last non-point
        visited right before the first point.

       -- then the 8-neighbor clockwise scan pattern for the next point begins (while storing the last non-point)
          - the current cell becomes the most recently found point and the offset that it starts with is the offset needed
        between the current cell and the previous non-point from the last 8-neighbor scan.
          - NOTE: one edge case is that the region is on an image boundary.  the offset needs to be allowed to go outside
          of the bounds by 1 pixel in order to reach the next region pixel sometimes.
         */

        if (points.isEmpty()) {
            return new PairIntArray();
        }

        if (imgWidth < 1) {
            throw new IllegalArgumentException("imgWidth must be > 0,  received " + imgWidth);
        }

        if (imgHeight < 1) {
            throw new IllegalArgumentException("imgHeight must be > 0,  received " + imgHeight);
        }

        if (points.size() == 1) {
            PairIntArray p = new PairIntArray();
            int idx = points.iterator().next();
            p.add(idx % imgWidth, idx / imgWidth);
            return p;
        }

        // previous (c, r) for the first point
        int[] beforeCROffset = Arrays.copyOf(offsets8[0], 2);
        int[] startCR = new int[]{-10, -10};

        // current boundary r, c and its last non-labeled point
        int[] cr = new int[]{-10, -10};
        int[] beforeCR = Arrays.copyOf(offsets8[0], 2);

        Map<PairInt, Integer> pCounts = new HashMap<>();

        PairIntArray p = new PairIntArray();

        // find the starter point.  O(w * h)
        int idx;
        for (int c = 0; c < imgWidth; ++c) {
            for (int r = 0; r < imgHeight; ++r) {
                idx = r * imgWidth + c;
                if (points.contains(idx)) {
                    startCR[0] = c;
                    startCR[1] = r;
                    p.add(c, r);
                    pCounts.put(new PairInt(c, r), 1);
                    if (r == 0) {
                        beforeCR[0] = c;
                        beforeCR[1] = r - 1;
                    }
                    break;
                } else {
                    beforeCR[0] = c;
                    beforeCR[1] = r;
                }
            }
            if (p.getN() > 0) {
                break;
            }
        }

        // we start with cr assigned to last point
        System.arraycopy(startCR, 0, cr, 0, 2);
        beforeCROffset[0] = beforeCR[0] - startCR[0];
        beforeCROffset[1] = beforeCR[1] - startCR[1];

        // stopping criteria are that we end up back at the starting point AND the previous non-points point of it is the same
        // as the stopping previous non-points point also
        int[] cr2 = new int[2];
        int idx2;
        boolean done = false;
        int nIter = 0;
        final int nMaxIter = 2 * imgHeight * imgWidth;
        while (!done) {
            if (nIter == nMaxIter) {
                System.out.printf("  ERROR:  exceeded nMaxIter.  cr=(%s), beforeCROffset=(%s), beforeCR=(%s)\n",
                    Arrays.toString(cr), Arrays.toString(beforeCROffset), Arrays.toString(beforeCR));
                break;
            }
            ++nIter;
            // nOffIter is the count of 8 neighbors visited on this scan
            int nOffIter = 0;
            //System.out.printf("current (c=%d,r=%d)\n", cr[0], cr[1]);
            // pt is used to find the offset relative to the new number
            PairInt pt = new PairInt(beforeCR[0] - cr[0], beforeCR[1] - cr[1]);
            {//DEBUG
                if (!offsetIndexMap.containsKey(pt)) {
                    System.out.printf("  ERROR:  cr=(%s), beforeCR=(%s), beforeCROffset=(%s), => offset(x,y=%s)\n",
                        Arrays.toString(cr), Arrays.toString(beforeCR), Arrays.toString(beforeCROffset), pt.toString());
                }
            }

            int offI = offsetIndexMap.get(pt);
            while (nOffIter < 8) {
                int[] offset = offsets16[offI];
                ++offI;
                ++nOffIter;
                for (int j = 0; j < 2; ++j) {
                    cr2[j] = cr[j] + offset[j];
                }
                if (cr2[0] < 0 || cr2[1] < 0 || cr2[0] == imgWidth || cr2[1] == imgHeight) {
                    System.arraycopy(cr2, 0, beforeCR, 0, 2);
                    continue;
                }
                idx2 = cr2[1] * imgWidth + cr2[0];
                //System.out.printf("  (c2=%d,r2=%d)\n", cr2[0], cr2[1]);
                if (!points.contains(idx2)) {
                    System.arraycopy(cr2, 0, beforeCR, 0, 2);
                    continue;
                }
                PairInt c = new PairInt(cr2[0], cr2[1]);

                // check for stopping criteria
                if (cr2[0] == startCR[0] && cr2[1] == startCR[1]) {
                    //System.out.printf("   offI=%d\n", offI);
                    if (offI > 0) {
                        //System.out.printf("   %s - %s = (offset=%d,%d) (start off = (%d, %d)\n", Arrays.toString(offsets16[offI-1]), Arrays.toString(offsets16[offI]),
                        //        offsets16[offI-1][0] - offsets16[offI][0], offsets16[offI-1][1] - offsets16[offI][1], beforeCROffset[0], beforeCROffset[1]);

                        if (offsets16[offI-1][0] - offsets16[offI][0] == beforeCROffset[0] && offsets16[offI-1][1] - offsets16[offI][1] == beforeCROffset[1]) {
                            done = true;
                        }
                    } else if (offset[0] == beforeCROffset[0] && offset[1] == beforeCROffset[1]){
                        done = true;
                    }
                } else if (pCounts.getOrDefault(c, 0) > 0) {
                    // check if previous point == first point, and this point == 2nd point,
                    // then we've reached the start again, but approached it from another direction
                    if (p.getX(1) == cr2[0] && p.getY(1) == cr2[1]
                        && p.getX(p.getN() - 1) == startCR[0] && p.getY(p.getN() - 1) == startCR[1]) {
                        p.removeRange(p.getN() - 1, p.getN() - 1);
                        done = true;
                    }
                }
                if (done) {
                    break;
                }
                //System.out.printf("  * add:(c2=%d, r2=%d), prev:(c=%d, r=%d)\n", cr2[0], cr2[1], beforeCR[0], beforeCR[1]);
                // this is the next point.
                p.add(cr2[0], cr2[1]);
                pCounts.put(c, pCounts.getOrDefault(c, 0) + 1);

                // the next round should start with cr2
                System.arraycopy(cr2, 0, cr, 0, 2);
                break;
            }
        }
        return p;
    }

    /**
     * create a map of all points for each label
     *  given a labeled array whose indexes are image indexes
     *      * where idx = row * imageWidth + col.
     * @param labels labels of each pixel in image
     * @return a map of all points for each label
     */
    public static Map<Integer, Set<Integer>> getLabeledPointsMap(int[] labels) {

        int n = labels.length;

        Map<Integer, Set<Integer>> map = new HashMap<>();

        for (int i = 0; i < n; ++i) {
            map.putIfAbsent(labels[i], new HashSet<>());
            map.get(labels[i]).add(i);
        }

        return map;
    }

    /**
     * create a label adjacency map given a labeled array whose indexes are image indexes
     * where idx = row * imageWidth + col.
     * @param labels labels of each pixel in image
     * @param imgWidth width of image
     * @param imgHeight height of image
     * @return the adjacency map of the labels
     */
    public static Map<Integer, Set<Integer>> getLabelAdjMap(int[] labels, int imgWidth, int imgHeight) {

        int n = labels.length;

        Map<Integer, Set<Integer>> map = new HashMap<>();

        int r2, c2, label1, label2;
        for (int r = 0; r < imgHeight; ++r) {
            for (int c = 0; c < imgWidth; ++c) {
                label1 = labels[r * imgWidth + c];
                for (int[] offset : offsets8) {
                    r2 = r + offset[0];
                    c2 = c + offset[1];
                    if (r2 < 0 || c2 < 0 || r2 == imgHeight || c2 == imgWidth) {
                        continue;
                    }
                    label2 = labels[r2 * imgWidth + c2];
                    if (label1 == label2) {
                        continue;
                    }
                    map.putIfAbsent(label1, new HashSet<>());
                    map.putIfAbsent(label2, new HashSet<>());
                    map.get(label1).add(label2);
                    map.get(label2).add(label1);
                }
            }
        }
        return map;
    }
}
