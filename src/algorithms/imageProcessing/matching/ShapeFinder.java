package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.ORB;
import static algorithms.imageProcessing.matching.ORBMatcher.maxNumberOfGaps;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.TwoDFloatArray;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTree;

/**
 * uses PartialShapeMatcher and search patterns to
 * find the best fitting match between a
 * group of adjacent segmented cells
 * to a template shape.  Any additional information
 * that helps limit the search should be done before
 * this stage and included or excluded in the
 * adjacency map.
 *
 * Also note that the class uses internal cache which requires single
 * threaded use of this class.
 *
 * @author nichole
 */
public class ShapeFinder {

    // partial shape matcher sampling distance
    private final int dp = 1;
    private final boolean useSameNSampl = false;

    private final PairIntArray bounds1;
    private final int[][] ch1;
    private final float scale1;
    private final float sz1;
    private final List<Set<PairInt>> listOfSets2;
    private final List<TwoDIntArray> listOfCH2s;
    private final float scale2;
    private final Map<OneDIntArray, PairIntArray> keyBoundsMap;
    private final float intersectionLimit;

    private final TObjectIntMap<PairInt> pointIndexes2Map;
    private final TIntObjectMap<TIntSet> adj2Map;
    private final List<PairInt> xyCen2List;
    private final List<QuadInt> xyMinMax2List;

    private final int xMax1;
    private final int yMax1;
    private final int xMax2;
    private final int yMax2;

    private SIGMA sigma = SIGMA.ZEROPOINTFIVE;

    /**
     * NOT READY FOR USE.
     *
     * note that this is the data of one image in pyramid1 compared to the data
     * in one image of pyramid2 and it is assumed that one such pair will have
     * template and true match of nearly the same scale.
     * On other words, this instance of ShapeFinder might not hold the true
     * match at same scale, but another instance (holding other scles of the
     * pyramids) should if the true match is
     * present.
     *
     * @param bounds1 clockwise ordered boundaries of template object.
     * @param ch1 1-D color histogram of template object
     * @param scale1 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid1.
     * @param listOfSets2 labeled point sets of dataset2 to be searched for template
     * object.listOfCH2s the color histograms for listOfSets2
     * @param scale2 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid2.
     * @param keyIndexMap an in/out variable for use in combination wtth
     * indexBoundsMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = sorted indexes of listOfSets2, value = the lookup key for
     * indexBoundsMap,
     * @param indexBoundsMap an in/out variable for use in combination wtth
     * keyIndexMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = the lookup key from keyIndexMap, value = clockwise ordered
     * bounding points of the aggregated sets of labels from the key in
     * indexBoundsMap,
     */
    public ShapeFinder(PairIntArray bounds1, int[][] ch1, float scale1, float sz1,
        int xMax1, int yMax1,
        List<Set<PairInt>> listOfSets2, List<TwoDIntArray> listOfCH2s,
        float scale2, 
        Map<OneDIntArray, PairIntArray> keyBoundsMap,
        int xMax2, int yMax2,
        float intersectionLimit) {

        this.bounds1 = bounds1;
        this.ch1 = ch1;
        this.scale1 = scale1;
        this.sz1 = sz1;
        this.listOfSets2 = listOfSets2;
        this.listOfCH2s = listOfCH2s;
        this.scale2 = scale2;
        this.keyBoundsMap = keyBoundsMap;
        this.intersectionLimit = intersectionLimit;
        this.xMax1 = xMax1;
        this.yMax1 = yMax1;
        this.xMax2 = xMax2;
        this.yMax2 = yMax2;

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        this.xyMinMax2List = new ArrayList<QuadInt>();
        this.xyCen2List = new ArrayList<PairInt>();
        this.pointIndexes2Map = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < listOfSets2.size(); ++i) {
            Set<PairInt> set = listOfSets2.get(i);
            for (PairInt p : set) {
                pointIndexes2Map.put(p, i);
            }
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            xyCen2List.add(new PairInt(
                (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1])));
            int[] minMaxXY = MiscMath.findMinMaxXY(set);
            xyMinMax2List.add(new QuadInt(minMaxXY[0], minMaxXY[1], minMaxXY[2],
                minMaxXY[3]));
        }

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        this.adj2Map = new TIntObjectHashMap<TIntSet>();
        for (int idx1 = 0; idx1 < listOfSets2.size(); ++idx1) {
            Set<PairInt> set = listOfSets2.get(idx1);
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                // find adjacent labels that are not idx1
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    int idx2 = pointIndexes2Map.get(p2);
                    if (idx2 == idx1) {
                        continue;
                    }
                    TIntSet set2 = adj2Map.get(idx1);
                    if (set2 == null) {
                        set2 = new TIntHashSet();
                        adj2Map.put(idx1, set2);
                    }
                    set2.add(idx2);
                }
            }
        }
    }

// DEBUG
public static TwoDFloatArray pyr1 = null;
public static TwoDFloatArray pyr2 = null;
public static String lbl = "";
public static int oct1 = -1;
public static int oct2 = -1;

    /**
     * NOT READY FOR USE.
     * uses PartialShapeMatcher and search patterns to
     * find the best fitting match between a
     * group of adjacent segmented cells
     * to a template shape.
     *
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @return
     */
    public ShapeFinderResult findAggregated() {

        /*
        some notes that will be updated when return to the search
        task.

        basically, have a global grid search,
        wrapping local FloydWarshall searches where search
            is the best fitting bounds of aggregated labels sets and
            their color histograms

        step (1) O(N_cells) + O(N_cells_in_bin * log_2(N_cells_in_bin)):
           - 2D bins, that is x and y of size sz1.
           - sort the cells within each bin by x centroid and ties by y centroid
        step (2) local search pattern
            The bigger picture is a sliding window of the 2D bin
            where the core of the search and aggregation is done
            within the 2D bin, but the aggregation can continue
            into the next bin to the right, next bin above,
            or above and to the right.

           this is the part that is changing.
           -- made an impl of Floyd Warshals all pairs.
           -- will impl a dijkstra's for every source
               but with the path constrained to be within
               an association distance of the template object size

        */

        //ShapeFinderResult sr = globalGridSearch();

        ShapeFinderResult sr = multiSourceBFS();

        if (sr != null) {
            print(sr, Long.toString(MiscDebug.getCurrentTimeFormatted()));
        }
        
        return sr;
    }

    private ShapeFinderResult globalGridSearch() {

        int w = xMax2 + 1;
        int h = yMax2 + 1;

        List<ShapeFinderResult> results = new ArrayList<ShapeFinderResult>();

        int nPPB = 1;

        int nX = nPPB * ((int)Math.floor(w/sz1) + 1);
        int nY = nPPB * ((int)Math.floor(h/sz1) + 1);
        int binSz = (int)Math.ceil(sz1);

        QuadTree<Integer, Integer> centroidQT = new QuadTree<Integer, Integer>();

        // populated bin numbers, x, y
        Set<PairInt> binNumbers = new HashSet<PairInt>();

        // store the segmented cell centroids into quadtree
        for (int i = 0; i < xyCen2List.size(); ++i) {
            PairInt xy = xyCen2List.get(i);

            centroidQT.insert(xy.getX(), xy.getY(), Integer.valueOf(i));

            int xBin = xy.getX()/binSz;
            int yBin = xy.getY()/binSz;
            binNumbers.add(new PairInt(xBin, yBin));
        }

        // visit the bins by smallest x,y first,
        //   and start each search with the labeled set which has the
        //   smallest x,y in that bin.
        //   any search can agregate from an adjacent bin, but not beyond
        //   that so a restricted list of set indexes is created
        //   for each search start.

        // sort so that all x for same y are listed, then next x for larger y
        List<PairInt> sortedBinNumbers = new ArrayList<PairInt>(binNumbers);
        Collections.sort(sortedBinNumbers, new XYSort());

        for (int i = 0; i < sortedBinNumbers.size(); ++i) {
            PairInt xyBin = sortedBinNumbers.get(i);
            int startX = xyBin.getX() * binSz;
            int startY = xyBin.getY() * binSz;

            int stopX = startX + binSz;
            if (stopX > (w - 1)) {
                stopX = w - 1;
            }

            int stopY = startY + binSz;
            if (stopY > (h - 1)) {
                stopY = h - 1;
            }

            Interval<Integer> intX = new Interval<Integer>(startX, stopX);
            Interval<Integer> intY = new Interval<Integer>(startY, stopY);
            Interval2D<Integer> rect = new Interval2D<Integer>(intX, intY);

            // centroid list indexes (these are same indexes are for label2 sets)
            List<Integer> indexes = centroidQT.query2D(rect);

            // aggregate adjacent indexes of larger x or larger y
            TIntSet idxs = new TIntHashSet();

            assert(!indexes.isEmpty());
            for (int idx2 : indexes) {
                QuadInt minMax2 =  xyMinMax2List.get(idx2);
                int maxDim = minMaxXY(minMax2);
                if (maxDim <= 1.15 * sz1) {
                    idxs.add(idx2);
                }
            }

            PairInt nextBin = new PairInt(xyBin.getX() + 1, xyBin.getY());
            if (binNumbers.contains(nextBin)) {
                Interval<Integer> intX2 = new Interval<Integer>(startX + binSz,
                    stopX + binSz);
                Interval<Integer> intY2 = new Interval<Integer>(startY, stopY);
                Interval2D<Integer> rect2 = new Interval2D<Integer>(intX2, intY2);
                // centroid list indexes
                List<Integer> indexes2 = centroidQT.query2D(rect2);
                if (indexes2 != null) {
                    for (int idx2 : indexes2) {
                        QuadInt minMax2 = xyMinMax2List.get(idx2);
                        int maxDim = minMaxXY(minMax2);
                        if (maxDim <= 1.15 * sz1) {
                            idxs.add(idx2);
                        }
                    }
                }
            }
            nextBin = new PairInt(xyBin.getX() + 1, xyBin.getY() + 1);
            if (binNumbers.contains(nextBin)) {
                Interval<Integer> intX2 = new Interval<Integer>(startX + binSz,
                    stopX + binSz);
                Interval<Integer> intY2 = new Interval<Integer>(startY + binSz,
                    stopY + binSz);
                Interval2D<Integer> rect2 = new Interval2D<Integer>(intX2, intY2);
                // centroid list indexes
                List<Integer> indexes2 = centroidQT.query2D(rect2);
                if (indexes2 != null) {
                    for (int idx2 : indexes2) {
                        QuadInt minMax2 = xyMinMax2List.get(idx2);
                        int maxDim = minMaxXY(minMax2);
                        if (maxDim <= 1.15 * sz1) {
                            idxs.add(idx2);
                        }
                    }
                }
            }
            nextBin = new PairInt(xyBin.getX(), xyBin.getY() + 1);
            if (binNumbers.contains(nextBin)) {
                Interval<Integer> intX2 = new Interval<Integer>(startX, stopX);
                Interval<Integer> intY2 = new Interval<Integer>(startY + binSz,
                    stopY + binSz);
                Interval2D<Integer> rect2 = new Interval2D<Integer>(intX2, intY2);
                // centroid list indexes
                List<Integer> indexes2 = centroidQT.query2D(rect2);
                if (indexes2 != null) {
                    for (int idx2 : indexes2) {
                        QuadInt minMax2 = xyMinMax2List.get(idx2);
                        int maxDim = minMaxXY(minMax2);
                        if (maxDim <= 1.15 * sz1) {
                            idxs.add(idx2);
                        }
                    }
                }
            }

            // TODO: invoke the local search with restricted label range
            ShapeFinderResult r = searchUsingFloydWarshall(idxs);

            //TODO: when replace aspectj w/ another AOP library, assert that
            //   all list 2 set labels are covered in adjIdxs and srchIdx,
            //   after all FW invocations.

            if (r == null) {
                continue;
            }

            int sz2 = ORBMatcher.calculateObjectSize(r.bounds2);

            //if ((sz1 > sz2 && Math.abs(sz1 / sz2) > 1.4) ||
            //    (sz2 > sz1 && Math.abs(sz2 / sz1) > 1.4)) {
            if ((sz1 > sz2 && Math.abs(sz1 / sz2) > 1.15) ||
                (sz2 > sz1 && Math.abs(sz2 / sz1) > 1.15)) {
                continue;
            }

            results.add(r);

            //if (debug) {
            /*{
                String strI = Integer.toString(startX);
                while (strI.length() < 3) {
                    strI = "0" + strI;
                }
                String strJ = Integer.toString(startY);
                while (strJ.length() < 3) {
                    strJ = "0" + strJ;
                }
                String str = "bin_" + strI + "_" + strJ;
                print(r, str); 
            }*/
        }

        if (results.isEmpty()) {
            return null;
        } else if (results.size() == 1) {
            return results.get(0);
        }

        // once though to find max diff chord sum and max dist
        double maxAvgDiffChord = Double.MIN_VALUE;
        double maxAvgDist = Double.MIN_VALUE;
        for (int i = 0; i < results.size(); ++i) {
            ShapeFinderResult r = results.get(i);
            double avgCD = r.getChordDiffSum() / (double) r.getNumberOfMatches();
            if (avgCD > maxAvgDiffChord) {
                maxAvgDiffChord = avgCD;
            }
            double avgDist = r.getDistSum() / (double) r.getNumberOfMatches();
            if (avgDist > maxAvgDist) {
                maxAvgDist = avgDist;
            }
        }

        double minCost = Double.MAX_VALUE;
        ShapeFinderResult minCostR = null;

        for (int i = 0; i < results.size(); ++i) {

            ShapeFinderResult r = results.get(i);

            int nb1 = Math.round((float) r.bounds1.getN() / (float) dp);

            float np = r.getNumberOfMatches();

            float countComp = 1.0F - (np / (float) nb1);
            double chordComp = ((float) r.getChordDiffSum() / np) / maxAvgDiffChord;
            double avgDist = r.getDistSum() / np;
            double distComp = avgDist / maxAvgDist;

            int lGap = maxNumberOfGaps(r.bounds1, r)/dp;
            float gCountComp = (float)lGap/(float)nb1;

            double sd = chordComp + countComp + gCountComp + distComp
                + r.intersection;

            if (sd < minCost) {
                minCost = sd;
                minCostR = r;
            }
        }

        return minCostR;
    }

    private ShapeFinderResult searchUsingFloydWarshall(TIntSet idxs) {

        /*
        the Floyd-Warshal all oairs pattern is adapted from the pseudocode in
        Cormen et al. "Intro to Algorithms".
        the runtime is O(V^3).
        */

        // -- each bounds can be filtered for size near sz1

        // -- need to assert contiguous before partial shape matcher,
        //       so, can only compare if adjacent to a set member

        // -- the dynamic programming array can be made smaller by making
        //    a second assignment of labels here and the lookup maps

        TIntIntMap i2ToOrigIndexMap = new TIntIntHashMap();
        TIntIntMap origIndexToI2Map = new TIntIntHashMap();

        TIntIterator iter = idxs.iterator();
        while (iter.hasNext()) {
            int idx0 = iter.next();
            int idx2 = origIndexToI2Map.size();
            origIndexToI2Map.put(idx0, idx2);
            i2ToOrigIndexMap.put(idx2, idx0);
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        ColorHistogram cHist = new ColorHistogram();

        int nb1 = Math.round((float) bounds1.getN() / (float) dp);

        // for the cost calculations, need the maximum average chord diff sums
        //   and dist sums
        double maxAvgDiffChord = Double.MIN_VALUE;
        double maxAvgDist = Double.MIN_VALUE;

        // ----- initialize the arrays ----------

        int n = origIndexToI2Map.size();

        ShapeFinderResult[][] results = new ShapeFinderResult[n][];

        for (int i = 0; i < n; ++i) {

            results[i] = new ShapeFinderResult[n];

            int idxI0 = i2ToOrigIndexMap.get(i);

            // note, the size filter before the invocation of this method
            // assures that each set is within size sz1

            int[][] ch2 = listOfCH2s.get(idxI0).a;
            float intersection = cHist.intersection(ch1, ch2);
            if (intersection < intersectionLimit) {
                continue;
            }

            // ----- calculate the cost of the shape match as init "w"s -----
            OneDIntArray keysI = new OneDIntArray(new int[]{idxI0});
            PairIntArray boundsI;
            if (keyBoundsMap.containsKey(keysI)) {
                boundsI = keyBoundsMap.get(keysI);
            } else {
                Set<PairInt> set1 = this.listOfSets2.get(idxI0);
                boundsI = imageProcessor.extractSmoothedOrderedBoundary(
                    new HashSet(set1), sigma, xMax2 + 1, yMax2 + 1);
                keyBoundsMap.put(keysI, boundsI);
            }

            int sz2 = ORBMatcher.calculateObjectSize(boundsI);

            if (sz2 > 1.15 * sz1) {
                continue;
            }

            PartialShapeMatcher matcher = new PartialShapeMatcher();
            matcher._overrideToThreshhold(0.2f);
            matcher.overrideSamplingDistance(dp);
            matcher.setToRemoveOutliers();
            //matcher.setToDebug();
            PartialShapeMatcher.Result r = matcher.match(bounds1, boundsI);
            if (r == null) {
                continue;
            }
            double avgCD = r.getChordDiffSum() / (double) r.getNumberOfMatches();
            if (avgCD > maxAvgDiffChord) {
                maxAvgDiffChord = avgCD;
            }
            double avgDist = r.getDistSum() / (double) r.getNumberOfMatches();
            if (avgDist > maxAvgDist) {
                maxAvgDist = avgDist;
            }
            ShapeFinderResult sr = new ShapeFinderResult(r, bounds1, boundsI,
                keysI.a);
            sr.intersection = intersection;

            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    // replacing default 0
                    results[i][j] = sr;
                    continue;
                }
                int idxJ0 = i2ToOrigIndexMap.get(j);
                if (adj2Map.containsKey(idxI0) &&
                    adj2Map.get(idxI0).contains(idxJ0)) {
                    results[i][j] = sr;
                }
            }
        }

        // ----- populate the arrays ----------

        // W is an nxn matrix of w[i][j]
        //     where each entry is the edge weight if any b etween i and j.
        // W has values 0 if i==j and inf where there is no connection.
        // D is an nxn matrix of d[i][j]
        //     where each entry is the weight of the shortest path between i and j.
        // d[i][j]_k = w[i][j]                                   if k = 0
        //           = min( d[i][k]_(k-1) + d[k][j]_(k-1) )      if k >= 1

        for (int k = 0; k < n; k++) {
            int idxK0 = i2ToOrigIndexMap.get(k);
            for (int i = 0; i < n; i++) {
                int idxI0 = i2ToOrigIndexMap.get(i);
                for (int j = 0; j < n; j++) {
                    int idxJ0 = i2ToOrigIndexMap.get(j);
                    if (i == j) {
                        //d[i][j] = 0;
                        continue;
                    }
                    ShapeFinderResult rIJ = results[i][j];
                    ShapeFinderResult rIK = results[i][k];
                    ShapeFinderResult rKJ = results[k][j];

                    if (rIK == null && rKJ == null) {
                        // i,j remains as is
                        continue;
                    }

                    ShapeFinderResult rIKPlusKJ = aggregateAndMatch(rIK, rKJ);

                    if (rIKPlusKJ == null) {
                        // i,j remains as is
                        continue;
                    }

                    if (rIJ == null) {
                        results[i][j] = rIKPlusKJ;
                        continue;
                    }

                    //TODO: consider whether to re-do the cost everytime
                    // with new max.  should be fine to re-calc as needed and
                    // then at the end.
                    double avgCD = rIJ.getChordDiffSum() / (double) rIJ.getNumberOfMatches();
                    if (avgCD > maxAvgDiffChord) {
                        maxAvgDiffChord = avgCD;
                    }
                    double avgDist = rIJ.getDistSum() / (double) rIJ.getNumberOfMatches();
                    if (avgDist > maxAvgDist) {
                        maxAvgDist = avgDist;
                    }
                    avgCD = rIKPlusKJ.getChordDiffSum() /
                        (double) rIKPlusKJ.getNumberOfMatches();
                    if (avgCD > maxAvgDiffChord) {
                        maxAvgDiffChord = avgCD;
                    }
                    avgDist = rIKPlusKJ.getDistSum() /
                        (double) rIKPlusKJ.getNumberOfMatches();
                    if (avgDist > maxAvgDist) {
                        maxAvgDist = avgDist;
                    }

                    float np = rIJ.getNumberOfMatches();
                    float countComp = 1.0F - (np / (float) nb1);
                    double chordComp = ((float) rIJ.getChordDiffSum() / np)
                        / maxAvgDiffChord;
                    avgDist = rIJ.getDistSum() / np;
                    double distComp = avgDist / maxAvgDist;
                    int lGap = ORBMatcher.maxNumberOfGaps(bounds1, rIJ)/dp;
                    float gCountComp = (float)lGap/(float)nb1;
                    double sd = chordComp + countComp + gCountComp + distComp
                        + rIJ.intersection;
                    double s0 = sd;

                    np = rIKPlusKJ.getNumberOfMatches();
                    countComp = 1.0F - (np / (float) nb1);
                    chordComp = ((float) rIKPlusKJ.getChordDiffSum() / np)
                        / maxAvgDiffChord;
                    avgDist = rIKPlusKJ.getDistSum() / np;
                    distComp = avgDist / maxAvgDist;
                    lGap = ORBMatcher.maxNumberOfGaps(bounds1, rIKPlusKJ)/dp;
                    gCountComp = (float)lGap/(float)nb1;
                    sd = chordComp + countComp + gCountComp + distComp
                        + rIKPlusKJ.intersection;
                    double s1 = sd;

                    if (s0 <= s1) {
                       // d[i][j] = s0;
                    } else {
                        //d[i][j] = s1;
                        // TODO: revisit this
                        //prev[i][j] = prev[k][j];
                        results[i][j] = rIKPlusKJ;
                    }
                }
            }
        }

        // --- review results, updating cost and return the mincost result

        double minCost = Double.MAX_VALUE;
        ShapeFinderResult minCostSR = null;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                ShapeFinderResult sr = results[i][j];
                if (sr == null) {
                    continue;
                }
                float np = sr.getNumberOfMatches();
                float countComp = 1.0F - (np / (float) nb1);
                double chordComp = ((float) sr.getChordDiffSum() / np)
                    / maxAvgDiffChord;
                double avgDist = sr.getDistSum() / np;
                double distComp = avgDist / maxAvgDist;
                int lGap = ORBMatcher.maxNumberOfGaps(bounds1, sr) / dp;
                float gCountComp = (float) lGap / (float) nb1;
                double sd = chordComp + countComp + gCountComp + distComp
                    + sr.intersection;

                if (sd < minCost) {
                    minCost = sd;
                    minCostSR = sr;
                }
            }
        }

        if (minCostSR != null) {
            minCostSR.data = new Object[]{Double.valueOf(minCost)};
        }

        return minCostSR;
    }

    private int findSmallestXYCentroid(List<Integer> indexes,
        List<PairInt> xyCenList) {

        int minX = Integer.MAX_VALUE;
        int minX_Y = Integer.MAX_VALUE;
        int minIdx = -1;

        for (int i = 0; i < indexes.size(); ++i) {
            int idx = indexes.get(i);
            int x = xyCenList.get(idx).getX();
            int y = xyCenList.get(idx).getY();
            if (x < minX) {
                minX = x;
                minX_Y = y;
                minIdx = idx;
            } else if (x == minX) {
                if (y < minX_Y) {
                    minX_Y = y;
                    minIdx = idx;
                }
            }
        }

        assert(minIdx > -1);

        return minIdx;
    }

    private int minMaxXY(QuadInt srchXYMinMax, QuadInt otherXYMinMax) {
        int diffX = Math.max(
            Math.abs(srchXYMinMax.getA() - otherXYMinMax.getB()),
            Math.abs(srchXYMinMax.getB() - otherXYMinMax.getA()));

        int diffY = Math.max(
            Math.abs(srchXYMinMax.getC() - otherXYMinMax.getD()),
            Math.abs(srchXYMinMax.getD() - otherXYMinMax.getC()));

        // consider taking the diagonal here
        return Math.max(diffX, diffY);
    }

    private int minMaxXY(QuadInt xyMinMax) {

        int diffX = Math.abs(xyMinMax.getA() - xyMinMax.getB());

        int diffY = Math.abs(xyMinMax.getC() - xyMinMax.getD());

        // consider taking the diagonal here
        return Math.max(diffX, diffY);
    }

    private ShapeFinderResult aggregateAndMatch(ShapeFinderResult rIK,
        ShapeFinderResult rKJ) {

        // -- either might be null but not both
        if (rIK == null) {
            return rKJ;
        } else if (rKJ == null) {
            return rIK;
        }

        // if both have same keys, assert same results and return one of them
        int[] l1 = rIK.labels2;
        int[] l2 = rKJ.labels2;
        if (Arrays.equals(l1, l2)) {
            // TODO: assert same content
            return rIK;
        }

        // -- check that the two have at least one adjacent set between them
        boolean foundAdj = false;
        for (int i = 0; i < l1.length; ++i) {
            int idx1 = l1[i];
            for (int j = 0; j < l2.length; ++j) {
                int idx2 = l2[j];
                if (this.adj2Map.containsKey(idx1) && adj2Map.get(idx1).contains(idx2)) {
                    foundAdj = true;
                    break;
                }
            }
            if (foundAdj) {
                break;
            }
        }
        if (!foundAdj) {
            return null;
        }

        TIntSet combIdxs = new TIntHashSet(l1.length + l2.length);
        combIdxs.addAll(l1);
        combIdxs.addAll(l2);

        ColorHistogram cHist = new ColorHistogram();

        int[][] ch2 = cHist.createWithDefaultSize();
        TIntIterator iter = combIdxs.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            cHist.add2To1(ch2, listOfCH2s.get(idx).a);
        }

        float intersection = cHist.intersection(ch1, ch2);
        if (intersection < intersectionLimit) {
            return null;
        }

        OneDIntArray keysI = new OneDIntArray(
            combIdxs.toArray(new int[combIdxs.size()]));
        Arrays.sort(keysI.a);

        PairIntArray boundsI;
        if (keyBoundsMap.containsKey(keysI)) {
            boundsI = keyBoundsMap.get(keysI);
        } else {
            Set<PairInt> combSet = new HashSet<PairInt>();
            iter = combIdxs.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                combSet.addAll(this.listOfSets2.get(idx));
            }
            ImageProcessor imageProcessor = new ImageProcessor();
            boundsI = imageProcessor.extractSmoothedOrderedBoundary(
                new HashSet(combSet), sigma, xMax2 + 1, yMax2 + 1);
            keyBoundsMap.put(keysI, boundsI);
        }

        int sz2 = ORBMatcher.calculateObjectSize(boundsI);

        if (sz2 > 1.15 * sz1) {
            return null;
        }

        PartialShapeMatcher matcher = new PartialShapeMatcher();
        matcher.overrideSamplingDistance(dp);
        //matcher.setToDebug();
        //matcher.setToUseSameNumberOfPoints();
        PartialShapeMatcher.Result r = matcher.match(bounds1, boundsI);

        if (r == null) {
            return null;
        }

        ShapeFinderResult sr = new ShapeFinderResult(r, bounds1, boundsI,
            keysI.a);
        sr.intersection = intersection;

        return sr;
    }

    private ShapeFinderResult multiSourceBFS() {

        // -- make a cache of search results
        //    to resuse
        // -- launch a BFS search from each label2 set.

        Map<OneDIntArray, ShapeFinderResult> cacheResults
            = new HashMap<OneDIntArray, ShapeFinderResult>();

        TIntSet idxs = new TIntHashSet();

        double[] maxChordAndDiff = matchSingly(cacheResults, idxs);
        
        TIntIterator iter = idxs.iterator();

        List<ShapeFinderResult> results = new ArrayList<ShapeFinderResult>();

        while (iter.hasNext()) {

            int i = iter.next();

            ShapeFinderResult sr = bfs(i, cacheResults, idxs,
                maxChordAndDiff);

            if (sr != null) {
                results.add(sr);               
            }
        }

        ShapeFinderResult minCostR = null;
        double minCost = Double.MAX_VALUE;
        
        for (int i = 0; i < results.size(); ++i) {
            
            ShapeFinderResult r = results.get(i);
            
            if (r == null) {
                continue;
            }
            double sd = calcCost(r, maxChordAndDiff[0], maxChordAndDiff[1]);
        
            if (sd < minCost) {
                minCost = sd;
                minCostR = r;
            }
        }

        return minCostR;
    }

    // match the label2 sets that are near size sz1 and cache the
    // results and return the maximum chord diff avg and maximum
    // dist avg
    private double[] matchSingly(Map<OneDIntArray, ShapeFinderResult> 
        cacheResults, TIntSet outIdxs) {

        // once though to find max diff chord sum and max dist
        double maxAvgDiffChord = Double.MIN_VALUE;
        double maxAvgDist = Double.MIN_VALUE;

        ImageProcessor imageProcessor = new ImageProcessor();
        ColorHistogram cHist = new ColorHistogram();
        
        for (int i = 0; i < listOfSets2.size(); ++i) {

            int[][] ch2 = listOfCH2s.get(i).a;
            float intersection = cHist.intersection(ch1, ch2);
            if (intersection < intersectionLimit) {                
                continue;
            }

            OneDIntArray keysI = new OneDIntArray(new int[]{i});
            PairIntArray boundsI;
            if (keyBoundsMap.containsKey(keysI)) {
                boundsI = keyBoundsMap.get(keysI);
            } else {
                Set<PairInt> set1 = this.listOfSets2.get(i);
                boundsI = imageProcessor.extractSmoothedOrderedBoundary(
                    new HashSet(set1), sigma, xMax2 + 1, yMax2 + 1);
                keyBoundsMap.put(keysI, boundsI);
            }

            int sz2 = ORBMatcher.calculateObjectSize(boundsI);
       
            if (sz2 > 1.15 * sz1) {                
                continue;
            }

            PartialShapeMatcher matcher = new PartialShapeMatcher();
            matcher.overrideSamplingDistance(dp);
            //matcher.setToDebug();
            //matcher.setToUseSameNumberOfPoints();
            PartialShapeMatcher.Result r = matcher.match(bounds1, boundsI);
            if (r == null) {                
                continue;
            }
            
            outIdxs.add(i);

            double avgCD = r.getChordDiffSum() / (double) r.getNumberOfMatches();
            if (avgCD > maxAvgDiffChord) {
                maxAvgDiffChord = avgCD;
            }
            double avgDist = r.getDistSum() / (double) r.getNumberOfMatches();
            if (avgDist > maxAvgDist) {
                maxAvgDist = avgDist;
            }
            ShapeFinderResult sr = new ShapeFinderResult(r, bounds1, boundsI,
                keysI.a);
            sr.intersection = intersection;
            
            cacheResults.put(keysI, sr);
        }

        return new double[]{maxAvgDiffChord, maxAvgDist};
    }

    private ShapeFinderResult bfs(int srcIdx, 
        Map<OneDIntArray, ShapeFinderResult> cacheResults, TIntSet idxs,
        double[] maxChordAndDiff) {

        OneDIntArray srcKeys = new OneDIntArray(new int[]{srcIdx});
        if (!cacheResults.containsKey(srcKeys)) {
            return null;
        }

        // only searching within radius of sz1

        // TODO: re-map the indexes to 0 to idxs.size
        
        int n = listOfSets2.size();
        double[] dist = new double[n];
        ShapeFinderResult[] results = new ShapeFinderResult[n];
        // 0 = white, 1 = grey, 2 = visited
        int[] visited = new int[n];
        Arrays.fill(dist, Double.MAX_VALUE);

        visited[srcIdx] = 1;
        results[srcIdx] = cacheResults.get(srcKeys);
        dist[srcIdx] = calcCost(results[srcIdx], 
            maxChordAndDiff[0],
            maxChordAndDiff[1]);

        // make a FIFO queue.
        ArrayDeque<Integer> queue = new ArrayDeque<Integer>();

        queue.add(Integer.valueOf(srcIdx));

        double mc = maxChordAndDiff[0];
        double md = maxChordAndDiff[1];
        
        QuadInt srcMinMax2 = xyMinMax2List.get(srcIdx);
        
        while (!queue.isEmpty()) {

            Integer uIndex = queue.pop();
            int uIdx = uIndex.intValue();
           
            TIntSet adjSet = adj2Map.get(uIdx);
            
            if (adjSet == null) {
                continue;
            }

            TIntIterator iter = adjSet.iterator();
            while (iter.hasNext()) {
                int vIdx = iter.next();
           
                if (!idxs.contains(vIdx)) {
                    continue;
                }
                if (visited[vIdx] == 2) {
                    continue;
                }
                if (results[vIdx] != null && 
                    (Arrays.binarySearch(results[vIdx].labels2,
                        uIdx) > -1)) {
                    continue;
                }
           
                assert(results[uIdx] != null);
                
                ShapeFinderResult rV = results[vIdx];
                double sdV = dist[vIdx];
                if (rV == null) {
                    OneDIntArray vKeys = new OneDIntArray(new int[]{vIdx});
                    rV = cacheResults.get(vKeys);
                }
               
                ShapeFinderResult uPlusV = aggregateAndMatch(results[uIdx], rV);
           
                if (uPlusV == null) {
                    continue;
                }
                
                // store results in cache
                OneDIntArray uvKeys = new OneDIntArray(
                    Arrays.copyOf(uPlusV.labels2, uPlusV.labels2.length));
                cacheResults.put(uvKeys, uPlusV);
                
                int sz2 = ORBMatcher.calculateObjectSize(uPlusV.bounds2);

                if (sz2 > 1.15 * sz1) {
                    continue;
                }

                double sdUPlusV = calcCost(uPlusV, maxChordAndDiff[0],
                    maxChordAndDiff[1]);
        
                //visited[vIdx] = 1;
                if (sdUPlusV < sdV) {
                    
                    dist[vIdx] = sdUPlusV;
                    results[vIdx] = uPlusV;
                    queue.add(vIdx);
                    
                    // update the maxes 
                    double avgCD = uPlusV.getChordDiffSum() / 
                        (double) uPlusV.getNumberOfMatches();
                    if (avgCD > mc) {
                        mc = avgCD;
                    }
                    double avgDist = uPlusV.getDistSum() /
                        (double) uPlusV.getNumberOfMatches();
                    if (avgDist > md) {
                        md = avgDist;
                    }
                }
            }
            
            visited[uIdx] = 2;
        }
        
        // re-do costs, and update max vars
        maxChordAndDiff[0] = mc;
        maxChordAndDiff[1] = md;

        ShapeFinderResult minCostR = null;
        double minCost = Double.MAX_VALUE;
        
        for (int i = 0; i < results.length; ++i) {
            
            ShapeFinderResult r = results[i];
           
            if (r == null) {
                continue;
            }
            double sd = calcCost(r, maxChordAndDiff[0], maxChordAndDiff[1]);
       
            if (sd < minCost) {
                minCost = sd;
                minCostR = r;
            }
        }
{
ShapeFinderResult r = minCostR;
int nb1 = Math.round((float) r.bounds1.getN() / (float) dp);
float np = r.getNumberOfMatches();
int lGap = maxNumberOfGaps(r.bounds1, r)/dp;
float gCountComp = (float)lGap/(float)nb1;

System.out.println("min: cost=" + minCost + " count=" + np + 
" chordAvg=" + ((float) r.getChordDiffSum() / np) +
" distAvg=" + (float)(r.getDistSum() / np) +
" gapAvg=" + gCountComp);
}
        return minCostR;
    }

    private double calcCost(ShapeFinderResult r, double maxAvgDiffChord, 
        double maxAvgDist) {
        
        int nb1 = Math.round((float) r.bounds1.getN() / (float) dp);

        float np = r.getNumberOfMatches();

        float countComp = 1.0F - (np / (float) nb1);
        double chordComp = ((float) r.getChordDiffSum() / np) / maxAvgDiffChord;
        double avgDist = r.getDistSum() / np;
        double distComp = avgDist / maxAvgDist;

        int lGap = maxNumberOfGaps(r.bounds1, r)/dp;
        float gCountComp = (float)lGap/(float)nb1;

        // these should be squared
        double sd = chordComp + countComp + gCountComp + distComp
            + r.intersection;
        
        return sd;
    }

    private void print(ShapeFinderResult r, String lbl2) {
        
        Image img1 = ORB.convertToImage(pyr1);
        Image img2 = ORB.convertToImage(pyr2);
        try {
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                //r.bounds1, r.bounds2);
                img1, img2);
            for (int ii = 0; ii < r.getNumberOfMatches(); ++ii) {
                int idx1 = r.getIdx1(ii);
                int idx2 = r.getIdx2(ii);
                int x1 = r.bounds1.getX(idx1);
                int y1 = r.bounds1.getY(idx1);
                int x2 = r.bounds2.getX(idx2);
                int y2 = r.bounds2.getY(idx2);
                if ((ii % 4) == 0) {
                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                }
            }
            String str = lbl + "bin_" + lbl2;
            String filePath = plotter.writeImage("_shape_" + str);
        } catch (Throwable t) {
        }
    }

    private void print2(ShapeFinderResult r, String lbl2) {
        
        Image img1 = ORB.convertToImage(pyr1);
        Image img2 = ORB.convertToImage(pyr2);
        try {
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                r.bounds1, r.bounds2);
                //img1, img2);
            for (int ii = 0; ii < r.getNumberOfMatches(); ++ii) {
                int idx1 = r.getIdx1(ii);
                int idx2 = r.getIdx2(ii);
                int x1 = r.bounds1.getX(idx1);
                int y1 = r.bounds1.getY(idx1);
                int x2 = r.bounds2.getX(idx2);
                int y2 = r.bounds2.getY(idx2);
                if ((ii % 4) == 0) {
                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                }
            }
            String str = lbl + "bin_" + lbl2;
            String filePath = plotter.writeImage("_shape_" + str);
        } catch (Throwable t) {
        }
    }

    public static class ShapeFinderResult extends Result {

        protected float intersection;
        protected final PairIntArray bounds1;
        protected final PairIntArray bounds2;
        protected final int[] labels2;

        public ShapeFinderResult(Result r, PairIntArray bounds1,
            PairIntArray bounds2, int[] keys) {

            super(r.n1, r.n2, r.getOriginalN1());

            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
            this.chordDiffSum = r.chordDiffSum;
            this.distSum = r.distSum;
            this.chordsNeedUpdates = r.chordsNeedUpdates;
            this.idx1s = r.idx1s;
            this.idx2s = r.idx2s;
            this.params = r.params;
            this.data = r.data;
            this.labels2 = keys;
        }

        public ShapeFinderResult(int n1, int n2, int offset, PairIntArray bounds1,
            PairIntArray bounds2, int[] keys) {

            super(n1, n2, offset);

            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
            this.labels2 = keys;
        }

        public ShapeFinderResult(int n1, int n2, int nOriginal,
            int offset, PairIntArray bounds1,
            PairIntArray bounds2, int[] keys) {

            super(n1, n2, nOriginal);

            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
            this.labels2 = keys;
        }
    }

    // sorted so that all X of same Y are first, then all X of larger Y
    private class XYSort implements Comparator<PairInt> {

        @Override
        public int compare(PairInt o1, PairInt o2) {

            int comp = Integer.compare(o1.getY(), o2.getY());
            if (comp != 0) {
                return comp;
            }
            comp = Integer.compare(o1.getX(), o2.getX());
            return comp;
        }

    }
}
