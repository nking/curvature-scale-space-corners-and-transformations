package algorithms.imageProcessing.matching;

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
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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
public class ShapeFinder2 {

    // partial shape matcher sampling distance
    private final int dp = 1;
    private final boolean useSameNSampl = false;

    private final PairIntArray bounds1;
    private final float scale1;
    private final float sz1;
    private final TIntObjectMap<Set<PairInt>> mapOfSets2;
    private final float scale2;
    private final Map<OneDIntArray, PairIntArray> keyBoundsMap;

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
     * @param mapOfSets2 labeled point sets of dataset2 to be searched for template
     * @param scale2 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid2.
     * @param keyIndexMap an in/out variable for use in combination wtth
     * indexBoundsMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = sorted keys of mapOfSets2, value = the lookup key for
     * indexBoundsMap,
     * @param indexBoundsMap an in/out variable for use in combination wtth
     * keyIndexMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = the lookup key from keyIndexMap, value = clockwise ordered
     * bounding points of the aggregated sets of labels from the key in
     * indexBoundsMap,
     */
    public ShapeFinder2(PairIntArray bounds1, float scale1, float sz1,
        int xMax1, int yMax1,
        TIntObjectMap<Set<PairInt>> mapOfSets2,
        float scale2, 
        Map<OneDIntArray, PairIntArray> keyBoundsMap,
        int xMax2, int yMax2) {

        this.bounds1 = bounds1;
        this.scale1 = scale1;
        this.sz1 = sz1;
        this.mapOfSets2 = mapOfSets2;
        this.scale2 = scale2;
        this.keyBoundsMap = keyBoundsMap;
        this.xMax1 = xMax1;
        this.yMax1 = yMax1;
        this.xMax2 = xMax2;
        this.yMax2 = yMax2;

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        this.xyMinMax2List = new ArrayList<QuadInt>();
        this.xyCen2List = new ArrayList<PairInt>();
        this.pointIndexes2Map = new TObjectIntHashMap<PairInt>();
        
        TIntObjectIterator<Set<PairInt>> iter = mapOfSets2.iterator();
        for (int i = 0; i < mapOfSets2.size(); ++i) {
            iter.advance();
            int segIdx = iter.key();
            Set<PairInt> set = iter.value();
            for (PairInt p : set) {
                pointIndexes2Map.put(p, segIdx);
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
        
        iter = mapOfSets2.iterator();
        for (int i = 0; i < mapOfSets2.size(); ++i) {
            iter.advance();
            int segIdx = iter.key();
            Set<PairInt> set = iter.value();
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                // find adjacent labels that are not idx1
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    int segIdx2 = pointIndexes2Map.get(p2);
                    if (segIdx2 == segIdx) {
                        continue;
                    }
                    TIntSet set2 = adj2Map.get(segIdx);
                    if (set2 == null) {
                        set2 = new TIntHashSet();
                        adj2Map.put(segIdx, set2);
                    }
                    set2.add(segIdx2);
                }
            }
        }
    }

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

        ShapeFinderResult sr = multiSourceBFS();

        return sr;
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

        OneDIntArray keysI = new OneDIntArray(
            combIdxs.toArray(new int[combIdxs.size()]));
        Arrays.sort(keysI.a);

        PairIntArray boundsI;
        if (keyBoundsMap.containsKey(keysI)) {
            boundsI = keyBoundsMap.get(keysI);
        } else {
            Set<PairInt> combSet = new HashSet<PairInt>();
            TIntIterator iter = combIdxs.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                combSet.addAll(this.mapOfSets2.get(idx));
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
        
        TIntObjectIterator<Set<PairInt>> iter = mapOfSets2.iterator();
        for (int i = 0; i < mapOfSets2.size(); ++i) {
            iter.advance();
            int segIdx = iter.key();

            OneDIntArray keysI = new OneDIntArray(new int[]{segIdx});
            PairIntArray boundsI;
            if (keyBoundsMap.containsKey(keysI)) {
                boundsI = keyBoundsMap.get(keysI);
            } else {
                Set<PairInt> set1 = iter.value();
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
        
        int n = mapOfSets2.size();
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
        double sd = chordComp + countComp + gCountComp + distComp;
        
        return sd;
    }

    public static class ShapeFinderResult extends Result {

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
