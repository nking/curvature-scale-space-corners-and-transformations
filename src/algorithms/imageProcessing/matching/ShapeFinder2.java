package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.misc.Misc;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
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
 * NOTE: not ready for use.   was an experiment before segmentation
 * was improved with MSER and before started using HOGs.
 * This will be re-done with the better segmentation and then
 * the better segmentation + HOGs, HCPT, color histograms, etc.
 * The current MSER matcher uses everything except shape.
 * 
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
     * match at same scale, but another instance (holding other scales of the
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

        this.pointIndexes2Map = new TObjectIntHashMap<PairInt>();
        
        TIntObjectIterator<Set<PairInt>> iter = mapOfSets2.iterator();
        for (int i = 0; i < mapOfSets2.size(); ++i) {
            iter.advance();
            int segIdx = iter.key();
            Set<PairInt> set = iter.value();
            for (PairInt p : set) {
                pointIndexes2Map.put(p, segIdx);
            }
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
     * Editing this now... the code no longer checks for a restricted 
     * distance or size in aggregating the labeled regions - it's up
     * to the user to filter the input sets given to this instance
     * themselves for those purposes.
     * 
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @return
     */
    public ShapeFinderResult findAggregated() {
       
        ShapeFinderResult sr = multiSourceBFS();

        return sr;
    }

    /**
     * NOT READY FOR USE.
     * uses PartialShapeMatcher and search patterns to
     * find the best fitting match between a
     * group of adjacent segmented cells
     * to a template shape.
     *
     * Editing this now... the code no longer checks for a restricted 
     * distance or size in aggregating the labeled regions - it's up
     * to the user to filter the input sets given to this instance
     * themselves for those purposes.
     * 
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @return
     */
    public ShapeFinderResult findAggregated(int srcIdx) {
       
        System.out.println("searching " + mapOfSets2.size() + " labeled regions");
        
        Map<OneDIntArray, ShapeFinderResult> cacheResults
            = new HashMap<OneDIntArray, ShapeFinderResult>();

        TIntSet idxs = new TIntHashSet();

        long t0 = System.currentTimeMillis();
        
        double[] maxChordAvgDiff = new double[]{matchSingly(cacheResults, idxs)};
        
        ShapeFinderResult sr = bfs(srcIdx, cacheResults, idxs, maxChordAvgDiff);

        long t1 = System.currentTimeMillis();
        
        System.out.println("a single bfs took " + ((t1 - t0)/1000) + " seconds");
        
        return sr;
    }

    @SuppressWarnings({"unchecked"})
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
                new HashSet<PairInt>(combSet), sigma, xMax2 + 1, yMax2 + 1);
            keyBoundsMap.put(keysI, boundsI);
        }

        PartialShapeMatcher matcher = new PartialShapeMatcher();
        matcher.overrideSamplingDistance(dp);
        matcher._overrideToThreshhold(0.2f);
        matcher.setToRemoveOutliers();
        //matcher.setToDebug();
        PartialShapeMatcher.Result r = matcher.match(bounds1, boundsI);

        if (r == null) {
            return null;
        }

        ShapeFinderResult sr = new ShapeFinderResult(r, bounds1, boundsI,
            keysI.a);

        return sr;
    }

    private ShapeFinderResult multiSourceBFS() {

        // -- make a cache of search results to resuse
        // -- launch a BFS search from each label2 set.

        System.out.println("searching " + mapOfSets2.size() + " labeled regions");
        
        Map<OneDIntArray, ShapeFinderResult> cacheResults
            = new HashMap<OneDIntArray, ShapeFinderResult>();

        TIntSet idxs = new TIntHashSet();

        long t0 = System.currentTimeMillis();
        
        double[] maxChordAvgDiff = new double[]{matchSingly(cacheResults, idxs)};
        
        long t1 = System.currentTimeMillis();
        
        System.out.println("matchSingly took " + ((t1 - t0)/1000) + " seconds");
        
        TIntIterator iter = idxs.iterator();

        List<ShapeFinderResult> results = new ArrayList<ShapeFinderResult>();

        while (iter.hasNext()) {

            int label = iter.next();

            t0 = System.currentTimeMillis();
            
            ShapeFinderResult sr = bfs(label, cacheResults, idxs, maxChordAvgDiff);

            t1 = System.currentTimeMillis();
        
            System.out.println("a single bfs took " + ((t1 - t0)/1000) + " seconds");
        
            if (sr != null) {
                results.add(sr);               
            }
        }

        ShapeFinderResult minCostR = null;
        double minCost = Double.MAX_VALUE;
       
        for (ShapeFinderResult r : results) {
                        
            if (r == null) {
                continue;
            }
            double sd = calcCost(r, maxChordAvgDiff[0]);
        
            if (sd < minCost) {
                minCost = sd;
                minCostR = r;
            }
        }

        return minCostR;
    }

    // run partial shape matcher on individual sets, store the labels and 
    // results and return the maximum chord diff avg 
    @SuppressWarnings({"unchecked"})
    private double matchSingly(Map<OneDIntArray, ShapeFinderResult> 
        cacheResults, TIntSet outIdxs) {

        double maxAvgDiffChord = Double.MIN_VALUE;

        ImageProcessor imageProcessor = new ImageProcessor();
        
        TIntObjectIterator<Set<PairInt>> iter = mapOfSets2.iterator();
        for (int i = 0; i < mapOfSets2.size(); ++i) {
            
            iter.advance();
            
            int label = iter.key();

            OneDIntArray keysI = new OneDIntArray(new int[]{label});
            PairIntArray boundsI;
            if (keyBoundsMap.containsKey(keysI)) {
                boundsI = keyBoundsMap.get(keysI);
            } else {
                Set<PairInt> set1 = iter.value();
                boundsI = imageProcessor.extractSmoothedOrderedBoundary(
                    new HashSet<PairInt>(set1), sigma, xMax2 + 1, yMax2 + 1);
                keyBoundsMap.put(keysI, boundsI);
            }
            
            PartialShapeMatcher matcher = new PartialShapeMatcher();
            matcher.overrideSamplingDistance(dp);
            matcher._overrideToThreshhold(0.2f);
            matcher.setToRemoveOutliers();
            //matcher.setToDebug();
            PartialShapeMatcher.Result r = matcher.match(bounds1, boundsI);
            if (r == null) {                
                continue;
            }
            
            outIdxs.add(label);

            double avgCD = r.getChordDiffSum() / (double) r.getNumberOfMatches();
            if (avgCD > maxAvgDiffChord) {
                maxAvgDiffChord = avgCD;
            }
            
            ShapeFinderResult sr = new ShapeFinderResult(r, bounds1, boundsI,
                keysI.a);
            
            cacheResults.put(keysI, sr);
        }

        return maxAvgDiffChord;
    }

    private ShapeFinderResult bfs(int srcIdx, Map<OneDIntArray, 
        ShapeFinderResult> cacheResults, TIntSet idxs, 
        double[] maxChordAvgDiff) {

        OneDIntArray srcKeys = new OneDIntArray(new int[]{srcIdx});
        if (!cacheResults.containsKey(srcKeys)) {
            return null;
        }

        int n = mapOfSets2.size();
        
        // key = label (a.k.a. segIdx), value = cost
        TIntDoubleMap dist = new TIntDoubleHashMap(n);
        TIntObjectIterator<Set<PairInt>> iter0 = mapOfSets2.iterator();
        for (int i = 0; i < n; ++i) {
            iter0.advance();;
            dist.put(iter0.key(), Double.MAX_VALUE);
        }
        
        // key = label (a.k.a. segIdx), value = state of visit
        // 0 = white, 1 = grey, 2 = visited
        TIntIntMap visited = new TIntIntHashMap(n);
        visited.put(srcIdx, 1);
        
        TIntObjectMap<ShapeFinderResult> results = new 
            TIntObjectHashMap<ShapeFinderResult>(n);
        
        results.put(srcIdx, cacheResults.get(srcKeys));
        dist.put(srcIdx, calcCost(results.get(srcIdx), maxChordAvgDiff[0]));
           
        // make a FIFO queue.
        ArrayDeque<Integer> queue = new ArrayDeque<Integer>();

        queue.add(Integer.valueOf(srcIdx));

        double mc = maxChordAvgDiff[0];

double sumDeltaT = 0;
int nDeltaT = 0;
long t = System.currentTimeMillis();
double sumDeltaT2 = 0;
int nDeltaT2 = 0;

        while (!queue.isEmpty()) {

            Integer uIndex = queue.pop();
            int uIdx = uIndex.intValue();
           
            TIntSet adjSet = adj2Map.get(uIdx);
            
            if (adjSet == null) {
                continue;
            }
            
            System.out.println("  => pop " + uIdx + " which has key=" +
                Arrays.toString(results.get(uIdx).labels2));

            TIntIterator iter = adjSet.iterator();
            while (iter.hasNext()) {
                int vIdx = iter.next();
           
                if (!idxs.contains(vIdx)) {
                    continue;
                }
                if (visited.containsKey(vIdx) && visited.get(vIdx) == 2) {
                    continue;
                }
                if (results.containsKey(vIdx) && 
                    (Arrays.binarySearch(results.get(vIdx).labels2, uIdx) > -1)
                    ) {
                    continue;
                }
           
                assert(results.containsKey(uIdx));
                
                ShapeFinderResult rV = results.get(vIdx);
                double sdV = dist.get(vIdx);
                if (Double.isNaN(sdV)) {
                    continue;
                }
                
                if (rV == null) {
                    OneDIntArray vKeys = new OneDIntArray(new int[]{vIdx});
                    rV = cacheResults.get(vKeys);
                }
               
                System.out.println("  => agg  u with v=" + vIdx 
                    + " which has key=" +
                    Arrays.toString(rV.labels2));
               
long t2 = System.currentTimeMillis();                
                
                ShapeFinderResult uPlusV = aggregateAndMatch(results.get(uIdx), rV);
                
long t2_2 = System.currentTimeMillis();
nDeltaT2++;
sumDeltaT2 += ((double)t2_2 - (double)t2) / 1000.;
t = t2_2;
                
                if (uPlusV == null) {
                    continue;
                }
                
                // store results in cache
                OneDIntArray uvKeys = new OneDIntArray(
                    Arrays.copyOf(uPlusV.labels2, uPlusV.labels2.length));
                cacheResults.put(uvKeys, uPlusV);
                
                double sdUPlusV = calcCost(uPlusV, maxChordAvgDiff[0]);
        
                if (Double.isNaN(sdUPlusV)) {
                    continue;
                }
                
                System.out.println("  => sdUPlusV=" + sdUPlusV + " sdV=" + sdV);
                    
                //visited[vIdx] = 1;
                if (sdUPlusV < sdV) {
                    
                    dist.put(vIdx, sdUPlusV);
                    results.put(vIdx, uPlusV);
                    queue.add(Integer.valueOf(vIdx));
                    
                    // update the maxes 
                    double avgCD = uPlusV.getChordDiffSum() / 
                        (double) uPlusV.getNumberOfMatches();
                    if (avgCD > mc) {
                        mc = avgCD;
                    }
                }
            }
        
            visited.put(uIdx, 2);
            
long t_2 = System.currentTimeMillis();
nDeltaT++;
sumDeltaT += ((double)t_2 - (double)t) / 1000.;
t = t_2;
        }
        
System.out.println("avg loop iteration = " + (sumDeltaT/nDeltaT) + " sec");

System.out.println("avg agg srch within loop iteration = " + (sumDeltaT2/nDeltaT2) + " sec");

        // re-do costs, and update max var
        maxChordAvgDiff[0] = mc;

        ShapeFinderResult minCostR = null;
        double minCost = Double.MAX_VALUE;
        
        TIntObjectIterator<ShapeFinderResult> iter2 = results.iterator();
        for (int i = 0; i < results.size(); ++i) {
            iter2.advance();
            ShapeFinderResult r = iter2.value();
           
            if (r == null) {
                continue;
            }
            double sd = calcCost(r, maxChordAvgDiff[0]);
       
            if (Double.isNaN(sd)) {
                continue;
            }
            
            if (sd < minCost) {
                minCost = sd;
                minCostR = r;
            }
        }
if (minCostR != null) {
ShapeFinderResult r = minCostR;
int nb1 = Math.round((float) r.bounds1.getN() / (float) dp);
float np = r.getNumberOfMatches();
int lGap = maxNumberOfGaps(r.bounds1, r)/dp;
float gCountComp = (float)lGap/(float)nb1;

System.out.println("min: cost=" + minCost + " count=" + np + 
" chordAvg=" + ((float) r.getChordDiffSum() / np) +
" gapAvg=" + gCountComp);
}
        return minCostR;
    }

    private double calcCost(ShapeFinderResult r, double maxAvgDiffChord) {
        
        int nb1 = Math.round((float) r.bounds1.getN() / (float) dp);

        float np = r.getNumberOfMatches();

        float countComp = 1.0F - (np / (float) nb1);
        double chordComp = ((float) r.getChordDiffSum() / np) / maxAvgDiffChord;

        int lGap = maxNumberOfGaps(r.bounds1, r)/dp;
        float gCountComp = (float)lGap/(float)nb1;

        // these should be squared
        double sd = chordComp*chordComp + countComp*countComp + gCountComp*gCountComp;
        
        if (Double.isNaN(sd)) {
            System.out.println("NaN: nb1=" + nb1 + " np=" + np
                + " chordsum=" + r.getChordDiffSum() + 
                " countComp=" + countComp + " gCountComp=" + gCountComp);
        }
        
        r.contextCost = sd;
        
        return sd;
    }

    public static class ShapeFinderResult extends Result {

        protected final PairIntArray bounds1;
        protected final PairIntArray bounds2;
        
        /**
         * cost calculated in the result's context of max chord diff avg and
         * max number of matchable.
         * Note that this value may be the square sum and not
         * the square root of the square sums (that's the case for
         * use throughout ShapeFinder2.java)
         */
        protected double contextCost = Double.MAX_VALUE;
        
        /**
         * a list of sorted labels, bounds2 is the boundary of the combined 
         * regions in labels2
         */
        protected final int[] labels2;

        public ShapeFinderResult(Result r, PairIntArray bounds1,
            PairIntArray bounds2, int[] keys) {

            super(r.n1, r.n2, r.getOriginalN1());

            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
            this.chordDiffSum = r.chordDiffSum;
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
        
        public double getContextCost() {
            return contextCost;
        }

        @Override
        public String toString() {
            String str = super.toString();
            str = str + " contextCost=" + contextCost;
            return str;
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
    
    private int maxNumberOfGaps(PairIntArray bounds,
        PartialShapeMatcher.Result r) {

        TIntSet mIdxs = new TIntHashSet(r.getNumberOfMatches());
        for (int i = 0; i < r.getNumberOfMatches(); ++i) {
            mIdxs.add(r.getIdx1(i));
        }

        int maxGapStartIdx = -1;
        int maxGap = 0;
        int cStartIdx = -1;
        int cGap = 0;

        // handling for startIdx of 0 to check for wraparound
        // of gap at end of block
        int gap0 = 0;

        for (int i = 0; i < bounds.getN(); ++i) {
            if (!mIdxs.contains(i)) {
                // is a gap
                if (cStartIdx == -1) {
                    cStartIdx = i;
                }
                cGap++;
                if (i == (bounds.getN() - 1)) {
                    if (gap0 > 0) {
                        // 0 1 2 3 4 5
                        // g g     g g
                        // gap0=2
                        // cGap=2 cStartIdx=4
                        if (cStartIdx > (gap0 - 1)) {
                            gap0 += cGap;
                        }
                    }
                    if (cGap > maxGap) {
                        maxGap = cGap;
                        maxGapStartIdx = cStartIdx;
                    }
                    if (gap0 > maxGap) {
                        maxGap = gap0;
                        maxGapStartIdx = 0;
                    }
                }
            } else {
                // is not a gap
                if (cStartIdx > -1) {
                    if (cGap > maxGap) {
                        maxGap = cGap;
                        maxGapStartIdx = cStartIdx;
                    }
                    if (cStartIdx == 0) {
                        gap0 = cGap;
                    }
                    cStartIdx = -1;
                    cGap = 0;
                }
            }
        }

        return maxGap;
    }

}
