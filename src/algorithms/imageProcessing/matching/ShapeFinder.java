package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.imageProcessing.transform.EuclideanEvaluator;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.util.PairInt;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * uses PartialShapeMatcher and search patterns to
 * find the best fitting match between a
 * group of adjacent segmented cells
 * to a template shape.  Any additional information
 * that helps limit the search should be done before
 * this stage and included or excluded in the
 * adjacency map.  Note that this class could change
 * in the future to use weights in the adjacency map,
 * but for now it's a binary relationship map.
 *
 * @author nichole
 */
public class ShapeFinder {
    
    private double maxDiffChordSum = Float.MAX_VALUE;
    
    private double maxDistTransformSum = Float.MAX_VALUE;

    private final float areaFactor = 2.f;

    // partial shape matcher sampling distance
    private final int dp = 2;

    //DEBUG
    static TIntSet expectedIndexes = null;
    static {
        expectedIndexes = new TIntHashSet();
        expectedIndexes.add(18);
        expectedIndexes.add(30);
        expectedIndexes.add(33);
        expectedIndexes.add(53);
        expectedIndexes.add(54);
        expectedIndexes.add(64);
        expectedIndexes.add(93);
        expectedIndexes.add(97);
        expectedIndexes.add(124);
    }

    /**
     * NOT READY FOR USE.
     * uses PartialShapeMatcher and search patterns to
     * find the best fitting match between a
     * group of adjacent segmented cells
     * to a template shape.  Any additional information
     * that helps limit the search should be done before
     * this stage and included or excluded in the
     * adjacency map.  Note that this class could change
     * in the future to use weights in the adjacency map,
     * but for now it's a binary relationship map.
     * NOTE also that methods to return the top k results
     * or the equivalent best within a tolerance might be
     * made in the future.
     *
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @param orderedBoundaries
     * @param pointsList
     * @param adjacencyMap
     * @return
     */
    public Result findMatchingCells(
        List<PairIntArray> orderedBoundaries,
        List<Set<PairInt>> pointsList,
        TIntObjectMap<TIntSet> adjacencyMap, PairIntArray template) {

        /*
        step (1) O(N_cells) + O(N_cells_in_bin * log_2(N_cells_in_bin)):
           the bin dimensions will be max of template dimensions.
           actually, should probably be the diagonal for max dist
           between any two points but that could be done more
           accurately.
           - 2D bins, that is x and y.
           - sort the cells by x centroid and ties by y centroid
        step (2) search pattern, 3 are presented.
            The bigger picture is a sliding window of the 2D bin
            where the core of the search and aggregation is done
            within the 2D bin, but the aggregation can continue
            into the next bin to the right, next bin above,
            or above and to the right.
            The search of the cells could proceed in 3 different
            ways.  will implement all 3 and keep the fastest of
            the robust of them.
            (2a)
               try all combinations in a bin, then all combinations
               that can be aggregated on to them that are in
               the adjacent "next" bins.
               this method is thorough brute force, but is
               combinatorial within bounds of area 3*binsize (or smaller
               if adj map is sparse) in runtime,
               so a faster but still accurate method would be better.
            (2b and c)
               over the whole image, use the partial shape matcher
               with cell boundaries that have area less than
               two or so times the dimension of the template.
               store the top k of those matches rather than just
               the best matching presumably.  this is the first
               step of 2 different searches, 2b and 2c.
               (2b)
                   visit each bin determining whether a match
                   to the whole template is findable within the
                   top k search results of each, aggregating
                   adjacent where it is matching or where it is in between
                   such matches (that is, embedded) ...
                   this is solving the problem on a cell level and combining
                   results to make a better solution.
                   difficult to see this as robust, but
                   if significant fraction of matching cells were
                   within a small top k of 5 or less of individul cell shape matching,
                   this might be feasible and it would be polynomial
                   search, bounded by the segmentation fragmentation
                   template size, and the choice of k.
               (2c)
                   for each bin, choose as the start of a search
                   the best fitting cell, then proceed with a path
                   like aggregation search (dijkstra's, for example).
                   more robust would be to let each cell with
                   a reasonable matching cost, be a start seed in
                   a separate path search.
                   if the later is done, need to keep track of
                   already searched groups across path solutions
                   and reuse those result.
                   Dijkstra's runtime for a single path search is O(E lg V).
                   **NOTE that this pattern of a dijkstra's search for each
                   feasible start could be done over the whole image's cells
                   instead of a bin's cells and could presumably be done
                   over a smaller number of starts than all cells by sorting
                   the start list by the partial match cost (score)
                   and reuse of the aggregated cell partial results.

               there are then 4 search patterns to implement here.
        */

        if (true) {
            return wideFWSearch(orderedBoundaries, pointsList, adjacencyMap,
                template);
        }

        throw new UnsupportedOperationException("not yet implemented");
    }

    private Result wideFWSearch(List<PairIntArray> orderedBoundaries,
        List<Set<PairInt>> pointsList,
        TIntObjectMap<TIntSet> adjacencyMap, PairIntArray template) {

        /*
        NOTE: if this method is to be used, it needs some improvements.
        -- possibly needs improvements to the evaluation.
        -- current test shows that it might need additional information used
           to filter the adjacency map given to this class.
           -- the main goal is to match the shape independent of color
              changes and texture so it could match the outline of
              the gingerbread man in different lighting and poses.
              also to have at some point, ability to match mountain
              silhouettes such as half dome, for example.

              the android statues 02 test might need additional information
              that distinguishes the gbman from falsely matched regions.
              considering if can identify foreground or background regions
              in the segmentation to filter out of the adjacency maps,
              for example, (or, conversely, could consider a separate match
              for those to find the inverse profile too...)
        */

        int n = orderedBoundaries.size();

        List<Result> outputSortedResults = new ArrayList<Result>(n);
        List<Integer> outputSortedIndexes = new ArrayList<Integer>(n);

        Map<VeryLongBitString, PairIntArray> aggregatedBoundaries =
            new HashMap<VeryLongBitString, PairIntArray>();

        Map<VeryLongBitString, Result> aggregatedResultMap =
            new HashMap<VeryLongBitString, Result>();

        // this can be derived from values in aggregatedResultMap,
        // but is currently kept to make debugging easier
        Map<VeryLongBitString, Double> aggregatedCostMap =
            new HashMap<VeryLongBitString, Double>();
        
        // maxDiffChordSum is calculated in here
        matchAndOrderByIncrCost(orderedBoundaries, template,
            outputSortedResults, outputSortedIndexes,
            aggregatedBoundaries, aggregatedResultMap,
            aggregatedCostMap);

        assert(aggregatedBoundaries.size() == aggregatedResultMap.size());
        assert(aggregatedCostMap.size() == aggregatedResultMap.size());

        List<PairInt> orderedBoundaryCentroids = calculateCentroids(
            orderedBoundaries);

        double minCost = Double.MAX_VALUE;
        VeryLongBitString minBitString = null;

        for (int i = 0; i < outputSortedIndexes.size(); ++i) {
        //for (int i = 11; i < 12; ++i) {

            Integer index = outputSortedIndexes.get(i);

            if (index.intValue() != 54) {
                continue;
            }
            //184,160  debugging
            PairInt pCen = orderedBoundaryCentroids.get(index);
            System.out.println("pCen=" + pCen + " idx=" + index + " i=" + i);


            // the in-out variables store reusable calculations and
            // also the resulting cost of this best search result
            VeryLongBitString bitString = minCostAggregationFW(
                orderedBoundaries, pointsList, adjacencyMap, template, index,
                orderedBoundaryCentroids,
                aggregatedBoundaries,
                aggregatedResultMap, aggregatedCostMap);

            if (bitString == null) {
                continue;
            }

            double cost = aggregatedCostMap.get(bitString).doubleValue();
            if (cost < minCost) {
                minCost = cost;
                minBitString = bitString;
            }
        }

        if (minBitString == null) {
            return null;
        }

        Result r = aggregatedResultMap.get(minBitString);
        r.setData(aggregatedBoundaries.get(minBitString));

        return r;
    }

    /**
     * has the side effect of maxDiffChordSum of maxDistTransformSum
     * 
     * @param orderedBoundaries
     * @param template
     * @param outputSortedResults - output variable
     * @param outputSortedIndexes - output variable
     * @param aggregatedBoundaries - output variable
     * @param aggregatedResultMap - output variable
     * @param aggregatedCostMap  - output variable
     */
    private void matchAndOrderByIncrCost(List<PairIntArray> orderedBoundaries,
        PairIntArray template, List<Result> outputSortedResults,
        List<Integer> outputSortedIndexes,
        Map<VeryLongBitString, PairIntArray> aggregatedBoundaries,
        Map<VeryLongBitString, Result> aggregatedResultMap,
        Map<VeryLongBitString, Double> aggregatedCostMap) {

        // key=result, value=index of orderedBoundaries item
        TObjectIntMap<Result> cellMatchResults = matchIndividually(
            orderedBoundaries, template);

        int n = cellMatchResults.size();
        TObjectIntIterator<Result> iter = cellMatchResults.iterator();
        
        double maxChord = Double.MIN_VALUE;
        double maxDist = Double.MIN_VALUE;
        for (int i = 0; i < n; ++i) {
            iter.advance();
            Result r = iter.key();
            int idx = iter.value();
            double d = r.getChordDiffSum();
            if (d > maxChord) {
                maxChord = d;
            }
            d = calcTransformationDistanceSum(r, orderedBoundaries.get(idx),
                template);
            if (d > maxDist) {
                maxDist = d;
            }
        }
        maxDiffChordSum = maxChord;
        maxDistTransformSum = maxDist;

        int nB = orderedBoundaries.size();
        
        int nT1 = template.getN();

        float[] costs = new float[n];
        int[] indexes = new int[n];
        Result[] results = new Result[n];

        iter = cellMatchResults.iterator();
        for (int i = 0; i < n; ++i) {
            iter.advance();

            Result r = iter.key();
            int idx = iter.value();
            
            double dist = calcTransformationDistanceSum(r, 
                orderedBoundaries.get(idx),
                template)/maxDistTransformSum;

            // calculating salukwdze dist^2 to use same reference, nT1
            int nI = r.getNumberOfMatches();
            float f = 1.f - ((float)nI/(float)nT1);
            double d = r.getChordDiffSum()/maxChord;
            float s = (float)(f * f + d * d + dist * dist);

            costs[i] = s;
            indexes[i] = i;
            results[i] = r;

            VeryLongBitString bs = new VeryLongBitString(nB);
            bs.setBit(idx);

            PairIntArray put = aggregatedBoundaries.put(bs, 
                orderedBoundaries.get(idx));
            assert(put == null);
            aggregatedResultMap.put(bs, r);
            aggregatedCostMap.put(bs, Double.valueOf(s));
        }

        QuickSort.sortBy1stArg(costs, indexes);

        for (int i = 0; i < n; ++i) {
            int index = indexes[i];
            Result r = results[index];

            outputSortedResults.add(r);

            int listIndex = cellMatchResults.get(r);
            outputSortedIndexes.add(Integer.valueOf(listIndex));
        }

    }

    private TObjectIntMap<Result> matchIndividually(
        List<PairIntArray> orderedBoundaries, PairIntArray template) {

        // key=result, value=index of orderedBoundaries item
        TObjectIntMap<Result> output = new TObjectIntHashMap<Result>();

        int[] dimensionsT = calcDimensions(template);

        int areaT = (int)Math.round(Math.sqrt(dimensionsT[0] * dimensionsT[0] +
            dimensionsT[1] * dimensionsT[1]));

        for (int i = 0; i < orderedBoundaries.size(); ++i) {

            PairIntArray p = orderedBoundaries.get(i);

            if (p.getN() < 6) {
                continue;
            }

            int[] dimensions = calcDimensions(template);

            int area = (int)Math.round(Math.sqrt(dimensions[0] * dimensions[0] +
                dimensions[1] * dimensions[1]));

            if (area > (areaFactor * areaT)) {
                // NOTE: this assumes the segmentation of the object from
                // a large background or foreground is complete separation.
                // If that isn't true, then a search pattern which is
                // more like a skyline partial matching pattern must be
                // used instead.
                // TODO: need a way to recognize that or at least mention
                // it in the documentation. the skyline search pattern will
                // be implemented soon.
                continue;
            }

            if (area < (areaT/areaFactor)) {
                //TODO: may need to revise this limit
                continue;
            }

            PartialShapeMatcher matcher = new PartialShapeMatcher();
            matcher.overrideSamplingDistance(dp);
            Result r = matcher.match(p, template);
            if (r == null) {
                continue;
            }

            output.put(r, i);
        }

        return output;
    }

    private int[] calcDimensions(PairIntArray a) {

        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;

        for (int i = 0; i < a.getN(); ++i) {
            int x = a.getX(i);
            int y = a.getY(i);
            if (x < minX) {
                minX = x;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
        }

        return new int[]{maxX - minX + 1, maxY - minY + 1};
    }

    private List<PairInt> calculateCentroids(List<PairIntArray> orderedBoundaries) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<PairInt> centroids = new ArrayList<PairInt>();

        for (int i = 0; i < orderedBoundaries.size(); ++i) {

            double[] xyCen = curveHelper.calculateXYCentroids(
                orderedBoundaries.get(i));

            PairInt xy = new PairInt((int)Math.round(xyCen[0]),
                (int)Math.round(xyCen[1]));

            centroids.add(xy);
        }

        return centroids;
    }

    /**
     * a Floyd Warshall search to find min-cost aggregation of segmented
     * cells within a limited distance of adjacency.
     *
     * @param orderedBoundaries
     * @param adjacencyMap
     * @param template
     * @param index
     * @param orderedBoundaryCentroids
     * @param maxDiffChordSum
     * @param aggregatedBoundaries - in out variable
     * @param aggregatedResultMap - in out variable
     * @param aggregatedCostMap - in out variable
     * @return
     */
    private VeryLongBitString minCostAggregationFW(
        List<PairIntArray> orderedBoundaries,
        List<Set<PairInt>> pointsList,
        TIntObjectMap<TIntSet> adjacencyMap,
        PairIntArray template, int index,
        List<PairInt> orderedBoundaryCentroids,
        Map<VeryLongBitString, PairIntArray> aggregatedBoundaries,
        Map<VeryLongBitString, Result> aggregatedResultMap,
        Map<VeryLongBitString, Double> aggregatedCostMap) {

        int[] dimensionsT = calcDimensions(template);

        int areaT = (int)Math.round(Math.sqrt(dimensionsT[0] * dimensionsT[0] +
            dimensionsT[1] * dimensionsT[1]));

        int n = orderedBoundaries.size();

        // NOTE: in the global image wide method, this will
        // be replaced with int[][] dist and int[][] prev.
        // This method for a specific index is asserting the logic in detail first.
        TObjectDoubleMap<PairInt> distMap = new TObjectDoubleHashMap<PairInt>();
        Map<PairInt, VeryLongBitString> indexesMap = new HashMap<PairInt,
            VeryLongBitString>();

        int nB = orderedBoundaries.size();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        TIntList indexes = new TIntArrayList();

        PairInt indexXY = orderedBoundaryCentroids.get(index);

        int nT1 = template.getN();

        float factor = 2.f;
        float maxDist = factor * (float)Math.max(dimensionsT[0], dimensionsT[1]);

        // ------- initialize the local maps and store partial results in
        //            output vars too

        // insert items within distance of index centroid and areaFactor times
        // the dimensionsT
        for (int i = 0; i < n; ++i) {

            if (orderedBoundaries.get(i).getN() < 6) {
                continue;
            }

            PairInt xy = orderedBoundaryCentroids.get(i);
            double dist = distance(indexXY, xy);
            if (dist > maxDist) {
                continue;
            }

            int[] dimensions = calcDimensions(orderedBoundaries.get(i));
            int area = (int)Math.round(Math.sqrt(dimensions[0] * dimensions[0] +
            dimensions[1] * dimensions[1]));
            if (area > areaFactor*areaT || (dimensions[0] > dimensionsT[0]) ||
                (dimensions[1] > dimensionsT[1])) {
                continue;
            }

            indexes.add(i);
        }

        System.out.println("indexes=" + indexes);

        for (int i = 0; i < indexes.size(); ++i) {

            //TODO: when convert to global method,
            // consider re-using the aggregate map information
            // IF it is present for key and IF all members
            // of the bitstring for the key are present in indexes
            // (the check for presence should help limit the result to adjacent
            // within a size constraint).

            int idx1 = indexes.get(i);

            if (!adjacencyMap.containsKey(idx1)) {
                continue;
            }

            VeryLongBitString bs1 = new VeryLongBitString(nB);
            bs1.setBit(idx1);

            Double cost1 = aggregatedCostMap.get(bs1);

            PairIntArray boundary1 = orderedBoundaries.get(idx1);

            for (int j = 0; j < indexes.size(); ++j) {

                int idx2 = indexes.get(j);

                if (!adjacencyMap.containsKey(idx1) ||
                    !adjacencyMap.get(idx1).contains(idx2)) {
                    continue;
                }

                if (idx1 == idx2) {
                    PairInt key = new PairInt(idx1, idx1);
                    distMap.put(key, cost1);
                    indexesMap.put(key, bs1);
                    continue;
                }

                VeryLongBitString bs2 = bs1.copy();
                bs2.setBit(idx2);

                PairInt key2;
                key2 = new PairInt(idx1, idx2);
                indexesMap.put(key2, bs2);

                PairIntArray boundary12 = aggregatedBoundaries.get(bs2);
                Double cost12 = aggregatedCostMap.get(bs2);

                if (cost12 == null) {

                    if (boundary12 == null) {
                        boundary12 = mergeAdjacentOrderedBorders(bs2,
                            orderedBoundaries, pointsList);
                    }

                    cost12 = calcAndStoreMatchCost(boundary12, template, bs2,
                        aggregatedBoundaries, aggregatedResultMap, 
                        aggregatedCostMap);
                }

                if (cost12 != null) {

                    distMap.put(key2, cost12);

                    //DEBUG
                    if (expectedIndexes.contains(idx1) ||
                        expectedIndexes.contains(idx2)) {
                        System.out.println("contains an expected: " +
                            " key2=" + key2 + " cost=" + cost12);
                    }
                }
            }
        }
        // ---- end initializing local maps

        PairInt minCostIdx = null;
        double minCost = Double.MAX_VALUE;

        /*
        s0 = costIJ;
        s1 = costIK + costKJ;
        */

        for (int i0 = 0; i0 < indexes.size(); ++i0) {
            int k = indexes.get(i0);
            for (int i1 = 0; i1 < indexes.size(); ++i1) {
                int i = indexes.get(i1);
                for (int i2 = 0; i2 < indexes.size(); ++i2) {
                    int j = indexes.get(i2);

                    // set these by adjacency.  false if not adjacent
                    boolean tIJ = true;
                    boolean tIK = true;
                    boolean tKJ = true;

                    if (!adjacencyMap.containsKey(i) ||
                        !adjacencyMap.get(i).contains(j)) {
                        tIJ = false;
                    }
                    if (!adjacencyMap.containsKey(i) ||
                        !adjacencyMap.get(i).contains(k)) {
                        tIK = false;
                    }
                    if (!adjacencyMap.containsKey(k) ||
                        !adjacencyMap.get(k).contains(j)) {
                        tKJ = false;
                    }

                    if (!tIJ && (!tIK || !tKJ)) {
                        continue;
                    }

                    boolean setPrev = true;
                    if (i == j) {
                        setPrev = false;
                    }

                    PairInt keyIJ;
                    keyIJ = new PairInt(i, j);

                    if (i == j) {
                        // this was set in initialization block above
                        if (tIJ) {
                            assert(distMap.containsKey(keyIJ));
                            assert(indexesMap.containsKey(keyIJ));
                            if (distMap.get(keyIJ) < minCost) {
                                minCost = distMap.get(keyIJ);
                                minCostIdx = keyIJ;
                                assert(indexesMap.get(minCostIdx) != null);
                            }
                        }
                        continue;
                    }

                    PairInt keyIK, keyKJ;
                    keyIK = new PairInt(i, k);
                    keyKJ = new PairInt(j, k);

                    Double s0 = null;
                    Double s1 = null;
                    //s0 = costIJ;
                    //s1 = costIK + costKJ;

                    VeryLongBitString bs0 = null;
                    VeryLongBitString bs1 = null;

                    if (tIJ) {
                        if (distMap.containsKey(keyIJ)) {
                            s0 = distMap.get(keyIJ);
                            bs0 = indexesMap.get(keyIJ);
                        } else {
                            PairInt p = new PairInt(keyIJ.getX(), keyIJ.getX());
                            if (indexesMap.get(p) == null) {
                                bs0 = new VeryLongBitString(nB);
                                bs0.setBit(keyIJ.getX());
                                //indexesMap.put(p, bs0.copy());
                            } else {
                                bs0 = indexesMap.get(p).copy();
                            }
                            bs0.setBit(keyIJ.getY());
                            if (!aggregatedCostMap.containsKey(bs0)) {
                                PairIntArray boundary = orderedBoundaries.get(
                                    keyIJ.getX());
                                if (keyIJ.getX() != keyIJ.getY()) {
                                    PerimeterFinder2 perFinder2 = new
                                       PerimeterFinder2();
                                    boundary = mergeAdjacentOrderedBorders(bs0,
                                        orderedBoundaries, pointsList);
                                }
                                Double s = calcAndStoreMatchCost(boundary, template, bs0,
                                    aggregatedBoundaries,
                                    aggregatedResultMap, aggregatedCostMap);
                                if (s == null) {
                                    tIJ = false;
                                }
                            }
                            if (tIJ) {
                                s0 = aggregatedCostMap.get(bs0);
                            }
                        }
                    }

                    if (tIK && tKJ) {
                        //s1 = costIK + costKJ;
                        VeryLongBitString bsIK = indexesMap.get(keyIK);
                        VeryLongBitString bsKJ = indexesMap.get(keyKJ);

                        if (bsIK == null) {
                            bsIK = new VeryLongBitString(nB);
                            bsIK.setBit(keyIK.getX());
                            bsIK.setBit(keyIK.getY());
                            //indexesMap.put(keyIK, bsIK);
                        }
                        if (bsKJ == null) {
                            bsKJ = new VeryLongBitString(nB);
                            bsKJ.setBit(keyKJ.getX());
                            bsKJ.setBit(keyKJ.getY());
                            //indexesMap.put(keyKJ, bsKJ);
                        }

                        bs1 = bsIK.or(bsKJ);
                        PairIntArray bIKKJ = aggregatedBoundaries.get(bs1);
                        s1 = aggregatedCostMap.get(bs1);
                        if (s1 == null) {
                            bIKKJ = mergeAdjacentOrderedBorders(bs1,
                                orderedBoundaries, pointsList);
                            Result r = null;
                            if (bIKKJ.getN() > 6) {
                                PartialShapeMatcher matcher = new PartialShapeMatcher();
                                matcher.overrideSamplingDistance(dp);
                                r = matcher.match(bIKKJ, template);
                            }

                            if (r == null) {
                                tIK = false;
                                tKJ = false;
                            } else {
                                
                                double dist = calcTransformationDistanceSum(r, 
                                    bIKKJ, template)/maxDistTransformSum;
                                
                                int nI = r.getNumberOfMatches();
                                float f = 1.f - ((float)nI/(float)nT1);
                                double d = r.getChordDiffSum()/maxDiffChordSum;
                                s1 = Double.valueOf(f * f + d * d + dist * dist);

                                aggregatedBoundaries.put(bs1, bIKKJ);
                                aggregatedCostMap.put(bs1, s1);
                                aggregatedResultMap.put(bs1, r);
                            }
                        }
                    }

                    if (tIJ && s0 != null && s1 != null &&
                        ((s0.doubleValue() <= s1.doubleValue()) || (!tIK || !tKJ))) {
                        assert(tIJ);
                        
if ((distMap.containsKey(keyIJ) && (distMap.get(keyIJ) > s0))
|| !distMap.containsKey(keyIJ)) {
                        distMap.put(keyIJ, s0);
                        indexesMap.put(keyIJ, bs0);
}


                        //DEBUG
                        if (expectedIndexes.contains(i) ||
                            expectedIndexes.contains(j)) {
                            if (keyIJ.getX() <= keyIJ.getY()) {
                                System.out.println("contains an expected: " +
                                    " keyIJ=" + keyIJ.getX() + "," 
                                    + keyIJ.getY() + " cost=" + s0);
                            } else {
                                System.out.println("contains an expected: " +
                                    " keyIJ=" + keyIJ.getY() + "," 
                                    + keyIJ.getX() + " *cost=" + s0);
                            }
                        }

                    } else if (tIK && tKJ && (s1 != null)) {
                        assert(tIK);
                        assert(tKJ);
                        assert(bs1 != null);
                        
if ((distMap.containsKey(keyIJ) && (distMap.get(keyIJ) > s1))
|| !distMap.containsKey(keyIJ)) {
                        distMap.put(keyIJ, s1);
                        indexesMap.put(keyIJ, bs1);
}
                        //DEBUG
                        if (expectedIndexes.contains(i) ||
                            expectedIndexes.contains(j)
                            || expectedIndexes.contains(j)) {
                            if (keyIJ.getX() <= keyIJ.getY()) {
                                System.out.println("contains an expected: " +
                                    " keyIJ=" + keyIJ.getX() + "," 
                                    + keyIJ.getY() + " ik+kh, keyIJ cost=" + s1);
                            } else {
                                System.out.println("contains an expected: " +
                                    " keyIJ=" + keyIJ.getY() + "," 
                                    + keyIJ.getX() + " *ik+kh, keyIJcost=" + s1);
                            }
                        }

                    } else {
                        continue;
                    }

 if (indexesMap.containsKey(keyIJ)) {
                    assert(indexesMap.get(keyIJ) != null);
                    if (distMap.get(keyIJ) < minCost) {
                        minCost = distMap.get(keyIJ);
                        minCostIdx = keyIJ;
                    }
}
                }
            }
        }

        if (minCostIdx == null) {
            return null;
        }

        {
            // expected:
            VeryLongBitString tBS = new VeryLongBitString(nB);
            tBS.setBit(18); tBS.setBit(30); tBS.setBit(33); tBS.setBit(53);
            tBS.setBit(54); tBS.setBit(64); tBS.setBit(93); tBS.setBit(97);
            tBS.setBit(124);
            Double cost = aggregatedCostMap.get(tBS);
            if (cost == null) {
                PairIntArray b =  mergeAdjacentOrderedBorders(tBS,
                    orderedBoundaries, pointsList);
                cost = calcAndStoreMatchCost(b, template, tBS,
                    aggregatedBoundaries,
                    aggregatedResultMap, aggregatedCostMap);
            }
            System.out.println("expected true answer cost=" + cost);
            for (Entry<PairInt, VeryLongBitString> entry : indexesMap.entrySet()) {
                if (entry.getValue().equals(tBS)) {
                    System.out.println("Expected found in distMap is " +
                        distMap.get(entry.getKey()));
                }
            }

            //found
            int[] minidxs = indexesMap.get(minCostIdx).getSetBits();
            for (int idx : minidxs) {
                PairIntArray pts = orderedBoundaries.get(idx);
                System.out.println("CEN=" +
                    Arrays.toString
                    (curveHelper.calculateXYCentroids(pts))
                    + " n=" + pts.getN()
                );
            }
        }


        VeryLongBitString minBitString = indexesMap.get(minCostIdx);

        System.out.println("minCostIdx=" + minCostIdx +
            " minCost=" + minCost + " bs=" +
            Arrays.toString(minBitString.getSetBits()));

        return minBitString;
    }

    private double distance(PairInt xy1, PairInt xy2) {

        int diffX = xy1.getX() - xy2.getX();
        int diffY = xy1.getY() - xy2.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }

    private Double calcAndStoreMatchCost(PairIntArray boundary,
        PairIntArray template, VeryLongBitString bs,
        Map<VeryLongBitString, PairIntArray> aggregatedBoundaries,
        Map<VeryLongBitString, Result> aggregatedResultMap,
        Map<VeryLongBitString, Double> aggregatedCostMap) {

        aggregatedBoundaries.put(bs, boundary);

        if (boundary.getN() < 6) {
            return null;
        }

        PartialShapeMatcher matcher = new PartialShapeMatcher();
        matcher.overrideSamplingDistance(dp);
        Result result12 = matcher.match(boundary, template);
        if (result12 == null) {
            return null;
        }
        aggregatedResultMap.put(bs, result12);

        double dist = calcTransformationDistanceSum(result12, 
            boundary, template)/maxDistTransformSum;
                       
        int nT1 = template.getN();
        // NOTE: this may need to be revised.  using the max diff chord
        // sum from only the single segmented cell matches as
        // the normalization for the Salukwdze distance^2.
        // May need to re-examine the bounds values and adjust the
        // normalizations.
        int nI = result12.getNumberOfMatches();
        float f = 1.f - ((float) nI / (float) nT1);
        double d = result12.getChordDiffSum() / maxDiffChordSum;
        double s = (f * f + d * d + dist * dist);
        Double cost12 = Double.valueOf(s);

        aggregatedCostMap.put(bs, cost12);

        return cost12;
    }

    private PairIntArray mergeAdjacentOrderedBorders(VeryLongBitString bs,
        List<PairIntArray> orderedBorders, List<Set<PairInt>> pointsList) {

        int[] setBits = bs.getSetBits();
        Set<PairInt> allPoints = new HashSet<PairInt>();
        for (int idx : setBits) {
            Set<PairInt> set = pointsList.get(idx);
            allPoints.addAll(set);
        }

        System.out.println("allPoints.size=" + allPoints.size()
         + " setBits=" + Arrays.toString(setBits));

        PerimeterFinder2 f2 = new PerimeterFinder2();

        return f2.extractOrderedBorder(allPoints);
    }

    private double calcTransformationDistanceSum(Result r, PairIntArray p,
        PairIntArray q) {
        
        if (r.getTransformationParameters() == null) {
            return Double.MAX_VALUE;
        }
        
        PairIntArray left = new PairIntArray(r.getNumberOfMatches());
        PairIntArray right = new PairIntArray(r.getNumberOfMatches());
        for (int i = 0; i < r.getNumberOfMatches(); ++i) {
            int idx = r.getIdx1(i);
            left.add(p.getX(idx), p.getY(idx));
            idx = r.getIdx2(i);
            right.add(q.getX(idx), q.getY(idx));
        }
        
        int tolerance = 5;
        
        EuclideanEvaluator evaluator = new EuclideanEvaluator();
        EuclideanTransformationFit fit = evaluator.evaluate(left,
            right, r.getTransformationParameters(), tolerance);
        
        List<Double> distances = fit.getErrors();
        
        double sum = 0;
        for (Double d : distances) {
            sum += d;
        }
        
        return sum;
    }
}
