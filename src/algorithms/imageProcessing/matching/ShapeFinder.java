package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.util.PairInt;
import algorithms.util.VeryLargeNumber;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;

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
    
    private final float areaFactor = 2.f;
    
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
     * @param adjacencyMap
     * @return 
     */
    public Result findMatchingCells(
        List<PairIntArray> orderedBoundaries,
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
            return wideDijkstraSearch(orderedBoundaries, adjacencyMap,
                template);
        }
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private Result wideDijkstraSearch(List<PairIntArray> orderedBoundaries, 
        TIntObjectMap<TIntSet> adjacencyMap, PairIntArray template) {
        
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
        
        double[] maxDiffChordSum = new double[1];
        
        matchAndOrderByIncrCost(orderedBoundaries, template, 
            outputSortedResults, outputSortedIndexes,
            aggregatedBoundaries, aggregatedResultMap, 
            aggregatedCostMap, maxDiffChordSum);
        
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
                orderedBoundaries, adjacencyMap, template, index,
                orderedBoundaryCentroids, maxDiffChordSum,
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
     * 
     * @param orderedBoundaries
     * @param template
     * @param outputSortedResults - output variable
     * @param outputSortedIndexes - output variable
     * @param aggregatedBoundaries - output variable
     * @param aggregatedResultMap - output variable
     * @param aggregatedCostMap  - output variable
     * @param maxDiffChordSum - output variable
     */
    private void matchAndOrderByIncrCost(List<PairIntArray> orderedBoundaries, 
        PairIntArray template, List<Result> outputSortedResults,
        List<Integer> outputSortedIndexes,
        Map<VeryLongBitString, PairIntArray> aggregatedBoundaries,
        Map<VeryLongBitString, Result> aggregatedResultMap,
        Map<VeryLongBitString, Double> aggregatedCostMap,
        double[] maxDiffChordSum) {

        TObjectIntMap<Result> cellMatchResults = matchIndividually(
            orderedBoundaries, template);
        
        int n = cellMatchResults.size();
        TObjectIntIterator<Result> iter = cellMatchResults.iterator();
        
        double maxChord = Double.MIN_VALUE;
        for (int i = 0; i < n; ++i) {
            iter.advance();
            Result r = iter.key();
            double d = r.getChordDiffSum();
            if (d > maxChord) {
                maxChord = d;
            }
        }
        maxDiffChordSum[0] = maxChord;
        
        int nT1 = template.getN();
        
        float[] costs = new float[n];
        int[] indexes = new int[n];
        Result[] results = new Result[n];
        
        iter = cellMatchResults.iterator();
        for (int i = 0; i < n; ++i) {
            iter.advance();            
            Result r = iter.key();
            
            // calculating salukwdze dist to use same reference, nT1
            int nI = r.getNumberOfMatches();
            float f = 1.f - ((float)nI/(float)nT1);            
            double d = r.getChordDiffSum()/maxChord;
            float s = (float)Math.sqrt(f * f + d * d);
            
            costs[i] = s;
            indexes[i] = i;
            results[i] = r;
            
            VeryLongBitString bs = new VeryLongBitString(n);
            bs.setBit(i);
            
            PairIntArray put = aggregatedBoundaries.put(bs, orderedBoundaries.get(i));
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
     * a Dijkstra's search to find min-cost aggregation of segmented
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
        List<PairIntArray> orderedBoundaries, TIntObjectMap<TIntSet> adjacencyMap, 
        PairIntArray template, int index, 
        List<PairInt> orderedBoundaryCentroids,
        double[] maxDiffChordSum,
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
        TObjectIntMap<PairInt> prevMap = new TObjectIntHashMap<PairInt>();
        Map<PairInt, VeryLongBitString> indexesMap = new HashMap<PairInt,
            VeryLongBitString>();
                
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        PerimeterFinder2 pFinder2 = new PerimeterFinder2();
        
        TIntList indexes = new TIntArrayList();
        
        PairInt indexXY = orderedBoundaryCentroids.get(index);
        
        int nT1 = template.getN();
        
        float factor = 2.f;
        float maxDist = factor * (float)Math.max(dimensionsT[0], dimensionsT[1]);
        
        // ------- initialize the local maps and store partial results in
        //            output vars too
        
        // to keep the bitstring smaller:
        int maxIdx = Integer.MIN_VALUE;
        
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
            
            if (i > maxIdx) {
                maxIdx = i;
            }
        }
        
        for (int i = 0; i < indexes.size(); ++i) {
            
            //TODO: when convert to global method,
            // consider re-using the aggregate map information
            // IF it is present for key and IF all members
            // of the bitstring for the key are present in indexes
            // (the check for presence should help limit the result to adjacent
            // within a size constraint).
            
            int idx1 = indexes.get(i);
            
            PairInt key1 = new PairInt(idx1, idx1);
            
            Double cost1 = aggregatedCostMap.get(key1);
            
            VeryLongBitString bs1 = new VeryLongBitString(maxIdx);
            bs1.setBit(idx1);
            
            PairIntArray boundary1 = orderedBoundaries.get(idx1);
            
            for (int j = 0; j < indexes.size(); ++j) {
                
                if (i == j) {
                    distMap.put(key1, cost1);
                    prevMap.put(key1, Integer.valueOf(-1));
                    indexesMap.put(key1, bs1);
                    continue;
                }
                
                int idx2 = indexes.get(j);
                
                if (!adjacencyMap.containsKey(idx1) || 
                    !adjacencyMap.get(idx1).contains(idx2)) {
                    continue;
                }
             
                VeryLongBitString bs2 = bs1.copy();
                bs2.setBit(idx2);
                
                PairInt key2;
                if (idx1 < idx2) {
                    key2 = new PairInt(idx1, idx2);
                    prevMap.put(key2, Integer.valueOf(idx1));
                } else {
                    key2 = new PairInt(idx2, idx1);
                    prevMap.put(key2, Integer.valueOf(idx2));
                }
                indexesMap.put(key2, bs2);
                
                PairIntArray boundary12 = aggregatedBoundaries.get(bs2);
                Result result12 = aggregatedResultMap.get(bs2);
                Double cost12 = aggregatedCostMap.get(bs2);
                VeryLongBitString bs12 = new VeryLongBitString(maxIdx);
                bs12.setBit(idx1);  
                bs12.setBit(idx2);
                if (boundary12 == null) {
                    PairIntArray boundary2 = orderedBoundaries.get(idx2);
                    boundary12 = pFinder2.mergeAdjacentOrderedBorders(
                        boundary1, boundary2);
                    
                    bs12 = bs2;
                    
                    aggregatedBoundaries.put(bs12, boundary12);
                                        
                    if (boundary12.getN() < 6) {
                        continue;
                    }
                                        
                    PartialShapeMatcher matcher = new PartialShapeMatcher();
                    result12 = matcher.match(boundary12, template);
                    if (result12 == null) {
                        continue;
                    }
                    aggregatedResultMap.put(bs12, result12);
                
                    // NOTE: this may need to be revised.  using the max diff chord
                    // sum from only the single segmented cell matches as
                    // the normalization for the Salukwdze distance.
                    // May need to re-examine the bounds values and adjust the
                    // normalizations.
                    int nI = result12.getNumberOfMatches();
                    float f = 1.f - ((float)nI/(float)nT1);            
                    double d = result12.getChordDiffSum()/maxDiffChordSum[0];
                    float s = (float)Math.sqrt(f * f + d * d);
                    cost12 = Double.valueOf(s);
                
                    aggregatedCostMap.put(bs12, cost12);                    
                }
                
                distMap.put(key2, cost12);
            }
        }
        // ---- end initializing local maps
        
        PairInt minCostIdx = null;
        double minCost = Double.MAX_VALUE;
         
        /*
        s0 = costIJ;
        s1 = costIK + costKJ;
        ------
        indexes = [0, 3, 5]
        fills in distMap and prevMap for
            0 0
            0 3
            0 5
            3 3
            3 5
            5 5
        
        k=[0,3,5]
          i=[0,3,5]
            j=[0,3,5]
              k=0,i=0,j=0
                
        for k=0, using already calculated pairs of segments:
        cost[0,0] = cost(bs0+0) already calc'ed
        cost[0,3] = cost(bs0+3) already calc'ed
        
        for k=3,i=0,j=0: s0=cost[0,0], s1=cost[0,3] 
        k=3,i=0,j=3: s0=cost[0,0], s1=cost[0,3] + cost[3,3]
                                   in this case, need to calc combined
                                   results.
                                    - get bitstrings from both to make combined bs
                                    - create combined boundary
                                    - cost = template match
        
        */
        
        for (int i0 = 0; i0 < indexes.size(); ++i0) {
            int k = indexes.get(i0);
            for (int i1 = 0; i1 < indexes.size(); ++i1) {
                int i = indexes.get(i1);
                for (int i2 = 0; i2 < indexes.size(); ++i2) {
                    int j = indexes.get(i2);
                    
                    boolean setPrev = true;                 
                    if (i == j) {
                        setPrev = false;
                    }
                    
                    PairInt keyIJ, keyIK, keyKJ;
                    if (i <= j) {
                        keyIJ = new PairInt(i, j);
                    } else {
                        keyIJ = new PairInt(j, i);
                    }
                    if (i <= k) {
                        keyIK = new PairInt(i, k);
                    } else {
                        keyIK = new PairInt(k, i);
                    }
                    if (k <= j) {
                        keyKJ = new PairInt(k, j);
                    } else {
                        keyKJ = new PairInt(j, k);
                    }
                    
                    // set these by adjacency.  false if not adjacent
                    boolean tIJ = false;
                    boolean tIK = false; 
                    boolean tKJ = false;
                                        
                    double s0 = Double.MAX_VALUE;
                    double s1 = Double.MAX_VALUE;
                    //s0 = costIJ;
                    //s1 = costIK + costKJ;
                    
                    VeryLongBitString bs0 = null;
                    VeryLongBitString bs1 = null;
                    
                    if (adjacencyMap.containsKey(i) && adjacencyMap.get(i).contains(j)) {
                        tIJ = true;
                        if (distMap.containsKey(keyIJ)) {
                            s0 = distMap.get(keyIJ);
                            bs0 = indexesMap.get(keyIJ);
                        } else {
                            bs0 = indexesMap.get(
                                new PairInt(keyIJ.getX(), keyIJ.getX())).copy();
                            bs0.setBit(keyIJ.getY());
                            assert(aggregatedCostMap.containsKey(bs0));
                            s0 = aggregatedCostMap.get(bs0);
                        }
                    }
                    
                    if (adjacencyMap.containsKey(i) && adjacencyMap.get(i).contains(k)) {
                        tIK = true;
                    }
                    if (adjacencyMap.containsKey(k) && adjacencyMap.get(k).contains(j)) {
                        tKJ = true;
                    }
                    if (tIK && tKJ) {
                        //s1 = costIK + costKJ;
                        VeryLongBitString bsIK = indexesMap.get(keyIK);
                        VeryLongBitString bsKJ = indexesMap.get(keyKJ);
                        VeryLongBitString bsIKKJ = bsIK.or(bsKJ);
                        PairIntArray bIKKJ = aggregatedBoundaries.get(bsIKKJ);
                        s1 = aggregatedCostMap.get(bsIKKJ);
                        
                        if (bIKKJ == null) {
                            bIKKJ = createAggregatedBoundary(orderedBoundaries,
                                bsIKKJ);
                            PartialShapeMatcher matcher = new PartialShapeMatcher();
                            Result r = matcher.match(bIKKJ, template);
                            
                            int nI = r.getNumberOfMatches();
                            float f = 1.f - ((float)nI/(float)nT1);            
                            double d = r.getChordDiffSum()/maxDiffChordSum[0];
                            s1 = (float)Math.sqrt(f * f + d * d);
                           
                            aggregatedBoundaries.put(bsIKKJ, bIKKJ);
                            aggregatedCostMap.put(bsIKKJ, Double.valueOf(s1));
                            aggregatedResultMap.put(bsIKKJ, r);
                        }
                    }
                    
                    if (i == j) {
                        // this is already set to the shape match... should not be 0
                        //dist[i][j] = 0;
                    } else if ((s0 <= s1) || (!tIK || !tKJ)) {
                        assert(tIJ);
                        distMap.put(keyIJ, Double.valueOf(s0));
                        indexesMap.put(keyIJ, bs0);
                    } else {
                        assert(tIK);
                        assert(tKJ);
                        distMap.put(keyIJ, Double.valueOf(s1));
                        if (setPrev) {
                            prevMap.put(keyIJ, Integer.valueOf(prevMap.get(tKJ)));
                        }
                        indexesMap.put(keyIJ, bs1);
                    }
                    
                    if (distMap.get(keyIJ) < minCost) {
                        minCost = distMap.get(keyIJ);
                        minCostIdx = keyIJ;
                    }
                }
            }
        }
        
        if (minCostIdx == null) {
            return null;
        }
        
        VeryLongBitString bs = indexesMap.get(minCostIdx);
        
        System.out.println("minCostIdx=" + minCostIdx + 
            " minCost=" + minCost + " bs=" +
            Arrays.toString(bs.getSetBits()));
        
        return bs;
    }
    
    private double distance(PairInt xy1, PairInt xy2) {
        
        int diffX = xy1.getX() - xy2.getX();
        int diffY = xy1.getY() - xy2.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return dist;
    }

    private PairIntArray createAggregatedBoundary(
        List<PairIntArray> orderedBoundaries, VeryLongBitString bsIKKJ) {
        
        PerimeterFinder2 perFinder2 = new PerimeterFinder2();
        
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
