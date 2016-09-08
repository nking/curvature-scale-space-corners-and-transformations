package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TObjectIntIterator;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

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
        TIntIntMap adjacencyMap, PairIntArray template) {
        
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
        TIntIntMap adjacencyMap, PairIntArray template) {
    
        List<Result> outputSortedResults = new ArrayList<Result>();
        List<Integer> outputSortedIndexes = new ArrayList<Integer>();
     
        matchAndOrderByIncrCost(orderedBoundaries, template, 
            outputSortedResults, outputSortedIndexes);
       
        int[] dimensionsT = calcDimensions(template);
        
        int areaT = (int)Math.round(Math.sqrt(dimensionsT[0] * dimensionsT[0] +
            dimensionsT[1] * dimensionsT[1]));
        
        Map<VeryLongBitString, PairIntArray> aggregatedBoundaries =
            new HashMap<VeryLongBitString, PairIntArray>();
        
        Map<VeryLongBitString, Result> aggregatedResultMap =
            new HashMap<VeryLongBitString, Result>();
        
        for (int i = 0; i < outputSortedIndexes.size(); ++i) {
            
            Integer index = outputSortedIndexes.get(i);
            
            /* Dijkstra's search:
            source is index, no dest, but search only extends to
                 cells within distance to 2 or so times template
                 dimensions from index centroid.  that is done
                 when the heap is populated.
            
            VeryLongBitString[] aggregatedKeys
            double[] costsFromS
            int[] prev
            while (!heap.isEmpty()) {
                u = heap.pop
                for (int v : adjacentToU) {
                    // calc cost from current aggregated u + v bounds
                    // after check for existing.
                    altCost = costsFromS[u] + cost[u->v]
                    if (altCost < costsFromS[v]) {
                        costsFromS[v] = altCost
                        heap.decreaseKey(v, altCost);
                        prev[v] = u
                        VeryLongBitString vKey = aKey.copy() + 
                            set bit for v
                        store vKey in AggregAtedKeys
                        store results in aggregatedBoundaries and 
                            aggregatedResultMap
            */
        }
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private void matchAndOrderByIncrCost(List<PairIntArray> orderedBoundaries, 
        PairIntArray template, List<Result> outputSortedResults,
        List<Integer> outputSortedIndexes) {

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
    
    private TObjectIntMap<PartialShapeMatcher.Result> matchIndividually(
        List<PairIntArray> orderedBoundaries, PairIntArray template) {
    
        // key=result, value=index of orderedBoundaries item
        TObjectIntMap<PartialShapeMatcher.Result> output = 
            new TObjectIntHashMap<PartialShapeMatcher.Result>();
        
        int[] dimensionsT = calcDimensions(template);
        
        int areaT = (int)Math.round(Math.sqrt(dimensionsT[0] * dimensionsT[0] +
            dimensionsT[1] * dimensionsT[1]));
        
        for (int i = 0; i < orderedBoundaries.size(); ++i) {
            
            PairIntArray p = orderedBoundaries.get(i);
            
            int[] dimensions = calcDimensions(template);
        
            int area = (int)Math.round(Math.sqrt(dimensions[0] * dimensions[0] +
                dimensions[1] * dimensions[1]));
            
            if (area > (3 * areaT)) {
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
            
            if (area < (areaT/10)) {
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
        
        return new int[]{minX, maxX, minY, maxY};
    }
    
}
