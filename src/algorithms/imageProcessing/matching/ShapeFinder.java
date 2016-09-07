package algorithms.imageProcessing.matching;

import algorithms.util.PairIntArray;
import gnu.trove.map.TIntIntMap;
import java.util.List;

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
    public PartialShapeMatcher.Result findMatchingCells(
        List<PairIntArray> orderedBoundaries,
        TIntIntMap adjacencyMap) {
        
        /*
        (1) O(N_cells) + O(N_cells_in_bin * log_2(N_cells_in_bin)):
           the bin dimensions will be max of template dimensions.
           actually, should probably be the diagonal for max dist
           between any two points but that could be done more
           accurately.
           - 2D bins, that is x and y.
           - sort the cells by x centroid and ties by y centroid
        (2) 
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
        */
        
        throw new UnsupportedOperationException("not yet implemented");
    }
}
