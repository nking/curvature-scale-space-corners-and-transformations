package algorithms.mst;

import algorithms.SubsetChooser;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import java.util.Arrays;

/**
 class to handle path sum and edits of
 a tour, used by TSPPrimsNST.
  
 NOTE, the method and class
 are meant for single threaded use only due to
 time and memory saving internal variable use.
 
 @author nichole
 */
public class TourHandler {
   
    /**
     * array w/ indexes being the tour order and
     * values being the original graph
     * vertexes.
     * NOTE: not editing last vertex because it
     * it a repeat of the first, unless the
     * first vertex is changed.
    */
    private final int[] tour;
    
    /**
     adjCostMap key = vertex index1, 
     value = map with key = index2 and value = 
     cost for edge index1 to index2. 
    */
    private final TIntObjectMap<TIntIntMap> 
        adjCostMap;
    
    /**
    map with key = vertex indexes, that is the
    values of tour[i] values which are also the indexes
    of the vertexes given as a graph to the TSP class.
    map value = index of tour array.
    */
    private final TIntIntMap tourValueToIndexMap; 
    
    private int pathSum = 0;

    /**
     * reserved as a cache for use in
     * method findNonIntersectingBestSwap
    */
    private int[] cache0 = new int[4];
    private int[] cache1 = new int[4];
    
    /**
     * constructor for tour handler.  Note that theTour
     * is not copied, and is modified by changePaths(). 
     * @param theTour
     * @param theAdjacencyCostMap 
     */
    public TourHandler(int[] theTour, 
        TIntObjectMap<TIntIntMap> theAdjacencyCostMap) {
        
        tour = theTour;
        
        adjCostMap = theAdjacencyCostMap;
        
        tourValueToIndexMap = new TIntIntHashMap();
    
        pathSum = 0;
    }
    
    private void init() {
        
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx = tour[i];
            tourValueToIndexMap.put(cIdx, i);
        }
        
        pathSum = 0;
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx1 = tour[i];
            int cIdx2 = tour[i + 1];
            int cost = adjCostMap.get(cIdx1).get(cIdx2);
            pathSum += cost;
        }
    }
    
    /**
     * given first indexes of edge1 and edge2 in tour but
     * in the reference frame of vertex indexes, 
     * find which swap among the 2 pairs
     * produces the longest path and return the pathSum
     * and populate outputIndexes with the suggested 
     * swap indexes.  NOTE, the method and class
     * are meant for single threaded use only due to
     * time and memory savings internal variable use.
     * @param idxEdgeAVertex1 edge A vertex 1 index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     * @param idxEdgeBVertex1 edge B vertex 1 index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     * @param outputVertexIdxs array of size 4
     * holding the combination of vertexes that result
     * in best path sum.
     * @return best path sum for a combination of swapping
     * edge vertexes between the 2 edges given else -1
     * if a better combination than current was not found.
     */
    public int findNonIntersectingBestSwap(
        int idxEdgeAVertex1, int idxEdgeBVertex1,
        int[] outputVertexIdxs) {
        
        int tIdxA1 = getTourIndex(idxEdgeAVertex1);
        int tIdxA2 = getNextTourIndex(tIdxA1);
        int tIdxPrevA1 = getPrevTourIndex(tIdxA1);
        int tIdxNextA2 = getNextTourIndex(tIdxA2);
            
        int tIdxB1 = getTourIndex(idxEdgeBVertex1);
        int tIdxB2 = getTourIndex(tIdxB1);
        int tIdxPrevB1 = getPrevTourIndex(tIdxB1);
        int tIdxNextB2 = getNextTourIndex(tIdxB2);
        
        int minSum = Integer.MAX_VALUE;
        
        cache0[0] = tIdxA1;
        cache0[1] = tIdxA2;
        cache0[2] = tIdxB1;
        cache0[3] = tIdxB2;
        
        //6 permutations to try
        SubsetChooser sc = new SubsetChooser(4, 2);
        for (int i = 0; i < 6; ++i) {
            
            sc.getNextSubset(cache1);
            
            //TODO: place hard-wired skip of same permuation
            // here and remove Arrays.equals when have
            // the i index
            
            for (int j = 0; j < cache1.length; ++j) {
                cache1[j] = cache0[cache1[j]];
            }
            if (Arrays.equals(cache0, cache1)) {
                // TODO: capture this and hardwire the
                // skip for this i
                System.out.println("combination index where "
                    + "they are the same=" + i);
                continue;
            }
            boolean isValid = validateEdges(tIdxPrevA1,
                tIdxNextA2, tIdxPrevB1, tIdxNextB2,
                cache1);
            
            if (!isValid) {
                continue;
            }
           
            int sum = peekSumPathChanges(tIdxPrevA1,
                tIdxA1, tIdxA2, tIdxNextA2, 
                tIdxPrevB1, tIdxB1, tIdxB2, tIdxNextB2,
                cache1);
            
            if (sum < minSum) {
                sum = minSum;
                System.arraycopy(cache1, 0, outputVertexIdxs, 
                    0, cache1.length);
            }
        }
        
        return (minSum == Integer.MAX_VALUE) ? -1 : minSum;
    }
    
    /**
     * given the vertex index from the original graph,
     * return the vertex index the follows it in the
     * tour.
     * @param vertexIdx
     * @return 
     */
    public int getNextVertexIndex(int vertexIdx) {
        int tIdx = getTourIndex(vertexIdx);
        return getNextTourIndex(tIdx);
    }
    
    private int getNextTourIndex(int tIdx) {
        // last index is same as first.  always
        // return first index.
        if (tIdx == (tour.length - 2)) {
            return 0;
        }
        return tIdx + 1;
    }
    
    private int getPrevTourIndex(int tIdx) {
        // last index is same as first.  always
        // return first index.
        if (tIdx == 0) {
            return tour.length - 2;
        }
        return tIdx - 1;
    }
    
    /**
     * given tour indexes A1, A2, B1, B2 and the
     * vertexes to change those to, calculate what
     * the path sum would be with those changes,
     * but do not apply the changes.
     * The tour indexes before and after edgeA and 
     * edgeB are passed in to avoid recalculating.
     * The invoking method has the responsibility of 
     * passing valid arguments to this method.
     @param tIdxPrevA1 is the tour index before that
     * containing vertexIdxs0[0].
     @param tIdxA1 tour index for first point of edgeA
     @param tIdxA2 tour index for last point of edgeA
     @param tIdxNextA2 is the tour index after that
       containing vertexIdxs0[1]
     @param tIdxPrevB1 is the tour index before that
     * containing vertexIdxs0[2]
     @param tIdxB1 tour index for first point of edgeB
     @param tIdxB2 tour index for last point of edgeB
     @param tIdxNextB2 is the tour index after that
       containing vertexIdxs0[3]
     @param vertexIdxs0 contains the indexes
     * for edge A vertex 1 and 2, then edge B vertex 
     * 1 and 2 that are currently in the tour.
     @param vertexIdxs1 contains the candidate 
     * change indexes
     * for edge A vertex 1 and 2, then edge B vertex 
     * 1 and 2.
     * @return the calculated path sum for the path
     * changed from vertexes vertexIdxs0 to vertexIdxs1.
     */
    private int peekSumPathChanges(int tIdxPrevA1,
        int tIdxA1, int tIdxA2, int tIdxNextA2, 
        int tIdxPrevB1, int tIdxB1, int tIdxB2, int tIdxNextB2,
        int[] vertexIdxs1) {
        
        int sum = pathSum;
        
        int minusA = 
            adjCostMap.get(tIdxPrevA1).get(tIdxA1) +
            adjCostMap.get(tIdxA1).get(tIdxA2) +
            adjCostMap.get(tIdxA2).get(tIdxNextA2);
        
        int minusB = 
            adjCostMap.get(tIdxPrevB1).get(tIdxB1) +
            adjCostMap.get(tIdxB1).get(tIdxB2) +
            adjCostMap.get(tIdxB2).get(tIdxNextB2);
        
        int plusA =
            adjCostMap.get(tIdxPrevA1)
                .get(getTourIndex(vertexIdxs1[0])) +
            adjCostMap.get(getTourIndex(vertexIdxs1[0]))
                .get(getTourIndex(vertexIdxs1[1])) +
            adjCostMap.get(getTourIndex(vertexIdxs1[1]))
                .get(tIdxNextA2);
        
        int plusB =
            adjCostMap.get(tIdxPrevB1)
                .get(getTourIndex(vertexIdxs1[2])) +
            adjCostMap.get(getTourIndex(vertexIdxs1[2]))
                .get(getTourIndex(vertexIdxs1[3])) +
            adjCostMap.get(getTourIndex(vertexIdxs1[3]))
                .get(tIdxNextB2);
        
        sum -= minusA;
        sum -= minusA;
        sum += plusA;
        sum += plusB;
        
        return sum;
    }
    
    /**
     * given tour values (which are the original graph
     * vertexes), change to vertexIdxs1 and return the
     * updated pathSum.
     @param idxEdgeAVertex1 edge A vertex 1 index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     @param idxEdgeBVertex1 edge B vertex 1 index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     @param vertexIdxs the vertexes to change edge A
     * and edge B to in thr tour
     @return the updated path sum after the changes
     * are applied.
    */
    public int changePaths(int idxEdgeAVertex1, 
        int idxEdgeBVertex1, int[] vertexIdxs) {
        
        int tIdxA1 = getTourIndex(idxEdgeAVertex1);
        int tIdxA2 = getNextTourIndex(tIdxA1);
        int tIdxPrevA1 = getPrevTourIndex(tIdxA1);
        int tIdxNextA2 = getNextTourIndex(tIdxA2);
            
        int tIdxB1 = getTourIndex(idxEdgeBVertex1);
        int tIdxB2 = getTourIndex(tIdxB1);
        int tIdxPrevB1 = getPrevTourIndex(tIdxB1);
        int tIdxNextB2 = getNextTourIndex(tIdxB2);
        
        int sum = peekSumPathChanges(tIdxPrevA1,
            tIdxA1, tIdxA2, tIdxNextA2, 
            tIdxPrevB1, tIdxB1, tIdxB2, tIdxNextB2, 
            vertexIdxs);
        
        tour[tIdxA1] = vertexIdxs[0];
        tour[tIdxA2] = vertexIdxs[1];
        tour[tIdxB1] = vertexIdxs[2];
        tour[tIdxB2] = vertexIdxs[3];
        
        return sum;
    }
  
    public int getPathSum() {
        return pathSum;
    }
    
    public int[] getTour() {
        return tour;
    }
    
    /**
     * given the original graph vertex index,
     * return the index of that point in the 
     * tour.  Note that if the result is 0,
     * the vertex is also present as the last
     * item of the tour, but that index is not returned.
     * @param vertexIdx
     * @return 
     */
    public int getTourIndex(int vertexIdx) {
        if (vertexIdx < 0 || (vertexIdx > (tour.length - 2))) {
            throw new IllegalArgumentException(
            "vertexIdx is out of bounds");
        }
        int tIdx = tourValueToIndexMap.get(vertexIdx);
        return tIdx;
    }
    
    public int getVertexIndex(int tourIdx) {
        if (tourIdx < 0 || (tourIdx > (tour.length - 1))) {
            throw new IllegalArgumentException(
            "tourIdx is out of bounds");
        }
        return tour[tourIdx];
    }

    private boolean validateEdges(int tIdxPrevA1, 
        int tIdxNextA2, int tIdxPrevB1, int tIdxNextB2,
        int[] tIndexes) {
        
        // validate that edge exists
        for (int j = 0; j < tIndexes.length; ++j) {
            int edge1, edge2;
            if (j == 0) {
                edge1 = tIdxPrevA1;
                edge2 = tIndexes[j];
            } else if (j == 1) {
                edge1 = tIndexes[j];
                edge2 = tIdxNextA2;
            } else if (j == 2) {
                edge1 = tIdxPrevB1;
                edge2 = tIndexes[j];
            } else {
                edge1 = tIndexes[j];
                edge2 = tIdxNextB2;
            }
            edge1 = getVertexIndex(edge1);
            edge2 = getVertexIndex(edge2);
            if (!(adjCostMap.containsKey(edge1)
                && adjCostMap.get(edge1)
                .containsKey(edge2))) {
                return false;
            }
            if (j == (tIndexes.length - 1)) {
                // try the combination if passed other checks
                edge1 = getVertexIndex(tIndexes[0]);
                edge2 = getVertexIndex(tIndexes[1]);
                if (!(adjCostMap.containsKey(edge1)
                    && adjCostMap.get(edge1)
                    .containsKey(edge2))) {
                    return false;
                }
                edge1 = getVertexIndex(tIndexes[2]);
                edge2 = getVertexIndex(tIndexes[3]);
                if (!(adjCostMap.containsKey(edge1)
                    && adjCostMap.get(edge1)
                    .containsKey(edge2))) {
                    return false;
                }
            }
        }
        return true;
    }
}
