package algorithms.mst;

import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;

/**
 * class to handle path sum and edits of
 * a tour, used by TSPPrimsNST.
 * 
 * @author nichole
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
     * adjCostMap key = vertex index1, 
     *   value = map with key = index2 and value = 
     *   cost for edge index1 to index2. 
     */
    private final TIntObjectMap<TIntIntMap> 
        adjCostMap;
    
    /**
     * map with key = vertex indexes, that is the
       values of tour[i] values which are also the indexes
       of the vertexes given as a graph to the TSP class.
       map value = index of tour array.
     */
    private final TIntIntMap tourValueToIndexMap; 
    
    private int pathSum = 0;

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
     * given indexes0, find which swap among the 2 pairs
     * produces the longest path and return the pathSum
     * and populate outputIndexes with the suggested 
     * swap indexes.
     * @param vertexIdxs0
     * @param outputVertexIdxs
     * @return 
     */
    public int findNonIntersectingBestSwap(
        int[] vertexIdxs0,
        int[] outputVertexIdxs) {
        
        // NOTE: need to make sure the directions are fine
        // with adjCostMap
        // and the best for shortes path sum
        
        throw new UnsupportedOperationException(
            "not yet implemented");
    }
    
    /**
     * given tour values (which are the original graph
     * vertexes), calculate what the pathSum would be
     * if those were changed to vertexIdxs1.
     * @param vertexIdxs0
     * @param vertexIdxs1
     * @return 
     */
    protected int peekSumPathChanges(int[] vertexIdxs0, 
        int[] vertexIdxs1) {
        
        throw new UnsupportedOperationException(
            "not yet implemented");
    }
    
    /**
     * given tour values (which are the original graph
     * vertexes), change to vertexIdxs1 and return the
     * updated pathSum.
     * @param vertexIdxs0
     * @param vertexIdxs1
     * @return 
     */
    public int changePaths(int[] vertexIdxs0, 
        int[] vertexIdxs1) {
        
        throw new UnsupportedOperationException(
            "not yet implemented");
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
}
