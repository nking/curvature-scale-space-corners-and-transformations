package algorithms.mst;

import algorithms.SubsetChooser;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairInt;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.Arrays;
import java.util.List;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTreeInterval2D;

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
 
    private final PairInt[] coordinates;
    
    private final QuadTreeInterval2D<Integer, PairInt> qt;
        
    // bounding boxes of tour[i] to tour[i+1]
    // stored by keys that are the first value,
    // tour[i] which is the vertex index
    private final TIntObjectMap<Interval2D<Integer>> indexEdgeBounds;
        
    private final TObjectIntMap<Interval2D<Integer>> edgeIndexBounds;
        
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
        TIntObjectMap<TIntIntMap> theAdjacencyCostMap,
        PairInt[] theCoordinates) {
        
        tour = theTour;
        
        adjCostMap = theAdjacencyCostMap;
        
        tourValueToIndexMap = new TIntIntHashMap();
    
        pathSum = 0;
        
        coordinates = theCoordinates;
        
        qt = new QuadTreeInterval2D<Integer, PairInt>();

        indexEdgeBounds
            = new TIntObjectHashMap<Interval2D<Integer>>();

        edgeIndexBounds
            = new TObjectIntHashMap<Interval2D<Integer>>();

        init();
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
            
            insertEdgeBox(cIdx1, cIdx2);
        }   
    }
    
    /**
     * given first index n the reference frame of 
     * vertex indexes,
     * find the edges in thr tour that it intersects
     * with and return the be improvement if any.
     * The result will be a -1 if the edge did not 
     * intersect.
     * 
     * @param idxEdgeAVertex1 edge A vertex 1 index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     * @param outputVertexIdxs array of size 4
     * holding the combination of vertexes that result
     * in best path sum.
     * @param outputIdxEdgeBVertex1 an array of size 1
     * to hold the vertex b of best edge to swap to
     * uncross intersecting edges.
     * @return best path sum for a combination of swapping
     * edge vertexes between the 2 edges given else -1
     * if a better combination than current was not found.
     */
    public int findNonIntersectingBestSwap(
        int idxEdgeAVertex1, 
        int[] outputIdxEdgeBVertex1, 
        int[] outputVertexIdxs) {
        
        int tIdxA1 = getTourIndex(idxEdgeAVertex1);
        int tIdxA2 = getNextTourIndex(tIdxA1);
        int tIdxPrevA1 = getPrevTourIndex(tIdxA1);
        int tIdxNextA2 = getNextTourIndex(tIdxA2);
 
        int idxEdgeAVertex2 = getVertexIndex(tIdxA2);

        Interval2D<Integer> box12 = 
            indexEdgeBounds.get(idxEdgeAVertex1);

        // find intersection boxes and look for 
        // edges that intersect with these

        List<Interval2D<Integer>> list = qt.query2D(box12);

        if (list.size() < 2) {
            return -1;
        }
        
        int[] tmp = new int[4];
        int minSum = Integer.MAX_VALUE;
        
        int x1 = coordinates[idxEdgeAVertex1].getX();
        int y1 = coordinates[idxEdgeAVertex1].getY();

        int x2 = coordinates[idxEdgeAVertex2].getX();
        int y2 = coordinates[idxEdgeAVertex2].getY();

        for (int listIdx = 0; listIdx < list.size();
            ++listIdx) {

            Interval2D<Integer> box34 = list.get(listIdx);

            if (box34.equals(box12)) {
                continue;
            }

            int idxEdgeBVertex1 = edgeIndexBounds.get(box34);
            int tIdxB1 = getTourIndex(idxEdgeBVertex1);
            int tIdxB2 = getNextTourIndex(tIdxB1);
            int idxEdgeBVertex2 = getVertexIndex(tIdxB2);
            
            int x3 = coordinates[idxEdgeBVertex1].getX();
            int y3 = coordinates[idxEdgeBVertex1].getY();

            int x4 = coordinates[idxEdgeBVertex2].getX();
            int y4 = coordinates[idxEdgeBVertex2].getY();

            // determine if p12 p34 intersect
            if (!LinesAndAngles.linesIntersect(
                x1, y1, x2, y2, x3, y3, x4, y4)) {
                continue;
            }
            
            int sum = findNonIntersectingBestSwap(
                idxEdgeAVertex1, idxEdgeAVertex2,
                idxEdgeBVertex1, idxEdgeBVertex2, tmp);
         
            if (sum < minSum) {
                sum = minSum;
                System.arraycopy(tmp, 0, outputVertexIdxs, 
                    0, tmp.length);
                outputIdxEdgeBVertex1[0] = idxEdgeBVertex1;
            }
        }
        
        return (minSum == Integer.MAX_VALUE) ? -1 : minSum;
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
     * @param idxEdgeAVertex2 edge A vertex 3 index,
     * the second edge index that follows 
     * idxEdgeAVertex1 in the tour array.
     * @param idxEdgeBVertex1 edge B vertex 1 index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     * @param idxEdgeBVertex2 edge B vertex 3 index,
     * the second edge index that follows 
     * idxEdgeBVertex1 in the tour array.
     * @param outputVertexIdxs array of size 4
     * holding the combination of vertexes that result
     * in best path sum.
     * @return best path sum for a combination of swapping
     * edge vertexes between the 2 edges given else -1
     * if a better combination than current was not found.
     */
    protected int findNonIntersectingBestSwap(
        int idxEdgeAVertex1, int idxEdgeAVertex2,
        int idxEdgeBVertex1, int idxEdgeBVertex2,
        int[] outputVertexIdxs) {

        int tIdxA1 = getTourIndex(idxEdgeAVertex1);
        int tIdxA2 = getTourIndex(idxEdgeAVertex2);
        assert(tIdxA2 == getNextTourIndex(tIdxA1));
        int tIdxPrevA1 = getPrevTourIndex(tIdxA1);
        int tIdxNextA2 = getNextTourIndex(tIdxA2);

        int tIdxB1 = getTourIndex(idxEdgeBVertex1);
        int tIdxB2 = getTourIndex(idxEdgeBVertex2);
        assert(tIdxB2 == getNextTourIndex(tIdxB1));
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
        return tour[getNextTourIndex(tIdx)];
    }
    
    /**
     * given the vertex index from the original graph,
     * return the vertex index the precedes it in the
     * tour.
     * @param vertexIdx
     * @return 
     */
    public int getPrevVertexIndex(int vertexIdx) {
        int tIdx = getTourIndex(vertexIdx);
        return tour[getPrevTourIndex(tIdx)];
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
    
    private void removeEdgeBoxes(int tIdxPrev1, int tIdx1, 
        int tIdx2, int tIdxNext2) {

        if (true) {
            throw new UnsupportedOperationException(
                "not yet implemented");
        }
        int cIdxPrev1 = getVertexIndex(tIdxPrev1);
        int cIdx1 = getVertexIndex(tIdx1);
        int cIdx2 = getVertexIndex(tIdx2);
             
        Interval2D<Integer> box;
        
        box = indexEdgeBounds.remove(cIdxPrev1);
        assert(box != null);
        edgeIndexBounds.remove(box);
        qt.remove(box);

        box = indexEdgeBounds.remove(cIdx1);
        assert(box != null);
        edgeIndexBounds.remove(box);
        qt.remove(box);
        
        box = indexEdgeBounds.remove(cIdx2);
        assert(box != null);
        edgeIndexBounds.remove(box);
        qt.remove(box);        
    }
    
    private void insertEdgeBox(int idxEdgeVertex1,
        int idxEdgeVertex2) {

        Interval2D<Integer> box;
        int x1, y1, x2, y2;
        Interval<Integer> xi1, yi1;
        PairInt desc;
        
        x1 = coordinates[idxEdgeVertex1].getX();
        y1 = coordinates[idxEdgeVertex1].getY();
        x2 = coordinates[idxEdgeVertex2].getX();
        y2 = coordinates[idxEdgeVertex2].getY();
        xi1 = new Interval<Integer>(x1, x2);
        yi1 = new Interval<Integer>(y1, y2);
        box = new Interval2D<Integer>(xi1, yi1);
        desc = new PairInt(idxEdgeVertex1, idxEdgeVertex2);
        qt.insert(box, desc);
        edgeIndexBounds.put(box, idxEdgeVertex1);
        indexEdgeBounds.put(idxEdgeVertex1, box);
        
    }
    
    private void insertEdgeBoxes(int idxEdgeVertex1,
        int idxEdgeVertex2,
        int idxEdgePrevVertex1, int idxEdgeNextVertex2) {

        int cIdx1 = idxEdgeVertex1;
        int cIdx2 = idxEdgeVertex2;
        int cIdxPrev1 = idxEdgePrevVertex1;
        int cIdxNext2 = idxEdgeNextVertex2;

        insertEdgeBox(cIdxPrev1, cIdx1);
        insertEdgeBox(cIdx1, cIdx2);
        insertEdgeBox(cIdx2, cIdxNext2);
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
        
        // ---- update the edge bounding boxes -----
        removeEdgeBoxes(tIdxPrevA1, tIdxA1, tIdxA2, 
            tIdxNextA2);
        removeEdgeBoxes(tIdxPrevB1, tIdxB1, tIdxB2, 
            tIdxNextB2);
        
        insertEdgeBoxes(vertexIdxs[0], vertexIdxs[1],
            tIdxPrevA1, tIdxNextA2);
        insertEdgeBoxes(vertexIdxs[2], vertexIdxs[3],
            tIdxPrevB1, tIdxNextB2);
        
        // ---- update the tour -----
        
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
