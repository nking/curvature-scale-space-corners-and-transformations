package algorithms.mst;

import algorithms.SubsetChooser;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
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
    private final int[] cache0 = new int[3];
    private final int[][] cache1 = new int[6][3];
    private final int[] cache2 = new int[2];
    
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
        
        for (int i = 0; i < 6; ++i) {
            cache1[i] = new int[3];
        }
        
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx = tour[i];
            tourValueToIndexMap.put(cIdx, i);
        }
        assert(tourValueToIndexMap.size() == tour.length - 1);
        
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
     * edge vertexes between the 2 edges given else Integer.MAX_VALUE;
     * if a better combination than current was not found.
     */
    public int findNonIntersectingBestSwap(
        final int idxEdgeAVertex1, 
        int[] outputIdxEdgeBVertex1, 
        int[] outputVertexIdxs) {
        
        int tIdxA1 = getTourIndex(idxEdgeAVertex1);
        int tIdxA2 = getNextTourIndex(tIdxA1);
 
        final int idxEdgeAVertex2 = getVertexIndex(tIdxA2);

        Interval2D<Integer> box12 = 
            indexEdgeBounds.get(idxEdgeAVertex1);
        assert(box12 != null);
        
        // find intersection boxes and look for 
        // edges that intersect with these

        List<Interval2D<Integer>> list = qt.query2D(box12);

        if (list.size() < 2) {
            return Integer.MAX_VALUE;
        }
        
        //vertex indexes returned for min path sum
        int[] tmp = new int[4];
        int minSum = pathSum;
        
        int x1 = coordinates[idxEdgeAVertex1].getX();
        int y1 = coordinates[idxEdgeAVertex1].getY();

        int x2 = coordinates[idxEdgeAVertex2].getX();
        int y2 = coordinates[idxEdgeAVertex2].getY();

        for (int listIdx = 0; listIdx < list.size();
            ++listIdx) {

            Interval2D<Integer> box34 = list.get(listIdx);
            assert(box34 != null);
            
            if (box34.equals(box12)) {
                continue;
            }

            final int idxEdgeBVertex1 = edgeIndexBounds.get(box34);
            int tIdxB1 = getTourIndex(idxEdgeBVertex1);
            int tIdxB2 = getNextTourIndex(tIdxB1);
            final int idxEdgeBVertex2 = getVertexIndex(tIdxB2);
            assert(tour[tIdxB1] == idxEdgeBVertex1);
            assert(tour[tIdxB2] == getNextVertexIndex(idxEdgeBVertex1));            
            assert(tour[tIdxB2] == idxEdgeBVertex2);
            assert(tIdxB1 == getTourIndex(idxEdgeBVertex1));
       
            // vertex disjoint
            if ((tIdxA1 == tIdxB1) || (tIdxA1 == tIdxB2) ||
                (tIdxA2 == tIdxB1) || (tIdxA2 == tIdxB2)) {
                continue;
            } 
            
            int x3 = coordinates[idxEdgeBVertex1].getX();
            int y3 = coordinates[idxEdgeBVertex1].getY();

            int x4 = coordinates[idxEdgeBVertex2].getX();
            int y4 = coordinates[idxEdgeBVertex2].getY();

            // determine if p12 p34 intersect
            if (!LinesAndAngles.linesIntersect(
                x1, y1, x2, y2, x3, y3, x4, y4)) {
                continue;
            }

            int sum = findMinValidBestSwap(
                idxEdgeAVertex1, idxEdgeAVertex2,
                idxEdgeBVertex1, idxEdgeBVertex2, tmp);
         
            if (sum < minSum) {
                
                minSum = sum;
                
                // tmp are graph vertex indexes
                
                System.arraycopy(tmp, 0, outputVertexIdxs, 
                    0, tmp.length);
 
                outputIdxEdgeBVertex1[0] = idxEdgeBVertex1;
                
 assert(assertSameSets(
 getVertexIndex(tIdxA1), 
 getVertexIndex(tIdxA2), 
 getVertexIndex(tIdxB1), 
 getVertexIndex(tIdxB2),
 outputVertexIdxs[0], 
 outputVertexIdxs[1], 
 outputVertexIdxs[2],
 outputVertexIdxs[3]));
 assert(assertSameSets(
 tIdxA1, tIdxA2, tIdxB1, tIdxB2,
 getTourIndex(outputVertexIdxs[0]), 
 getTourIndex(outputVertexIdxs[1]), 
 getTourIndex(outputVertexIdxs[2]),
 getTourIndex(outputVertexIdxs[3])));
 
            }
        }
        
        if (minSum < pathSum) {
            return minSum;
        }
        
        // ---- adding logic to handle swapping edge test while
        //      still have information that edge intersects others
        
        int tIdxPrevA1 = getPrevTourIndex(tIdxA1);
        int tIdxNextA2 = getNextTourIndex(tIdxA2);

        // else, try to swap leading edge alone
        int sum = peekSumPathChangesReverseEdge(tIdxPrevA1,
            tIdxA1, tIdxA2, tIdxNextA2);
        if (sum < minSum) {
            // not returning sum to allow existing logic to 
            // remain, it acts upon data from 2 edges
            sum = changePathsToReverseEdge(tIdxA1);
        }
        
        return (minSum == pathSum) ? Integer.MAX_VALUE : minSum;
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
    protected int findMinValidBestSwap(
        int idxEdgeAVertex1, int idxEdgeAVertex2,
        final int idxEdgeBVertex1, final int idxEdgeBVertex2,
        int[] outputVertexIdxs) {

        assert(outputVertexIdxs.length == 4);
        
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

        int minSum = pathSum;

        // fixing first integer and permuting the other 3
        cache0[0] = tIdxA2;
        cache0[1] = tIdxB1;
        cache0[2] = tIdxB2;
        
        cache2[0] = 0;
        permutation(cache0, cache1, 3, cache2);
        
        for (int i = 0; i < 6; ++i) {

            boolean isValid = validateEdges(tIdxPrevA1,
                tIdxNextA2, tIdxPrevB1, tIdxNextB2,
                tIdxA1, cache1[i][0], cache1[i][1], cache1[i][2]);

            if (!isValid) {
                continue;
            }
            
            int sum = peekSumPathChanges(tIdxPrevA1,
                tIdxA1, tIdxA2, tIdxNextA2,
                tIdxPrevB1, tIdxB1, tIdxB2, tIdxNextB2,
                tIdxA1, cache1[i][0], cache1[i][1], cache1[i][2]);

            if (sum < minSum) {
                
                minSum = sum;
                
                // outputVertexIdxs is graph vertex indexes
                // instead of tour indexes.
                // cache1[i] are tour indexes.
                
                outputVertexIdxs[0] = getVertexIndex(tIdxA1);
                for (int j = 0; j < cache1[i].length; ++j) {
                    outputVertexIdxs[j + 1] = 
                        getVertexIndex(cache1[i][j]);
                }

assert(assertSameSets(
 getVertexIndex(tIdxA1), 
 getVertexIndex(tIdxA2), 
 getVertexIndex(tIdxB1), 
 getVertexIndex(tIdxB2),
 outputVertexIdxs[0], 
 outputVertexIdxs[1], 
 outputVertexIdxs[2],
 outputVertexIdxs[3]));
 assert(assertSameSets(
 tIdxA1, tIdxA2, tIdxB1, tIdxB2,
 getTourIndex(outputVertexIdxs[0]), 
 getTourIndex(outputVertexIdxs[1]), 
 getTourIndex(outputVertexIdxs[2]),
 getTourIndex(outputVertexIdxs[3])));
 
            }
        }
        
        return (minSum == pathSum) ? Integer.MAX_VALUE : minSum;
    }

    void permutation(int a[], int[][] out, int size, int[] count) {
        if (size == 1) {
            System.arraycopy(a, 0, out[count[0]], 0, a.length);
            count[0]++;
            return;
        }
        for (int i = 0; i < size; i++) {
            permutation(a, out, size - 1, count);
            int swap = a[size - 1];
            if (size % 2 == 1) {
                a[size - 1] = a[0];
                a[0] = swap;
            } else {
                a[size - 1] = a[i];
                a[i] = swap;
            }
        }
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
     @param tIdx1 candidate change edge A vertex1 tour index
     @param tIdx2 candidate change edge A vertex2 tour index
     @param tIdx3 candidate change edge B vertex1 tour index
     @param tIdx4 candidate change edge B vertex2 tour index
     * @return the calculated path sum for the path
     * changed from vertexes vertexIdxs0 to vertexIdxs1.
     */
    protected int peekSumPathChanges(int tIdxPrevA1,
        int tIdxA1, int tIdxA2, int tIdxNextA2, 
        int tIdxPrevB1, int tIdxB1, int tIdxB2, int tIdxNextB2,
        int tIdx1, int tIdx2, int tIdx3, int tIdx4) {
        
        int sum = pathSum;
        
        int a1 = adjCostMap.get(getVertexIndex(tIdxPrevA1))
            .get(getVertexIndex(tIdxA1));
        int a2 = adjCostMap.get(getVertexIndex(tIdxA1))
            .get(getVertexIndex(tIdxA2));
        int a3 = adjCostMap.get(getVertexIndex(tIdxA2))
            .get(getVertexIndex(tIdxNextA2));
        
        int minusA = a1 + a2 + a3;
         
        int b1 = adjCostMap.get(getVertexIndex(tIdxPrevB1))
            .get(getVertexIndex(tIdxB1));
        int b2 = adjCostMap.get(getVertexIndex(tIdxB1))
            .get(getVertexIndex(tIdxB2));
        int b3 = adjCostMap.get(getVertexIndex(tIdxB2))
            .get(getVertexIndex(tIdxNextB2));
        
        int minusB = b1 + b2 + b3;
        
        int plusA1 = adjCostMap.get(getVertexIndex(tIdxPrevA1))
            .get(getVertexIndex(tIdx1));
        int plusA2 = adjCostMap.get(getVertexIndex(tIdx1))
            .get(getVertexIndex(tIdx2));
        int plusA3 = adjCostMap.get(getVertexIndex(tIdx2))
            .get(getVertexIndex(tIdxNextA2));
        
        int plusA = plusA1 + plusA2 + plusA3;
        
        int plusB1 = adjCostMap.get(getVertexIndex(tIdxPrevB1))
            .get(getVertexIndex(tIdx3));
        int plusB2 = adjCostMap.get(getVertexIndex(tIdx3))
            .get(getVertexIndex(tIdx4));
        int plusB3 = adjCostMap.get(getVertexIndex(tIdx4))
            .get(getVertexIndex(tIdxNextB2));
        
        int plusB = plusB1 + plusB2 + plusB3;
        
        sum -= minusA;
        sum -= minusB;
        sum += plusA;
        sum += plusB;
        
        return sum;
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
     
     * @return the calculated path sum for the path
     * changed from vertexes vertexIdxs0 to vertexIdxs1.
     */
    protected int peekSumPathChangesReverseEdge(int tIdxPrevA1,
        int tIdxA1, int tIdxA2, int tIdxNextA2) {
        
        int sum = pathSum;
        
        int a1 = adjCostMap.get(getVertexIndex(tIdxPrevA1))
            .get(getVertexIndex(tIdxA1));
        int a2 = adjCostMap.get(getVertexIndex(tIdxA1))
            .get(getVertexIndex(tIdxA2));
        int a3 = adjCostMap.get(getVertexIndex(tIdxA2))
            .get(getVertexIndex(tIdxNextA2));
        
        int minusA = a1 + a2 + a3;
                
        int plusA1 = adjCostMap.get(getVertexIndex(tIdxPrevA1))
            .get(getVertexIndex(tIdxA2));
        int plusA2 = adjCostMap.get(getVertexIndex(tIdxA2))
            .get(getVertexIndex(tIdxA1));
        int plusA3 = adjCostMap.get(getVertexIndex(tIdxA1))
            .get(getVertexIndex(tIdxNextA2));
        
        int plusA = plusA1 + plusA2 + plusA3;
                
        sum -= minusA;
        sum += plusA;
        
        return sum;
    }
    
    private void removeEdgeBoxes(int tIdxPrev1, int tIdx1, 
        int tIdx2, int tIdxNext2) {

        int cIdxPrev1 = getVertexIndex(tIdxPrev1);
        int cIdx1 = getVertexIndex(tIdx1);
        int cIdx2 = getVertexIndex(tIdx2);
             
        Interval2D<Integer> box;
        
        box = indexEdgeBounds.remove(cIdxPrev1);
        if (box != null) {
            edgeIndexBounds.remove(box);
            qt.remove(box);
        }

        box = indexEdgeBounds.remove(cIdx1);
        if (box != null) {
            edgeIndexBounds.remove(box);
            qt.remove(box);
        }
        
        box = indexEdgeBounds.remove(cIdx2);
        if (box != null) {
            edgeIndexBounds.remove(box);
            qt.remove(box);
        }
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
        if (x1 <= x2) {
            xi1 = new Interval<Integer>(x1, x2);
        } else  {
            xi1 = new Interval<Integer>(x2, x1);
        }
        if (y1 <= y2) {
            yi1 = new Interval<Integer>(y1, y2);
        } else {
            yi1 = new Interval<Integer>(y2, y1);
        }
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
        int tIdxB2 = getNextTourIndex(tIdxB1);
        int tIdxPrevB1 = getPrevTourIndex(tIdxB1);
        int tIdxNextB2 = getNextTourIndex(tIdxB2);
        
        // ---- update the edge bounding boxes -----
        removeEdgeBoxes(tIdxPrevA1, tIdxA1, tIdxA2, 
            tIdxNextA2);
        removeEdgeBoxes(tIdxPrevB1, tIdxB1, tIdxB2, 
            tIdxNextB2);
        
        insertEdgeBoxes(vertexIdxs[0], vertexIdxs[1],
            getVertexIndex(tIdxPrevA1), 
            getVertexIndex(tIdxNextA2));
        insertEdgeBoxes(vertexIdxs[2], vertexIdxs[3],
            getVertexIndex(tIdxPrevB1), 
            getVertexIndex(tIdxNextB2));
        
        // ---- update the tour -----
    
        pathSum = peekSumPathChanges(tIdxPrevA1,
            tIdxA1, tIdxA2, tIdxNextA2, 
            tIdxPrevB1, tIdxB1, tIdxB2, tIdxNextB2, 
            getTourIndex(vertexIdxs[0]), 
            getTourIndex(vertexIdxs[1]), 
            getTourIndex(vertexIdxs[2]),
            getTourIndex(vertexIdxs[3]));
 
 assert(assertSameSets(
 getVertexIndex(tIdxA1), 
 getVertexIndex(tIdxA2), 
 getVertexIndex(tIdxB1), 
 getVertexIndex(tIdxB2),
 vertexIdxs[0], 
 vertexIdxs[1], 
 vertexIdxs[2],
 vertexIdxs[3]));
 assert(assertSameSets(
 tIdxA1, tIdxA2, tIdxB1, tIdxB2,
 getTourIndex(vertexIdxs[0]), 
 getTourIndex(vertexIdxs[1]), 
 getTourIndex(vertexIdxs[2]),
 getTourIndex(vertexIdxs[3])));
 
 assert(assertTourData());
 System.out.println(String.format("change tour[%d] = %d to %d", 
 tIdxA1, tour[tIdxA1], vertexIdxs[0]));
 System.out.println(String.format("change tour[%d] = %d to %d", 
 tIdxA2, tour[tIdxA2], vertexIdxs[1]));
 System.out.println(String.format("change tour[%d] = %d to %d", 
 tIdxB1, tour[tIdxB1], vertexIdxs[2]));
 System.out.println(String.format("change tour[%d] = %d to %d", 
 tIdxB2, tour[tIdxB2], vertexIdxs[3]));
 
        tour[tIdxA1] = vertexIdxs[0];
        tour[tIdxA2] = vertexIdxs[1];
        tour[tIdxB1] = vertexIdxs[2];
        tour[tIdxB2] = vertexIdxs[3];
 
        if (tIdxA2 == 0) {
            tour[tour.length - 1] = vertexIdxs[1];
        } else if (tIdxB1 == 0) {
            tour[tour.length - 1] = vertexIdxs[2];
        } else if (tIdxB2 == 0) {
            tour[tour.length - 1] = vertexIdxs[3];
        }
        
        // ---- update the values of tour map
        // since the changes are swaps, all changed items
        // are present here in the updates
        tourValueToIndexMap.put(vertexIdxs[0], tIdxA1);
        tourValueToIndexMap.put(vertexIdxs[1], tIdxA2);
        tourValueToIndexMap.put(vertexIdxs[2], tIdxB1);
        tourValueToIndexMap.put(vertexIdxs[3], tIdxB2);
         
        assert(assertTourData());
       
System.out.println("new pathSum==" + pathSum + 
" tour=" + Arrays.toString(tour));

        return pathSum;
    }
  
    /**
     * given tour values (which are the original graph
     * vertexes), change to vertexIdxs1 and return the
     * updated pathSum.
     @param tIdxA1 edge A vertex 1 tour index,
     * the second edge index is implicitly the one
     * that follows this in the tour array.
     @return the updated path sum after the changes
     * are applied.
    */
    public int changePathsToReverseEdge(int tIdxA1) {
    
        int tIdxA2 = getNextTourIndex(tIdxA1);
        int tIdxPrevA1 = getPrevTourIndex(tIdxA1);
        int tIdxNextA2 = getNextTourIndex(tIdxA2);
        
        // ---- update the edge bounding boxes -----
        removeEdgeBoxes(tIdxPrevA1, tIdxA1, tIdxA2, 
            tIdxNextA2);
        
        insertEdgeBoxes(
            getVertexIndex(tIdxA1),
            getVertexIndex(tIdxA2),
            getVertexIndex(tIdxPrevA1), 
            getVertexIndex(tIdxNextA2));
        
        // ---- update the tour -----
    
        pathSum = peekSumPathChangesReverseEdge(
            tIdxPrevA1, tIdxA1, tIdxA2, tIdxNextA2);
 
        tour[tIdxA1] = getVertexIndex(tIdxA2);
        tour[tIdxA2] = getVertexIndex(tIdxA1);
 
        if (tIdxA2 == 0) {
            tour[tour.length - 1] = tour[tIdxA2];
        }
        
        // ---- update the values of tour map
        // since the changes are swaps, all changed items
        // are present here in the updates
        tourValueToIndexMap.put(tour[tIdxA1], tIdxA1);
        tourValueToIndexMap.put(tour[tIdxA2], tIdxA2);
         
        assert(assertTourData());
       
System.out.println("new pathSum==" + pathSum + 
" tour=" + Arrays.toString(tour));

        return pathSum;
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

    /**
     * 
     * @param tIdxPrevA1
     * @param tIdxNextA2
     * @param tIdxPrevB1
     * @param tIdxNextB2
     @param tIdx1 candidate change edge A 1 tour vertex
     @param tIdx2 candidate change edge A 2 tour vertex
     @param tIdx3 candidate change edge B 1 tour vertex
     @param tIdx4 candidate change edge B 2 tour vertex
     * @return 
     */
    private boolean validateEdges(int tIdxPrevA1, 
        int tIdxNextA2, int tIdxPrevB1, int tIdxNextB2,
        int tIdx1, int tIdx2, int tIdx3, int tIdx4) {
        
        // validate that edge exists
        for (int j = 0; j < 4; ++j) {
            int edge1, edge2;
            if (j == 0) {
                edge1 = getVertexIndex(tIdxPrevA1);
                edge2 = getVertexIndex(tIdx1);
            } else if (j == 1) {
                edge1 = getVertexIndex(tIdx2);
                edge2 = getVertexIndex(tIdxNextA2);
            } else if (j == 2) {
                edge1 = getVertexIndex(tIdxPrevB1);
                edge2 = getVertexIndex(tIdx3);
            } else {
                edge1 = getVertexIndex(tIdx4);
                edge2 = getVertexIndex(tIdxNextB2);
            }
           
            if (!(adjCostMap.containsKey(edge1)
                && adjCostMap.get(edge1)
                .containsKey(edge2))) {
                return false;
            }
            if (j == 3) {
                // try the combination if passed other checks
                edge1 = getVertexIndex(tIdx1);
                edge2 = getVertexIndex(tIdx2);
                if (!(adjCostMap.containsKey(edge1)
                    && adjCostMap.get(edge1)
                    .containsKey(edge2))) {
                    return false;
                }
                edge1 = getVertexIndex(tIdx3);
                edge2 = getVertexIndex(tIdx4);
                if (!(adjCostMap.containsKey(edge1)
                    && adjCostMap.get(edge1)
                    .containsKey(edge2))) {
                    return false;
                }
            }
        }
        return true;
    }
    
    private boolean assertTourData() {
        
        assert(tourValueToIndexMap.size() == tour.length - 1);
        
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx = tour[i];
            if (!tourValueToIndexMap.containsKey(cIdx)) {
                return false;
            }
            int tIdx = tourValueToIndexMap.get(cIdx);
            if (tIdx != i) {
                return false;
            }
        }
        
        return true;
    }
    
    private boolean assertSameSets(int t1, int t2, int t3, int t4,
        int s1, int s2, int s3, int s4) {
        TIntSet set1 = new TIntHashSet();
        TIntSet set2 = new TIntHashSet();
        
        set1.add(t1);
        set1.add(t2);
        set1.add(t3);
        set1.add(t4);
        
        set2.add(s1);
        set2.add(s2);
        set2.add(s3);
        set2.add(s4);
        
        if (set1.size() != set2.size()) {
            return false;
        }
        TIntIterator iter = set1.iterator();
        while (iter.hasNext()) {
            int t = iter.next();
            if (!set2.contains(t)) {
                return false;
            }
        }
        return true;
    }
}
