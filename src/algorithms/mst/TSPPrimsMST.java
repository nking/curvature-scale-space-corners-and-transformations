package algorithms.mst;

import algorithms.QuickSort;
import algorithms.Rotate;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.Arrays;

/**
 * adapted from Cormen et al. Introduction to Algorithms
 *
 * An approximate solution to the traveling salesman problem.
 *
 * Input:
 *    -- G = (V, E) is a complete undirected graph
 *    -- each edge (u, v) in E has a positive floating point 
 *       number cost c(u, v)
 * 
 * Find a Hamiltonian cycle (= a tour of G, that is each node visited exactly
 * once) w/ minimum cost.  Below, using
 * notation that a subset of the edges called A will be associated with c(A) cost.
  
 * Uses a minimum spanning tree:
 * Minimum spanning tree is the minimal network that spans all nodes in a tree
 * and has the smallest cost (sum of edges).
 *
 * The MST is implemented from pseudo code in Cormen et al. Introduction to Algorithms and
 * from http://en.wikipedia.org/wiki/Prim's_algorithm.
   Useful also was
  http://www.i-programmer.info/projects/61-algorithms/534-minimum-spanning-tree.html?start=1
  
       
 * @author nichole
 */
public class TSPPrimsMST {
       
    /**
     * NOT READY FOR USe YET
     * 
     * The approximate TSP tour calculated from given vertexes
     * and adjacency map is refined to remove crossing 
     * edges.  
     * 
     * NOTE: user must ensure that the range of keys
     * in the adjCostMap is between 0 and nVertexes - 1.
     *  
     * @param coordinates
     * @param adjCostMap key = vertex index1, 
     *   value = map with key = index2 and value = 
     *   cost for edge index1 to index2. 
     */
    public int[] approxTSPTour(
        PairInt[] coordinates,
        final TIntObjectMap<TIntIntMap> adjCostMap, 
        boolean doSort) {
        
        final int nVertexes = coordinates.length;
        
        if (!doSort) {
            return approxTSPTour(coordinates, adjCostMap);
        }
        
        PairInt[] coordinates2 = Arrays.copyOf(coordinates,
            coordinates.length);
        int[] indexes = new int[coordinates.length];
        for (int i = 0; i < coordinates.length; ++i) {
            indexes[i] = i;
        }
        QuickSort.sortByDecrYThenIncrX(coordinates2, indexes);
      
        int[] revIndexes = new int[indexes.length];
        for (int i = 0; i < indexes.length; ++i) {
            int oIdx = indexes[i];
            revIndexes[oIdx] = i;
        }
        
        TIntObjectMap<TIntIntMap> adjCostMap2 = 
            new TIntObjectHashMap<TIntIntMap>();
        TIntObjectIterator<TIntIntMap> iter = adjCostMap.iterator();
        for (int i = 0; i < adjCostMap.size(); ++i) {
            iter.advance();
            int idx1 = iter.key();
            TIntIntMap map1 = iter.value();
            TIntIntIterator iter2 = map1.iterator();
            
            int idx2 = revIndexes[idx1];
            TIntIntMap map2 = new TIntIntHashMap();
            adjCostMap2.put(idx2, map2);
            for (int j = 0; j < map1.size(); ++j) {
                iter2.advance();
                int oIdx2 = iter2.key();
                int cost = iter2.value();
                map2.put(revIndexes[oIdx2], cost);
            }
        }
        
        int[] tour = approxTSPTour(coordinates2, 
            adjCostMap2);
        
        // transform the indexes back to original reference frame
        for (int i = 0; i < tour.length; ++i) {
            int vIdx = tour[i];
            tour[i] = indexes[vIdx];
        }
        
        // rotate 0 back to position 0
        int zIdx = -1;
        for (int i = 0; i < tour.length; ++i) {
            if (tour[i] == 0) {
                zIdx = i;
                break;
            }
        }
        assert(zIdx > -1);
        
        Rotate r = new Rotate();
        r.rotate(tour, zIdx);
        
        int lastVOffset = tour.length - 1 - zIdx;
        // shift items below vOffset up by one
        for (int i = (lastVOffset + 1); i < tour.length; ++i) {
            tour[i - 1] = tour[i];
        }
        tour[tour.length - 1] = tour[0];
        
        System.out.println("refined tour (in orig ref frame)=" 
            + Arrays.toString(tour));
        
        return tour;
    }
    
    /**
     * 
     * @param nVertexes
     * @param coordinates
     * @param adjCostMap
     * @param doSort
     * @return 
     */
    private int[] approxTSPTour(PairInt[] coordinates,
        final TIntObjectMap<TIntIntMap> adjCostMap) {
        
        final int nVertexes = coordinates.length;
        
        int[] tour = approxTSPTour(nVertexes, adjCostMap);
     
        System.out.println("approx tour=" + Arrays.toString(tour));
    
        //TODO: because sorting and MST are N log N which is
        // fast for this algorithm, could experiment
        // with different sorting methods quickly and return
        // the best results.
        // OR experiment with modifying the pre and post order
        // MST merge of walks to also use a nearest neighbors to avoid large
        // gaps.
        
        // uncross edges where feasible
        
        TourHandler tourHandler = new TourHandler(
            tour, adjCostMap, coordinates);
        
        tourHandler.modifyTourIntersectingEdges();
        
        return tour;
    }

    /**
     * 
     * NOT READY FOr USe YET
     * The approximate TSP tour created from Prim's Minimum Spanning Tree is
     * returned. Note that the result may contain crossing edges, hence not
     * optimal.
     *
     * NOTE: user must ensure that the range of keys in the adjCostMap is
     * between 0 and nVertexes - 1.
     *
     * @param nVertexes
     * @param adjCostMap key = vertex index1, 
     *   value = map with key = index2 and value = 
     *   cost for edge index1 to index2. 
     */
    public int[] approxTSPTour(
        final int nVertexes,
        final TIntObjectMap<TIntIntMap> adjCostMap) {
        
        /* Approx TSP-Tour(G, c) {
         *     -- select a vertex r in V[G] as the 'root' vertex
         *     -- compute a minimum spanning tree T for G from root r using MST-PRIM(G, c, r)
         *     -- let L be the list of vertices visited in a preorder tree walk of T
         *     -- return the hamiltonian cycle H that visits the vertices in the order L 
         */
        
        PrimsMST prims = new PrimsMST();
        
        prims.calculateMinimumSpanningTree(nVertexes, 
            adjCostMap);

        int[] walk = prims.getPreOrderPostOrderWalk();
        
        int[] tour = Arrays.copyOf(walk, walk.length + 1);
        
        tour[tour.length - 1] = tour[0];
    
        return tour;
    }

    protected int distance(int x1, int y1, int x2, int y2) {
        int diffX = x1 - x2;
        int diffY = y1 - y2;
        return diffX * diffX + diffY * diffY;
    }    
}
