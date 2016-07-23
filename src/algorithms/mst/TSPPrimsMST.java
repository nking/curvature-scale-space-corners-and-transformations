package algorithms.mst;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairInt;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.Arrays;
import java.util.List;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTreeInterval2D;

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
     * The approximate TSP tour calculated from given vertexes
     * and adjacency map is refined to remove crossing 
     * edges.  
     * 
     * NOTE: user must ensure that the range of keys
     * in the adjCostMap is between 0 and nVertexes - 1.
     *  
     * @param nVertexes
     * @param coordinates
     * @param adjCostMap key = vertex index1, 
     *   value = map with key = index2 and value = 
     *   cost for edge index1 to index2. 
     */
    public int[] approxTSPTour(
        final int nVertexes, PairInt[] coordinates,
        final TIntObjectMap<TIntIntMap> adjCostMap) {
        
        int[] tour = approxTSPTour(nVertexes, adjCostMap);
     
        // uncross edges where feasible
        
        TourHandler tourHandler = new TourHandler(
            tour, adjCostMap, coordinates);
                
        int nIter = 0;
        int nMaxIter = 10;
        int nChanged = 0;
        
        do {
            
            nChanged = 0;
            
            int bestVertexIdxA = -1;
            int bestVertexIdxB = -1;
            int[] bestVertexIdxs1 = new int[4];
            int minPathSum = tourHandler.getPathSum();
            int[] outputVertexIdxs = new int[4];
            
            int[] outputVertexB = new int[1];
            
            //TODO: consider on nIter=0, making a list
            // of vertex indexes which have intersecting
            // lines, and then after first iteration,
            // only iterate over those as cIdx1.
            
            for (int cIdx1 = 0; cIdx1 < coordinates.length;
                ++cIdx1) {

                int sum = 
                    tourHandler.findNonIntersectingBestSwap(
                    cIdx1, outputVertexB, 
                    outputVertexIdxs);

                if (sum < minPathSum) {
                    minPathSum = sum;
                    bestVertexIdxA = cIdx1;
                    bestVertexIdxB = outputVertexB[0];
                    System.arraycopy(outputVertexIdxs, 
                        0, bestVertexIdxs1, 0, 
                        outputVertexIdxs.length);
                }
            }
            
            if (minPathSum < Integer.MAX_VALUE) {
            //if (minPathSum < tourHandler.getPathSum()) {
                
                int sum = tourHandler.changePaths(
                    bestVertexIdxA, bestVertexIdxB, 
                    bestVertexIdxs1);
                System.out.println("sum=" + sum + " min=" + minPathSum);
                assert(sum == minPathSum);
            
                //TODO: 
                // - handle changes to quadtree
                // - handle changes to edge maps
                
                nChanged++;
            }
        
            nIter++;
            
        } while ((nChanged != 0) && (nIter < nMaxIter));
        
        return tour;
    }

    /**
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
                
        int[] preorderWalk = prims.getPreOrderWalkOfTree();
        
        int[] tour = Arrays.copyOf(preorderWalk,
            preorderWalk.length + 1);
        
        tour[tour.length - 1] = tour[0];
    
        return tour;
    }

    protected int distance(int x1, int y1, int x2, int y2) {
        int diffX = x1 - x2;
        int diffY = y1 - y2;
        return diffX * diffX + diffY * diffY;
    }    
}
