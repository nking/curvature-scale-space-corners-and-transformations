package algorithms.mst;

import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
     * NOTE: user must ensure that the range of keys
     * in the adjCostMap is between 0 and nVertexes - 1.
     * 
     * @param nVertexes
     * @param adjCostMap key=vertex index, 
     *   value=set of pairint where each pairint has 
     *   x = adjacent vertex and y = cost of edge.
     */
    public int[] approxTSPTour(
        final int nVertexes, final TIntObjectMap<Set<PairInt>>
            adjCostMap) {
        
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
       
        //TODO: here is where need to optimize the tour,
        // making the edges non-crossing.
        
        /*
        crossed edges:
        
        The alternative connecting pairs from two crossing 
        segments have shorter separation than the 
        crossed.
        
        need a way to find them that is faster
        than O(N^2) (for dense graphs) though...
        can probably use a 2d tree such as range
        tree for closer to O(N*log_2(N))
        */
        
        tour[tour.length - 1] = tour[0];
    
        return tour;
    }
    
}
