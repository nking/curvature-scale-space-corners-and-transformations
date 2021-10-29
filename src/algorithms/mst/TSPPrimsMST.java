package algorithms.mst;

import gnu.trove.list.TIntList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;

/**
 * adapted from Cormen et al. Introduction to Algorithms
 *
 * An approximate solution to the traveling salesman problem.
 * The algorithms provides solution which is at most, twice the cost of an
 * optimal tour.
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
 * Closer to optimal in:
 * https://en.wikipedia.org/wiki/Travelling_salesman_problem
 * 
 * For some datasets, especially large and clustered, it's likely a common
 * practice to solve TSP within the clusters, then between those clusters...?
 
 * Exact dynamic programming TSP can be found in the shared project
 * src/test/java/algorithms/exponential/TSPDynamic.java
 * 
 * @author nichole
 */
public class TSPPrimsMST {
     
    /**
     * 
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
        
        int maxCost = PrimsMST.maxEdgeCost(adjCostMap);
        
        PrimsMST prims = new PrimsMST();
        
        prims.calculateMinimumSpanningTree(adjCostMap, maxCost);

        TIntList treeWalk = prims.getPreorderIndexes();
        
        // the tour is the hamiltonian cycle H that visits the vertices in the order L.
        // where Hamiltonian, is a simple cycle that includes all vertices.
        treeWalk.add(treeWalk.get(0));
        
        int[] tour = treeWalk.toArray();
    
        return tour;
    }

    protected int distance(int x1, int y1, int x2, int y2) {
        int diffX = x1 - x2;
        int diffY = y1 - y2;
        return diffX * diffX + diffY * diffY;
    }    
    
}
