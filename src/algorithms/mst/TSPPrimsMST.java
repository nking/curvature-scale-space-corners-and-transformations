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
     
        /*
        crossed edges:
        
        The alternative connecting pairs from two crossing 
        segments have shorter separation than the 
        crossed.
        
        need a way to find them that is faster
        than O(N^2) (for dense graphs) so have
        added a quadtree for 2d intervals to the
        project.
        */
        
        TourHandler tourHandler = new TourHandler(
            tour, adjCostMap);
        
        QuadTreeInterval2D<Integer, PairInt> qt =
            new QuadTreeInterval2D<Integer, PairInt>();
        
        // bounding boxes of tour[i] to tour[i+1]
        // stored by keys that are the first value,
        // tour[i] which is the vertex index
        TIntObjectMap<Interval2D<Integer>> indexEdgeBounds =
            new TIntObjectHashMap<Interval2D<Integer>>();
        
        TObjectIntMap<Interval2D<Integer>> edgeIndexBounds =
            new TObjectIntHashMap<Interval2D<Integer>>();
        
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx1 = tourHandler.getVertexIndex(i);
            int x1 = coordinates[i].getX();
            int y1 = coordinates[i].getY();
            int x2, y2, cIdx2;
            if (i == (tour.length - 2)) {
                x2 = coordinates[0].getX();
                y2 = coordinates[0].getY();
                cIdx2 = tourHandler.getVertexIndex(0);
            } else {
                x2 = coordinates[i + 1].getX();
                y2 = coordinates[i + 1].getY();
                cIdx2 = tourHandler.getVertexIndex(i + 1);
            }
            Interval<Integer> xi1 = new Interval<Integer>(x1, x2);
            Interval<Integer> yi1 = new Interval<Integer>(y1, y2);
            Interval2D<Integer> box = new Interval2D<Integer>(xi1, yi1);

            PairInt desc = new PairInt(cIdx1, cIdx2);
            qt.insert(box, desc);
            
            edgeIndexBounds.put(box, cIdx1);
            indexEdgeBounds.put(cIdx1, box);
        }
        
        //TODO: implement a "delete" in qt
        
        int nIter = 0;
        int nMaxIter = 10;
        int nChanged = 0;
        
        do {
            
            nChanged = 0;
            
            int bestVertexIdxA = -1;
            int bestVertexIdxB = -1;
            int[] bestVertexIdxs1 = new int[4];
            int minPathSum = Integer.MAX_VALUE;
            int[] outputVertexIdxs = new int[4];
            
            //TODO: consider on nIter=0, making a list
            // of vertex indexes which have intersecting
            // lines, and then after first iteration,
            // only iterate over those as cIdx1.
            
            for (int cIdx1 = 0; cIdx1 < coordinates.length;
                ++cIdx1) {

                int tIdx1 = tourHandler.getTourIndex(cIdx1);
                int tIdx2 = (tIdx1 < (coordinates.length - 1)) ?
                    tIdx1 + 1 : 0;

                int cIdx2 = tourHandler.getVertexIndex(tIdx2);

                Interval2D<Integer> box12 = indexEdgeBounds.get(cIdx1);

                // find intersection boxes and look for 
                // edges that intersect with these

                List<Interval2D<Integer>> list = 
                    qt.query2D(box12);

                if (list.size() < 2) {
                    continue;
                }

                int x1 = coordinates[cIdx1].getX();
                int y1 = coordinates[cIdx1].getY();

                int x2 = coordinates[cIdx2].getX();
                int y2 = coordinates[cIdx2].getY();

                for (int listIdx = 0; listIdx < list.size();
                    ++listIdx) {

                    Interval2D<Integer> box34 = list.get(listIdx);

                    if (box34.equals(box12)) {
                        continue;
                    }

                    int cIdx3 = edgeIndexBounds.get(box34);

                    int tIdx3 = tourHandler.getTourIndex(cIdx3);

                    int tIdx4 = (tIdx3 < (coordinates.length - 1)) ?
                        tIdx3 + 1 : 0;

                    int cIdx4 = tourHandler.getVertexIndex(tIdx4);

                    int x3 = coordinates[cIdx3].getX();
                    int y3 = coordinates[cIdx3].getY();

                    int x4 = coordinates[cIdx4].getX();
                    int y4 = coordinates[cIdx4].getY();

                    // determine if p12 p34 intersect
                    if (!LinesAndAngles.linesIntersect(
                        x1, y1, x2, y2, x3, y3, x4, y4)) {
                        continue;
                    }

                    int sum = tourHandler.findNonIntersectingBestSwap(
                        cIdx1, cIdx3, outputVertexIdxs);
                    
                    if (sum < minPathSum) {
                        minPathSum = sum;
                        bestVertexIdxA = cIdx1;
                        bestVertexIdxB = cIdx3;
                        System.arraycopy(outputVertexIdxs, 
                            0, bestVertexIdxs1, 0, 
                            outputVertexIdxs.length);
                    }
                }
            }
            
            if (minPathSum < Integer.MAX_VALUE) {
                
                int sum = tourHandler.changePaths(
                    bestVertexIdxA, bestVertexIdxB, 
                    bestVertexIdxs1);
                
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
