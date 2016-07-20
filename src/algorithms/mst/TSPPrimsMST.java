package algorithms.mst;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairFloat;
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
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
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
     * @param adjCostMap key=vertex index, 
     *   value=set of pairint where each pairint has 
     *   x = adjacent vertex and y = cost of edge.
     */
    public int[] approxTSPTour(
        final int nVertexes, PairInt[] coordinates,
        final TIntObjectMap<Set<PairInt>> adjCostMap) {
        
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
        
        QuadTreeInterval2D<Integer, PairInt> qt =
            new QuadTreeInterval2D<Integer, PairInt>();
        
        // bounding boxes of tour[i] to tour[i+1]
        // stored by keys that are the first value,
        // tour[i] which is the vertex index
        TIntObjectMap<Interval2D<Integer>> indexEdgeBounds =
            new TIntObjectHashMap<Interval2D<Integer>>();
        
        TObjectIntMap<Interval2D<Integer>> edgeIndexBounds =
            new TObjectIntHashMap<Interval2D<Integer>>();
        
        // map with key = vertex indexes, that is the
        // values of tour[i] which are also the indexes
        // of the coordinates array.
        // map value = index of tour array.
        TIntIntMap tourValueToIndexMap = 
            new TIntIntHashMap();
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx = tour[i];
            tourValueToIndexMap.put(cIdx, i);
        }
        
        for (int i = 0; i < (tour.length - 1); ++i) {
            int cIdx1 = tour[i];
            int x1 = coordinates[i].getX();
            int y1 = coordinates[i].getY();
            int x2, y2, cIdx2;
            if (i == (tour.length - 2)) {
                x2 = coordinates[0].getX();
                y2 = coordinates[0].getY();
                cIdx2 = tour[0];
            } else {
                x2 = coordinates[i + 1].getX();
                y2 = coordinates[i + 1].getY();
                cIdx2 = tour[i + 1];
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
                
        Stack<Integer> stack = new Stack<Integer>();
        for (int i = 0; i < coordinates.length; ++i) {
            stack.add(Integer.valueOf(i));
        }
        
        TIntSet visited = new TIntHashSet();
        
        while (!stack.isEmpty()) {
            Integer cIndex1 = stack.pop();
            int cIdx1 = cIndex1.intValue();
            if (visited.contains(cIdx1)) {
                continue;
            }
            int tIdx1 = tourValueToIndexMap.get(cIdx1);
            
            int tIdx2 = (tIdx1 < (coordinates.length - 1)) ?
                tIdx1 + 1 : 0;
            
            int cIdx2 = tour[tIdx2];
            
            Interval2D<Integer> box12 = indexEdgeBounds.get(cIdx1);
            
            // find intersection boxes and look for 
            // edges that intersect with these
            
            List<Interval2D<Integer>> list = 
                qt.query2D(box12);
            
            if (list.size() < 2) {
                continue;
            }
            
            boolean found = false;
            
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
                
                int tIdx3 = tourValueToIndexMap.get(cIdx3);
            
                int tIdx4 = (tIdx3 < (coordinates.length - 1)) ?
                    tIdx3 + 1 : 0;
            
                int cIdx4 = tour[tIdx4];
                
                int x3 = coordinates[cIdx3].getX();
                int y3 = coordinates[cIdx3].getY();
            
                int x4 = coordinates[cIdx4].getX();
                int y4 = coordinates[cIdx4].getY();
                
                // determine if p12 p34 intersect
                if (!LinesAndAngles.linesIntersect(
                    x1, y1, x2, y2, x3, y3, x4, y4)) {
                    continue;
                }
                found = true;
                
                // TODO:
                // swap tour positions
                // and update associated data
                // --> need to implement qt.delete too
        
                
                /*    
                a to b and c to d is current edge pair
                
                find the shorter total and uncrossed
                among ac and bd, or ad and bc.
                
                - swap positions in tour        
                  use tour value map to find tour index
                  then swap tour values and update map
                - remove boxAB and boxCD from qr
                  remove boxAB from edgemap by cIdx1
                  remove boxCD from edgemap by cIdx3
                - create box AC and BD and add to 
                  qt and add to edgemap
                */
               
                if (true) {
                    throw new UnsupportedOperationException(
                    "not yet implemented");
                }
                
                // add the leading cIdxs back unto stack
                
                break;
            }
            
            if (!found) {
                visited.add(cIdx1);
            }
        }
        
        tour[tour.length - 1] = tour[0];

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
     * @param adjCostMap key=vertex index, value=set of pairint where each
     * pairint has x = adjacent vertex and y = cost of edge.
     */
    public int[] approxTSPTour(
        final int nVertexes, final TIntObjectMap<Set<PairInt>> adjCostMap) {

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
    
}
