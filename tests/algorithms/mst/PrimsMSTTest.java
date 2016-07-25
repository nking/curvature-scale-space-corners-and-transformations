package algorithms.mst;

import algorithms.compGeometry.convexHull.PolarAngleQuickSort;
import algorithms.util.PairInt;
import algorithms.util.ScatterPointPlotterPNG;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PrimsMSTTest extends TestCase {
    
    public PrimsMSTTest() {
    }
    
    /**
     * test created from Cormen et al. Chap 24, MST, Fig 23.4.
     * @throws Exception 
     */
    public void est0() throws Exception {
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4    |          2     \         |    9
         *  [a]       11      [i]         4     14      [e]
         *       8    |     7                    |   10
         *         \ [h] / -- 1  --[g]-- 2 -- \ [f]/
         *   
         * letter  index        
         *    a      0
         *    b      1
         *    c      2
         *    d      3
         *    e      4
         *    f      5
         *    g      6
         *    h      7
         *    i      8
         */

        int nVertexes = 9;
        TIntObjectMap<TIntIntMap> adjCostMap 
            = new TIntObjectHashMap<TIntIntMap>();
        
        TIntIntMap map = new TIntIntHashMap();
        adjCostMap.put(0, map);
        map.put(1, 4);
        map.put(7, 8);
        
        map = new TIntIntHashMap();
        adjCostMap.put(1, map);
        map.put(0, 4);
        map.put(7, 11);
        
        map = new TIntIntHashMap();
        adjCostMap.put(2, map);
        map.put(3, 7);
        map.put(5, 4);
        map.put(8, 2);
        
        map = new TIntIntHashMap();
        adjCostMap.put(3, map);
        map.put(2, 7);
        map.put(4, 9);
        map.put(5, 14);
        
        map = new TIntIntHashMap();
        adjCostMap.put(4, map);
        map.put(5, 10);
        map.put(3, 9);
        
        map = new TIntIntHashMap();
        adjCostMap.put(5, map);
        map.put(4, 10);
        map.put(2, 4);
        map.put(6, 2);
        
        map = new TIntIntHashMap();
        adjCostMap.put(6, map);
        map.put(5, 2);
        map.put(7, 1);
        
        map = new TIntIntHashMap();
        adjCostMap.put(7, map);
        map.put(0, 8);        
        map.put(1, 11);
        map.put(8, 7);
        map.put(6, 1);
      
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(
            nVertexes, adjCostMap);                
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4               2     \              9
         *  [a]              [i]          4            [e]
         *       8                                  
         *         \ [h]   -- 1  --[g]-- 2 -- \ [f]
         */
        /*
        0 1  
        2 3  
        2 5  
        2 8  
        3 4  
        5 6  
        6 7  
        7 0  
        */
        
        int[] predecessorArray = prims.getPrecessorArray();
        assertTrue(Arrays.equals(
            new int[]{-1, 0, 5, 2, 3, 6, 7, 0, 2}, 
            predecessorArray));
        
        int[] treeWalk = prims.getPreOrderWalkOfTree();
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4               2     \              9
         *  [a]              [i]          4            [e]
         *       8                                  
         *         \ [h]   -- 1  --[g]-- 2 -- \ [f]
        */
        
        // traverse the tree to find these from root=0 
        /*
        0 1, 0 7
        7 6
        6 5
        5 2
        2 8, 2 3
        3 4 
        */
        
    }

    public void est1() throws Exception {
        
        /*         / [b]           /[c]-- 7  -- [d]\
         *       4    |          2     \         |    9
         *  [a]       11      [i]         4     14      [e]
         *       8    |     7                    |   10
         *         \ [h] / -- 1  --[g]-- 2 -- \ [f]/
         *   
         * letter  index        
         *    a      0
         *    b      1
         *    c      2
         *    d      3
         *    e      4
         *    f      5
         *    g      6
         *    h      7
         *    i      8
        
         *         / [1]           /[2]-- 7  -- [3]\
         *       4    |          2     \         |    9
         *  [0]       11      [8]         4     14      [4]
         *       8    |     7                    |   10
         *         \ [7] / -- 1  --[6]-- 2 -- \ [5]/
         */

        int nVertexes = 9;
        TIntObjectMap<TIntIntMap> adjCostMap 
            = new TIntObjectHashMap<TIntIntMap>();
        
        TIntIntMap map = new TIntIntHashMap();
        adjCostMap.put(0, map);
        map.put(1, 4);
        map.put(7, 8);
        
        map = new TIntIntHashMap();
        adjCostMap.put(2, map);
        map.put(3, 7);
        map.put(8, 2);
        
        map = new TIntIntHashMap();
        adjCostMap.put(3, map);
        map.put(4, 9);
        
        map = new TIntIntHashMap();
        adjCostMap.put(5, map);
        map.put(2, 4);
        
        map = new TIntIntHashMap();
        adjCostMap.put(6, map);
        map.put(5, 2);
        
        map = new TIntIntHashMap();
        adjCostMap.put(7, map);
        map.put(6, 1);
      
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(
            nVertexes, adjCostMap);                
        
        /*         / [1]           /[2]-- 7  -- [3]\
         *       4               2     \              9
         *  [0]              [8]          4            [4]
         *       8
         *         \ [7]   -- 1  --[6]-- 2 -- \ [5]
        
        [-1, 0, 5, 2, 3, 6, 7, 0, 2]
             *  *  *  *  *  *  *  
         */
     
        int[] predecessorArray = prims.getPrecessorArray();
       
        assertTrue(Arrays.equals(
            new int[]{-1, 0, 5, 2, 3, 6, 7, 0, 2}, 
            predecessorArray));
                
    }
    
    public void est2() throws Exception {
        
        /*
        testing a pre-order, post-order merged traversal.        
        */
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        sr.setSeed(seed);
        
        List<PairInt> points = new ArrayList<PairInt>();
        Map<PairInt, Integer> pointMap = new HashMap<PairInt, Integer>();
        
        populateWithTest2Data(points, pointMap, true);
    
        // ------ default order ----
        TIntObjectMap<TIntIntMap> adjCostMap = 
            createAdjacencyMap(points);
        
        PrimsMST mst = new PrimsMST();
        mst.calculateMinimumSpanningTree(points.size(), 
            adjCostMap);
        
        int[] walk = mst.getPreOrderPostOrderWalk();
         
        print(walk, points, pointMap);
    
        // ------ shuffle the order excepting node 0 --
        for (int i = 0; i < 5; ++i) {
            
            shuffle(sr, points);
                        
            adjCostMap = createAdjacencyMap(points);
        
            mst = new PrimsMST();
            mst.calculateMinimumSpanningTree(points.size(), 
                adjCostMap);
        
            walk = mst.getPreOrderPostOrderWalk();
        
            print(walk, points, pointMap);
        }        
    }
    
    public void test3() throws Exception {
        
        /*
        testing a pre-order, post-order merged traversal.
        
        testing whether point order changes need to be made.
        */
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        sr.setSeed(seed);
        
        List<PairInt> points = new ArrayList<PairInt>();
        Map<PairInt, Integer> pointMap = new HashMap<PairInt, Integer>();
        
        populateWithRandomData(sr, points, pointMap, true);
    
        // ------ default order ----
        TIntObjectMap<TIntIntMap> adjCostMap = 
            createAdjacencyMap(points);
        
        PrimsMST mst = new PrimsMST();
        mst.calculateMinimumSpanningTree(points.size(), 
            adjCostMap);
        
        int[] walk = mst.getPreOrderPostOrderWalk();
        
        print2(walk, points, pointMap, 0);
       
        // ------ shuffle the order excepting node 0 --
        for (int i = 0; i < 5; ++i) {
            
            shuffle(sr, points);
        
            Collections.sort(points, new PointComparator());
            
            adjCostMap = createAdjacencyMap(points);
        
            mst = new PrimsMST();
            mst.calculateMinimumSpanningTree(points.size(), 
                adjCostMap);
        
            walk = mst.getPreOrderPostOrderWalk();
        
            print2(walk, points, pointMap, i + 1);
        }        
    }
    
    private void populateWithTest2Data(List<PairInt> points,
        Map<PairInt, Integer> pointMap, boolean doSort) {
        
        /*
         * From Cormen et al. chap 35 Approximation Algorithms, 35.2 Traveling-Salesman problem, Fig 35.2
         * 
         * 
         *   6  |---|---|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   5  |---|--[A]0-|--[D]3-|---|---|----
         *      |   |   |   |   |   |   |   |
         *   4  |---|---|---|---|--[E]4-|---|----
         *      |   |   |   |   |   |   |   |
         *   3  |---|--[B]1-|--[F]5-|--[G]6-|----
         *      |   |   |   |   |   |   |   |
         *   2  |--[C]2-|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   1  |---|---|--[H]7-|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   0  |---|---|---|---|---|-------|----
         *      0   1   2   3   4   5   6   7
        
        Prims for default point ordering:
                    0
               1        3
            2         4
                    5  6
                  7
        */
       
        points.add(new PairInt(2, 5));      
        points.add(new PairInt(2, 3));        
        points.add(new PairInt(1, 2));        
        points.add(new PairInt(4, 5));        
        points.add(new PairInt(5, 4));        
        points.add(new PairInt(4, 3));        
        points.add(new PairInt(6, 3));        
        points.add(new PairInt(3, 1));

        if (doSort) {
            Collections.sort(points,
                new PointComparator());
        }
        
        populateMap(points, pointMap);
    }
    
    private void populateWithRandomData(
        SecureRandom sr, List<PairInt> points,
        Map<PairInt, Integer> pointMap, boolean doSort) 
        throws Exception {
                
        int min = 0;
        int max = 100;
        int n = 100;
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        // set first point 
        {
            PairInt p = new PairInt(0, 100);
            points.add(p);
            added.add(p);
        }
        
        for (int i = 1; i < n; ++i) {
            int x = sr.nextInt(max);
            int y = sr.nextInt(max);
            PairInt p = new PairInt(x, y);
            while (added.contains(p)) {
                x = sr.nextInt(max);
                y = sr.nextInt(max);
                p = new PairInt(x, y);
            }
            points.add(p);
            added.add(p);
        }
        
        if (doSort) {
            Collections.sort(points, new PointComparator());
        }
        
        populateMap(points, pointMap);
        
    }
    
    private void populateMap(
        List<PairInt> points,
        Map<PairInt, Integer> pointMap) {
        
        for (int i = 0; i < points.size(); ++i) {
            PairInt p = points.get(i);
            pointMap.put(p, Integer.valueOf(i));
        }
    }
    
    private TIntObjectMap<TIntIntMap> 
        createAdjacencyMap(List<PairInt> points) {
        
        TIntObjectMap<TIntIntMap>
            adjCostMap = new TIntObjectHashMap<TIntIntMap>();
        
        for (int i = 0; i < points.size(); ++i) {
            int x1 = points.get(i).getX();
            int y1 = points.get(i).getY();
            TIntIntMap map1 = adjCostMap.get(i);            
            if (map1 == null) {
                map1 = new TIntIntHashMap();
                adjCostMap.put(i, map1);
            }
            
            for (int j = 0; j < points.size(); ++j) {
                if (i == j) {
                    continue;
                }
                int x2 = points.get(j).getX();
                int y2 = points.get(j).getY();                
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                int dist = diffX * diffX + diffY * diffY;
                map1.put(j, dist);
            }
        }
        
        return adjCostMap;
    }
        
    public void shuffle(SecureRandom sr, 
        List<PairInt> points) throws NoSuchAlgorithmException {
        
        // shuffle all except first point
                
        for (int i = (points.size() - 1); i > 1; --i) {
            int j = sr.nextInt(i + 1);
            while (j == 0) {
                // avoiding change of points[0]
                j = sr.nextInt(i + 1);
            }
            PairInt swap = points.get(j);
            points.set(j, points.get(i));
            points.set(i, swap);
        }
    }

    private void print(int[] walk, List<PairInt> points,
        Map<PairInt, Integer> pointMap) {

        StringBuilder sb = new StringBuilder("");
        for (int v : walk) {
            PairInt p = points.get(v);
            String r = pointMap.get(p).toString();
            sb.append(r).append(",");
        }
        System.out.println(sb.toString());
    }

    private void print2(int[] walk, List<PairInt> points, 
        Map<PairInt, Integer> pointMap, int plotNumber) 
        throws Exception {
      
        ScatterPointPlotterPNG plotter = new
            ScatterPointPlotterPNG();
        
        float[] xPoints = new float[points.size()];
        float[] yPoints = new float[points.size()];
        for (int i = 0; i < points.size(); ++i) {
            PairInt p = points.get(i);
            int idx = pointMap.get(p);
            xPoints[i] = points.get(idx).getX();
            yPoints[i] = points.get(idx).getY();
        }
        
        plotter.plotLabeledPoints(0, 100, 0, 100, 
            xPoints, yPoints, 
            Integer.toString(plotNumber), 
            "X", "Y");
    
        plotter.writeFile(plotNumber);
        
        print(walk, points, pointMap);
    }
    
    private void sort(PairInt p0, List<PairInt> points) {
        PolarAngleQuickSort.sort2(p0, points);
    }
    
    public class PointComparator implements Comparator<PairInt> {
        @Override
        public int compare(PairInt o1, PairInt o2) {
            // sort by decr y then incr x
            int c = Integer.compare(o2.getY(), o1.getY());
            if (c != 0) {
                return c;
            }
            return Integer.compare(o1.getX(), o2.getX());
        }
    }
    
}
