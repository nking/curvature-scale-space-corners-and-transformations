package algorithms.mst;

import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPPrimsMSTTest extends TestCase {
   
    public void testTourHandler() {
            
        /* Test of traveling salesman approximate tour. 
         * 
         * From Cormen et al. chap 35 Approximation Algorithms, 35.2 Traveling-Salesman problem, Fig 35.2
         * 
         * 
         *   6  |---|---|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   5  |---|--[A]--|--[D]--|---|---|----
         *      |   |   |   |   |   |   |   |
         *   4  |---|---|---|---|--[E]--|---|----
         *      |   |   |   |   |   |   |   |
         *   3  |---|--[B]--|--[F]--|--[G]--|----
         *      |   |   |   |   |   |   |   |
         *   2  |--[C]--|---|---|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   1  |---|---|--[H]--|---|---|---|----
         *      |   |   |   |   |   |   |   |
         *   0  |---|---|---|---|---|-------|----
         *      0   1   2   3   4   5   6   7
        
         TSP-approx:
        
          L is pre-order walk of mst prims tree
        H is the tour of those visits in that order

        an MST is
                  A0
            B1         D3
           C2 H7       E4
                      F5 G6
        
        Prims here is giving:
                    0
               1        3
            2         4
                    5  6
                  7

        pre-order is
        root, left subtree, right subtree
        // given the top node as the starter
        //level 0            [0]
        //level 1     [1]           [4]
        //level 2   [2] [3]       [5] [6]
        // process node sees 0,1,2,3,4,5,6

        then pre-order trversal as a tour:
        a, b, c, b, h, b, a, d, e, f, e, g, e, d, a
        removing vertexes already encountered:
        a, b, c,  h, d, e, f, g,

        a pre-order for prims mst here is:
                   0
               1        3
            2         4
                    5  6
                  7
        
         0, 1, 2, 3, 4, 5, 7, 6
        
        In contrast, an optimal cost, non-crossing:
        a, b, c,  h, f, g, e, d
        0  1  2   7  5  6  4  3
        */
                
        PairInt[] points = new PairInt[8];
        points[0] = new PairInt(2, 5);
        points[1] = new PairInt(2, 3);
        points[2] = new PairInt(1, 2);
        points[3] = new PairInt(4, 5);
        points[4] = new PairInt(5, 4);
        points[5] = new PairInt(4, 3);
        points[6] = new PairInt(6, 3);
        points[7] = new PairInt(3, 1);
        
        int[] expected = new int[9];
        expected[0] = 0;
        expected[1] = 1;
        expected[2] = 2;
        expected[3] = 7;
        expected[4] = 3;
        expected[5] = 4;
        expected[6] = 5;
        expected[7] = 6;
        expected[8] = 0;
        /*A B C D E F G H
          0 1 2 3 4 5 6 7
        */
         
        TIntObjectMap<TIntIntMap>
            adjCostMap = new TIntObjectHashMap<TIntIntMap>();
        
        for (int i = 0; i < points.length; ++i) {
            int x1 = points[i].getX();
            int y1 = points[i].getY();
            TIntIntMap map1 = adjCostMap.get(i);            
            if (map1 == null) {
                map1 = new TIntIntHashMap();
                adjCostMap.put(i, map1);
            }
            
            for (int j = 0; j < points.length; ++j) {
                if (i == j) {
                    continue;
                }
                int x2 = points[j].getX();
                int y2 = points[j].getY();
                
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                int dist = diffX * diffX + diffY * diffY;
//System.out.println("i=" + i + " j=" + j + "  dist=" + dist);                
                map1.put(j, dist);
            }
        }
        
        /* a pre-order for prims mst here is:
                   0
               1        3
            2         4
                    5  6
                  7
        
         0, 1, 2, 3, 4, 5, 7, 6
        
               2, 7        6  3  
        */
        
        int[] tour = new int[]{0, 1, 2, 3, 4, 5, 7, 6, 0};
        
    }
    
    public void test0() {
        
        /* Test of traveling salesman approximate tour. 
         * 
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
        
         TSP-approx:
        
          L is pre-order walk of mst prims tree
        H is the tour of those visits in that order

        an MST is
                  A0
            B1         D3
           C2 H7       E4
                      F5 G6
        
        Prims here is giving:
                    0
               1        3
            2         4
                    5  6
                  7

        pre-order is
        root, left subtree, right subtree
        // given the top node as the starter
        //level 0            [0]
        //level 1     [1]           [4]
        //level 2   [2] [3]       [5] [6]
        // process node sees 0,1,2,3,4,5,6

        then pre-order trversal as a tour:
        a, b, c, b, h, b, a, d, e, f, e, g, e, d, a
        removing vertexes already encountered:
        a, b, c,  h, d, e, f, g,

        a pre-order for prims mst here is:
                   0
               1        3
            2         4
                    5  6
                  7
        
         0, 1, 2, 3, 4, 5, 7, 6
        
        In contrast, an optimal cost, non-crossing:
        a, b, c,  h, f, g, e, d
        0  1  2   7  5  6  4  3
        */
                
        PairInt[] points = new PairInt[8];
        points[0] = new PairInt(2, 5);
        points[1] = new PairInt(2, 3);
        points[2] = new PairInt(1, 2);
        points[3] = new PairInt(4, 5);
        points[4] = new PairInt(5, 4);
        points[5] = new PairInt(4, 3);
        points[6] = new PairInt(6, 3);
        points[7] = new PairInt(3, 1);
        
        int[] expected = new int[9];
        expected[0] = 0;
        expected[1] = 1;
        expected[2] = 2;
        expected[3] = 7;
        expected[4] = 5;
        expected[5] = 6;
        expected[6] = 4;
        expected[7] = 3;
        expected[8] = 0;
        /*A B C D E F G H
          0 1 2 3 4 5 6 7
        */
         
        TIntObjectMap<TIntIntMap>
            adjCostMap = new TIntObjectHashMap<TIntIntMap>();
        
        for (int i = 0; i < points.length; ++i) {
            int x1 = points[i].getX();
            int y1 = points[i].getY();
            TIntIntMap map1 = adjCostMap.get(i);            
            if (map1 == null) {
                map1 = new TIntIntHashMap();
                adjCostMap.put(i, map1);
            }
            
            for (int j = 0; j < points.length; ++j) {
                if (i == j) {
                    continue;
                }
                int x2 = points[j].getX();
                int y2 = points[j].getY();
                
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                int dist = diffX * diffX + diffY * diffY;
//System.out.println("i=" + i + " j=" + j + "  dist=" + dist);                
                map1.put(j, dist);
            }
        }
        
        
        int maxCost = PrimsMST.maxEdgeCost(adjCostMap);
        PrimsMST prims = new PrimsMST();
        prims.calculateMinimumSpanningTree(adjCostMap, maxCost);
        int[] walk = prims.getPreorderIndexes().toArray();
        System.out.println("prims walk=" +
            Arrays.toString(walk));
        assertTrue(Arrays.equals(
            new int[]{0, 1, 2, 3, 4, 5, 7, 6}, walk));
       
    }
    
    public void testATT48() throws Exception {
        
        List<PairInt> pointList = new ArrayList<PairInt>();
        TIntList expectedBestTour = new TIntArrayList();
        TIntObjectMap<TIntIntMap>
            adjCostMap = new TIntObjectHashMap<TIntIntMap>();
        
        populateATT48(pointList, expectedBestTour,
            adjCostMap);
        
        PairInt[] points = pointList.toArray(new PairInt[pointList.size()]);
        
        TSPPrimsMST tsp = new TSPPrimsMST();
        int[] tour = tsp.approxTSPTour(adjCostMap.size(), adjCostMap);
        
        /*ScatterPointPlotterPNG plotter = new
            ScatterPointPlotterPNG();
        */
        float[] xPoints = new float[pointList.size() + 1];
        float[] yPoints = new float[pointList.size() + 1];
        float[] xPointsE = new float[pointList.size() + 1];
        float[] yPointsE = new float[pointList.size() + 1];
        for (int i = 0; i < tour.length; ++i) {
            int vIdx = tour[i];
            PairInt p = pointList.get(vIdx);
            xPoints[i] = p.getX();
            yPoints[i] = p.getY();
            
            vIdx = expectedBestTour.get(i);
            p = pointList.get(vIdx);
            xPointsE[i] = p.getX();
            yPointsE[i] = p.getY();
        }
        
        /*
        plotter.plotLabeledPoints(0, 8000, 0, 5500, 
            xPoints, yPoints, 
            Integer.toString(100), "X", "Y");
    
        plotter.writeFile(100);
        */
        PolygonAndPointPlotter plotter2 = 
            new PolygonAndPointPlotter();
        plotter2.addPlot(xPoints, yPoints, 
            xPoints, yPoints, "tour Att48");
        plotter2.addPlot(xPoints, yPoints, 
            xPointsE, yPointsE, "BEST tour Att48");
        plotter2.writeFile(101);
        
        /*assertEquals(expectedBestTour.size(), tour.length);
        for (int i = 0; i < tour.length; ++i) {
            int tourIdx = tour[i];
            assertEquals(expectedBestTour.get(i), tourIdx);
        }*/
    }

    private void populateATT48(List<PairInt> pointList, 
        TIntList expectedBestTour, 
        TIntObjectMap<TIntIntMap> adjCostMap) 
        throws Exception {

        String dir = 
            ResourceFinder.findTestResourcesDirectory()
            + "/att48";
            
        readCoords(dir + "/att48_xy.txt", pointList);
    
        readBestTour(dir + "/att48_s.txt", 
            expectedBestTour);
        
        readAdjacency(dir + "/att48_d.txt", adjCostMap);
    }

    private void readCoords(String filePath, 
        List<PairInt> pointList) throws Exception {
        
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            while (line != null) {
                
                String item1 = line.substring(0, 4).trim();
                
                String item2 = line.substring(5).trim();
                
                Integer x =
                    Integer.parseInt(item1);
                
                Integer y =
                    Integer.parseInt(item2);
                
                PairInt p = new PairInt(x.intValue(), y.intValue());
                
                pointList.add(p);
                
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }   
    }

    private void readBestTour(String filePath, 
        TIntList expectedBestTour) throws Exception {
        
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            while (line != null) {
                
                Integer index = Integer.parseInt(line.trim());
                
                // change to zero based indexes
                expectedBestTour.add(
                    index.intValue() - 1);
                
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
    }

    private void readAdjacency(String filePath, 
        TIntObjectMap<TIntIntMap> adjCostMap) 
        throws Exception {

        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine();
            
            int i = 0;
            
            while (line != null) {
            
                TIntIntMap idx2CostMap = new TIntIntHashMap();
                adjCostMap.put(i, idx2CostMap);
            
                String[] items = line.split("\\s+");
                if (items[0] == null || items[0].equals("")) {
                    items = Arrays.copyOfRange(items, 1, 
                        items.length);
                }
                
                for (int j = 0; j < items.length; ++j) {
                    Integer c = Integer.parseInt(items[j].trim());
                    if (!c.equals(Integer.valueOf(0))) {
                        idx2CostMap.put(j, c.intValue());
                    }
                }
                
                line = in.readLine();
                ++i;
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
    }
}
