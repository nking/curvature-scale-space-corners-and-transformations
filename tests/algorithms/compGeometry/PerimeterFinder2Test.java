package algorithms.compGeometry;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PerimeterFinder2Test extends TestCase {
    
    public void testExtractBorder() {
        
        /*
         4 5 6 7 8 9 
           @ @ @ @ %  5
         % @ @ @ @    4
           % % %      3
        */
        
        Set<PairInt> data = new HashSet<PairInt>();
        for (int i = 5; i < 10; ++i) {
            data.add(new PairInt(i, 5));
        }
        for (int i = 4; i < 9; ++i) {
            data.add(new PairInt(i, 4));
        }
        for (int i = 5; i < 8; ++i) {
            data.add(new PairInt(i, 3));
        }
        
        Set<PairInt> expected = new HashSet<PairInt>(data);
        expected.remove(new PairInt(6, 4));
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        Set<PairInt> boundary = finder.extractBorder(data);
        
        for (PairInt p : boundary) {
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
    }
    
    public void test1() {
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        PairIntArray expected = getSet1Boundaries();
        Set<PairInt> contiguousPoints = getSet1();
        
        Set<PairInt> expectedSet = Misc.convert(expected);
        Set<PairInt> resultUnordered = finder.extractBorder(
            contiguousPoints);
        assertEquals(expectedSet.size(), resultUnordered.size());
        for (PairInt p : resultUnordered) {
            assertTrue(expectedSet.remove(p));
        }
        assertEquals(0, expectedSet.size());

        // --- test medial azis points -- 
        expectedSet = new HashSet<PairInt>();
        expectedSet.add(new PairInt(1, 7));
        expectedSet.add(new PairInt(1, 8));
        expectedSet.add(new PairInt(2, 8));
        expectedSet.add(new PairInt(8, 7));
        expectedSet.add(new PairInt(8, 8));
        MedialAxis1 medAxis1 = new MedialAxis1(contiguousPoints, 
            resultUnordered);
        medAxis1.findMedialAxis();
        Set<PairInt> medAxisPts = medAxis1.getMedialAxisPoints();
        assertEquals(expectedSet.size(), medAxisPts.size());
        for (PairInt p : medAxisPts) {
            assertTrue(expectedSet.remove(p));
        }
        assertEquals(0, expectedSet.size());
        // ----
        
        contiguousPoints = getSet1();
        finder = new PerimeterFinder2();
        PairIntArray results = finder.extractOrderedBorder(
            contiguousPoints);
        
        System.out.println("results=" + results.toString());
        
        expected = new PairIntArray();
        expected.add(0, 3); expected.add(0, 4); expected.add(0, 5); 
        expected.add(0, 6); expected.add(0, 7); expected.add(0, 8);
        expected.add(0, 9); expected.add(1, 9); expected.add(2, 9);
        expected.add(3, 8); expected.add(4, 8); expected.add(5, 8); 
        expected.add(6, 9); expected.add(7, 9); expected.add(8, 9);
        expected.add(9, 8); expected.add(9, 7); 
        expected.add(9, 6); expected.add(9, 5); 
        expected.add(9, 4); expected.add(9, 3); expected.add(9, 2);
        expected.add(9, 1); expected.add(9, 0); 
        expected.add(8, 6); expected.add(7, 7); expected.add(7, 8); 
        expected.add(3, 7); expected.add(2, 6); expected.add(1, 6);

        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //    expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
        
    }
    
    /**
     * returns lists of sequential CW ordered contiguous
     * boundary points with concave and convex turns
     * and single pixel width spurs that lead to
     * pixels that only have one neighbor and
     * 2 pixel width regions (in which the.
     * 
     * @return 
     */
    private PairIntArray getSet1Boundaries() {
        
        /*
        //medial axis points: 1,7  1,8  2,8
        
            0 1 2 3 4 5 6 7 8 9
        12                @
        11              @ @
        10            @ @
         9  @ @ @ @ @ @ @ @ @ @
         8  @ * * @ @ @ @ @ * @
         7  @ * @ @       @ * @
         6  @ @ @         @ @ @
         5  @ @             @ @
         4  @               @ @
         3  @               @ @
         2                  @ @
         1                    @
         0                    @
        
            0 1 2 3 4 5 6 7 8 9
        */
        PairIntArray list0 = new PairIntArray();
            
        for (int i = 3; i <= 9; ++i) {
            list0.add(0, i);
        }
        for (int i = 1; i <= 4; ++i) {
            list0.add(i, 9);
        }
        list0.add(5, 9);
        list0.add(5, 10);
        list0.add(6, 11);
        list0.add(7, 12);
        list0.add(7, 11);
        list0.add(6, 10);
        for (int i = 6; i <= 9; ++i) {
            list0.add(i, 9);
        }
        for (int i = 8; i >= 0; --i) {
            list0.add(9, i);
        }
        for (int i = 2; i <= 6; ++i) {
            list0.add(8, i);
        }
        for (int i = 6; i <= 8; ++i) {
            list0.add(7, i);
        }
        for (int i = 6; i >= 3; --i) {
            list0.add(i, 8);
        }
        list0.add(3, 7);
        list0.add(2, 7);
        list0.add(2, 6);
        list0.add(1, 6);
        list0.add(1, 5);
                
        return list0;
    }
   
    /*
    //medial axis points: 1,7  1,8  2,8

        0 1 2 3 4 5 6 7 8 9
    12                @
    11              @ @
    10            @ @
     9  @ @ @ @ @ @ @ @ @ @
     8  @ * * @ @ @ @ @ * @ 
     7  @ * @ @       @ * @ 
     6  @ @ @         @ @ @
     5  @ @             @ @
     4  @               @ @
     3  @               @ @
     2                  @ @
     1                    @
     0                    @

        0 1 2 3 4 5 6 7 8 9
    */
    private Set<PairInt> getSet1() {
         
        Set<PairInt> set1 = Misc.convert(getSet1Boundaries());
        set1.add(new PairInt(1, 7));
        set1.add(new PairInt(1, 8));
        set1.add(new PairInt(2, 8));
        set1.add(new PairInt(8, 7));
        set1.add(new PairInt(8, 8));
        
        return set1;
    }
   
    public void testFindMinXY() {
        
        Set<PairInt> set = new HashSet<PairInt>();
        set.add(new PairInt(10, 10));
        set.add(new PairInt(10, 3));
        set.add(new PairInt(12, 8));
        set.add(new PairInt(18, 0));
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairInt xy = finder.findMinXY(set);
        
        assertEquals(10, xy.getX());
        assertEquals(3, xy.getY());
        
    }
    
    public void testOrdered0() {
        /*
        4
        3  *     *
        2     *  
        1  *     *
        0  1  2  3  4
        */
        
        Set<PairInt> medialAxis = new HashSet<PairInt>();
        medialAxis.add(new PairInt(1, 1));
        medialAxis.add(new PairInt(2, 2));
        medialAxis.add(new PairInt(1, 3));
        medialAxis.add(new PairInt(3, 1));
        medialAxis.add(new PairInt(3, 3));
        
        Set<PairInt> boundary = new HashSet<PairInt>();
        for (int i = 3; i >= 0; --i) {
            boundary.add(new PairInt(4, i));
        }
        for (int i = 3; i >= 1; --i) {
            boundary.add(new PairInt(i, 0));
        }
        for (int i = 0; i <= 4; ++i) {
            boundary.add(new PairInt(0, i));
        }
        for (int i = 1; i <= 4; ++i) {
            boundary.add(new PairInt(i, 4));
        }
                
        
        PairIntArray expected = new PairIntArray(boundary.size());
        for (int i = 0; i <= 4; ++i) {
            expected.add(0, i);
        }
        for (int i = 1; i <= 4; ++i) {
            expected.add(i, 4);
        }
        for (int i = 3; i >= 0; --i) {
            expected.add(4, i);
        }
        for (int i = 3; i >= 1; --i) {
            expected.add(i, 0);
        }
       /*
        4
        3  *     *
        2     *  
        1  *     *
        0  1  2  3  4
        
        20  21  22  23  24     4
        15  16* 17  18* 19     3
        10  11  12* 13  14     2
        5   6*  7   8*  9      1
        0   1   2   3   4      0 
        */
       Set<PairInt> shapePoints = new HashSet<PairInt>();
       for (int i = 0; i < 5; ++i) {
           for (int j = 0; j < 5; ++j) {
               shapePoints.add(new PairInt(i, j));
           }
       }
       
       PerimeterFinder2 finder = new PerimeterFinder2();
        
        PairIntArray results = finder.extractOrderedBorder(
            boundary, medialAxis, shapePoints);
        
        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //    expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
    }
    
    public void testOrdered2() {
        Set<PairInt> shapePoints = new HashSet<PairInt>();
        for (int i = 1; i < 7; ++i) {
            for (int j = 1; j < 4; ++j) {
                shapePoints.add(new PairInt(i, j));
            }
        }
        for (int i = 1; i < 4; ++i) {
            for (int j = 4; j < 7; ++j) {
                shapePoints.add(new PairInt(i, j));
            }
        }
        
        /*
        test w/ a concave section
        
        6  *  *  *
        5  *  .  *
        4  *  .  *
        3  *  .  *  *  *  *
        2  *  .  .  .  .  *
        1  *  *  *  *  *  *
        0  1  2  3  4  5  6
        

         5    35 36 3738 39 40 41
         4    28 29 3031 32 33 34
         3    21 22 2324 25 26 27
         2    14 15 1617 18 19 20 
         1    7  8  9 10 11 12 13
         0    0  1  2  3  4  5  6
        */
        
        Set<PairInt> medialAxis = new HashSet<PairInt>();
        for (int i = 5; i >= 2; --i) {
            medialAxis.add(new PairInt(2, i));
        }
        for (int i = 3; i <= 5; ++i) {
            medialAxis.add(new PairInt(i, 2));
        }
        
        Set<PairInt> boundary = new HashSet<PairInt>();
        for (int i = 1; i <= 3; ++i) {
            for (int j = 1; j <= 6; ++j) {
                boundary.add(new PairInt(i, j));
            }
        }
        for (int i = 4; i <= 6; ++i) {
            for (int j = 1; j <= 3; ++j) {
                boundary.add(new PairInt(i, j));
            }
        }
        boundary.removeAll(medialAxis);
        
        PairIntArray expected = new PairIntArray(boundary.size());
        for (int i = 1; i <= 6; ++i) {
            expected.add(1, i);
        }
        expected.add(2, 6);
        for (int i = 6; i >= 3; --i) {
            expected.add(3, i);
        }
        for (int i = 4; i <= 6; ++i) {
            expected.add(i, 3);
        }
        expected.add(6, 2);
        for (int i = 6; i >= 2; --i) {
            expected.add(i, 1);
        }
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        PairIntArray results = finder.extractOrderedBorder(
            boundary, medialAxis, shapePoints);
        
        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
    }

    private TIntObjectMap<Set<PairInt>> createCostAdjacencyMap(
        List<PairInt> points) {
        
        TIntObjectMap<Set<PairInt>> map = 
            new TIntObjectHashMap<Set<PairInt>>();
        
        TObjectIntMap<PairInt> pointIndexMap
             = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < points.size(); ++i) {
            PairInt p = points.get(i);
            pointIndexMap.put(p, i);
        }
                
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            int idx1 = pointIndexMap.get(p);
            PairInt p1 = new PairInt(idx1, 1);
            
            Set<PairInt> set1 = map.get(idx1);
            if (set1 == null) {
                set1 = new HashSet<PairInt>();
                map.put(idx1, set1);
            }
            
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (!pointIndexMap.containsKey(p2)) {
                    continue;
                }
                int idx2 = pointIndexMap.get(p2);
                
                Set<PairInt> set2 = map.get(idx2);
                if (set2 == null) {
                    set2 = new HashSet<PairInt>();
                    map.put(idx2, set2);
                }
                
                PairInt p3 = new PairInt(idx2, 1);
                
                set1.add(p3);
                set2.add(p1);
            }
        }
        
        assertEquals(points.size(), map.size());
        
        return map;
    }
    
    public void testMergeAdjacentOrderedBorders() {
        
        /*
        6  *  *  *
        5  *  .  *
        4  *  .  *
        3  *  .  * 
        2  *  .  * 
        1  *  *  *
        0  1  2  3  4  5  6
        */
        PairIntArray a1 = new PairIntArray();
        for (int i = 1; i <= 6; ++i) {
            a1.add(1, i);
        }
        a1.add(2, 6);
        for (int i = 6; i >= 1; --i) {
            a1.add(3, i);
        }
        a1.add(2, 1);
        
        /*
        
        3           *  *  *
        2           *  .  *
        1           *  *  *
        0  1  2  3  4  5  6
        */
        PairIntArray a2 = new PairIntArray();
        for (int i = 1; i <= 3; ++i) {
            a2.add(4, i);
        }
        a2.add(5, 3);
        for (int i = 3; i >= 1; --i) {
            a2.add(6, i);
        }
        a2.add(5, 1);
       
        /*
        6  *  *  *
        5  *  .  *
        4  *  .  *9
        3  *  .  *  *2 *  *
        2  *  .  .  .  .  *
        1 0*  *12*  *8 *  *
        0  1  2  3  4  5  6
        */
        PairIntArray expected = new PairIntArray();
        for (int i = 1; i <= 6; ++i) {
            expected.add(1, i);
        }
        expected.add(2, 6);
        for (int i = 6; i >= 4; --i) {
            expected.add(3, i);
        }
        for (int i = 4; i <= 6; ++i) {
            expected.add(i, 3);
        }
        expected.add(6, 2);
        for (int i = 6; i >= 2; --i) {
            expected.add(i, 1);
        }
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray r = finder.mergeAdjacentOrderedBorders(a1, a2);
        
        assertEquals(expected.getN(), r.getN());
        
        for (int i = 0; i < r.getN(); ++i) {
            assertEquals(expected.getX(i), r.getX(i));
            assertEquals(expected.getY(i), r.getY(i));
        }
    }
   
}
