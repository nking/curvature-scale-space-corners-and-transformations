package algorithms.compGeometry;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.HashSet;
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
    
    public void test1() throws Exception {
        
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
        
        //System.out.println("UNORDERED n=" + resultUnordered.size() 
        //    + " => " + Misc.convertWithoutOrder(resultUnordered));
        // ----
        
        contiguousPoints = getSet1();
        finder = new PerimeterFinder2();
        PairIntArray results = finder.extractOrderedBorder(
            contiguousPoints);
        
        /*
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
        System.out.println("results=" + results.toString());
        
        expected = new PairIntArray();
        expected.add(0, 3); expected.add(0, 4); expected.add(0, 5); 
        expected.add(0, 6); expected.add(0, 7); expected.add(0, 8);
        expected.add(0, 9); expected.add(1, 9); expected.add(2, 9);
        expected.add(3, 9); expected.add(4, 9); 
        expected.add(5, 10); expected.add(6, 11); expected.add(7, 12); 
        expected.add(7, 11); expected.add(6, 10);
        expected.add(7, 9); expected.add(8, 9); expected.add(9, 9);
        expected.add(9, 8); expected.add(9, 7); 
        expected.add(9, 6); expected.add(9, 5); 
        expected.add(9, 4); expected.add(9, 3); expected.add(9, 2);
        expected.add(9, 1); expected.add(9, 0); 
        expected.add(8, 2); expected.add(8, 3); expected.add(8, 4);
        expected.add(8, 5);
        expected.add(7, 6); expected.add(7, 7); 
        expected.add(6, 8); expected.add(5, 8); expected.add(4, 8); 
        expected.add(3, 7); expected.add(2, 6); expected.add(1, 5);

        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //    expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
        
        // --- same test with pixel indexes
        int w = 30;
        int h = 30;
        TIntSet cSets = new TIntHashSet();
        Set<PairInt> set = getSet1();
        for (PairInt p : set) {
            int pixIdx = (p.getY() * w) + p.getX();
            cSets.add(pixIdx);
        }
        results = finder.extractOrderedBorder(cSets, w, h);
        
        /*
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
        System.out.println("results=" + results.toString());
        
        expected = new PairIntArray();
        expected.add(0, 3); expected.add(0, 4); expected.add(0, 5); 
        expected.add(0, 6); expected.add(0, 7); expected.add(0, 8);
        expected.add(0, 9); expected.add(1, 9); expected.add(2, 9);
        expected.add(3, 9); expected.add(4, 9); 
        expected.add(5, 10); expected.add(6, 11); expected.add(7, 12); 
        expected.add(7, 11); expected.add(6, 10);
        expected.add(7, 9); expected.add(8, 9); expected.add(9, 9);
        expected.add(9, 8); expected.add(9, 7); 
        expected.add(9, 6); expected.add(9, 5); 
        expected.add(9, 4); expected.add(9, 3); expected.add(9, 2);
        expected.add(9, 1); expected.add(9, 0); 
        expected.add(8, 2); expected.add(8, 3); expected.add(8, 4);
        expected.add(8, 5);
        expected.add(7, 6); expected.add(7, 7); 
        expected.add(6, 8); expected.add(5, 8); expected.add(4, 8); 
        expected.add(3, 7); expected.add(2, 6); expected.add(1, 5);

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
            shapePoints);
        
        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //    expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
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
        
        /*
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray r = finder.mergeAdjacentOrderedBorders(a1, a2);
        
        assertEquals(expected.getN(), r.getN());
        
        for (int i = 0; i < r.getN(); ++i) {
            assertEquals(expected.getX(i), r.getX(i));
            assertEquals(expected.getY(i), r.getY(i));
        }
        */
    }
   
}
