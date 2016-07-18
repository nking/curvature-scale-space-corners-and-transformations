package algorithms.compGeometry;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
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
        
        Set<PairInt> data = getSet1();
        
        Set<PairInt> expected = new HashSet<PairInt>();
        List<List<PairInt>> orderedBounds = getSet1Boundaries();
        for (List<PairInt> list : orderedBounds) {
            expected.addAll(list);
        }
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        Set<PairInt> boundary = finder.extractBorder(data);
        
        for (PairInt p : boundary) {
            //System.out.println("boundary p=" + p);
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
    }
    
    /**
     * returns lists of sequential CCW ordered contiguous
     * boundary points.
     * 
     * @return 
     */
    private List<List<PairInt>> getSet1Boundaries() {
        /*
        points with gaps, similar to st. louis arch image gaps:
            0 1 2 3 4 5 6 7 8 9
         0  @ @ @ @ @ @ @ @ @ @
         1  @ @ @ @ @ @ @ @ @ @
         2  @ @ @ @       @ @ @
         3  @ @ @     @   @ @ @
         4  @ @     @ @     @ @
         5  @     @ @ @     @ @
         6  @   @ @ @ @ @   @ @
         7      @ @ @ @ @   @ @
         8      @ @ @ @ @   @ @
         9      @ @ @ @ @   @ @
        
        Can see from point (0,6) and the disconnected contiguous
        points that the boundaries would have to be considered
        sequential contiguous points whose endpoints
        may have a gap due to resolution of single pixel width
        extremity.
        
        */
        List<PairInt> list0 = new ArrayList<PairInt>();
        
        // ---- HERE is complex part of ordering the border ----
        list0.add(new PairInt(0, 6));
    
        for (int i = 5; i >= 1; --i) {
            list0.add(new PairInt(0, i));
        }
        for (int i = 0; i <= 9; ++i) {
            list0.add(new PairInt(i, 0));
        }
        for (int i = 1; i <= 9; ++i) {
            list0.add(new PairInt(9, i));
        }
        for (int i = 9; i >= 3; --i) {
            list0.add(new PairInt(8, i));
        }
        for (int i = 3; i >= 1; --i) {
            list0.add(new PairInt(7, i));
        }
        for (int i = 6; i >= 3; --i) {
            list0.add(new PairInt(i, 1));
        }
        list0.add(new PairInt(3, 2));
        list0.add(new PairInt(2, 2));
        list0.add(new PairInt(2, 3));
        list0.add(new PairInt(1, 3));
        list0.add(new PairInt(1, 4));
        
        // ---- second contiguous region ----
        
        List<PairInt> list1 = new ArrayList<PairInt>();
        for (int i = 6; i >= 2; --i) {
            list1.add(new PairInt(i, 9));
        }
        for (int i = 8; i >= 6; --i) {
            list1.add(new PairInt(2, i));
        }
        int y = 6;
        for (int i = 3; i <= 4; ++i) {
            list1.add(new PairInt(i, y));
            list1.add(new PairInt(i, y - 1));
            y--;
        }
        for (int i = 3; i <= 6; ++i) {
            list1.add(new PairInt(5, i));
        }
        for (int i = 6; i <= 8; ++i) {
            list1.add(new PairInt(6, i));
        }
        
        List<List<PairInt>> output = new ArrayList<List<PairInt>>();
        output.add(list0);
        output.add(list1);
        
        return output;
    }
 
    /*
    points with gaps, similar to st. louis arch image gaps:
        0 1 2 3 4 5 6 7 8 9
     0  @ @ @ @ @ @ @ @ @ @
     1  @ @ @ @ @ @ @ @ @ @
     2  @ @ @ @       @ @ @
     3  @ @ @     @   @ @ @
     4  @ @     @ @     @ @
     5  @     @ @ @     @ @
     6  @   @ @ @ @ @   @ @
     7      @ @ @ @ @   @ @
     8      @ @ @ @ @   @ @
     9      @ @ @ @ @   @ @
    */
    private Set<PairInt> getSet1() {
        Set<PairInt> points = getSet(10, 10);
        points.remove(new PairInt(4, 2));
        points.remove(new PairInt(5, 2));
        points.remove(new PairInt(6, 2));
        points.remove(new PairInt(3, 3));
        points.remove(new PairInt(4, 3));
        points.remove(new PairInt(6, 3));
        points.remove(new PairInt(2, 4));
        points.remove(new PairInt(3, 4));
        points.remove(new PairInt(6, 4));
        points.remove(new PairInt(7, 4));
        points.remove(new PairInt(1, 5));
        points.remove(new PairInt(2, 5));
        points.remove(new PairInt(6, 5));
        points.remove(new PairInt(7, 5));
        points.remove(new PairInt(1, 6));
        points.remove(new PairInt(7, 6));
        points.remove(new PairInt(0, 7));
        points.remove(new PairInt(1, 7));
        points.remove(new PairInt(7, 7));
        points.remove(new PairInt(0, 8));
        points.remove(new PairInt(1, 8));
        points.remove(new PairInt(7, 8));
        points.remove(new PairInt(0, 9));
        points.remove(new PairInt(1, 9));
        points.remove(new PairInt(7, 9));
        return points;
    }
    
    private Set<PairInt> getSet(int nCols, int nRows) {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                points.add(new PairInt(col, row));
            }
        }
        
        return points;
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
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        PairIntArray results = finder.extractOrderedBorder(
            boundary, medialAxis);
        
        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //    expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
    }
    
    public void testOrdered2() {
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
            boundary, medialAxis);
        
        assertEquals(expected.getN(), results.getN());
        
        for (int i = 0; i < expected.getN(); ++i) {
            //System.out.println("i=" + i + " " +
            //expected.getX(i) + " " + expected.getY(i));
            assertEquals(expected.getX(i), results.getX(i));
            assertEquals(expected.getY(i), results.getY(i));
        }
    }
}
