package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.BlobMedialAxes;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
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
            System.out.println("boundary p=" + p);
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
}
