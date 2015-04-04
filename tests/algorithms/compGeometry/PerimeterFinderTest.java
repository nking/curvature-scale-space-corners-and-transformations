package algorithms.compGeometry;

import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;
import org.junit.After;
import org.junit.Before;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PerimeterFinderTest extends TestCase {
    
    public PerimeterFinderTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    public void testFind() throws Exception {
        
        PerimeterFinder finder = new PerimeterFinder();
        
        /*
                   0  1 2 3 4 5 6
         row 0:    @  . . @ . . @ . .
         row 1:    .  . . . . . . . .
         row 2:    .  . . . . . . . .
         row 3:    .  . . @ . . @ . .                    
        */
        
        int[] rowMinMax = new int[2];
        
        Set<PairInt> points = new HashSet<PairInt>();
        points.add(new PairInt(0, 0));
        points.add(new PairInt(3, 0));
        points.add(new PairInt(6, 0));
        points.add(new PairInt(3, 3));
        points.add(new PairInt(6, 3));
        
        rowMinMax = new int[2];
        Map<Integer, PairInt> rowColBounds = finder.find(points, rowMinMax);
        
        assertTrue(rowColBounds.size() == 4);
        
        assertTrue(rowMinMax[0] == 0);
        assertTrue(rowMinMax[1] == 3);
        
        PairInt cb = rowColBounds.get(Integer.valueOf(0));
        assertNotNull(cb);
        assertTrue(cb.getX() == 0);
        assertTrue(cb.getY() == 6);
        
        cb = rowColBounds.get(Integer.valueOf(3));
        assertNotNull(cb);
        assertTrue(cb.getX() == 3);
        assertTrue(cb.getY() == 6);
        
        cb = rowColBounds.get(Integer.valueOf(1));
        assertNotNull(cb);
        assertTrue(cb.getX() == 1);
        assertTrue(cb.getY() == 6);
        
        cb = rowColBounds.get(Integer.valueOf(2));
        assertNotNull(cb);
        assertTrue(cb.getX() == 2);
        assertTrue(cb.getY() == 6);
    }
    
    public void testFind2() throws Exception {
        
        PerimeterFinder finder = new PerimeterFinder();
        
        /*
                   0  1 2 3 4 5 6 7 8 9
         row 0:    .  . . @ . . @ @ @ @
         row 1:    .  . . . . . . . . .
         row 2:    .  @ . . @ . . . . .
         row 3:    @  . . @ . . @ . . .                  
        */
       
        int[] rowMinMax = new int[2];
                
        Set<PairInt> points2 = new HashSet<PairInt>();
        points2.add(new PairInt(3, 0));
        points2.add(new PairInt(6, 0));
        points2.add(new PairInt(7, 0));
        points2.add(new PairInt(8, 0));
        points2.add(new PairInt(9, 0));
        
        points2.add(new PairInt(1, 2));
        points2.add(new PairInt(4, 2));
        
        points2.add(new PairInt(0, 3));
        points2.add(new PairInt(3, 3));
        points2.add(new PairInt(6, 3));
        
        rowMinMax = new int[2];
        Map<Integer, PairInt> rowColBounds = finder.find(points2, rowMinMax);
        
        assertTrue(rowColBounds.size() == 4);
        
        assertTrue(rowMinMax[0] == 0);
        assertTrue(rowMinMax[1] == 3);
        
        PairInt cb = rowColBounds.get(Integer.valueOf(0));
        assertNotNull(cb);
        assertTrue(cb.getX() == 3);
        assertTrue(cb.getY() == 9);
        
        cb = rowColBounds.get(Integer.valueOf(3));
        assertNotNull(cb);
        assertTrue(cb.getX() == 0);
        assertTrue(cb.getY() == 6);
        
        cb = rowColBounds.get(Integer.valueOf(1));
        assertNotNull(cb);
        assertTrue(cb.getX() == 2);
        assertTrue(cb.getY() == 7);
        
        cb = rowColBounds.get(Integer.valueOf(2));
        assertNotNull(cb);
        assertTrue(cb.getX() == 1);
        assertTrue(cb.getY() == 4);
        
    }
    
    public void testFind3() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources("test_mask_0.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);

        finder.findGroups(0);
                
        int n = finder.getNumberOfGroups();
        
        assertTrue(n == 1);
        List<Set<PairInt> > list = new ArrayList<Set<PairInt> >();        
        for (int gIdx = 0; gIdx < n; gIdx++) {
            Set<PairInt> set = new HashSet<PairInt>();            
            finder.getXY(gIdx, set);
            list.add(set);
        }        
        Collections.sort(list, new ListSizeComparator());
        
        Set<PairInt> skyPoints = list.get(list.size() - 1);
        
        int h = img.getHeight();
        int w = img.getWidth();
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        int[] outputRowMinMax = new int[2];
        Map<Integer, PairInt> rowColRanges = perimeterFinder.find(
            skyPoints, outputRowMinMax);
            
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 240);
        
        PairInt cRange = rowColRanges.get(Integer.valueOf(0));
        assertTrue(cRange.getX() == 0);
        assertTrue(cRange.getY() == (w - 1));
        
        cRange = rowColRanges.get(Integer.valueOf(5));
        assertTrue(cRange.getX() == 3);
        assertTrue(cRange.getY() == (w - 1));
        
        cRange = rowColRanges.get(Integer.valueOf(9));
        assertTrue(cRange.getX() == 7);
        assertTrue(cRange.getY() == (w - 1));
        
        cRange = rowColRanges.get(Integer.valueOf(141));
        assertTrue(cRange.getX() == 0);
        assertTrue(cRange.getY() == (w - 1));
        
        cRange = rowColRanges.get(Integer.valueOf(178));
        assertTrue(cRange.getX() == 55);
        assertTrue(cRange.getY() == 268);
        
        cRange = rowColRanges.get(Integer.valueOf(240));
        assertTrue(cRange.getX() == 115);
        assertTrue(cRange.getY() == 128);
    }
    
    public void testFindRowColRanges0() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet0();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 3, 0, 1);
        assertTrue(rowColRanges.size() == 2);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
    }
    
    public void testFindRowColRanges1() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet1();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 7, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 3);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 0);
        
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 2);
        assertTrue(colRange.getY() == 3);
        
        colRange = colRanges.get(2);
        assertTrue(colRange.getX() == 6);
        assertTrue(colRange.getY() == 7);
    }
    
    public void testFindRowColRanges2() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet2();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 6, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 3);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 0);
        
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 2);
        assertTrue(colRange.getY() == 3);
        
        colRange = colRanges.get(2);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 5);
        
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
    }
    
    public void testFindRowColRanges3() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet3();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 6, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 0);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
    }
    
    public void testFindRowColRanges4() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet4();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 6, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 2);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
    }
    
    public void testFindRowColRanges5() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet5();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 6, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 2);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 0);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
    }
    
    public void testBoundedByPointsInLowerRows0() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet0();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 3, 0, 1);
        assertTrue(rowColRanges.size() == 2);
        
        int row = 0;
        boolean bounded = perimeterFinder.boundedByPointsInHigherRows(row, 
            1, 2, 1, rowColRanges);
        assertTrue(bounded);
        
        row = 1;
        bounded = perimeterFinder.boundedByPointsInHigherRows(row, 
            1, 2, 1, rowColRanges);
        assertFalse(bounded);
        
        
        row = 1;
        int minRow = 0;
        bounded = perimeterFinder.boundedByPointsInLowerRows(row, 
            1, 2, minRow, rowColRanges);
        assertTrue(bounded);
        
        row = 0;
        bounded = perimeterFinder.boundedByPointsInLowerRows(row, 
            1, 2, minRow, rowColRanges);
        assertFalse(bounded);
        
    }
    
    public void testBoundedByPointsInLowerRows1() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet1();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 7, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        int row = 1;
        int maxRow = 2;
        boolean bounded = perimeterFinder.boundedByPointsInHigherRows(row, 
            1, 1, maxRow, rowColRanges);
        assertTrue(bounded);
        
        row = 1;
        int minRow = 0;
        bounded = perimeterFinder.boundedByPointsInLowerRows(row, 
            1, 1, minRow, rowColRanges);
        assertTrue(bounded);
        
        
        row = 1;
        maxRow = 2;
        bounded = perimeterFinder.boundedByPointsInHigherRows(row, 
            4, 5, maxRow, rowColRanges);
        assertTrue(bounded);
        
        row = 1;
        minRow = 0;
        bounded = perimeterFinder.boundedByPointsInLowerRows(row, 
            4, 5, minRow, rowColRanges);
        assertTrue(bounded);
    }
    
    public void testBoundedByPointsInLowerRows4() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet4();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 7, 0, 2);
        assertTrue(rowColRanges.size() == 3);
        
        int row = 1;
        int maxRow = 2;
        boolean bounded = perimeterFinder.boundedByPointsInHigherRows(row, 
            4, 4, maxRow, rowColRanges);
        assertFalse(bounded);
        
        row = 1;
        int minRow = 0;
        bounded = perimeterFinder.boundedByPointsInLowerRows(row, 
            4, 4, minRow, rowColRanges);
        assertFalse(bounded);
       
    }
    
    public void testFind2_0() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet0();

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find2(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 2);

        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 1);

        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);

        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
    }
    
    public void testFind2_1() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet1();

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find2(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 2);
        
        /*
        1: group of points with points all along the perimeter
        0 1 2 3 4 5 6 7
        @ @ @ @ @ @ @ @
        @   @ @     @ @
        @ @ @ @ @ @ @ @
        */
        assertTrue(outputEmbeddedGapPoints.size() == 3);
        Set<PairInt> expectedEmbedded = new HashSet<PairInt>();
        expectedEmbedded.add(new PairInt(1, 1));
        expectedEmbedded.add(new PairInt(4, 1));
        expectedEmbedded.add(new PairInt(5, 1));
        for (PairInt embedded : outputEmbeddedGapPoints) {
            assertTrue(expectedEmbedded.remove(embedded));
        }
        assertTrue(expectedEmbedded.isEmpty());
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
    }
    
    public void testFind2_2() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet2();

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find2(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 2);
        
        
        assertTrue(outputEmbeddedGapPoints.size() == 2);
        Set<PairInt> expectedEmbedded = new HashSet<PairInt>();
        expectedEmbedded.add(new PairInt(1, 1));
        expectedEmbedded.add(new PairInt(4, 1));
        for (PairInt embedded : outputEmbeddedGapPoints) {
            assertTrue(expectedEmbedded.remove(embedded));
        }
        assertTrue(expectedEmbedded.isEmpty());
        
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 5);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
    }
    
    public void testFind2_3() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet3();

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find2(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 2);
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 0);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
    }
    
    public void testFind2_4() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet4();

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find2(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 2);
        
        assertTrue(outputEmbeddedGapPoints.isEmpty());
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 2);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(1);
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(2);
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
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
    public void testFind2_6() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet6();

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find2(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 10);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 9);
        
        assertTrue(outputEmbeddedGapPoints.isEmpty());
        
        List<PairInt> colRanges = rowColRanges.get(0);
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 9);
        
        colRanges = rowColRanges.get(9);
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 2);
        assertTrue(colRange.getY() == 6);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 8);
        assertTrue(colRange.getY() == 9);
    }
    
    /*
    0: dense group of points w/o gaps
        @ @ @ @
        @ @ @ @
    1: group of points with points all along the perimeter
        @ @ @ @ @ @ @ @
        @   @ @     @ @
        @ @ @ @ @ @ @ @
    2: group of points with points not all along the perimeter
        @ @ @ @ @ @ @
        @   @ @   @
        @ @ @ @ @ @ @
    3: group of points missing rows in some places
        @ @ @ @ @ @ @

        @ @ @ @ @ @ @
    4: group of points missing columns in some places
        @ @ @ @   @ @
        @ @ @ @   @ @
        @ @ @ @   @ @
    5: group of points missing intersecting column and row in some places
        @ @ @ @   @ @

        @ @ @ @   @ @
    */
    private Set<PairInt> getSet0() {
        return getSet(4, 2);        
    }
    private Set<PairInt> getSet1() {
        Set<PairInt> points = getSet(8, 3);
        points.remove(new PairInt(1, 1));
        points.remove(new PairInt(4, 1));
        points.remove(new PairInt(5, 1));
        return points;
    }
    private Set<PairInt> getSet2() {
        Set<PairInt> points = getSet(7, 3);
        points.remove(new PairInt(1, 1));
        points.remove(new PairInt(4, 1));
        points.remove(new PairInt(6, 1));
        return points;
    }
    private Set<PairInt> getSet3() {
        Set<PairInt> points = getSet(7, 3);
        for (int col = 0; col < 7; col++) {
            points.remove(new PairInt(col, 1));
        }
        return points;
    }
    private Set<PairInt> getSet4() {
        Set<PairInt> points = getSet(7, 3);
        for (int row = 0; row < 3; row++) {
            points.remove(new PairInt(4, row));
        }
        return points;
    }
    private Set<PairInt> getSet5() {
        Set<PairInt> points = getSet4();
        for (int col = 0; col < 7; col++) {
            points.remove(new PairInt(col, 1));
        }
        return points;
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
    private Set<PairInt> getSet6() {
        Set<PairInt> points = new HashSet<PairInt>();
        for (int row = 0; row < 2; row++) {
            for (int col = 0; col < 10; col++) {
                points.add(new PairInt(col, row));
            }
        }
        points.add(new PairInt(0, 2));
        points.add(new PairInt(1, 2));
        points.add(new PairInt(2, 2));
        points.add(new PairInt(3, 2));
        points.add(new PairInt(7, 2));
        points.add(new PairInt(8, 2));
        points.add(new PairInt(9, 2));
        
        points.add(new PairInt(0, 3));
        points.add(new PairInt(1, 3));
        points.add(new PairInt(2, 3));
        points.add(new PairInt(5, 3));
        points.add(new PairInt(7, 3));
        points.add(new PairInt(8, 3));
        points.add(new PairInt(9, 3));
        
        points.add(new PairInt(0, 4));
        points.add(new PairInt(1, 4));
        points.add(new PairInt(4, 4));
        points.add(new PairInt(5, 4));
        points.add(new PairInt(8, 4));
        points.add(new PairInt(9, 4));
        
        points.add(new PairInt(0, 5));
        points.add(new PairInt(3, 5));
        points.add(new PairInt(4, 5));
        points.add(new PairInt(5, 5));
        points.add(new PairInt(8, 5));
        points.add(new PairInt(9, 5));
        
        points.add(new PairInt(0, 6));
        points.add(new PairInt(2, 6));
        points.add(new PairInt(3, 6));
        points.add(new PairInt(4, 6));
        points.add(new PairInt(5, 6));
        points.add(new PairInt(6, 6));
        points.add(new PairInt(8, 6));
        points.add(new PairInt(9, 6));
        
        for (int row = 7; row < 10; row++) {
            points.add(new PairInt(2, row));
            points.add(new PairInt(3, row));
            points.add(new PairInt(4, row));
            points.add(new PairInt(5, row));
            points.add(new PairInt(6, row));
            points.add(new PairInt(8, row));
            points.add(new PairInt(9, row));
        }
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
    
    private class ListSizeComparator implements Comparator< Set<PairInt> > {

        @Override
        public int compare(Set<PairInt> o1, Set<PairInt> o2) {
            
            if (o1 == null && o2 != null) {
                return 1;
            } 
            if (o1 != null && o2 == null) {
                return -1;
            }
            
            return Integer.compare(o1.size(), o2.size());
        }
        
    }
}
