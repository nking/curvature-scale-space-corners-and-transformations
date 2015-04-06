package algorithms.compGeometry;

import algorithms.compGeometry.PerimeterFinder.Gap;
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
import java.util.Stack;
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
    
    public void testFindRowColRanges0() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> points = getSet0();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.findRowColRanges(
            points, 0, 3, 0, 1);
        assertTrue(rowColRanges.size() == 2);
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
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
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
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
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
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
        
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.isEmpty());
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 2);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 2);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 0);
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 2);

        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 1);

        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);

        colRanges = rowColRanges.get(Integer.valueOf(1));
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
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
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
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 7);
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
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
        
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 5);
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 2);
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 0);
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 2);
        
        assertTrue(outputEmbeddedGapPoints.isEmpty());
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 2);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(2));
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
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 10);
        
        assertTrue(outputRowMinMax[0] == 0);
        assertTrue(outputRowMinMax[1] == 9);
        
        assertTrue(outputEmbeddedGapPoints.isEmpty());
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(0));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 0);
        assertTrue(colRange.getY() == 9);
        
        colRanges = rowColRanges.get(Integer.valueOf(9));
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 2);
        assertTrue(colRange.getY() == 6);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 8);
        assertTrue(colRange.getY() == 9);
        
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
     6  @ # @ @ @ @ @   @ @
     7      @ @ @ @ @   @ @
     8      @ @ @ @ @   @ @
     9      @ @ @ @ @ # @ @
    */
    public void testFindContiguousGaps_6_0() throws Exception {

        Set<PairInt> points = getSet6();
        points.add(new PairInt(1, 6));
        points.add(new PairInt(7, 9));

        int minX = 0;
        int maxX = 9; 
        int minY = 0; 
        int maxY = 9;
    
        List<List<Gap>> expectedGapGroups = new ArrayList<List<Gap>>();
        List<Gap> contigGaps = new ArrayList<Gap>();
        contigGaps.add(new PerimeterFinder.Gap(2, 4, 6));
        contigGaps.add(new PerimeterFinder.Gap(3, 3, 4));
        contigGaps.add(new PerimeterFinder.Gap(3, 6, 6));
        contigGaps.add(new PerimeterFinder.Gap(4, 2, 3));
        contigGaps.add(new PerimeterFinder.Gap(4, 6, 7));
        contigGaps.add(new PerimeterFinder.Gap(5, 1, 2));
        contigGaps.add(new PerimeterFinder.Gap(5, 6, 7));
        contigGaps.add(new PerimeterFinder.Gap(6, 7, 7));
        contigGaps.add(new PerimeterFinder.Gap(7, 7, 7));
        contigGaps.add(new PerimeterFinder.Gap(8, 7, 7));
        expectedGapGroups.add(contigGaps);
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Map<Integer, List<PairInt>> rowColRanges = 
            perimeterFinder.findRowColRanges(points, minX, maxX, minY, maxY);
        assertTrue(rowColRanges.size() == 10);
        
        Stack<Gap> gaps = perimeterFinder.findGaps(rowColRanges, 
            minX, maxX, minY, maxY);
        
        Set<Gap> gaps2 = new HashSet<Gap>(gaps);
        
        assertTrue(gaps.size() == 10);
        
        //--- assert order of items in stack, and assert that all expected 
        //    are present
        
        Set<Gap> expected = new HashSet<Gap>(expectedGapGroups.get(0));
        int lastRow = Integer.MIN_VALUE;
        int lastCol = Integer.MIN_VALUE;
        while (!gaps.isEmpty()) {
            Gap gap = gaps.pop();
            
            int row = gap.getRow();
            int col = gap.getStart();
            
            if (row == lastRow) {
                assertTrue(col > lastCol);
            } else {
                assertTrue(row > lastRow);
            }
            
            assertTrue(expected.remove(gap));
            
            lastCol = gap.getStopInclusive();
            lastRow = row;
        }
        assertTrue(expected.isEmpty());
        
        
        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        Set<PairInt> expectedPoints = new HashSet<PairInt>();
        for (Gap gap : gaps2) {
            for (int i = gap.getStart(); i <= gap.getStopInclusive(); i++) {
                PairInt p = new PairInt(i, gap.getRow());
                expectedPoints.add(p);
            }
        }
        for (PairInt p : outputEmbeddedGapPoints) {
            assertTrue(expectedPoints.remove(p));
        }
        assertTrue(expectedPoints.isEmpty());
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
     6  @ # @ @ @ @ @   @ @
     7      @ @ @ @ @   @ @
     8      @ @ @ @ @   @ @
     9      @ @ @ @ @ # @ @
    */
    public void testFindContiguousGaps_6() throws Exception {

        Set<PairInt> points = getSet6();
        points.add(new PairInt(1, 6));
        points.add(new PairInt(7, 9));

        int minX = 0;
        int maxX = 9; 
        int minY = 0; 
        int maxY = 9;
    
        List<List<Gap>> expectedGapGroups = new ArrayList<List<Gap>>();
        List<Gap> contigGaps = new ArrayList<Gap>();
        contigGaps.add(new PerimeterFinder.Gap(2, 4, 6));
        contigGaps.add(new PerimeterFinder.Gap(3, 3, 4));
        contigGaps.add(new PerimeterFinder.Gap(3, 6, 6));
        contigGaps.add(new PerimeterFinder.Gap(4, 2, 3));
        contigGaps.add(new PerimeterFinder.Gap(4, 6, 7));
        contigGaps.add(new PerimeterFinder.Gap(5, 1, 2));
        contigGaps.add(new PerimeterFinder.Gap(5, 6, 7));
        contigGaps.add(new PerimeterFinder.Gap(6, 7, 7));
        contigGaps.add(new PerimeterFinder.Gap(7, 7, 7));
        contigGaps.add(new PerimeterFinder.Gap(8, 7, 7));
        expectedGapGroups.add(contigGaps);
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Map<Integer, List<PairInt>> rowColRanges = 
            perimeterFinder.findRowColRanges(points, minX, maxX, minY, maxY);
        assertTrue(rowColRanges.size() == 10);
        
        List<List<Gap>> contiguousGapGroups = perimeterFinder.findContiguousGaps(rowColRanges, 
            minX, maxX, minY, maxY);
        
        assertTrue(contiguousGapGroups.size() == 1);
        
        Set<Gap> expected = new HashSet<Gap>(expectedGapGroups.get(0));        
        for (Gap gap : contiguousGapGroups.get(0)) {
            assertTrue(expected.remove(gap));
        }
        assertTrue(expected.isEmpty());
        
        
        //-----------
        List<List<Gap>> boundedGaps = perimeterFinder.findBoundedGaps(
            contiguousGapGroups, minY, maxY, rowColRanges);
        
        expected = new HashSet<Gap>(expectedGapGroups.get(0)); 
        for (List<Gap> gapGroups : boundedGaps) {
            for (Gap gap : gapGroups) {
                assertTrue(expected.remove(gap));
            }
        }
        assertTrue(expected.isEmpty());
    }
    
    /*2: group of points with points not all along the perimeter
        0 1 2 3 4 5 6
        @ @ @ @ @ @ @
        @   @ @   @
        @ @ @ @ @ @ @
    */
    public void testFindContiguousGaps_2() throws Exception {

        Set<PairInt> points = getSet2();

        int minX = 0;
        int maxX = 6; 
        int minY = 0; 
        int maxY = 2;
    
        List<List<Gap>> expectedGapGroups = new ArrayList<List<Gap>>();
        List<Gap> contigGaps = new ArrayList<Gap>();
        contigGaps.add(new PerimeterFinder.Gap(1, 1, 1));
        expectedGapGroups.add(contigGaps);
        contigGaps = new ArrayList<Gap>();
        contigGaps.add(new PerimeterFinder.Gap(1, 4, 4));
        expectedGapGroups.add(contigGaps);
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Map<Integer, List<PairInt>> rowColRanges = 
            perimeterFinder.findRowColRanges(points, minX, maxX, minY, maxY);
        assertTrue(rowColRanges.size() == 3);
        
        List<List<Gap>> contiguousGapGroups = perimeterFinder.findContiguousGaps(rowColRanges, 
            minX, maxX, minY, maxY);
        
        assertTrue(contiguousGapGroups.size() == 2);
        
        Set<Gap> expected = new HashSet<Gap>(expectedGapGroups.get(0));
        expected.addAll(expectedGapGroups.get(1));
        for (List<Gap> contiguousGapGroup : contiguousGapGroups) {
            for (Gap gap : contiguousGapGroup) {
                assertTrue(expected.remove(gap));
            }
        }
        assertTrue(expected.isEmpty());
        
        
        //-----------
        List<List<Gap>> boundedGaps = perimeterFinder.findBoundedGaps(
            contiguousGapGroups, minY, maxY, rowColRanges);
        
        expected = new HashSet<Gap>(expectedGapGroups.get(0));
        expected.addAll(expectedGapGroups.get(1));
        for (List<Gap> contiguousGapGroup : contiguousGapGroups) {
            for (Gap gap : contiguousGapGroup) {
                assertTrue(expected.remove(gap));
            }
        }
        assertTrue(expected.isEmpty());
    }
    
    public void testUpdateRowColRangesForAddedPoints_0() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        Set<PairInt> points = getSet(4, 2, 5, 5);

        int[] outputRowMinMax = new int[2];

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 2);

        /*
         4 5 6 7 8 9 
           @ @ @ @ %  5
         % @ @ @ @    6
           %          7
        */
        Set<PairInt> add = new HashSet<PairInt>();
        add.add(new PairInt(4, 6));
        add.add(new PairInt(5, 7));
        add.add(new PairInt(9, 5));
        
        perimeterFinder.updateRowColRangesForAddedPoints(rowColRanges, 
            outputRowMinMax, add);
    
        assertTrue(outputRowMinMax[0] == 5);
        assertTrue(outputRowMinMax[1] == 7);
        
        assertTrue(rowColRanges.size() == 3);
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(5));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 9);

        colRanges = rowColRanges.get(Integer.valueOf(6));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 4);
        assertTrue(colRange.getY() == 8);
        
        colRanges = rowColRanges.get(Integer.valueOf(7));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 5);
        assertTrue(colRange.getY() == 5);
        
        /*
         4 5 6 7 8 9 
           @ @ @ @ %  5
         % @ @ @ @    6
           %          7
        */
        
        int imageMaxColumn = 9;
        int imageMaxRow = 10;
        
        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            rowColRanges, outputRowMinMax, imageMaxColumn, imageMaxRow);
        
        assertTrue(borderPixels.size() == 10);
        
        Set<PairInt> expected = new HashSet<PairInt>(points);
        expected.addAll(add);
        expected.remove(new PairInt(5, 6));
        
        for (PairInt p : borderPixels) {
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
    }
    
    public void testUpdateRowColRangesForAddedPoints_5() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        /*
            3 4 5 6 7 8 9
            @ @ @ @   @ @  2
                           3
            @ @ @ @   @ @  4
        */
        Set<PairInt> points = getSet(7, 3, 3, 2);
        for (int col = 3; col < 10; col++) {
            points.remove(new PairInt(col, 3));
        }
        points.remove(new PairInt(7, 2));
        points.remove(new PairInt(7, 4));

        int[] outputRowMinMax = new int[2];
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        /*
            3 4 5 6 7 8 9
            %              1
            @ @ @ @ % @ @  2
            %     %        3
            @ @ @ @ % @ @  4
        */
        Set<PairInt> add = new HashSet<PairInt>();
        add.add(new PairInt(3, 1));
        add.add(new PairInt(3, 3));
        add.add(new PairInt(6, 3));
        add.add(new PairInt(7, 2));
        add.add(new PairInt(7, 4));

        perimeterFinder.updateRowColRangesForAddedPoints(rowColRanges, 
            outputRowMinMax, add);
        
        assertTrue(outputRowMinMax[0] == 1);
        assertTrue(outputRowMinMax[1] == 4);
        
        assertTrue(rowColRanges.size() == 4);
        
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 3);

        colRanges = rowColRanges.get(Integer.valueOf(2));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 9);
        
        colRanges = rowColRanges.get(Integer.valueOf(3));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(4));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 9);
    }
    
    public void testUpdateForAddedPoints_5() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        /*
            3 4 5 6 7 8 9
            @ @ @ @   @ @  2
                           3
            @ @ @ @   @ @  4
        */
        Set<PairInt> points = getSet(7, 3, 3, 2);
        for (int col = 3; col < 10; col++) {
            points.remove(new PairInt(col, 3));
        }
        points.remove(new PairInt(7, 2));
        points.remove(new PairInt(7, 4));

        int[] outputRowMinMax = new int[2];
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        /*
            3 4 5 6 7 8 9
            %              1
            @ @ @ @ % @ @  2
            %     %        3
            @ @ @ @ % @ @  4
        */
        Set<PairInt> add = new HashSet<PairInt>();
        add.add(new PairInt(3, 1));
        add.add(new PairInt(3, 3));
        add.add(new PairInt(6, 3));
        add.add(new PairInt(7, 2));
        add.add(new PairInt(7, 4));
        
        // ------ this is the intermediate stage, not yet condensed ranges -----
        boolean allAreWithinExistingRange = 
            perimeterFinder.updateForAddedPoints(rowColRanges, 
            outputRowMinMax, add);
        
        assertFalse(allAreWithinExistingRange);
        
        assertTrue(outputRowMinMax[0] == 1);
        assertTrue(outputRowMinMax[1] == 4);
        
        assertTrue(rowColRanges.size() == 4);
        
        /*
            3 4 5 6 7 8 9
            %              1
            @ @ @ @ % @ @  2
            %     %        3
            @ @ @ @ % @ @  4
        */
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(1));
        assertTrue(colRanges.size() == 1);
        PairInt colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 3);

        colRanges = rowColRanges.get(Integer.valueOf(2));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 9);
        
        colRanges = rowColRanges.get(Integer.valueOf(3));
        assertTrue(colRanges.size() == 2);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 3);
        colRange = colRanges.get(1);
        assertTrue(colRange.getX() == 6);
        assertTrue(colRange.getY() == 6);
        
        colRanges = rowColRanges.get(Integer.valueOf(4));
        assertTrue(colRanges.size() == 1);
        colRange = colRanges.get(0);
        assertTrue(colRange.getX() == 3);
        assertTrue(colRange.getY() == 9);
    }
    
    public void testGetBorderPixels_5() throws Exception {

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        /*
        5: group of points missing intersecting column and row in some places
        0 1 2 3 4 5 6  
        @ @ @ @   @ @  0
          % %     % %  1
        @ @ @ @   @ @  2
        
        TODO:  make a test out of this and make sure that pixels in top 
        row in cols 0-2 are in borderPixels
        0 1 2 3 4 5 6  
        @ @ @ @ @ @ @   0
              % @   %   1
        @ @ @ @ @ @ @   2
        */        
        Set<PairInt> points = getSet5();
        points.add(new PairInt(1, 1));
        points.add(new PairInt(2, 1));
        points.add(new PairInt(5, 1));
        points.add(new PairInt(6, 1));

        int[] outputRowMinMax = new int[2];
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> rowColRanges = perimeterFinder.find(
            points, outputRowMinMax, outputEmbeddedGapPoints);
        assertTrue(rowColRanges.size() == 3);
        
        int imageMaxColumn = 10;
        int imageMaxRow = 10;
        
        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            rowColRanges, outputRowMinMax, imageMaxColumn, imageMaxRow);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(1, 1));
        expected.add(new PairInt(3, 0));
        expected.add(new PairInt(2, 1));
        expected.add(new PairInt(3, 2));
        expected.add(new PairInt(2, 2));
        expected.add(new PairInt(1, 2));
        expected.add(new PairInt(0, 2));
        expected.add(new PairInt(5, 0));
        expected.add(new PairInt(5, 1));
        expected.add(new PairInt(5, 2));
        expected.add(new PairInt(6, 0));
        expected.add(new PairInt(6, 1));
        expected.add(new PairInt(6, 2));
        
        for (PairInt p : borderPixels) {
            assertTrue(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
       
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
    
    private Set<PairInt> getSet(int nCols, int nRows, int xOffset, int yOffset) {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int row = 0; row < nRows; row++) {
            for (int col = 0; col < nCols; col++) {
                points.add(new PairInt(col + xOffset, row + yOffset));
            }
        }
        
        return points;
    }
    
}
