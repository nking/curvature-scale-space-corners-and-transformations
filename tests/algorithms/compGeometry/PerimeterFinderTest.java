package algorithms.compGeometry;

import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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
