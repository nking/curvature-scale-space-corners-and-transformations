package algorithms.compGeometry;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashSet;
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
        
        PairIntArray points = new PairIntArray();
        /*
                   0  1 2 3 4 5 6
         row 0:    @  . . @ . . @ . .
         row 1:    .  . . . . . . . .
         row 2:    .  . . . . . . . .
         row 3:    .  . . @ . . @ . .                    
        */
        points.add(0, 0);
        points.add(3, 0);
        points.add(6, 0);
        points.add(3, 3);
        points.add(6, 3);
        
        int[] rowMinMax = new int[2];
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
        
        //=============================================================
        Set<PairInt> points2 = new HashSet<PairInt>();
        points2.add(new PairInt(0, 0));
        points2.add(new PairInt(3, 0));
        points2.add(new PairInt(6, 0));
        points2.add(new PairInt(3, 3));
        points2.add(new PairInt(6, 3));
        
        rowMinMax = new int[2];
        rowColBounds = finder.find(points, rowMinMax);
        
        assertTrue(rowColBounds.size() == 4);
        
        assertTrue(rowMinMax[0] == 0);
        assertTrue(rowMinMax[1] == 3);
        
        cb = rowColBounds.get(Integer.valueOf(0));
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
        
        PairIntArray points = new PairIntArray();
        /*
                   0  1 2 3 4 5 6 7 8 9
         row 0:    .  . . @ . . @ @ @ @
         row 1:    .  . . . . . . . . .
         row 2:    .  @ . . @ . . . . .
         row 3:    @  . . @ . . @ . . .                  
        */
        points.add(3, 0);
        points.add(6, 0);
        points.add(7, 0);
        points.add(8, 0);
        points.add(9, 0);
        
        points.add(1, 2);
        points.add(4, 2);
        
        points.add(0, 3);
        points.add(3, 3);
        points.add(6, 3);
        
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColBounds = finder.find(points, rowMinMax);
        
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
        assertTrue(cb.getY() == 8);
        
        cb = rowColBounds.get(Integer.valueOf(2));
        assertNotNull(cb);
        assertTrue(cb.getX() == 1);
        assertTrue(cb.getY() == 7);
        
        //=============================================================
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
        rowColBounds = finder.find(points, rowMinMax);
        
        assertTrue(rowColBounds.size() == 4);
        
        assertTrue(rowMinMax[0] == 0);
        assertTrue(rowMinMax[1] == 3);
        
        cb = rowColBounds.get(Integer.valueOf(0));
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
        assertTrue(cb.getY() == 8);
        
        cb = rowColBounds.get(Integer.valueOf(2));
        assertNotNull(cb);
        assertTrue(cb.getX() == 1);
        assertTrue(cb.getY() == 7);
        
    }
}
