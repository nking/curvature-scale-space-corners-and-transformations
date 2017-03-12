package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HoughTransformTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public HoughTransformTest(String testName) {
        super(testName);
    }

    
    private double distance(CornerRegion cr, PairInt p) {
        int x1 = cr.getX()[cr.getKMaxIdx()];
        int y1 = cr.getY()[cr.getKMaxIdx()];

        int diffX = x1 - p.getX();
        int diffY = y1 - p.getY();
        
        return Math.sqrt(diffX * diffX + diffY * diffY);
    }
    
    private double distance(PairInt p1, PairInt p2) {

        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        
        return Math.sqrt(diffX * diffX + diffY * diffY);
    }

    // assuming cr was tested as further from endPoints than 2 pixels
    private boolean isInBetween(PairInt[] endPoints, double distBetweenEndPoints,
        CornerRegion cr) {
        
        double d0 = distance(cr, endPoints[0]);
        
        double d1 = distance(cr, endPoints[1]);
        
        return (d0 + d1) < distBetweenEndPoints;
    }

    public void testFindContiguousLines()  throws Exception {
        
        /*
        public List<Set<PairInt>> findContiguousLines(Set<PairInt> points, 
        GreyscaleImage theta360) {
        */
        /*
         *           Y
         *          90
         *     135   |    +45
         *           |
         *   180------------ 0   X
         *           |
         *    225    |   315
         *          270
         *
        
        2
        1 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
        0           @ @ @ @ @ @
        9         @             @
        8       @                 @
        7     @                     @     
        6   @                         @
        5 @                             @
        4   @                         @
        3     @                     @
        2       @                 @
        1         @             @
        0           @ @ @ @ @ @
          0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
        */
        /*Set<PairInt> points = new HashSet<PairInt>();
        
        points.add(new PairInt(5, 0)); points.add(new PairInt(6, 0)); points.add(new PairInt(7, 0)); 
        points.add(new PairInt(8, 0)); points.add(new PairInt(9, 0)); points.add(new PairInt(10, 0));
        points.add(new PairInt(4, 1)); points.add(new PairInt(11, 1));
        points.add(new PairInt(3, 2)); points.add(new PairInt(12, 2));
        points.add(new PairInt(2, 3)); points.add(new PairInt(13, 3));
        points.add(new PairInt(1, 4)); points.add(new PairInt(14, 4));
        points.add(new PairInt(0, 5)); points.add(new PairInt(15, 5));
        points.add(new PairInt(1, 6)); points.add(new PairInt(14, 6));
        points.add(new PairInt(2, 7)); points.add(new PairInt(13, 7));
        points.add(new PairInt(3, 8)); points.add(new PairInt(12, 8));
        points.add(new PairInt(4, 9)); points.add(new PairInt(11, 9));
        points.add(new PairInt(5, 10)); points.add(new PairInt(6, 10)); points.add(new PairInt(7, 10));
        points.add(new PairInt(8, 10)); points.add(new PairInt(9, 10)); points.add(new PairInt(10, 10));
        
        int nExpected = points.size();
        
        HoughTransform ht = new HoughTransform();
        Map<Set<PairInt>, PairInt> lines 
            = ht.findContiguousLines(points, 2);
        
        assertEquals(6, lines.size());
        
        int count = 0;
        for (Set<PairInt> line : lines.keySet()) {
            for (PairInt p : line) {
                boolean rm = points.remove(p);
                assertTrue(rm);
                count++;
            }
        }
        
        // allow for corners not being in the lines:
        assertTrue(points.isEmpty());
        */
    }
}
