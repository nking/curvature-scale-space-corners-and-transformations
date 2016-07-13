package algorithms.compGeometry;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MedialAxis1Test extends TestCase {
    
    public MedialAxis1Test() {
    }
    
    public void test0() {
        
        List<PairInt> border = new ArrayList<PairInt>();
        Set<PairInt> points = new HashSet<PairInt>();
        
        populateTestData(border, points);

        /*
        Figure 4 of 
        "Efficient Computation of A Simplified Medial Axis"
            by Foskey, Lin, and Manocha

                   ___        11,12
                   | |        10,11
          _________| |_______ 8,9
          |                 | 7,8
          |                 | 5,6
          |                 | 3,4
          |                 | 2,3
          |_________________| 0,1
          0123456789012345678
                    1
        */
        
        MedialAxis1 medAxis1 = new MedialAxis1(points, border);

        PairInt p = new PairInt(3, 3);

        Set<PairInt> nearestB = medAxis1.getNearestBoundaryPoints(p);        
        
        // assert nearestB are 3,0 and 0,3
        Set<PairInt> expectedNearest = new HashSet<PairInt>();
        expectedNearest.add(new PairInt(3, 0));
        expectedNearest.add(new PairInt(0, 3));
        for (PairInt npb : nearestB) {
            assertTrue(expectedNearest.remove(npb));
        }
        assertTrue(expectedNearest.isEmpty());
        
        List<MedialAxis1.MedialAxisPoint> output = new
            ArrayList<MedialAxis1.MedialAxisPoint>();
        //6,5  6,2
        medAxis1.intersectsMedialAxis(nearestB, p, output);
        
        int z = 1;
    }
    
    /*
    Figure 4 of 
    "Efficient Computation of A Simplified Medial Axis"
        by Foskey, Lin, and Manocha
    
               ___        11,12
               | |        9,10
      _________| |_______ 7,8
      |                 | 5,6
      |                 | 3,4   <---- y=4 for a long medial axis segment
      |                 | 2,3
      |_________________| 0,1
      0123456789012345678
                1
    */
    private void populateTestData(List<PairInt> border, 
        Set<PairInt> areaPoints) {
    
        // clockwise
        // 18,0 -> 0,0
        for (int i = 18; i > -1; --i) {
            border.add(new PairInt(i, 0));
        }
        // 0,1 -> 0,8
        for (int i = 1; i < 9; ++i) {
            border.add(new PairInt(0, i));
        }
        // (1,8)->(8,8)
        for (int i = 1; i < 9; ++i) {
            border.add(new PairInt(i, 8));
        }
        // (9,8)->(9,12)
        for (int i = 8; i <= 12; ++i) {
            border.add(new PairInt(9, i));
        }
        // 10,12 -> 11,12
        for (int i = 10; i <= 11; ++i) {
            border.add(new PairInt(i, 12));
        }
        // 11,11 -> 11,8
        for (int i = 11; i >= 8; --i) {
            border.add(new PairInt(11, i));
        }
        // 12,8 -> 18,8
        for (int i = 12; i < 19; ++i) {
            border.add(new PairInt(i, 8));
        }
        // 18,7 18,1
        for (int i = 7; i > 0; --i) {
            border.add(new PairInt(18, i));
        }
        // interior points.  this includes border, which is fine
        // because MedialAxis1 removes the border from them
        for (int i = 0; i < 19; ++i) {
            for (int j = 0; j <= 8; ++j) {
                areaPoints.add(new PairInt(i, j));
            }
        }
        for (int i = 9; i <= 11; ++i) {
            for (int j = 9; j <= 12; ++j) {
                areaPoints.add(new PairInt(i, j));
            }
        }
    }
}
