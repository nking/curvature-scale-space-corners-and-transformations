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
        
        populateTestData0(border, points);
 
        /*
        Figure 2 of Yang, Brock, and Moll
           _____________
           |           |
           |___________| height=11
           <-----24---->
                                  -              @ 10
                                           ?  @    9
                                           @       8
                                        @     ?    7
                       ?          *  @             6
            @@@@@@@@@@@?@@@@@@@@@@@@               5
                                     @             4
                                        ?  ?       3
                                     ?     @       2
                                              @    1
                                                 @ 0
          10 11 12 13 14 15 16 17 18 19 20 21 22 23
        */
        
        MedialAxis1 medAxis1 = new MedialAxis1(points, border);

        PairInt p = new PairInt(18, 6);

        // --- checking nearest neighbors ---
        Set<PairInt> nearestB = medAxis1.getNearestBoundaryPoints(p);        
        
        Set<PairInt> expectedNearest = new HashSet<PairInt>();
        expectedNearest.add(new PairInt(18, 10));
        for (PairInt npb : nearestB) {
            assertTrue(expectedNearest.remove(npb));
        }
        assertTrue(expectedNearest.isEmpty());
        
        Set<PairInt> nearestB2 = medAxis1.getNearestBoundaryPoints(new PairInt(22, 7));        
        expectedNearest = new HashSet<PairInt>();
        expectedNearest.add(new PairInt(23, 7));
        for (PairInt npb : nearestB2) {
            assertTrue(expectedNearest.remove(npb));
        }
        assertTrue(expectedNearest.isEmpty());
        //------
        
        List<MedialAxis1.MedialAxisPoint> output = new
            ArrayList<MedialAxis1.MedialAxisPoint>();
        
        medAxis1.intersectsMedialAxis(nearestB, p, output);
        
        // 23,9  within maxX and maxY of 23,10
        //207, 230   for 228
        int z = 1;
    }
    
    /*
    Figure 2 of Yang, Brock, and Moll

           _____________
           |           |
           |___________| height=11
           <-----24---->
    */
    private void populateTestData0(List<PairInt> border, 
        Set<PairInt> areaPoints) {
    
        // clockwise
        // 23,0 -> 0,0
        for (int i = 23; i > -1; --i) {
            border.add(new PairInt(i, 0));
        }
        // 0,1 -> 0,10
        for (int i = 1; i < 11; ++i) {
            border.add(new PairInt(0, i));
        }
        // (1,10)->(23,10)
        for (int i = 1; i < 24; ++i) {
            border.add(new PairInt(i, 10));
        }
        // 23,9 23,1
        for (int i = 9; i > 0; --i) {
            border.add(new PairInt(23, i));
        }
        // interior points.  this includes border, which is fine
        // because MedialAxis1 removes the border from them
        for (int i = 0; i < 24; ++i) {
            for (int j = 0; j < 11; ++j) {
                areaPoints.add(new PairInt(i, j));
            }
        }
    }
    
    public void est1() {
        
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
        //3,6
        medAxis1.intersectsMedialAxis(nearestB, p, output);
        
        int z = 1;
    }
    
    /*
    Figure 4 of 
    "Efficient Computation of A Simplified Medial Axis"
        by Foskey, Lin, and Manocha
    
               ___        19,20,21
               | |        16,17,18
      _________| |_______ 15
      |                 | 11
      |                 | 8   <---- y=7 for a long medial axis segment
      |                 | 5
      |_________________| 0-2
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
        // 0,1 -> 0,15
        for (int i = 1; i < 15; ++i) {
            border.add(new PairInt(0, i));
        }
        // (1,15)->(8,15)
        for (int i = 1; i < 9; ++i) {
            border.add(new PairInt(i, 15));
        }
        // (9,15)->(9,21)
        for (int i = 15; i <= 21; ++i) {
            border.add(new PairInt(9, i));
        }
        // 10,21 -> 11,21
        for (int i = 10; i <= 11; ++i) {
            border.add(new PairInt(i, 21));
        }
        // 11,21 -> 11,15
        for (int i = 11; i >= 15; --i) {
            border.add(new PairInt(11, i));
        }
        // 12,15 -> 18,15
        for (int i = 12; i < 19; ++i) {
            border.add(new PairInt(i, 15));
        }
        // 18,14 18,1
        for (int i = 14; i > 0; --i) {
            border.add(new PairInt(18, i));
        }
        // interior points.  this includes border, which is fine
        // because MedialAxis1 removes the border from them
        for (int i = 0; i < 19; ++i) {
            for (int j = 0; j <= 15; ++j) {
                areaPoints.add(new PairInt(i, j));
            }
        }
        for (int i = 9; i <= 11; ++i) {
            for (int j = 9; j <= 21; ++j) {
                areaPoints.add(new PairInt(i, j));
            }
        }
    }
}
