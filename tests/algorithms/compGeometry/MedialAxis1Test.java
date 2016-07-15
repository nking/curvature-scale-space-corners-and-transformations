package algorithms.compGeometry;

import algorithms.compGeometry.MedialAxis1.MedialAxisResults;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
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
        
        int status = medAxis1.intersectsMedialAxis(nearestB, p, output);
        
        assertEquals(1, status);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(22, 8));
        expected.add(new PairInt(14, 6));
        expected.add(new PairInt(20, 3));
        
        for (MedialAxis1.MedialAxisPoint mp : output) {
            MedialAxis1.PVector[] vs = mp.getVectors();
            for (MedialAxis1.PVector v : vs) {
                System.out.println("medial axis pt=" + v.getPoint());
                int x = v.getPoint().getX();
                int y = v.getPoint().getY();
                PairInt rm = null;
                for (PairInt p2 : expected) {
                    if ((Math.abs(x - p2.getX()) < 2) &&
                        (Math.abs(y - p2.getY()) < 2)) {
                        rm = p2;
                        break;
                    }
                }
                assertTrue(expected.remove(rm));
            }
        }
        assertTrue(expected.isEmpty());
   
        // --------- test status==-2 -----
        p = new PairInt(14, 2);

        nearestB = medAxis1.getNearestBoundaryPoints(p);        
        
        status = medAxis1.intersectsMedialAxis(nearestB, p, output);
        
        assertEquals(-2, status);
      
        MedialAxisResults results2 = medAxis1.findMedialAxesNearPoint(p);
        assertNotNull(results2);
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(18, 5));
        expected.add(new PairInt(10, 5));
        
        assertEquals(14, results2.centerAndDistance.p.getX());
        assertEquals(4, results2.centerAndDistance.p.getY());
        assertTrue(Math.abs(results2.centerAndDistance.delta - 4) < 0.1);
        output = results2.medialAxes;
        
        for (MedialAxis1.MedialAxisPoint mp : output) {
            MedialAxis1.PVector[] vs = mp.getVectors();
            for (MedialAxis1.PVector v : vs) {
                System.out.println("*medial axis pt=" + v.getPoint());
                int x = v.getPoint().getX();
                int y = v.getPoint().getY();
                PairInt rm = null;
                for (PairInt p2 : expected) {
                    if ((Math.abs(x - p2.getX()) < 2) &&
                        (Math.abs(y - p2.getY()) < 2)) {
                        rm = p2;
                        break;
                    }
                }
                assertTrue(expected.remove(rm));
            }
        }
        assertTrue(expected.isEmpty());
        
        // ---- test removing one of the results from the points ----
        Set<PairInt> removed = medAxis1.subtractFromPoints(
            results2.centerAndDistance.p, results2.centerAndDistance.delta);
        
        assertTrue(removed.size() > 1);
        
        // -- just running public method, until method is testable
        medAxis1 = new MedialAxis1(points, border);
        
        medAxis1.findMedialAxis();
        
        LinkedList<MedialAxis1.MedialAxisPoint>
            list = medAxis1.getMedAxisList();
        
        Set<PairInt> mAPs = new HashSet<PairInt>();
        for (MedialAxis1.MedialAxisPoint mp : list) {
            PairInt pp = mp.getVectors()[0].getPoint();
            System.out.println("**med axis pt = " + pp);
            assertFalse(mAPs.contains(pp));
            mAPs.add(pp);
        }
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
    
    public void test2() {
        
        List<PairInt> border = new ArrayList<PairInt>();
        Set<PairInt> points = new HashSet<PairInt>();
        
        populateTestData2(border, points);
 
        /*
        adding hole to data0 in center of object
           _____________
           |           |
           |___________| height=11
           <-----24---->
                                  -              @ 10
                                              @    9
                                  *        @       8
                         .      ?       ?          7
                      / -- \         @             6
                                                   5
                      \ __ /         @             4
                         .              @          3
                                           @       2
                                              @    1
                                                 @ 0
          10 11 12 13 14 15 16 17 18 19 20 21 22 23
        */
        
        MedialAxis1 medAxis1 = new MedialAxis1(points, border);

        PairInt p = new PairInt(18, 8);

        // --- checking nearest neighbors ---
        Set<PairInt> nearestB = medAxis1.getNearestBoundaryPoints(p);        
                
        List<MedialAxis1.MedialAxisPoint> output = new
            ArrayList<MedialAxis1.MedialAxisPoint>();
        
        medAxis1.intersectsMedialAxis(nearestB, p, output);
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(17, 7));
        expected.add(new PairInt(20, 7));
        
        for (MedialAxis1.MedialAxisPoint mp : output) {
            MedialAxis1.PVector[] vs = mp.getVectors();
            for (MedialAxis1.PVector v : vs) {
                System.out.println("medial axis pt=" + v.getPoint());
                int x = v.getPoint().getX();
                int y = v.getPoint().getY();
                PairInt rm = null;
                for (PairInt p2 : expected) {
                    if ((Math.abs(x - p2.getX()) < 2) &&
                        (Math.abs(y - p2.getY()) < 2)) {
                        rm = p2;
                        break;
                    }
                }
                assertTrue(expected.remove(rm));
            }
        }
        assertTrue(expected.isEmpty());
        
    }
    
    /*
    adding hole to data0 in center of object
       _____________
       |           |
       |___________| height=11
       <-----24---->
                                               10
                                               9
                                               8
                     .                         7
                  / -- \                       6
                 |      |                      5
                  \ __ /                       4
                     .                         3
                                               2
                                               1
                                               0
      10 11 12 13 14 15 16 17 18 19 20 21 22 23
    */
    private void populateTestData2(List<PairInt> border, 
        Set<PairInt> areaPoints) {
    
        populateTestData0(border, areaPoints);
        
        areaPoints.remove(new PairInt(14, 5));
        areaPoints.remove(new PairInt(15, 5));
        areaPoints.remove(new PairInt(16, 5));
        areaPoints.remove(new PairInt(15, 4));
        areaPoints.remove(new PairInt(15, 6));
    
        border.add(new PairInt(13, 5));
        border.add(new PairInt(14, 6));
        border.add(new PairInt(15, 7));
        border.add(new PairInt(16, 6));
        border.add(new PairInt(17, 5));
        border.add(new PairInt(16, 4));
        border.add(new PairInt(15, 3));
        border.add(new PairInt(14, 4));
    }
    
}
