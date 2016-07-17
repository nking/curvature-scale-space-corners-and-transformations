package algorithms.compGeometry;

import algorithms.compGeometry.MedialAxis1.MedialAxisResults;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
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
    
    public void test0() throws IOException {
        
        MedialAxis1 medAxis1 = null;
        Set<PairInt> mAPs = null;
        LinkedList<MedialAxis1.MedialAxisPoint> list = null;
        Set<PairInt> expected = null;
        
        Set<PairInt> border = new HashSet<PairInt>();
        Set<PairInt> points = new HashSet<PairInt>();
        
        populateTestData0(border, points);
 
        
        medAxis1 = new MedialAxis1(points, border);
        
        medAxis1.findMedialAxis();
        
        list = medAxis1.getMedAxisList();
        
        mAPs = new HashSet<PairInt>();
        for (MedialAxis1.MedialAxisPoint mp : list) {
            PairInt pp = mp.getCenter();
            //System.out.println("**med axis pt = " + pp);
            assertFalse(mAPs.contains(pp));
            mAPs.add(pp);
        }
        
        /*
        float[] x = new float[mAPs.size()];
        float[] y = new float[mAPs.size()];
        int count = 0;
        for (PairInt p2 : mAPs) {
            x[count] = p2.getX();
            y[count] = p2.getY();
            count++;
        }
        PolygonAndPointPlotter plotter 
            = new PolygonAndPointPlotter(0, 23, 0, 10);
        float[] xp = null;
        float[] yp = null;
        plotter.addPlot(x, y, xp, yp, "med axis");
        plotter.writeFile();
        */
            
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
        
        // add expected points
        expected = new HashSet<PairInt>();
        for (int i = 1; i < 6; ++i) {
            expected.add(new PairInt(i, i));
            expected.add(new PairInt(i, 10 - i));
            expected.add(new PairInt(17 + i, 4 + i));
            expected.add(new PairInt(23 - i, i));
        }
        for (int i = 6; i < 18; ++i) {
            expected.add(new PairInt(i, 5));
        }
        
        for (PairInt p2 : mAPs) {
            //System.out.println(" ->med axis pt = " + p2);
            assertTrue(expected.remove(p2));
        }
        assertTrue(expected.isEmpty());
        
        if (!(medAxis1.getSepAng() == Math.PI/9 && medAxis1.getNSampl()
            == 18)) {
            return;
        }
        
        expected.clear();
        
        medAxis1 = new MedialAxis1(points, border);

        PairInt p = new PairInt(18, 6);

        // --- checking nearest neighbors ---
        Set<PairInt> nearestB = medAxis1.getNearestBoundaryPoints(p);        
        
        Set<PairInt> expectedNearest = new HashSet<PairInt>();
        expectedNearest.add(new PairInt(18, 10));
        for (PairInt npb : nearestB) {
            assertTrue(expectedNearest.remove(npb));
        }
        assertTrue(expectedNearest.isEmpty());
        
        Set<PairInt> nearestB2 = 
            medAxis1.getNearestBoundaryPoints(new PairInt(22, 7));        
        expectedNearest = new HashSet<PairInt>();
        expectedNearest.add(new PairInt(23, 7));
        for (PairInt npb : nearestB2) {
            assertTrue(expectedNearest.remove(npb));
        }
        assertTrue(expectedNearest.isEmpty());
        //------
        
        List<MedialAxis1.MedialAxisPoint> output = new
            ArrayList<MedialAxis1.MedialAxisPoint>();
        
        nearestB2 = medAxis1.getNearestBoundaryPoints(p);
        int status = medAxis1.intersectsMedialAxis(
            medAxis1.findNearestBoundsAsArray(p.getX(), p.getY()), 
            p, output);
        
        assertEquals(1, status);
        
        //18,6 w/ r=4
        expected.clear();
        expected.add(new PairInt(22, 8));
        expected.add(new PairInt(14, 6));
        expected.add(new PairInt(20, 3));
        
        for (MedialAxis1.MedialAxisPoint mp : output) {
            PairInt center = mp.getCenter();
            //System.out.println("medial axis pt=" + center);
            int x1 = center.getX();
            int y1 = center.getY();
            PairInt rm = null;
            for (PairInt p2 : expected) {
                if ((Math.abs(x1 - p2.getX()) < 2) &&
                    (Math.abs(y1 - p2.getY()) < 2)) {
                    rm = p2;
                    break;
                }
            }
            assertTrue(expected.remove(rm));
        }
        assertTrue(expected.isEmpty());
   
        // --------- test status==-2 -----
        p = new PairInt(14, 2);

        nearestB = medAxis1.getNearestBoundaryPoints(p);        
        
        status = medAxis1.intersectsMedialAxis(
            medAxis1.findNearestBoundsAsArray(p.getX(), p.getY()),
            p, output);
        
        assertEquals(-2, status);
      
        MedialAxisResults results2 = medAxis1.findMedialAxesNearPoint(p);
        assertNotNull(results2);
        
        // depending upon changes to search pattern while debugging,
        // the point actually searched is changing
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(12, 6));
        expected.add(new PairInt(16, 5));
        expected.add(new PairInt(17, 6));
        
        assertEquals(14, results2.center.getX());
        assertEquals(7, results2.center.getY());
        assertTrue(Math.abs(results2.centerSrchR - 3) < 0.1);
        output = results2.medialAxes;
        
        for (MedialAxis1.MedialAxisPoint mp : output) {
            PairInt center = mp.getCenter();
            //System.out.println("*medial axis pt=" + center);
            int x1 = center.getX();
            int y1 = center.getY();
            PairInt rm = null;
            for (PairInt p2 : expected) {
                if ((Math.abs(x1 - p2.getX()) < 2) &&
                    (Math.abs(y1 - p2.getY()) < 2)) {
                    rm = p2;
                    break;
                }
            }
            assertTrue(expected.remove(rm));
        }
        assertTrue(expected.isEmpty());
        
        // ---- test removing one of the results from the points ----
        Set<PairInt> removed = medAxis1.subtractFromPoints(
            results2.center, results2.centerSrchR);
        
        assertTrue(removed.size() > 1);
        
        medAxis1 = new MedialAxis1(points, border);
        
        medAxis1.findMedialAxis();
        
        list = medAxis1.getMedAxisList();
        
        mAPs = new HashSet<PairInt>();
        for (MedialAxis1.MedialAxisPoint mp : list) {
            PairInt pp = mp.getCenter();
            //System.out.println("**med axis pt = " + pp);
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
    private void populateTestData0(Collection<PairInt> border, 
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
    
    public void test1() throws IOException {
                
        Set<PairInt> border = new HashSet<PairInt>();
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
        medAxis1.intersectsMedialAxis(
            medAxis1.findNearestBoundsAsArray(p.getX(), p.getY()), 
            p, output);

        // -- visualize the results ----        
        medAxis1 = new MedialAxis1(points, border);
        
        medAxis1.findMedialAxis();
        
        Set<PairInt> mSet = medAxis1.getMedialAxisPoints();
        
        float[] x = new float[mSet.size()];
        float[] y = new float[mSet.size()];
        int count = 0;
        for (PairInt p2 : mSet) {
            x[count] = p2.getX();
            y[count] = p2.getY();
            count++;
        }
        
        PolygonAndPointPlotter plotter 
            = new PolygonAndPointPlotter(0, 18, 0, 21);
        float[] xp = null;
        float[] yp = null;
        plotter.addPlot(x, y, xp, yp, "med axis");
        plotter.writeFile(20);
        
        // --- compare to erosion filter and line thinner ----
        border.clear();
        points.clear();
        populateTestData(border, points);
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(points, 0, 18, 0, 21);
        x = new float[points.size()];
        y = new float[points.size()];
        count = 0;
        for (PairInt p2 : points) {
            x[count] = p2.getX();
            y[count] = p2.getY();
            count++;
        }
        plotter = new PolygonAndPointPlotter(0, 18, 0, 21);
        xp = null;
        yp = null;
        plotter.addPlot(x, y, xp, yp, "med axis");
        plotter.writeFile(21);
        
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
    private void populateTestData(Collection<PairInt> border, 
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
    
    public void test2() throws IOException {
                
        MedialAxis1 medAxis1 = null;
        Set<PairInt> mAPs = null;
        LinkedList<MedialAxis1.MedialAxisPoint> list = null;
        Set<PairInt> expected = null;
        
        Set<PairInt> border = new HashSet<PairInt>();
        Set<PairInt> points = new HashSet<PairInt>();
        
        populateTestData2(border, points);
 
        
        medAxis1 = new MedialAxis1(points, border);
        
        medAxis1.findMedialAxis();
        
        Set<PairInt> mSet = medAxis1.getMedialAxisPoints();
        
        /*
        float[] x = new float[mSet.size()];
        float[] y = new float[mSet.size()];
        int count = 0;
        for (PairInt p2 : mSet) {
            x[count] = p2.getX();
            y[count] = p2.getY();
            count++;
        }
        
        PolygonAndPointPlotter plotter 
            = new PolygonAndPointPlotter(0, 23, 0, 10);
        float[] xp = null;
        float[] yp = null;
        plotter.addPlot(x, y, xp, yp, "med axis");
        plotter.writeFile(10);
        
        // --- compare to erosion filter and line thinner ----
        border.clear();
        points.clear();
        populateTestData2(border, points);
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(points, 0, 23, 0, 10);
        x = new float[points.size()];
        y = new float[points.size()];
        count = 0;
        for (PairInt p2 : points) {
            x[count] = p2.getX();
            y[count] = p2.getY();
            count++;
        }
        plotter = new PolygonAndPointPlotter(0, 23, 0, 10);
        xp = null;
        yp = null;
        plotter.addPlot(x, y, xp, yp, "med axis");
        plotter.writeFile(11);
        */    
        
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
        
        if (!(medAxis1.getSepAng() == Math.PI/9 && medAxis1.getNSampl()
            == 18)) {
            return;
        }
        
        medAxis1 = new MedialAxis1(points, border);

        PairInt p = new PairInt(18, 8);

        // --- checking nearest neighbors ---
        Set<PairInt> nearestB = medAxis1.getNearestBoundaryPoints(p);        
                
        List<MedialAxis1.MedialAxisPoint> output = new
            ArrayList<MedialAxis1.MedialAxisPoint>();
        
        medAxis1.intersectsMedialAxis(
            medAxis1.findNearestBoundsAsArray(p.getX(), p.getY()),
            p, output);
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(17, 7));
        expected.add(new PairInt(20, 7));
        
        for (MedialAxis1.MedialAxisPoint mp : output) {
            PairInt center = mp.getCenter();
            //System.out.println("medial axis pt=" + center);
            int x1 = center.getX();
            int y1 = center.getY();
            PairInt rm = null;
            for (PairInt p2 : expected) {
                if ((Math.abs(x1 - p2.getX()) < 2) &&
                    (Math.abs(y1 - p2.getY()) < 2)) {
                    rm = p2;
                    break;
                }
            }
            assertTrue(expected.remove(rm));
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
    private void populateTestData2(Collection<PairInt> border, 
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
