package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class EdgeExtractorTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EdgeExtractorTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void testMergeAdjacentEndPoints() throws Exception {
        
        System.out.println("testMergeAdjacentEndPoints");
           
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        
        PairIntArray edge0 = new PairIntArray();
        edge0.add(10, 10);
        edge0.add(11, 10);
        edge0.add(12, 10);
        
        for (int orientation = 0; orientation < 8; orientation ++) {
            
            PairIntArray edge1 = new PairIntArray();
            
            PairIntArray expectedOutput = new PairIntArray();
            
            switch(orientation) {
                case 0:
                    edge1.add(13, 10);
                    edge1.add(14, 10);
                    edge1.add(15, 10);
                    expectedOutput = edge0.copy();
                    expectedOutput.add(13, 10);
                    expectedOutput.add(14, 10);
                    expectedOutput.add(15, 10);
                    break;
                case 1:
                    edge1.add(13, 9);
                    edge1.add(14, 9);
                    edge1.add(15, 9);
                    expectedOutput = edge0.copy();
                    expectedOutput.add(13, 9);
                    expectedOutput.add(14, 9);
                    expectedOutput.add(15, 9);
                    break;
                case 2:
                    edge1.add(12, 9);
                    edge1.add(13, 9);
                    edge1.add(14, 9);
                    expectedOutput = edge0.copy();
                    expectedOutput.add(12, 9);
                    expectedOutput.add(13, 9);
                    expectedOutput.add(14, 9);
                    break;
                case 3:
                    edge1.add(11, 9);
                    edge1.add(12, 9);
                    edge1.add(13, 9);
                    expectedOutput = edge0.copy();
                    expectedOutput.add(11, 9);
                    expectedOutput.add(12, 9);
                    expectedOutput.add(13, 9);
                    break;
                case 4:
                    // reversed list to be inserted int front of edge0
                    edge1.add(9, 10);
                    edge1.add(8, 10);
                    edge1.add(7, 10);
                    expectedOutput.add(12, 10);
                    expectedOutput.add(11, 10);
                    expectedOutput.add(10, 10);
                    expectedOutput.add(9, 10);
                    expectedOutput.add(8, 10);
                    expectedOutput.add(7, 10);
                    expectedOutput.reverse();
                    break;
                case 5:
                    // forward list to be inserted int front of edge0
                    edge1.add(7, 11);
                    edge1.add(8, 11);
                    edge1.add(9, 11);
                    expectedOutput.add(7, 11);
                    expectedOutput.add(8, 11);
                    expectedOutput.add(9, 11);
                    expectedOutput.add(10, 10);
                    expectedOutput.add(11, 10);
                    expectedOutput.add(12, 10);
                    break;
                case 6:
                    // reversed list to be inserted int front of edge0
                    edge1.add(10, 11);
                    edge1.add(9, 11);
                    edge1.add(8, 11);
                    expectedOutput.add(12, 10);
                    expectedOutput.add(11, 10);
                    expectedOutput.add(10, 10);
                    expectedOutput.add(10, 11);
                    expectedOutput.add(9, 11);
                    expectedOutput.add(8, 11);
                    expectedOutput.reverse();
                    break;
                case 7:
                    // forward list to be inserted int front of edge0
                    edge1.add(9, 11);
                    edge1.add(10, 11);
                    edge1.add(11, 11);
                    expectedOutput.add(12, 10);
                    expectedOutput.add(11, 10);
                    expectedOutput.add(10, 10);
                    expectedOutput.add(9, 11);
                    expectedOutput.add(10, 11);
                    expectedOutput.add(11, 11);
                    expectedOutput.reverse();
                    break;
            }
            
            System.out.println("orientation case=" + orientation);
            
            List<PairIntArray> edges = new ArrayList<PairIntArray>();
            edges.add(edge0.copy());
            edges.add(edge1.copy());

            List<PairIntArray> output = contourExtractor.mergeAdjacentEndPoints(
                edges);

            assertTrue(output.size() == 1);
            assertTrue(output.get(0).getN() == (edge0.getN() + edge1.getN()));

            boolean[] found0 = new boolean[edge0.getN()];
            boolean[] found1 = new boolean[edge1.getN()];
            for (int i = 0; i < output.get(0).getN(); i++) {
                int x = output.get(0).getX(i);
                int y = output.get(0).getY(i);
                for (int ii = 0; ii < edge0.getN(); ii++) {
                    int xx = edge0.getX(ii);
                    int yy = edge0.getY(ii);
                    if ((xx == x) && (yy == y)) {
                        found0[ii] = true;
                        break;
                    }
                }
                for (int ii = 0; ii < edge1.getN(); ii++) {
                    int xx = edge1.getX(ii);
                    int yy = edge1.getY(ii);
                    if ((xx == x) && (yy == y)) {
                        found1[ii] = true;
                        break;
                    }
                }
            }
            for (int i = 0; i < found0.length; i++) {
                assertTrue(found0[i]);
            }
            for (int i = 0; i < found1.length; i++) {
                assertTrue(found1[i]);
            }
            
            assertTrue(expectedOutput.getN() == output.get(0).getN());
                        
            for (int i = 0; i < output.get(0).getN(); i++) {
                int x = output.get(0).getX(i);
                int y = output.get(0).getY(i);
                int xx = expectedOutput.getX(i);
                int yy = expectedOutput.getY(i);
                assertTrue((xx == x) && (yy == y));
            }
        }      
    }
    
    public void testFillInGaps() throws Exception {
        /*
        -- Test these gaps from the africa.png edges
        [junit] 0) SEE_4 (333,141) --- (149,43)
        [junit] 1) (308,338) --- (187,309)
        [junit] 2) (168,202) --- (35,150)
        [junit] 3) (43,114) --- (144,45)
        [junit] 4) (390,198) --- (335,143) GAP
        [junit] 5) (387,205) --- (336,253)
        [junit] 6) SEE_11 (197,288) --- (171,233)
        [junit] 7) (333,312) --- (313,327)
        [junit] 8) (179,211) --- (172,230)
        [junit] 9) SEE_13 (342,303) --- (336,262)
        [junit] 10) SEE_12 (45,135) --- (38,147)
        [junit] 11) GAP (197,290) --- (187,305)
        [junit] 12) GAP (45,133) --- (43,126)
        [junit] 13) GAP (342,305) --- (338,309)
        [junit] 14) (45,119) --- (45,123)
        */
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(446, 434));
        
        PairIntArray edge0 = new PairIntArray();
        edge0.add(333, 141);
        //...
        edge0.add(149, 43);
        
        PairIntArray edge1 = new PairIntArray();
        edge1.add(390, 198);
        //...
        edge1.add(335, 143);
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        List<PairIntArray> output = contourExtractor.fillInGaps(edges);
        
        assertTrue(output.size() == 1);
        assertTrue(output.get(0).getN() == (edge0.getN() + edge1.getN() + 1));

        //assert gap was filled:
        assertTrue(contourExtractor.getImage().getValue(334, 142) > 0);
        
        boolean[] found0 = new boolean[edge0.getN()];
        boolean[] found1 = new boolean[edge1.getN()];
        boolean foundFilledGap = false;
        for (int i = 0; i < output.get(0).getN(); i++) {
            int x = output.get(0).getX(i);
            int y = output.get(0).getY(i);
            for (int ii = 0; ii < edge0.getN(); ii++) {
                int xx = edge0.getX(ii);
                int yy = edge0.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found0[ii] = true;
                    break;
                }
            }
            for (int ii = 0; ii < edge1.getN(); ii++) {
                int xx = edge1.getX(ii);
                int yy = edge1.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found1[ii] = true;
                    break;
                }
            }
            if ((x == 334) && (y == 142)) {
                foundFilledGap = true;
            }
        }
        
        assertTrue(foundFilledGap);
        
        for (int i = 0; i < found0.length; i++) {
            assertTrue(found0[i]);
        }
        for (int i = 0; i < found1.length; i++) {
            assertTrue(found1[i]);
        }
        
        /*
        -- Test these gaps from the africa.png edges
        [junit] 6) SEE_11 (197,288) --- (171,233)
        [junit] 7) (333,312) --- (313,327)
        [junit] 8) (179,211) --- (172,230)
        [junit] 9) SEE_13 (342,303) --- (336,262)
        [junit] 10) SEE_12 (45,135) --- (38,147)
        [junit] 11) GAP (197,290) --- (187,305)
        [junit] 12) GAP (45,133) --- (43,126)
        [junit] 13) GAP (342,305) --- (338,309)
        [junit] 14) (45,119) --- (45,123)
        */
        edge0 = new PairIntArray();
        edge0.add(197, 288);
        //...
        edge0.add(171, 233);
        
        edge1 = new PairIntArray();
        edge1.add(197,290);
        //...
        edge1.add(187,305);
        
        edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        output = contourExtractor.fillInGaps(edges);
        
        assertTrue(output.size() == 1);
        assertTrue(output.get(0).getN() == (edge0.getN() + edge1.getN() + 1));

        //assert gap was filled:
        assertTrue(contourExtractor.getImage().getValue(197, 289) > 0);
        
        found0 = new boolean[edge0.getN()];
        found1 = new boolean[edge1.getN()];
        foundFilledGap = false;
        for (int i = 0; i < output.get(0).getN(); i++) {
            int x = output.get(0).getX(i);
            int y = output.get(0).getY(i);
            for (int ii = 0; ii < edge0.getN(); ii++) {
                int xx = edge0.getX(ii);
                int yy = edge0.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found0[ii] = true;
                    break;
                }
            }
            for (int ii = 0; ii < edge1.getN(); ii++) {
                int xx = edge1.getX(ii);
                int yy = edge1.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found1[ii] = true;
                    break;
                }
            }
            if ((x == 197) && (y == 289)) {
                foundFilledGap = true;
            }
        }
        
        assertTrue(foundFilledGap);
        
        for (int i = 0; i < found0.length; i++) {
            assertTrue(found0[i]);
        }
        for (int i = 0; i < found1.length; i++) {
            assertTrue(found1[i]);
        }
        
        /*
        -- Test these gaps from the africa.png edges
        [junit] 9) SEE_13 (342,303) --- (336,262)
        [junit] 10) SEE_12 (45,135) --- (38,147)
        [junit] 12) GAP (45,133) --- (43,126)
        [junit] 13) GAP (342,305) --- (338,309)
        [junit] 14) (45,119) --- (45,123)
        */
        edge0 = new PairIntArray();
        edge0.add(342, 303);
        //...
        edge0.add(336, 262);
        
        edge1 = new PairIntArray();
        edge1.add(342, 305);
        //...
        edge1.add(338, 309);
        
        edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        output = contourExtractor.fillInGaps(edges);
        
        assertTrue(output.size() == 1);
        assertTrue(output.get(0).getN() == (edge0.getN() + edge1.getN() + 1));

        //assert gap was filled:
        assertTrue(contourExtractor.getImage().getValue(342, 304) > 0);
        
        found0 = new boolean[edge0.getN()];
        found1 = new boolean[edge1.getN()];
        foundFilledGap = false;
        for (int i = 0; i < output.get(0).getN(); i++) {
            int x = output.get(0).getX(i);
            int y = output.get(0).getY(i);
            for (int ii = 0; ii < edge0.getN(); ii++) {
                int xx = edge0.getX(ii);
                int yy = edge0.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found0[ii] = true;
                    break;
                }
            }
            for (int ii = 0; ii < edge1.getN(); ii++) {
                int xx = edge1.getX(ii);
                int yy = edge1.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found1[ii] = true;
                    break;
                }
            }
            if ((x == 342) && (y == 304)) {
                foundFilledGap = true;
            }
        }
        
        assertTrue(foundFilledGap);
        
        for (int i = 0; i < found0.length; i++) {
            assertTrue(found0[i]);
        }
        for (int i = 0; i < found1.length; i++) {
            assertTrue(found1[i]);
        }
        
        /*
        -- Test these gaps from the africa.png edges
        [junit] 10) SEE_12 (45, 135) --- (38, 147)
        [junit] 12) GAP (45, 133) --- (43, 126)
        */
        edge0 = new PairIntArray();
        edge0.add(45, 135);
        //...
        edge0.add(38, 147);
        
        edge1 = new PairIntArray();
        edge1.add(45, 133);
        //...
        edge1.add(43, 126);
        
        edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        output = contourExtractor.fillInGaps(edges);
        
        assertTrue(output.size() == 1);
        assertTrue(output.get(0).getN() == (edge0.getN() + edge1.getN() + 1));

        //assert gap was filled:
        assertTrue(contourExtractor.getImage().getValue(45, 134) > 0);
        
        found0 = new boolean[edge0.getN()];
        found1 = new boolean[edge1.getN()];
        foundFilledGap = false;
        for (int i = 0; i < output.get(0).getN(); i++) {
            int x = output.get(0).getX(i);
            int y = output.get(0).getY(i);
            for (int ii = 0; ii < edge0.getN(); ii++) {
                int xx = edge0.getX(ii);
                int yy = edge0.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found0[ii] = true;
                    break;
                }
            }
            for (int ii = 0; ii < edge1.getN(); ii++) {
                int xx = edge1.getX(ii);
                int yy = edge1.getY(ii);
                if ((xx == x) && (yy == y)) {
                    found1[ii] = true;
                    break;
                }
            }
            if ((x == 45) && (y == 134)) {
                foundFilledGap = true;
            }
        }
        
        assertTrue(foundFilledGap);
        
        for (int i = 0; i < found0.length; i++) {
            assertTrue(found0[i]);
        }
        for (int i = 0; i < found1.length; i++) {
            assertTrue(found1[i]);
        }
    }
    
    public void testFindClosestPair() throws Exception {
        
        System.out.println("testFindClosestPair");
        
        PairIntArray curve0 = new PairIntArray();
        curve0.add(45, 133);
        curve0.add(45, 132);
        curve0.add(45, 131);
        curve0.add(45, 130);
        curve0.add(45, 129);
        curve0.add(45, 128);
        curve0.add(45, 127);
        curve0.add(44, 126);
        curve0.add(43, 125);
        curve0.add(42, 125);
        curve0.add(43, 124);
        curve0.add(44, 123);//11
        curve0.add(44, 124);
        curve0.add(44, 125);
        curve0.add(43, 126);
        curve0.add(45, 126);
        curve0.add(45, 125);
         
        PairIntArray curve1 = new PairIntArray();
        curve1.add(45, 123); // 0
        curve1.add(45, 122);
        curve1.add(45, 121);
        curve1.add(45, 120);
        curve1.add(44, 119);
        curve1.add(44, 118);
        curve1.add(43, 117);
        curve1.add(43, 116);
        curve1.add(44, 117);
        curve1.add(45, 118);
        curve1.add(45, 119);
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
        new GreyscaleImage(10, 10));
        
        int[] curve0Idx = new int[1];
        int[] curve1Idx = new int[1];
        
        double sep = contourExtractor.findClosestPair(curve0, curve1, 
            curve0Idx, curve1Idx);
        
        assertTrue(Math.abs(sep - 1) < 0.1);
        
        assertTrue(curve0Idx[0] == 11);
        
        assertTrue(curve1Idx[0] == 0);
    }
    
    public void testIsRangeConnected() throws Exception {
        
        PairIntArray edge0 = new PairIntArray();
        edge0.add(45, 133);
        edge0.add(45, 132);
        edge0.add(45, 131);
        edge0.add(45, 130);
        edge0.add(45, 129);
        edge0.add(45, 128);
        edge0.add(45, 127);
        edge0.add(44, 126);
        edge0.add(43, 125);
        edge0.add(42, 125);
        edge0.add(43, 124);
        edge0.add(44, 123); //11 closest point
        edge0.add(44, 124);
        edge0.add(44, 125);
        edge0.add(43, 126);
        
        PairIntArray edge1 = new PairIntArray();
        edge1.add(45, 119);
        edge1.add(45, 118);
        edge1.add(44, 117);
        edge1.add(43, 116);
        edge1.add(43, 117);
        edge1.add(44, 118);
        edge1.add(44, 119);
        edge1.add(45, 120);
        edge1.add(45, 121);
        edge1.add(45, 122);
        edge1.add(45, 123); //10 closest point
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(10, 10));
        
        assertTrue(contourExtractor.isRangeConnected(edge0, 0, 11));
        
        assertTrue(contourExtractor.isRangeConnected(edge1, 0, 10));
        
        edge0 = new PairIntArray();
        edge0.add(45, 133);
        edge0.add(45, 131);
        edge0.add(45, 130);
        edge0.add(45, 129);
        edge0.add(45, 128);
        
        assertFalse(contourExtractor.isRangeConnected(edge0, 0, 2));
    }
    
    public void estConnectClosestPointsIfCanTrim0() throws Exception {
        
        PairIntArray edge0 = new PairIntArray();
        for (int i = 0; i < 60; i++) {
            int y = 133 + 60 - i;
            int x = 45;
            edge0.add(x, y);
        }
        edge0.add(45, 133);
        edge0.add(45, 132);
        edge0.add(45, 131);
        edge0.add(45, 130);
        edge0.add(45, 129);
        edge0.add(45, 128);
        edge0.add(45, 127);
        edge0.add(44, 126);
        edge0.add(43, 125);
        edge0.add(42, 125);
        edge0.add(43, 124);
        edge0.add(44, 123); //<====closest point
        edge0.add(44, 124);
        edge0.add(44, 125);
        edge0.add(43, 126);
        
        PairInt[] expectedMissingPoints = new PairInt[3];
        expectedMissingPoints[0] = new PairInt(44, 124);
        expectedMissingPoints[1] = new PairInt(44, 125);
        expectedMissingPoints[2] = new PairInt(43, 126);
        
        PairIntArray edge1 = new PairIntArray();
        for (int i = 0; i < 30; i++) {
            int y = 119 - 30 + i;
            int x = 46;
            edge1.add(x, y);
        }
        edge1.add(45, 119);
        edge1.add(45, 118);
        edge1.add(44, 117);
        edge1.add(43, 116);
        edge1.add(43, 117);
        edge1.add(44, 118);
        edge1.add(44, 119);
        edge1.add(45, 120);
        edge1.add(45, 121);
        edge1.add(45, 122);
        edge1.add(45, 123); //<==== closest point
        
        int nTot = edge0.getN() + edge1.getN();
        int nExpected = nTot - expectedMissingPoints.length;
        
        PairInt edge0PointFirst = new PairInt(edge0.getX(0), edge0.getY(0));
        PairInt edge0JoinPoint = new PairInt(44, 123);
        PairInt edge1JoinPoint = new PairInt(45, 123);
        PairInt edge1PointLast = new PairInt(
            edge1.getX(edge1.getN() - 1), edge1.getY(edge1.getN() - 1));
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(200, 200));
        
        List<PairIntArray> output = 
            contourExtractor.connectClosestPointsIfCanTrim(edges);
        
        assertTrue(output.size() == 1);
        
        PairIntArray edge = output.get(0);
        
        assertTrue(edge.getN() == nExpected);
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            for (int j = 0; j < expectedMissingPoints.length; j++) {
                int xx = expectedMissingPoints[j].getX();
                int yy = expectedMissingPoints[j].getY();
                if ((x == xx) && (y == yy)) {
                    fail("points that should have been removed are present");
                }
            }
        }
        
        assertTrue((edge.getX(0) == edge0PointFirst.getX()) && 
            (edge.getY(0) == edge0PointFirst.getY()));
        
System.out.println(edge.getX(nExpected - 1) + ":" + edge1PointLast.getX() + " "
+ edge.getY(nExpected - 1) + ":" + edge1PointLast.getY());

        assertTrue((edge.getX(nExpected - 1) == edge1PointLast.getX()) && 
            (edge.getY(nExpected - 1) == edge1PointLast.getY()));

        boolean foundJP = false;
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            if (foundJP) {
                assertTrue((x == edge1JoinPoint.getX()) 
                    && (y == edge1JoinPoint.getY()));
                break;
            } else {
                if ((x == edge0JoinPoint.getX()) && (y == edge0JoinPoint.getY())) {
                    foundJP = true;
                }
            }
        }
        assertTrue(foundJP);
    }
    
    public void testConnectClosestPointsIfCanTrim1() throws Exception {
        PairIntArray edge0 = new PairIntArray();
        for (int i = 0; i < 60; i++) {
            int y = 133 + 60 - i;
            int x = 45;
            edge0.add(x, y);
        }
        edge0.add(45, 133);
        edge0.add(45, 132);
        edge0.add(45, 131);
        edge0.add(45, 130);
        edge0.add(45, 129);
        edge0.add(45, 128);
        edge0.add(45, 127);
        edge0.add(44, 126);
        edge0.add(43, 125);
        edge0.add(42, 125);
        edge0.add(43, 124);
        edge0.add(44, 123); //<====closest point
        edge0.add(44, 124);
        edge0.add(44, 125);
        edge0.add(43, 126);
        
        edge0.reverse();
        
        PairIntArray edge1 = new PairIntArray();
        for (int i = 0; i < 30; i++) {
            int y = 119 - 30 + i;
            int x = 46;
            edge1.add(x, y);
        }
        edge1.add(45, 119);
        edge1.add(45, 118);
        edge1.add(44, 117);
        edge1.add(43, 116);
        edge1.add(43, 117);
        edge1.add(44, 118);
        edge1.add(44, 119);
        edge1.add(45, 120);
        edge1.add(45, 121);
        edge1.add(45, 122);
        edge1.add(45, 123); //<==== closest point
        
        edge1.reverse();
        
        //========
        PairInt[] expectedMissingPoints = new PairInt[3];
        expectedMissingPoints[0] = new PairInt(44, 124);
        expectedMissingPoints[1] = new PairInt(44, 125);
        expectedMissingPoints[2] = new PairInt(43, 126);
        
        int nTot = edge0.getN() + edge1.getN();
        int nExpected = nTot - expectedMissingPoints.length;
        
        PairInt expectedFirstPoint = new PairInt(edge1.getX(edge1.getN() - 1), 
            edge1.getY(edge1.getN() - 1));
        PairInt expectedLastPoint = new PairInt(edge0.getX(edge0.getN() - 1), 
            edge0.getY(edge0.getN() - 1));
        PairInt edge0JoinPoint = new PairInt(44, 123);
        PairInt edge1JoinPoint = new PairInt(45, 123);
        
        //======
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(200, 200));
        
        List<PairIntArray> output = 
            contourExtractor.connectClosestPointsIfCanTrim(edges);
        
        assertTrue(output.size() == 1);
        
        PairIntArray edge = output.get(0);
        
        assertTrue(edge.getN() == nExpected);
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            for (int j = 0; j < expectedMissingPoints.length; j++) {
                int xx = expectedMissingPoints[j].getX();
                int yy = expectedMissingPoints[j].getY();
                if ((x == xx) && (y == yy)) {
                    fail("points that should have been removed are present");
                }
            }
        }
        
        assertTrue((edge.getX(0) == expectedFirstPoint.getX()) && 
            (edge.getY(0) == expectedFirstPoint.getY()));
        
        assertTrue((edge.getX(nExpected - 1) == expectedLastPoint.getX()) && 
            (edge.getY(nExpected - 1) == expectedLastPoint.getY()));

        boolean foundJP = false;
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            if (foundJP) {
                assertTrue((x == edge0JoinPoint.getX()) 
                    && (y == edge0JoinPoint.getY()));
                break;
            } else {
                if ((x == edge1JoinPoint.getX()) && (y == edge1JoinPoint.getY())) {
                    foundJP = true;
                }
            }
        }
        
        assertTrue(foundJP);
    }
    
    public void testConnectClosestPointsIfCanTrim2() throws Exception {
        
        PairIntArray edge0 = new PairIntArray();
        for (int i = 0; i < 60; i++) {
            int y = 133 + 60 - i;
            int x = 45;
            edge0.add(x, y);
        }
        edge0.add(45, 133);
        edge0.add(45, 132);
        edge0.add(45, 131);
        edge0.add(45, 130);
        edge0.add(45, 129);
        edge0.add(45, 128);
        edge0.add(45, 127);
        edge0.add(44, 126);
        edge0.add(43, 125);
        edge0.add(42, 125);
        edge0.add(43, 124);
        edge0.add(44, 123); //<====closest point
        edge0.add(44, 124);
        edge0.add(44, 125);
        edge0.add(43, 126);
        
        PairIntArray edge1 = new PairIntArray();
        for (int i = 0; i < 30; i++) {
            int y = 119 - 30 + i;
            int x = 46;
            edge1.add(x, y);
        }
        edge1.add(45, 119);
        edge1.add(45, 118);
        edge1.add(44, 117);
        edge1.add(43, 116);
        edge1.add(43, 117);
        edge1.add(44, 118);
        edge1.add(44, 119);
        edge1.add(45, 120);
        edge1.add(45, 121);
        edge1.add(45, 122);
        edge1.add(45, 123); //<==== closest point
        
        edge1.reverse();
        
        //========
        PairInt[] expectedMissingPoints = new PairInt[3];
        expectedMissingPoints[0] = new PairInt(44, 124);
        expectedMissingPoints[1] = new PairInt(44, 125);
        expectedMissingPoints[2] = new PairInt(43, 126);
        
        int nTot = edge0.getN() + edge1.getN();
        int nExpected = nTot - expectedMissingPoints.length;
        
        PairInt expectedFirstPoint = new PairInt(edge0.getX(0), 
            edge0.getY(0));
        PairInt expectedLastPoint = new PairInt(edge1.getX(edge1.getN() - 1), 
            edge1.getY(edge1.getN() - 1));
        PairInt edge0JoinPoint = new PairInt(44, 123);
        PairInt edge1JoinPoint = new PairInt(45, 123);
        
        //======
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(200, 200));
        
        List<PairIntArray> output = 
            contourExtractor.connectClosestPointsIfCanTrim(edges);
        
        assertTrue(output.size() == 1);
        
        PairIntArray edge = output.get(0);
        
        assertTrue(edge.getN() == nExpected);
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            for (int j = 0; j < expectedMissingPoints.length; j++) {
                int xx = expectedMissingPoints[j].getX();
                int yy = expectedMissingPoints[j].getY();
                if ((x == xx) && (y == yy)) {
                    fail("points that should have been removed are present");
                }
            }
        }
        
        assertTrue((edge.getX(0) == expectedFirstPoint.getX()) && 
            (edge.getY(0) == expectedFirstPoint.getY()));
        
        assertTrue((edge.getX(nExpected - 1) == expectedLastPoint.getX()) && 
            (edge.getY(nExpected - 1) == expectedLastPoint.getY()));

        boolean foundJP = false;
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            if (foundJP) {
                assertTrue((x == edge1JoinPoint.getX()) 
                    && (y == edge1JoinPoint.getY()));
                break;
            } else {
                if ((x == edge0JoinPoint.getX()) && (y == edge0JoinPoint.getY())) {
                    foundJP = true;
                }
            }
        }
        
        assertTrue(foundJP);
   
    }
    
    public void testConnectClosestPointsIfCanTrim3() throws Exception {
        
        PairIntArray edge0 = new PairIntArray();
        for (int i = 0; i < 60; i++) {
            int y = 133 + 60 - i;
            int x = 45;
            edge0.add(x, y);
        }
        edge0.add(45, 133);
        edge0.add(45, 132);
        edge0.add(45, 131);
        edge0.add(45, 130);
        edge0.add(45, 129);
        edge0.add(45, 128);
        edge0.add(45, 127);
        edge0.add(44, 126);
        edge0.add(43, 125);
        edge0.add(42, 125);
        edge0.add(43, 124);
        edge0.add(44, 123); //<====closest point
        edge0.add(44, 124);
        edge0.add(44, 125);
        edge0.add(43, 126);
        
        edge0.reverse();
        
        PairIntArray edge1 = new PairIntArray();
        for (int i = 0; i < 30; i++) {
            int y = 119 - 30 + i;
            int x = 46;
            edge1.add(x, y);
        }
        edge1.add(45, 119);
        edge1.add(45, 118);
        edge1.add(44, 117);
        edge1.add(43, 116);
        edge1.add(43, 117);
        edge1.add(44, 118);
        edge1.add(44, 119);
        edge1.add(45, 120);
        edge1.add(45, 121);
        edge1.add(45, 122);
        edge1.add(45, 123); //<==== closest point
                
        //========
        PairInt[] expectedMissingPoints = new PairInt[3];
        expectedMissingPoints[0] = new PairInt(44, 124);
        expectedMissingPoints[1] = new PairInt(44, 125);
        expectedMissingPoints[2] = new PairInt(43, 126);
        
        int nTot = edge0.getN() + edge1.getN();
        int nExpected = nTot - expectedMissingPoints.length;
        
        PairInt expectedFirstPoint = new PairInt(edge1.getX(0), 
            edge1.getY(0));
        PairInt expectedLastPoint = new PairInt(edge0.getX(edge0.getN() - 1), 
            edge0.getY(edge0.getN() - 1));
        PairInt edge0JoinPoint = new PairInt(44, 123);
        PairInt edge1JoinPoint = new PairInt(45, 123);
        
        //======
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(edge0.copy());
        edges.add(edge1.copy());
        
        EdgeExtractor contourExtractor = new EdgeExtractor(
            new GreyscaleImage(200, 200));
        
        List<PairIntArray> output = 
            contourExtractor.connectClosestPointsIfCanTrim(edges);
        
        assertTrue(output.size() == 1);
        
        PairIntArray edge = output.get(0);
        
        assertTrue(edge.getN() == nExpected);
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            for (int j = 0; j < expectedMissingPoints.length; j++) {
                int xx = expectedMissingPoints[j].getX();
                int yy = expectedMissingPoints[j].getY();
                if ((x == xx) && (y == yy)) {
                    fail("points that should have been removed are present");
                }
            }
        }
        
        assertTrue((edge.getX(0) == expectedFirstPoint.getX()) && 
            (edge.getY(0) == expectedFirstPoint.getY()));
        
        assertTrue((edge.getX(nExpected - 1) == expectedLastPoint.getX()) && 
            (edge.getY(nExpected - 1) == expectedLastPoint.getY()));

        boolean foundJP = false;
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            if (foundJP) {
                assertTrue((x == edge0JoinPoint.getX()) 
                    && (y == edge0JoinPoint.getY()));
                break;
            } else {
                if ((x == edge1JoinPoint.getX()) && (y == edge1JoinPoint.getY())) {
                    foundJP = true;
                }
            }
        }
        
        assertTrue(foundJP);
   
    }
    
    public void testFindEdges() throws Exception {
                
        String filePath = ResourceFinder.findFileInTestResources("africa.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        // get a line thinned image:
        CannyEdgeFilter edgeFilter = new CannyEdgeFilter();
        edgeFilter.applyFilter(img);
        
        EdgeExtractor contourExtractor = new EdgeExtractor(img);
        List<PairIntArray> edges = contourExtractor.findEdges();
        
        int clr = 0;
        Image img2 = new Image(img.getWidth(), img.getHeight());
        for (int i = 0; i < edges.size(); i++) {
            PairIntArray edge = edges.get(i);
            if (clr > 5) {
                clr = 0;
            }
            int c = Color.BLUE.getRGB();
            switch(clr) {
                case 1:
                    c = Color.PINK.getRGB();
                    break;
                case 2:
                    c = Color.GREEN.getRGB();
                    break;
                case 3:
                    c = Color.RED.getRGB();
                    break;
                case 4:
                    c = Color.CYAN.getRGB();
                    break;
                case 5:
                    c = Color.MAGENTA.getRGB();
                    break;
            }
            for (int ii = 0; ii < edge.getN(); ii++) {
                img2.setRGB(edge.getX(ii), edge.getY(ii), c);
            }
            clr++;
        }
        
        ImageDisplayer.displayImage("edge detected image", img2);
        
        assertTrue(!edges.isEmpty());
    }

    @Test
    public void testExtractEdges() {
        
        GreyscaleImage img = new GreyscaleImage(10, 10);
                        
        PairIntArray t = getSquare();
        for (int j = 0; j < t.getN(); j++) {                
            t.set(j, t.getX(j), t.getY(j));
            img.setValue(t.getX(j), t.getY(j), 250);                
        }        
        
        GreyscaleImage guide = img.copyImage();
        EdgeExtractor extractor = null;
        
        for (int k = 0; k < 2; k++) {
            
            if (k == 0) {
                extractor = new EdgeExtractor(img);
            } else {
                extractor = new EdgeExtractor(img, guide);
            }

            List<PairIntArray> result = extractor.findEdges();

            assertNotNull(result);
            assertTrue(result.size() == 1);

            // assert that they are all connected
            double sep = Math.sqrt(2);
            PairIntArray edge = result.get(0);
            for (int i = 1; i < (edge.getN() - 1); i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);

                double distPrev = Math.pow(edge.getX(i - 1) - x, 2) + 
                    Math.pow(edge.getY(i - 1) - y, 2);
                distPrev = Math.sqrt(distPrev);

                double distNext = Math.pow(edge.getX(i + 1) - x, 2) + 
                    Math.pow(edge.getY(i + 1) - y, 2);
                distNext = Math.sqrt(distNext);

                assertTrue(distPrev <= (sep + 0.01));
                assertTrue(distNext <= (sep + 0.01));
            }
        
            // assert that main lines are present though the corners may have
            //  been eroded
            /*
              6
              5   @ @ @ @ @ @ @ @ 
              4   @             @
              3   @             @  
              2   @             @  
              1   @ @ @ @ @ @ @ @
              0
                0 1 2 3 4 5 6 7 8 9
            */

            assertTrue(edge.getN() >= 18);

            // print result
            GreyscaleImage extracted = new GreyscaleImage(10, 10);
            for (int i = 0; i < edge.getN(); i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                extracted.setValue(x, y, 1);
            }

            /*
            for (int row = (extracted.getHeight() - 1); row > -1; row--) {
                StringBuilder sb = new StringBuilder();
                for (int col = 0; col < extracted.getWidth(); col++) {
                    int v = extracted.getValue(col, row);
                    if (v == 0) {
                        sb.append(" ");
                    } else {
                        sb.append(v);
                    }
                    sb.append(" ");
                }
                System.out.println(sb.toString());
            }
            System.out.println("\n");
            */

            // make a pair int array out of expected and remove the corners as found
            PairIntArray expected = getSquare();
            // remove the corners which may have been eroded
            for (int i = (expected.getN() - 1); i > -1; i--) {
                if ((expected.getX(i) == 1) && (expected.getY(i) == 1)) {
                    expected.removeRange(i, i);
                } else if ((expected.getX(i) == 1) && (expected.getY(i) == 5)) {
                    expected.removeRange(i, i);
                } else if ((expected.getX(i) == 8) && (expected.getY(i) == 1)) {
                    expected.removeRange(i, i);
                } else if ((expected.getX(i) == 8) && (expected.getY(i) == 5)) {
                    expected.removeRange(i, i);
                }
            }

            for (int i = 0; i < edge.getN(); i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                for (int j = 0; j < expected.getN(); j++) {
                    if ((expected.getX(j) == x) && (expected.getY(j) == y)) {
                        expected.removeRange(j, j);
                        break;
                    }
                }
            }

            assertTrue(expected.getN() == 0);
        }
    }
    
    @Test
    public void testRemoveEdgesShorterThan() {
        
        GreyscaleImage img = new GreyscaleImage(100, 100);
        
        List<PairIntArray> allEdges = new ArrayList<PairIntArray>();
        
        int minNumberOfPixelsInEdge = 10;
        PairIntArray t = getSquare();
        for (int j = 0; j < t.getN(); j++) {                
            t.set(j, t.getX(j) + 50, t.getY(j));
            img.setValue(t.getX(j), t.getY(j), 250);                
        }
        allEdges.add(t);
        
        t = new PairIntArray();
        for (int i = 0; i < (minNumberOfPixelsInEdge - 1); i++) {
            t.add(10, 10 + i);
            img.setValue(t.getX(i), t.getY(i), 250);
        }
        allEdges.add(t);
        
        t = getSquare();
        for (int j = 0; j < t.getN(); j++) {
            t.set(j, t.getX(j), t.getY(j) + 50);
            img.setValue(t.getX(j), t.getY(j), 250); 
        }
        allEdges.add(t);
        
        EdgeExtractor extractor = new EdgeExtractor(img);
        
        List<PairIntArray> result = extractor.findEdges();
        
        assertNotNull(result);
        assertTrue(result.size() == 2);
        
        // assert that they are all connected
        
        assertTrue(allEdges.size() == 3);
        extractor.removeEdgesShorterThan(allEdges, minNumberOfPixelsInEdge);
        
        assertNotNull(allEdges);
        assertTrue(allEdges.size() == 2);        
    }

    /**
      7    
      6
      5   @ @ @ @ @ @ @ @ 
      4   @             @
      3   @             @  
      2   @             @  
      1   @ @ @ @ @ @ @ @
      0
        0 1 2 3 4 5 6 7 8 9
    */
    private PairIntArray getSquare() {
        
        PairIntArray xy = new PairIntArray();
        for (int x = 1; x < 9; x++) {
            xy.add(x, 5);
        }
        for (int y = 4; y > 0; y--) {
            xy.add(8, y);
        }
        for (int x = 7; x > 0; x--) {
            xy.add(x, 1);
        }
        for (int y = 2; y < 5; y++) {
            xy.add(1, y);
        }
        
        return xy;
    }
    
}
