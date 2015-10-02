package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
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
    
    public List<PairIntArray> getJunctionEdges0() {
    
        /*
                        0      4   / 3
                     |-----||-----|
                                   \
                                  1 \     2
                                     \|-------
        */
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
                
        // add edge 0.   start=(15,10) end=(24,10)
        edges.add(new PairIntArray());
        for (int x = 15; x < 25; x++) {
            edges.get(edges.size() - 1).add(x, 10);
        }
        // add edge 1.  start=(35,11) end=(44,20)
        edges.add(new PairIntArray());
        for (int x = 35; x < 45; x++) {
            edges.get(edges.size() - 1).add(x, 11 + (x - 35));
        }
        // add edge 2.  start=(45,20) end=(54,20)
        edges.add(new PairIntArray());
        for (int x = 45; x < 55; x++) {
            edges.get(edges.size() - 1).add(x, 11 + (44 - 35)); 
        }
        // add edge 3.  start=(35,9) end=(37,7)
        edges.add(new PairIntArray());
        for (int x = 35; x < 38; x++) {
            edges.get(edges.size() - 1).add(x, 9 - (x - 35)); 
        }
        // add edge 4.   start=(25,10) end=(34,10)
        edges.add(new PairIntArray());
        for (int x = 25; x < 35; x++) {
            edges.get(edges.size() - 1).add(x, 10);
        }
        
        return edges;
    }
    
    public List<PairIntArray> getJunctionEdges1() {
    
        /*
                *\
             0 /  \ 1      result will be single closed edge 0:2:1
              /    \
             |------|
                 2
        */
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
                
        // add edge 0
        edges.add(new PairIntArray());
        // start = (25, 5), end = (23, 7)
        int y = 4;
        for (int x = 25; x >= 23; x--) {
            y++;
            edges.get(edges.size() - 1).add(x, y);
        }
        
        // add edge 1
        edges.add(new PairIntArray());
        // start = (28,7) end = (26, 5)
        y = 8;
        for (int x = 28; x >= 26; x--) {
            y--;
            edges.get(edges.size() - 1).add(x, y);
        }
        
        // add edge 2
        edges.add(new PairIntArray());
        y = 7;
        // start=(22,7)     end=(27,7)
        for (int x = 22; x < 28; x++) {
            edges.get(edges.size() - 1).add(x, y); 
        }
        
        return edges;
    }
    
    private void assertJunctionEdges0EndPoints(Map<PairInt, Integer> 
        endPointMap) {
        
        assertTrue(endPointMap.size() == 10);
        
        Integer index;
        
        // edge 0
        index = endPointMap.get(new PairInt(15, 10));
        assertNotNull(index);
        assertTrue(index.intValue() == 0);
        index = endPointMap.get(new PairInt(24, 10));
        assertNotNull(index);
        assertTrue(index.intValue() == 0);

        // edge 4
        index = endPointMap.get(new PairInt(25, 10));
        assertNotNull(index);
        assertTrue(index.intValue() == 4);
        index = endPointMap.get(new PairInt(34, 10));
        assertTrue(index.intValue() == 4);
        
        // edge 1
        index = endPointMap.get(new PairInt(35, 11));
        assertNotNull(index);
        assertTrue(index.intValue() == 1);
        index = endPointMap.get(new PairInt(44, 20));
        assertNotNull(index);
        assertTrue(index.intValue() == 1);
        
        // edge 2
        index = endPointMap.get(new PairInt(45, 20));
        assertNotNull(index);
        assertTrue(index.intValue() == 2);
        index = endPointMap.get(new PairInt(54, 20)); 
        assertNotNull(index);
        assertTrue(index.intValue() == 2);
        
        // edge 3
        index = endPointMap.get(new PairInt(35, 9));
        assertNotNull(index);
        assertTrue(index.intValue() == 3);
        index = endPointMap.get(new PairInt(37, 7));
        assertNotNull(index);
        assertTrue(index.intValue() == 3);
    }
    
    private void assertJunctionEdges0Output(List<PairIntArray> output) {
        
        assertTrue(output.size() == 2);
        
        List<PairIntArray> originalEdges = getJunctionEdges0();
        
        // assert output[0] is edge 0 : 4 : 1 : 2
        PairIntArray expectedOutput0 = new PairIntArray();
        expectedOutput0.addAll(originalEdges.get(0));
        expectedOutput0.addAll(originalEdges.get(4));
        expectedOutput0.addAll(originalEdges.get(1));
        expectedOutput0.addAll(originalEdges.get(2));
        
        assertTrue(expectedOutput0.getN() == output.get(0).getN());
        
        for (int i = 0; i < expectedOutput0.getN(); i++) {
            int eX = expectedOutput0.getX(i);
            int eY = expectedOutput0.getY(i);
            int x = output.get(0).getX(i);
            int y = output.get(0).getY(i);
            assertTrue(eX == x);
            assertTrue(eY == y);
        }
        
        // assert output[1] is edge 3
        PairIntArray expectedOutput1 = new PairIntArray();
        expectedOutput1.addAll(originalEdges.get(3));
        
        assertTrue(expectedOutput1.getN() == output.get(1).getN());
   
        // the order isn't specified, so assert that one ordering doesn't fail
        boolean foundAll = true;
        for (int i = 0; i < expectedOutput1.getN(); i++) {
            int eX = expectedOutput1.getX(i);
            int eY = expectedOutput1.getY(i);
            int x = output.get(1).getX(i);
            int y = output.get(1).getY(i);
            if (eX == x) {
                foundAll = false;
            }
            if (eY == y) {
                foundAll = false;
            }
        }
        if (!foundAll) {
            // reverse the pattern
            expectedOutput1.reverse();
            foundAll = true;
            for (int i = 0; i < expectedOutput1.getN(); i++) {
                int eX = expectedOutput1.getX(i);
                int eY = expectedOutput1.getY(i);
                int x = output.get(1).getX(i);
                int y = output.get(1).getY(i);
                if (eX == x) {
                    foundAll = false;
                }
                if (eY == y) {
                    foundAll = false;
                }
            }
        }
    }
    
    private void assertJunctionEdges1Output(List<PairIntArray> output) {
        
        assertTrue(output.size() == 1);
        
        List<PairIntArray> originalEdges = getJunctionEdges1();
        
        PairIntArray edge = output.get(0);
        
        // assert edge 0 : 2 : 1
        
        int count = 0;
        
        // assert edge 0
        for (int i = 0; i < originalEdges.get(0).getN(); i++) {
            int x = originalEdges.get(0).getX(i);
            int y = originalEdges.get(0).getY(i);
            int xt = edge.getX(count);
            int yt = edge.getY(count);
            assertTrue(xt == x);
            assertTrue(yt == y);
            count++;
        }
        
        //assert edge 2 is next
        for (int i = 0; i < originalEdges.get(2).getN(); i++) {
            int x = originalEdges.get(2).getX(i);
            int y = originalEdges.get(2).getY(i);
            int xt = edge.getX(count);
            int yt = edge.getY(count);
            assertTrue(xt == x);
            assertTrue(yt == y);
            count++;
        }
        
        // assert edge 1 is next
        for (int i = 0; i < originalEdges.get(1).getN(); i++) {
            int x = originalEdges.get(1).getX(i);
            int y = originalEdges.get(1).getY(i);
            int xt = edge.getX(count);
            int yt = edge.getY(count);
            assertTrue(xt == x);
            assertTrue(yt == y);
            count++;
        }
    }
    
    public void testGetOppositeEndPointOfEdge() throws Exception {
        
        List<PairIntArray> edges = getJunctionEdges0();
      
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));

        for (int i = 0; i < edges.size(); i++) {
            PairIntArray edge = edges.get(i);
            PairInt startPoint = new PairInt(edge.getX(0), edge.getY(0));
            int n = edge.getN();
            PairInt expectedEndPoint = new PairInt(edge.getX(n - 1), edge.getY(n - 1));
            
            PairInt endPoint = extractor.getOppositeEndPointOfEdge(startPoint, edge);
            
            assertTrue(expectedEndPoint.equals(endPoint));
        }
    }
    
    public void testFindStartingPoint() throws Exception {
        
        List<PairIntArray> edges = getJunctionEdges0();
        
        /*
                        0      4   / 3
                     |-----||-----|
                                   \
                                  1 \     2
                                     \|-------
        */
        
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        
        extractor.overrideEdgeSizeLowerLimit(1);
        
        Map<PairInt, Integer> endPointMap = extractor.createEndPointMap(
            edges);
        
        PairInt startPoint = new PairInt(edges.get(2).getX(0),
            edges.get(2).getY(0));
        
        PairInt earliestStartPoint = extractor.findStartingPoint(startPoint, 
            endPointMap, edges);
        
        assertTrue(earliestStartPoint.getX() == edges.get(0).getX(0));
        
        assertTrue(earliestStartPoint.getY() == edges.get(0).getY(0));
        
    }
    
    public void testFindStartingPoint1() throws Exception {
        
        List<PairIntArray> edges = getJunctionEdges1();
        
         /*
                *\
             0 /  \ 1      result will be single closed edge 0:2:1
              /    \
             |------|
                 2
        Edge 0: start = (25, 5), end = (23, 7)
        */
        
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        
        extractor.overrideEdgeSizeLowerLimit(1);
        
        Map<PairInt, Integer> endPointMap = extractor.createEndPointMap(
            edges);
        
        PairInt startPoint = new PairInt(edges.get(0).getX(0),
            edges.get(0).getY(0));
        
        PairInt earliestStartPoint = extractor.findStartingPoint(startPoint, 
            endPointMap, edges);
        
        assertTrue(earliestStartPoint.getX() == edges.get(0).getX(0));
        
        assertTrue(earliestStartPoint.getY() == edges.get(0).getY(0));
        
    }
    
    public void testCreateEndPointMap() throws Exception {
        
        List<PairIntArray> edges = getJunctionEdges0();
        
        /*
                        0      4   / 3
                     |-----||-----|
                                   \
                                  1 \     2
                                     \|-------
        */
        
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        
        extractor.overrideEdgeSizeLowerLimit(1);
        
        Map<PairInt, Integer> endPointMap = extractor.createEndPointMap(
            edges);
        
        assertJunctionEdges0EndPoints(endPointMap);
        
        // ===== reverse the 2nd item in edges ======
        edges = getJunctionEdges0();
        edges.get(1).reverse();
        
        endPointMap = extractor.createEndPointMap(edges);
        
        assertJunctionEdges0EndPoints(endPointMap);
        
        // ===== reverse the 3rd and 4th item in edges ======
        edges = getJunctionEdges0();
        edges.get(2).reverse();
        edges.get(3).reverse();
        
        endPointMap = extractor.createEndPointMap(edges);
        
        assertJunctionEdges0EndPoints(endPointMap);
        
        // ===== reverse edges 4 and 1 ======
        edges = getJunctionEdges0();
        edges.get(4).reverse();
        edges.get(1).reverse();
        
        endPointMap = extractor.createEndPointMap(edges);
        
        assertJunctionEdges0EndPoints(endPointMap);
        
    }
    
    public void testAppendToOutput() throws Exception {
        
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        
        List<PairIntArray> edges = getJunctionEdges1();
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        output.add(edges.get(0));
        
        // add expected rest of edges in order and assert them
        extractor.appendToOutput(output, edges.get(2));
        extractor.appendToOutput(output, edges.get(1));
        assertJunctionEdges1Output(output);
        
        // ===repeat test, but reverse edges 2 and 1 ====
        edges = getJunctionEdges1();
        output = new ArrayList<PairIntArray>();
        output.add(edges.get(0));
        edges.get(2).reverse();
        edges.get(1).reverse();
        extractor.appendToOutput(output, edges.get(2));
        extractor.appendToOutput(output, edges.get(1));
        assertJunctionEdges1Output(output);
    }
    
    public void testMergeAdjacentEndPoints2_1() throws Exception {
        
        List<PairIntArray> edges = getJunctionEdges1();
        
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        extractor.overrideEdgeSizeLowerLimit(1);
        
        List<PairIntArray> output = extractor.mergeAdjacentEndPoints(edges);
        
        assertJunctionEdges1Output(output);
        
        // ====== reverse edge 1  and test again =======        
        extractor = new EdgeExtractor(new GreyscaleImage(100, 100));
        extractor.overrideEdgeSizeLowerLimit(1);
        edges = getJunctionEdges1();
        edges.get(1).reverse();
        
        output = extractor.mergeAdjacentEndPoints(edges);
        
        assertJunctionEdges1Output(output);
       
        // ====== reverse edge 2  and test again =======
        extractor = new EdgeExtractor(new GreyscaleImage(100, 100));
        extractor.overrideEdgeSizeLowerLimit(1);
        edges = getJunctionEdges1();
        edges.get(2).reverse();
        
        output = extractor.mergeAdjacentEndPoints(edges);
        
        assertJunctionEdges1Output(output);
        
    }
    
    public void testMergeAdjacentEndPoints2_0() throws Exception {
                
        for (int nTest = 0; nTest < 4; nTest++) {
        
            List<PairIntArray> edges = getJunctionEdges0();
            
            switch(nTest) {
                case 1:
                    edges.get(4).reverse();
                    break;
                case 2:
                    edges.get(3).reverse();
                    edges.get(1).reverse();
                    break;
                case 3:
                    edges.get(2).reverse();
                    break;
                default:
                    break;
            }
        
            EdgeExtractor extractor = new EdgeExtractor(
                new GreyscaleImage(100, 100));
            extractor.overrideEdgeSizeLowerLimit(1);

            List<PairIntArray> output = extractor.mergeAdjacentEndPoints(edges);
            
            //System.out.println("nTest=" + nTest);
            assertJunctionEdges0Output(output);

        }
    }

    public void testMergeAdjacentEndPoints2_2() throws Exception {
               
        EdgeExtractor extractor = new EdgeExtractor(
            new GreyscaleImage(100, 100));
        extractor.overrideEdgeSizeLowerLimit(1);

        List<PairIntArray> output = extractor.mergeAdjacentEndPoints(
            new ArrayList<PairIntArray>());
        
        assertTrue(output.isEmpty());
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
       
        //TODO: redo tests for this method
        
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
                
        String fileName = "africa2.png";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        log.info("fileName=" + fileName);
        
        {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/test.png", img);
        }
        
        // get a line thinned image:
        CannyEdgeFilter edgeFilter = new CannyEdgeFilter();
        edgeFilter.applyFilter(img);
        
        IEdgeExtractor contourExtractor = new EdgeExtractorWithJunctions(img);
        List<PairIntArray> edges = contourExtractor.findEdges();
        
        log.info("edges.size()=" + edges.size());
        
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
        
        String dirPath = ResourceFinder.findDirectory("bin");
        String sep = System.getProperty("file.separator");
        ImageIOHelper.writeOutputImage(dirPath + sep + "africa_edges.png", img2);
        
        img.multiply(250);
        ImageIOHelper.writeOutputImage(dirPath + sep + "africa_thinned.png", 
            img);
        
        assertTrue(!edges.isEmpty());
    }
    
    public void testFindEdges2() throws Exception {
                
        //String fileName = "house.gif";
        //String fileName = "lab.gif";
        //String fileName = "susan-in.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "africa2.png";
        //String fileName = "valve_gaussian.png";
        //String fileName = "lena.jpg";
        String fileName = "middlebury_cones_im2.png";
        //String fileName = "venturi_mountain_j6_0001.png";
        
        log.info("fileName=" + fileName);
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        String dirPath = ResourceFinder.findDirectory("bin");
        String sep = System.getProperty("file.separator");
        
        int idx = fileName.lastIndexOf(".");
        String fileNameRoot = fileName.substring(0, idx);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        // to experiment w/ a color mapping instead of intensity:
        if (fileName.contains("cones") || fileName.contains("venturi")) {//cannot use on B&W images
            ImageExt clrImg = ImageIOHelper.readImageExt(filePath);
            //img = imageProcessor.createGreyscaleFromColorSegmentation(clrImg);
            img = imageSegmentation.applyUsingCIEXYPolarThetaThenKMPPThenHistEq(clrImg, 4, false);
            
            ImageIOHelper.writeOutputImage(dirPath + sep + fileNameRoot 
                + "_color_theta.png", img);
        }
        
        // get a line thinned image:
        CannyEdgeFilter edgeFilter = new CannyEdgeFilter();
        edgeFilter.applyFilter(img);
        
        IEdgeExtractor contourExtractor = new EdgeExtractorWithJunctions(img);
        List<PairIntArray> edges = contourExtractor.findEdges();
        
        log.info("edges.size()=" + edges.size());
        
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
        
        ImageIOHelper.writeOutputImage(dirPath + sep + fileNameRoot + 
            "_edges.png", img2);
        
        img.multiply(250);
        ImageIOHelper.writeOutputImage(dirPath + sep + fileNameRoot 
            + "_thinned.png", 
            img);
        
        //assertFalse(edges.isEmpty());
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
    
    public static void main(String[] args) {
        
        try {
            
            EdgeExtractorTest test = new EdgeExtractorTest();
            
            //test.testFindEdges2();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
