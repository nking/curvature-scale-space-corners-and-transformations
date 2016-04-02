package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class EdgeExtractorSimpleTest extends TestCase {
    
    public EdgeExtractorSimpleTest() {
    }
    
    public void test0() {
        
        /*
        7  #                #  
        6    #  #        #     #
        5          @  #        #
        4          #
        3          #
        2          #
        1          #
        0
          0  1  2  3  4  5  6  7
        */
        
        int[][] edgeImage = new int[8][];
        for (int i = 0; i < edgeImage.length; ++i) {
            edgeImage[i] = new int[8];
        }
        
        // [row][col]
        edgeImage[1][3] = 1;
        edgeImage[2][3] = 1;
        edgeImage[3][3] = 1;
        edgeImage[4][3] = 1;
        edgeImage[5][3] = 1;
        
        edgeImage[6][2] = 1;
        edgeImage[6][1] = 1;
        edgeImage[7][0] = 1;
        
        edgeImage[5][4] = 1;
        edgeImage[6][5] = 1;
        edgeImage[7][6] = 1;
        edgeImage[6][7] = 1;
        edgeImage[5][7] = 1;
        
        EdgeExtractorSimple edgeExtractor = new EdgeExtractorSimple(edgeImage);
        edgeExtractor.extractEdges();
        
        Set<PairInt> junctions = edgeExtractor.getJunctions();
        
        List<PairIntArray> edges = edgeExtractor.getEdges();
                
        Set<PairInt> expectedJunctions = new HashSet<PairInt>();
        expectedJunctions.add(new PairInt(4, 3));
        expectedJunctions.add(new PairInt(5, 3));
        expectedJunctions.add(new PairInt(5, 4));
        
        assertEquals(expectedJunctions.size(), junctions.size());
        Set<PairInt> rm = new HashSet<PairInt>();
        for (PairInt p : expectedJunctions) {
            assertTrue(expectedJunctions.contains(p));
            rm.add(p);
        }
        assertTrue(expectedJunctions.removeAll(rm));
        assertTrue(expectedJunctions.isEmpty());
        
        assertEquals(3, edges.size());
        
        Set<PairInt> line1 = new HashSet<PairInt>();
        line1.add(new PairInt(1, 3));
        line1.add(new PairInt(2, 3));
        line1.add(new PairInt(3, 3));
        
        Set<PairInt> line2 = new HashSet<PairInt>();
        line2.add(new PairInt(6, 2));
        line2.add(new PairInt(6, 1));
        line2.add(new PairInt(7, 0));
        
        Set<PairInt> line3 = new HashSet<PairInt>();
        line3.add(new PairInt(6, 5));
        line3.add(new PairInt(7, 6));
        line3.add(new PairInt(6, 7));
        line3.add(new PairInt(5, 7));
        
        //TODO: add a test for order of points in the extracted edges
        
        for (int i = 0; i < edges.size(); ++i) {
            
            PairIntArray edge = edges.get(i);
            
            assertTrue(edge.getN() > 0);
                        
            int currentList = 1;
            
            for (int idx = 0; idx < edge.getN(); ++idx) {
                
                PairInt p0 = new PairInt(edge.getX(idx), edge.getY(idx));
                
                int lineList = -1;
                if (line1.contains(p0)) {
                    lineList = 1;
                    line1.remove(p0);
                } else if (line2.contains(p0)) {
                    lineList = 2;
                    line2.remove(p0);
                } else if (line3.contains(p0)) {
                    lineList = 3;
                    line3.remove(p0);
                } else {
                    fail("could not find point " + p0 + " in an epected list");
                }
                
                if (idx == 0) {
                    currentList = lineList;
                } else {
                    if (currentList != lineList) {
                        fail("extracted line list doesn't contain expected");
                    }
                }
            }
            
        }
        assertTrue(line1.isEmpty());
        assertTrue(line2.isEmpty());
        assertTrue(line3.isEmpty());
    }
}
