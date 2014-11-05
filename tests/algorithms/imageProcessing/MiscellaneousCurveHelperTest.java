package algorithms.imageProcessing;

import Jama.Matrix;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class MiscellaneousCurveHelperTest extends TestCase {
    
    public MiscellaneousCurveHelperTest(String testName) {
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
    
    public void testCurveIsOrderedClockwise() throws Exception {
        
        PairIntArray cwCurve = new PairIntArray();
        
        int xc = 5;
        int yc = 5;
        
        cwCurve.add(xc + 4, yc);
        cwCurve.add(xc + 3, yc - 1);
        cwCurve.add(xc, yc - 2);
        cwCurve.add(xc - 3, yc - 2);
        cwCurve.add(xc - 2, yc - 1);
        cwCurve.add(xc - 3, yc);
        cwCurve.add(xc - 1, yc + 3);
        cwCurve.add(xc, yc + 4);
        cwCurve.add(xc + 2, yc + 2);
        
        MiscellaneousCurveHelper helper = new MiscellaneousCurveHelper();
        boolean isCW = helper.curveIsOrderedClockwise(cwCurve);
        
        assertTrue(isCW);
        
        cwCurve.reverse();
        
        isCW = helper.curveIsOrderedClockwise(cwCurve);
        
        assertFalse(isCW);
    }
   
    public void testCrossCorrelation() {
                
        PairIntArray curve0, curve1;
        int[] crossCorrelationOffset = new int[1];
        boolean isAdjacent;
        
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        
        //y = x^2
        for (int x = -4; x < 5; x++) {
            int y = x * x;
            curve0.add(x, y);
        }
        
        // create curve1, shorter and shifted in x by 1
        int expectedOffset = 2;
        for (int i = 0; i < 5; i++) {
            int idx = i + expectedOffset;
            curve1.add(curve0.getX(idx) + 1, curve0.getY(idx));
        }
        
        /*
        [junit]  -4 -3 -2 -1 0 1 2 3 4
        [junit]        -1  0 1 2 3
        */
                        
        MiscellaneousCurveHelper instance = new MiscellaneousCurveHelper();
        
        isAdjacent = instance.crossCorrelation(curve0, curve1, 
            crossCorrelationOffset);
        
        assertTrue(isAdjacent);
        
        assertTrue(crossCorrelationOffset[0] == expectedOffset);
        
        
        // ====== test case for curve1 is shorter than curve0 but starts
        // at a negative offset
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        //y = x^2
        for (int x = -4; x < 5; x++) {
            int y = x * x;
            curve0.add(x, y);
        }
        
        for (int x = -7; x < 3; x++) {
            int y = x * x;
            curve1.add(x - 1, y);
        }
        
        expectedOffset = 3;
        
        instance = new MiscellaneousCurveHelper();
        
        isAdjacent = instance.crossCorrelation(curve0, curve1, 
            crossCorrelationOffset);
        
        assertTrue(isAdjacent);
        //1
        assertTrue(crossCorrelationOffset[0] == expectedOffset);

        
        // ====== test case for curve0 and curve1 not overlapping
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        //y = x^2
        for (int x = -4; x < 5; x++) {
            int y = x * x;
            curve0.add(x, y);
        }
        for (int x = 20; x < 30; x++) {
            int y = x * x;
            curve1.add(x, y);
        }
        
        instance = new MiscellaneousCurveHelper();
        
        isAdjacent = instance.crossCorrelation(curve0, curve1, 
            crossCorrelationOffset);
        
        assertFalse(isAdjacent);
        
        
        // ====== test case for curve0 being shorter than curve1 and
        // having an offset w.r.t. curve1 that is such that curve0 end points
        // are not matched (further in offset space than curve1 extends to).
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        //y = x^2
        for (int x = -4; x < 5; x++) {
            int y = x * x;
            curve0.add(x, y);
        }
        for (int x = 3; x < 8; x++) {
            int y = x * x;
            curve1.add(x + 1, y);
        }

        //                          3
        //                          4  5  6  7  8
        // -4  -3  -2  -1  0  1  2  3
        
        expectedOffset = 7;
        
        instance = new MiscellaneousCurveHelper();
        
        isAdjacent = instance.crossCorrelation(curve0, curve1, 
            crossCorrelationOffset);
        
        assertTrue(isAdjacent);
        
        assertTrue(crossCorrelationOffset[0] == expectedOffset);

    }
   
    public void testProcessPair() throws Exception {
        
        /*
        given 2 edges, return true if they overlap. If they overlap
        curve0 is given the larger curve and any outlying check points.
        protected boolean processOverlappingPair(PairIntArrayWithColor curve0, 
            PairIntArrayWithColor curve1) {
        */
        
        PairIntArray curve0, curve1;
        
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        
        //y = x^2
        for (int x = -4; x < 5; x++) {
            int y = x * x;
            curve0.add(x, y);
        }
        
        // create curve1, shorter and shifted in x by 1
        int expectedOffset = 2;
        for (int i = 0; i < 5; i++) {
            int idx = i + expectedOffset;
            curve1.add(curve0.getX(idx) + 1, curve0.getY(idx));
        }
        
        PairIntArrayWithColor c0 = new PairIntArrayWithColor(curve0);
        PairIntArrayWithColor c1 = new PairIntArrayWithColor(curve1);
        
        /*
        [junit]  -4 -3 -2 -1 0 1 2 3 4
        [junit]        -1  0 1 2 3
        */
                        
        MiscellaneousCurveHelper instance = new MiscellaneousCurveHelper();
       
        boolean overlapped = instance.processOverlappingPair(c0, c1);
        
        assertTrue(overlapped);
        
        /*
        3                         3
        2          @  @  @        2               @
        1    @  @  @  @     ==>   1   @  @  @  @
        0                         0
          0  1  2  3  4  5         0  1  2  3  4  5 
        */
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        curve0.add(1, 1);
        curve0.add(2, 1);
        curve0.add(3, 1);
        curve0.add(4, 1);
        curve1.add(3, 2);
        curve1.add(4, 2);
        curve1.add(5, 2);
        c0 = new PairIntArrayWithColor(curve0);
        c1 = new PairIntArrayWithColor(curve1);
        overlapped = instance.processOverlappingPair(c0, c1);
        
        assertTrue(overlapped);
        assertTrue(c0.getN() == 5);
        
        /*
        6          @
        5       @
        4    @
        3                         
        2                  
        1    @  @  @  @    
        0                  
          0  1  2  3  4  5 
        */
        
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        curve0.add(1, 1);
        curve0.add(2, 1);
        curve0.add(3, 1);
        curve0.add(4, 1);
        curve1.add(1, 4);
        curve1.add(2, 5);
        curve1.add(3, 6);
        c0 = new PairIntArrayWithColor(curve0);
        c1 = new PairIntArrayWithColor(curve1);
        overlapped = instance.processOverlappingPair(c0, c1);
        
        assertFalse(overlapped);
    }
    
    public void testPruneIncludedAdjacentCurves() throws Exception {
        
        /*
        3                         3
        2          @  @  @        2               @
        1    @  @  @  @     ==>   1   @  @  @  @
        0                         0
          0  1  2  3  4  5         0  1  2  3  4  5 
        */
        PairIntArray curve0 = new PairIntArray();
        PairIntArray curve1 = new PairIntArray();
        curve0.add(1, 1);
        curve0.add(2, 1);
        curve0.add(3, 1);
        curve0.add(4, 1);
        curve1.add(3, 2);
        curve1.add(4, 2);
        curve1.add(5, 2);
        
        List<PairIntArray> edges = new ArrayList<PairIntArray>();
        edges.add(curve0);
        edges.add(curve1);
        
        MiscellaneousCurveHelper instance = new MiscellaneousCurveHelper();
        List<PairIntArray> output = instance.pruneAndIncludeAdjacentCurves(edges, 10);
       
        assertTrue(output.size() == 1);
        
        
        /*
        6          @
        5       @
        4    @
        3                         
        2                  
        1    @  @  @  @    
        0                  
          0  1  2  3  4  5 
        */
        
        curve0 = new PairIntArray();
        curve1 = new PairIntArray();
        curve0.add(1, 1);
        curve0.add(2, 1);
        curve0.add(3, 1);
        curve0.add(4, 1);
        curve1.add(1, 4);
        curve1.add(2, 5);
        curve1.add(3, 6);
        
        edges = new ArrayList<PairIntArray>();
        edges.add(curve0);
        edges.add(curve1);
        
        output = instance.pruneAndIncludeAdjacentCurves(edges, 10);
       
        assertTrue(output.size() == 2);
    }

    public void testFindMinIdx() {
        
        PairIntArray closedCurve = getSquare();
        
        MiscellaneousCurveHelper instance = new MiscellaneousCurveHelper();
        
        int result = instance.findMinIdx(closedCurve);
        
        assertTrue(18 == result);
    }
  
    public void testCalculateCentroid() throws Exception {
        
        MiscellaneousCurveHelper helper = new MiscellaneousCurveHelper();
        
        PairIntArray xy = new PairIntArray();
        xy.add(1, 1);
        xy.add(9, 1);
        xy.add(9, 5);
        xy.add(1, 5);
        float[] equalWeights = new float[]{1.f/xy.getN(), 1.f/xy.getN(),
            1.f/xy.getN(), 1.f/xy.getN()}; 
        
        double[] cenXY = helper.calculateXYCentroids(xy, equalWeights);
        assertTrue(cenXY[0] == 5);
        assertTrue(cenXY[1] == 3);
        
        cenXY = helper.calculateXYCentroids(xy);
        assertTrue(cenXY[0] == 5);
        assertTrue(cenXY[1] == 3);
        
        float[] unequalWeights = new float[]{0.5f/xy.getN(), 1.5f/xy.getN(),
            0.5f/xy.getN(), 1.5f/xy.getN()};
        cenXY = helper.calculateXYCentroids(xy, unequalWeights);
        assertTrue(cenXY[0] == 5);
        assertTrue(cenXY[1] == 3);
        
        double[][] xyM = new double[2][4];
        xyM[0] = new double[]{1, 9, 9, 1};
        xyM[1] = new double[]{1, 1, 5, 5};
        Matrix xyMatrix = new Matrix(xyM);
        
        cenXY = helper.calculateXYCentroids(xyMatrix);
        assertTrue(cenXY[0] == 5.0);
        assertTrue(cenXY[1] == 3.0);
        
        PairFloatArray xyf = new PairFloatArray();
        xyf.add(1, 1);
        xyf.add(9, 1);
        xyf.add(9, 5);
        xyf.add(1, 5);
        
        cenXY = helper.calculateXYCentroids(xyf);
        assertTrue(cenXY[0] == 5.0);
        assertTrue(cenXY[1] == 3.0);

        SearchableCurve searchableXY = new SearchableCurve(xy);
        
        cenXY = helper.calculateXYCentroids(searchableXY);
        assertTrue(cenXY[0] == 5.0);
        assertTrue(cenXY[1] == 3.0);
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
    
    public void testDistanceBetweenPointAnd2Lines() {
        
        /*
         5
         4
         3    #     @
         2       @
         1    @
         0
           0  1  2  3  4  5  6
        */
        float x0 = 1;
        float y0 = 1;
        float x1 = 3;
        float y1 = 3;
        
        float xP = 1;
        float yP = 3;
        double expectedDist = Math.sqrt(2);
                
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, 
            xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
        
        /*
         5
         4 #
         3   #      @
         2       @
         1    @
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 1;
        x1 = 3;
        y1 = 3;
        
        xP = 0;
        yP = 4;
        
        // side 1
        double side1 = Math.sqrt(2*2 + 2*2);
        expectedDist = Math.sqrt(4*4 - side1*side1);
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, 
            xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
        
        /*
         5       .  #
         4       .
         3       . 
         2    @  .
         1             @
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 2;
        x1 = 4;
        y1 = 1;
        
        xP = 3;
        yP = 5;
        expectedDist = 3.4785;
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
            
        
        /*
         5         #
         4       
         3        
         2    @        @
         1            
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 2;
        x1 = 4;
        y1 = 2;
        
        xP = 3;
        yP = 5;
        expectedDist = 3;
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
                
        /*
         5    @    
         4         #  
         3        
         2    @  
         1            
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 2;
        x1 = 1;
        y1 = 5;
        
        xP = 3;
        yP = 4;
        expectedDist = 2;
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
    }
   
    public void testDistanceBetweenPointAnd2Lines_2() {
        
        /*
         5
         4
         3    #     @
         2       @
         1    @
         0
           0  1  2  3  4  5  6
        */
        int x0 = 1;
        int y0 = 1;
        int x1 = 3;
        int y1 = 3;
        
        float xP = 1;
        float yP = 3;
        double expectedDist = Math.sqrt(2);
                
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, 
            xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
        
        /*
         5
         4 #
         3   #      @
         2       @
         1    @
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 1;
        x1 = 3;
        y1 = 3;
        
        xP = 0;
        yP = 4;
        
        // side 1
        double side1 = Math.sqrt(2*2 + 2*2);
        expectedDist = Math.sqrt(4*4 - side1*side1);
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, 
            xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
        
        /*
         5       .  #
         4       .
         3       . 
         2    @  .
         1             @
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 2;
        x1 = 4;
        y1 = 1;
        
        xP = 3;
        yP = 5;
        expectedDist = 3.4785;
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
            
        
        /*
         5         #
         4       
         3        
         2    @        @
         1            
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 2;
        x1 = 4;
        y1 = 2;
        
        xP = 3;
        yP = 5;
        expectedDist = 3;
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
                
        /*
         5    @    
         4         #  
         3        
         2    @  
         1            
         0
           0  1  2  3  4  5  6
        */
        x0 = 1;
        y0 = 2;
        x1 = 1;
        y1 = 5;
        
        xP = 3;
        yP = 4;
        expectedDist = 2;
                
        dist = curveHelper.distanceFromPointToALine(x0, y0, x1, y1, xP, yP);
        
        assertTrue(Math.abs(expectedDist - dist) < 0.01);
    }
   
}
