package algorithms.imageProcessing;

import Jama.Matrix;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class MiscellaneousCurveHelperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
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
    
    public void testCatalogNearlyStraightLineSegments() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        xy.add(0, 0);
        xy.add(1, 0);
        xy.add(2, 1);
        xy.add(3, 1);
        xy.add(4, 2);
        xy.add(5, 2);//5
        xy.add(6, 1);
        
        PairIntArray lineSegment;
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        lineSegment = curveHelper.findJaggedLineSegments(xy);
        assertTrue(lineSegment.getX(0) == 0);
        assertTrue(lineSegment.getY(0) == 5);
        
        int[] endIndex = curveHelper.validateJaggedLineSegment(xy,
        /*startIndex*/ 0, /*endIndex*/ 5, /*stepWidth*/ 2, /*dx*/ 1, /*dy*/ 1,
        /*widthIsAlongX*/ Boolean.TRUE);
        assertTrue(endIndex[0] == 5);
        
        endIndex = curveHelper.validateJaggedLineSegment(xy,
        /*startIndex*/ 0, /*endIndex*/ 6, /*stepWidth*/ 2, /*dx*/ 1, /*dy*/ 1,
        /*widthIsAlongX*/ Boolean.TRUE);
        assertTrue(endIndex[0] == 5);
        
        
        xy = new PairIntArray();
        xy.add(0, 0);
        xy.add(-1, 0);
        xy.add(-2, 1);
        xy.add(-3, 1);
        xy.add(-4, 2);
        xy.add(-5, 2);
        xy.add(-6, 1);
        lineSegment = 
            curveHelper.findJaggedLineSegments(xy);
        assertTrue(lineSegment.getX(0) == 0);
        assertTrue(lineSegment.getY(0) == 5);
        
        endIndex = curveHelper.validateJaggedLineSegment(xy,
        /*startIndex*/ 0, /*endIndex*/ 5, /*stepWidth*/ 2, /*dx*/ -1, /*dy*/ 1,
        /*widthIsAlongX*/ Boolean.TRUE);
        assertTrue(endIndex[0] == 5);
        
        endIndex = curveHelper.validateJaggedLineSegment(xy,
        /*startIndex*/ 0, /*endIndex*/ 6, /*stepWidth*/ 2, /*dx*/ -1, /*dy*/ 1,
        /*widthIsAlongX*/ Boolean.TRUE);
        assertTrue(endIndex[0] == 5);
    
    }
    
    public void testCatalogNearlyStraightLineSegments2() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        
        // ==== long jagged line whose avg step width should be 2 so this
        //      entire segment should be considered a line
        xy.add(147, 71); //0
        xy.add(146, 71);
        xy.add(145, 71);
        xy.add(144, 72);
        xy.add(143, 72);
        xy.add(142, 73);
        xy.add(141, 73);
        xy.add(140, 74);
        xy.add(139, 74);
        xy.add(138, 75);
        xy.add(137, 75); //10
        xy.add(136, 76);
        xy.add(135, 76);        
        xy.add(134, 77);
        xy.add(133, 77);
        xy.add(132, 77); //15        
        xy.add(131, 78); //16
        xy.add(130, 79); //17
        xy.add(129, 79);
        xy.add(128, 79);
        xy.add(127, 80); //20
        xy.add(126, 80);
        xy.add(125, 81);
        xy.add(124, 81);
        xy.add(123, 82);
        xy.add(122, 82);
        xy.add(121, 83);
        xy.add(120, 83);
        xy.add(119, 84);
        xy.add(118, 84);
        xy.add(117, 84);//30
        xy.add(116, 85);
        xy.add(115, 85);
        xy.add(114, 86);
        xy.add(113, 86);        
        xy.add(112, 87);
        xy.add(111, 87);
        xy.add(110, 88);
        xy.add(109, 88);
        xy.add(108, 89);
        xy.add(107, 89);//40
        xy.add(106, 89);
        xy.add(105, 90);
        xy.add(104, 90);
        xy.add(103, 91);//44
        
        // add a corner
        xy.add(102, 90); //45
        xy.add(101, 90);
        xy.add(100, 90);
        xy.add(99, 89);
        xy.add(98, 89);        
        xy.add(97, 88); //50
        xy.add(96, 88); 
        xy.add(95, 87);
        xy.add(94, 87);
        xy.add(93, 87);        
        xy.add(92, 86);
        xy.add(91, 86); 
        xy.add(90, 86);
        xy.add(89, 85); //58
     
        xy.add(89, 84); //59 <===
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        PairIntArray lineSegment = 
            curveHelper.findJaggedLineSegments(xy);
        
        
        assertTrue(lineSegment.getX(0) == 0);
        assertTrue(lineSegment.getY(0) == 43);
        
        assertTrue(lineSegment.getX(1) == 45);
        assertTrue(lineSegment.getY(1) == 57);
        
    }
    
    public void testCatalogNearlyStraightLineSegments3() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        
        for (int i = 71; i <= 94; i++) {
            xy.add(i, 62);
        }
        //24 points
        
        for (int i = 95; i <= 141; i++) {
            xy.add(i, 61);
        }
        //47 points
        
        for (int i = 142; i <= 166; i++) {
            xy.add(i, 60);
        }
        //25 points
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        PairIntArray lineSegment = 
            curveHelper.findJaggedLineSegments(xy);
        assertTrue(lineSegment.getX(0) == 0);
        assertTrue(lineSegment.getY(0) == (xy.getN() - 1));
        
    }
    
    public void testCatalogNearlyStraightLineSegments4() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        
        for (int j = 212; j <= 217; j++) {
            xy.add(187, j);
        }
        
        for (int j = 218; j <= 230; j++) {
            xy.add(186, j);
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        PairIntArray lineSegment = 
            curveHelper.findJaggedLineSegments(xy);
        assertTrue(lineSegment.getX(0) == 0);
        assertTrue(lineSegment.getY(0) == (xy.getN() - 1));
        
    }
    
    public void testJaggedLines6() throws Exception {
        PairIntArray xy = getEdge1();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        /*
         [junit] INFO: SEGMENT: 0 : 26
         [junit] INFO: SEGMENT: 37 : 70
         [junit] INFO: SEGMENT: 71 : 99
         [junit] INFO: SEGMENT: 109 : 144
         [junit] INFO: SEGMENT: 145 : 208
         [junit] INFO: SEGMENT: 215 : 256
         [junit] INFO: SEGMENT: 256 : 315
         [junit] INFO: SEGMENT: 318 : 348
         [junit] INFO: SEGMENT: 353 : 409
         [junit] INFO: SEGMENT: 410 : 459
         [junit] INFO: SEGMENT: 463 : 478
         [junit] INFO: SEGMENT: 478 : 498
         [junit] INFO: SEGMENT: 515 : 526
         [junit] INFO: SEGMENT: 539 : 547
         [junit] INFO: SEGMENT: 566 : 583
         [junit] INFO: SEGMENT: 596 : 614
        */
        
        PairIntArray segments = curveHelper.findJaggedLineSegments(xy);
        
        for (int i = 0; i < segments.getN(); i++) {
            log.info("SEGMENT: " + segments.getX(i) + " : " + segments.getY(i));
        }
    }
    
    public void testSortByX() throws Exception {
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
   
        PairIntArray xy = new PairIntArray();
        int nPoints = 100;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);

        Set<Integer> set = new HashSet<Integer>();
        
        for (int i = 0; i < nPoints; i++) {
            int start = sr.nextInt(1000000);
            while (set.contains(Integer.valueOf(start))) {
                start = sr.nextInt(1000000);
            }
            xy.add(start*100, (start*100 + 50));
        }
        
        curveHelper.sortByX(xy);
        
        int lastX = Integer.MIN_VALUE;
        for (int i = 0; i < nPoints; i++) {
            int x = xy.getX(i);
            assertTrue(x >= lastX);
            lastX = x;
        }
    }
    
    public PairIntArray getEdge0() {
        
        PairIntArray xy = new PairIntArray();
        
        xy.add(238, 203);
        xy.add(239, 204);
        xy.add(239, 205);
        xy.add(239, 206);
        xy.add(238, 207);
        xy.add(238, 208);
        xy.add(238, 209);
        xy.add(238, 210);
        xy.add(238, 211);
        xy.add(238, 212);
        xy.add(238, 213);
        xy.add(238, 214);
        xy.add(238, 215);
        xy.add(238, 216);
        xy.add(238, 217);
        xy.add(238, 218);
        xy.add(237, 219);
        xy.add(237, 220);
        xy.add(236, 221);
        xy.add(236, 222);
        xy.add(235, 223);
        xy.add(234, 224);
        xy.add(233, 225);
        xy.add(232, 225);
        xy.add(231, 226);
        xy.add(230, 226);
        xy.add(229, 227);
        xy.add(228, 227);
        xy.add(227, 228);
        xy.add(226, 228);
        xy.add(225, 229);
        xy.add(224, 230);
        xy.add(223, 230);
        xy.add(222, 231);
        xy.add(221, 232);
        xy.add(220, 232);
        xy.add(219, 233);
        xy.add(218, 233);
        xy.add(217, 234);
        xy.add(216, 235);
        xy.add(215, 235);
        xy.add(214, 236);
        xy.add(213, 236);
        xy.add(212, 237);
        xy.add(211, 238);
        xy.add(210, 238);
        xy.add(209, 239);
        xy.add(208, 239);
        xy.add(207, 240);
        xy.add(206, 241);
        xy.add(205, 241);
        xy.add(204, 241);
        xy.add(203, 241);
        xy.add(202, 241);
        xy.add(201, 240);
        xy.add(200, 240);
        xy.add(199, 239);
        xy.add(198, 238);
        xy.add(197, 238);
        xy.add(196, 237);
        xy.add(195, 237);
        xy.add(194, 236);
        xy.add(193, 235);
        xy.add(192, 235);
        xy.add(191, 234);
        xy.add(190, 234);
        xy.add(189, 233);
        xy.add(188, 232);
        xy.add(187, 231);
        xy.add(186, 230);
        xy.add(186, 229);
        xy.add(186, 228);
        xy.add(186, 227);
        xy.add(186, 226);
        xy.add(186, 225);
        xy.add(186, 224);
        xy.add(186, 223);
        xy.add(186, 222);
        xy.add(186, 221);
        xy.add(186, 220);
        xy.add(186, 219);
        xy.add(186, 218);
        xy.add(187, 217);
        xy.add(187, 216);
        xy.add(187, 215);
        xy.add(187, 214);
        xy.add(187, 213);
        xy.add(187, 212);
        xy.add(188, 211);
        xy.add(189, 211);
        xy.add(190, 212);
        xy.add(191, 213);
        xy.add(192, 214);
        xy.add(193, 214);
        xy.add(194, 215);
        xy.add(195, 216);
        xy.add(196, 217);
        xy.add(197, 217);
        xy.add(198, 218);
        xy.add(199, 218);
        xy.add(200, 219);
        xy.add(201, 220);
        xy.add(202, 220);
        xy.add(203, 221);
        xy.add(204, 221);
        xy.add(205, 221);
        xy.add(206, 221);
        xy.add(207, 221);
        xy.add(208, 221);
        xy.add(209, 220);
        xy.add(210, 219);
        xy.add(211, 219);
        xy.add(212, 218);
        xy.add(213, 218);
        xy.add(214, 217);
        xy.add(215, 216);
        xy.add(216, 216);
        xy.add(217, 215);
        xy.add(218, 215);
        xy.add(219, 214);
        xy.add(220, 213);
        xy.add(221, 213);
        xy.add(222, 212);
        xy.add(223, 212);
        xy.add(224, 211);
        xy.add(225, 210);
        xy.add(226, 210);
        xy.add(227, 209);
        xy.add(228, 209);
        xy.add(229, 208);
        xy.add(230, 208);
        xy.add(231, 207);
        xy.add(232, 206);
        xy.add(233, 206);
        xy.add(234, 205);
        xy.add(235, 204);
        xy.add(236, 204);
        xy.add(237, 203);
        xy.add(238, 202);
        xy.add(237, 201);
        xy.add(236, 200);
        xy.add(235, 199);
        xy.add(234, 199);
        xy.add(233, 198);
        xy.add(232, 198);
        xy.add(231, 197);
        xy.add(230, 197);
        xy.add(229, 196);
        xy.add(228, 196);
        xy.add(227, 195);
        xy.add(226, 195);
        xy.add(225, 194);
        xy.add(224, 194);
        xy.add(223, 193);
        xy.add(222, 193);
        xy.add(221, 193);
        xy.add(220, 193);
        xy.add(219, 193);
        xy.add(218, 193);
        xy.add(217, 194);
        xy.add(216, 194);
        xy.add(215, 195);
        xy.add(214, 196);
        xy.add(213, 196);
        xy.add(212, 196);
        xy.add(211, 197);
        xy.add(210, 198);
        xy.add(209, 198);
        xy.add(208, 199);
        xy.add(207, 199);
        xy.add(206, 200);
        xy.add(205, 201);
        xy.add(204, 201);
        xy.add(203, 201);
        xy.add(202, 202);
        xy.add(201, 203);
        xy.add(200, 203);
        xy.add(199, 204);
        xy.add(198, 204);
        xy.add(197, 205);
        xy.add(196, 206);
        xy.add(196, 207);

        return xy;
    }

    public PairIntArray getEdge1() {
        
        PairIntArray xy = new PairIntArray();
        xy.add(198, 145);
        xy.add(198, 146);
        xy.add(198, 147);
        xy.add(198, 148);
        xy.add(198, 149);
        xy.add(198, 150);
        xy.add(198, 151);
        xy.add(198, 152);
        xy.add(198, 153);
        xy.add(198, 154);
        xy.add(198, 155);
        xy.add(198, 156);
        xy.add(198, 157);
        xy.add(198, 158);
        xy.add(198, 159);
        xy.add(198, 160);
        xy.add(198, 161);
        xy.add(198, 162);
        xy.add(198, 163);
        xy.add(198, 164);
        xy.add(198, 165);
        xy.add(198, 166);
        xy.add(198, 167);
        xy.add(198, 168);
        xy.add(198, 169);
        xy.add(198, 170);
        xy.add(198, 171);
        xy.add(199, 172);
        xy.add(200, 173);
        xy.add(201, 172);
        xy.add(202, 171);
        xy.add(203, 170);
        xy.add(204, 169);
        xy.add(204, 168);
        xy.add(205, 167);
        xy.add(206, 166);
        xy.add(207, 165);
        xy.add(207, 164);//<=== 37
        xy.add(208, 163);
        xy.add(208, 162);
        xy.add(209, 161);
        xy.add(210, 160);
        xy.add(211, 159);
        xy.add(211, 158);
        xy.add(212, 157);
        xy.add(212, 156);
        xy.add(213, 155);
        xy.add(214, 154);
        xy.add(215, 153);
        xy.add(215, 152);
        xy.add(216, 151);
        xy.add(217, 150);
        xy.add(217, 149);
        xy.add(218, 148);
        xy.add(218, 147);
        xy.add(219, 146);
        xy.add(220, 145);
        xy.add(221, 144);
        xy.add(222, 143);
        xy.add(222, 142);
        xy.add(223, 141);
        xy.add(223, 140);
        xy.add(224, 139);
        xy.add(224, 138);
        xy.add(225, 137);
        xy.add(226, 136);
        xy.add(226, 135);
        xy.add(227, 134);
        xy.add(228, 133);
        xy.add(228, 132); //<=== 69
        xy.add(229, 131);
        xy.add(229, 130);
        xy.add(229, 129);
        xy.add(229, 128);
        xy.add(229, 127);
        xy.add(229, 126);
        xy.add(229, 125);
        xy.add(229, 124);
        xy.add(229, 123);
        xy.add(230, 122);
        xy.add(230, 121);
        xy.add(230, 120);
        xy.add(230, 119);
        xy.add(230, 118);
        xy.add(230, 117);
        xy.add(230, 116);
        xy.add(230, 115);
        xy.add(230, 114);
        xy.add(230, 113);
        xy.add(230, 112);
        xy.add(231, 111);
        xy.add(231, 110);
        xy.add(231, 109);
        xy.add(231, 108);
        xy.add(231, 107);
        xy.add(231, 106);
        xy.add(231, 105);
        xy.add(231, 104);
        xy.add(231, 103);
        xy.add(231, 102);
        xy.add(230, 101);
        xy.add(229, 101);
        xy.add(228, 101);
        xy.add(227, 102);
        xy.add(226, 103);
        xy.add(225, 104);
        xy.add(224, 105);
        xy.add(223, 106);
        xy.add(222, 107);
        xy.add(222, 108); //<=== 109
        xy.add(221, 109);
        xy.add(220, 110);
        xy.add(219, 111);
        xy.add(218, 112);
        xy.add(218, 113);
        xy.add(217, 114);
        xy.add(216, 115);
        xy.add(216, 116);
        xy.add(215, 117);
        xy.add(214, 118);
        xy.add(213, 119);
        xy.add(213, 120);
        xy.add(212, 121);
        xy.add(211, 122);
        xy.add(211, 123);
        xy.add(210, 124);
        xy.add(209, 125);
        xy.add(209, 126);
        xy.add(208, 127);
        xy.add(207, 128);
        xy.add(206, 129);
        xy.add(206, 130);
        xy.add(205, 131);
        xy.add(204, 132);
        xy.add(203, 133);
        xy.add(203, 134);
        xy.add(202, 135);
        xy.add(201, 136);
        xy.add(200, 137);
        xy.add(200, 138);
        xy.add(199, 139);
        xy.add(199, 140);
        xy.add(198, 141);
        xy.add(198, 142);
        xy.add(197, 143); //<=== 144
        xy.add(196, 143);
        xy.add(195, 143);
        xy.add(194, 143);
        xy.add(193, 143);
        xy.add(192, 143);
        xy.add(191, 143);
        xy.add(190, 143);
        xy.add(189, 143);
        xy.add(188, 143);
        xy.add(187, 144);
        xy.add(186, 144);
        xy.add(185, 144);
        xy.add(184, 144);
        xy.add(183, 144);
        xy.add(182, 144);
        xy.add(181, 144);
        xy.add(180, 144);
        xy.add(179, 144);
        xy.add(178, 144);
        xy.add(177, 144);
        xy.add(176, 144);
        xy.add(175, 144);
        xy.add(174, 144);
        xy.add(173, 144);
        xy.add(172, 144);
        xy.add(171, 144);
        xy.add(170, 144);
        xy.add(169, 144);
        xy.add(168, 144);
        xy.add(167, 144);
        xy.add(166, 144);
        xy.add(165, 144);
        xy.add(164, 144);
        xy.add(163, 144);
        xy.add(162, 144);
        xy.add(161, 144);
        xy.add(160, 144);
        xy.add(159, 144);
        xy.add(158, 145);
        xy.add(157, 145);
        xy.add(156, 145);
        xy.add(155, 145);
        xy.add(154, 145);
        xy.add(153, 145);
        xy.add(152, 145);
        xy.add(151, 145);
        xy.add(150, 145);
        xy.add(149, 145);
        xy.add(148, 145);
        xy.add(147, 145);
        xy.add(146, 145);
        xy.add(145, 145);
        xy.add(144, 145);
        xy.add(143, 145);
        xy.add(142, 145);
        xy.add(141, 145);
        xy.add(140, 145);
        xy.add(139, 145);
        xy.add(138, 145);
        xy.add(137, 145);
        xy.add(136, 145);
        xy.add(135, 145);
        xy.add(134, 145);
        xy.add(133, 145);
        xy.add(132, 146);
        xy.add(131, 146);
        xy.add(130, 145);
        xy.add(129, 145);
        xy.add(128, 145);
        xy.add(127, 146);
        xy.add(126, 146);
        xy.add(125, 146);
        xy.add(124, 146);
        xy.add(123, 146);
        xy.add(122, 146);
        xy.add(121, 146);
        xy.add(120, 146);
        xy.add(119, 146);
        xy.add(118, 146);
        xy.add(117, 146);
        xy.add(116, 146);
        xy.add(115, 146);
        xy.add(114, 146);
        xy.add(113, 146);
        xy.add(112, 146);
        xy.add(111, 146);
        xy.add(110, 146);
        xy.add(109, 146);
        xy.add(108, 146);
        xy.add(107, 146);
        xy.add(106, 146);
        xy.add(105, 146);
        xy.add(104, 146);
        xy.add(103, 146);
        xy.add(102, 147);
        xy.add(101, 147);
        xy.add(100, 147);
        xy.add(99, 147);
        xy.add(98, 147);
        xy.add(97, 147);
        xy.add(96, 147);
        xy.add(95, 147);
        xy.add(94, 147);
        xy.add(93, 147);
        xy.add(92, 147);
        xy.add(91, 147);
        xy.add(90, 147);
        xy.add(89, 147);
        xy.add(88, 147);
        xy.add(87, 147);
        xy.add(86, 147);
        xy.add(85, 147);
        xy.add(84, 146);
        xy.add(83, 145);
        xy.add(82, 144);
        xy.add(81, 144); //<---- 260
        xy.add(80, 143);
        xy.add(79, 142);
        xy.add(78, 141);
        xy.add(77, 141);
        xy.add(76, 140);
        xy.add(75, 139);
        xy.add(74, 138);
        xy.add(73, 138);
        xy.add(72, 137);
        xy.add(71, 136);
        xy.add(70, 135);
        xy.add(69, 135);
        xy.add(68, 134);
        xy.add(67, 133);
        xy.add(66, 133);
        xy.add(65, 132);
        xy.add(64, 131);
        xy.add(63, 130);
        xy.add(62, 130);
        xy.add(61, 129);
        xy.add(60, 128);
        xy.add(59, 127);
        xy.add(58, 127);
        xy.add(57, 126);
        xy.add(56, 125);
        xy.add(55, 124);
        xy.add(54, 124);
        xy.add(53, 123);
        xy.add(52, 122);
        xy.add(51, 121);
        xy.add(50, 120);
        xy.add(49, 120);
        xy.add(48, 119);
        xy.add(47, 118);
        xy.add(46, 118);
        xy.add(45, 117);
        xy.add(44, 116);
        xy.add(43, 115);
        xy.add(42, 115);
        xy.add(41, 114);
        xy.add(40, 113);
        xy.add(39, 113);
        xy.add(38, 112);
        xy.add(37, 111);
        xy.add(36, 111);
        xy.add(35, 110);
        xy.add(34, 109);
        xy.add(33, 108);
        xy.add(32, 107);
        xy.add(31, 107);
        xy.add(30, 106);
        xy.add(29, 105);
        xy.add(28, 105);
        xy.add(27, 104);
        xy.add(26, 103);  //<--- 315
        
        xy.add(25, 104);
        xy.add(24, 105);
        xy.add(24, 106);
        xy.add(24, 107);
        xy.add(24, 108);
        xy.add(24, 109);
        xy.add(24, 110);
        xy.add(24, 111);
        xy.add(24, 112);
        xy.add(24, 113);
        xy.add(24, 114);
        xy.add(24, 115);
        xy.add(24, 116);
        xy.add(24, 117);
        xy.add(24, 118);
        xy.add(24, 119);
        xy.add(24, 120);
        xy.add(24, 121);
        xy.add(24, 122);
        xy.add(24, 123);
        xy.add(24, 124);
        xy.add(24, 125);
        xy.add(24, 126);
        xy.add(24, 127);
        xy.add(24, 128);
        xy.add(24, 129);
        xy.add(24, 130);
        xy.add(24, 131);
        xy.add(24, 132);
        xy.add(24, 133);
        xy.add(24, 134);
        xy.add(24, 135);
        xy.add(24, 136);
        xy.add(25, 137);
        xy.add(25, 138);
        xy.add(26, 139);
        xy.add(27, 140);
        
        xy.add(28, 140);  // <--- 353
        xy.add(29, 141);
        xy.add(30, 142);
        xy.add(31, 143);
        xy.add(32, 144);
        xy.add(33, 144);
        xy.add(34, 145);
        xy.add(35, 146);
        xy.add(36, 147);
        xy.add(37, 147);
        xy.add(38, 148);
        xy.add(39, 149);
        xy.add(40, 150);
        xy.add(41, 151);
        xy.add(42, 151);
        xy.add(43, 152);
        xy.add(44, 153);
        xy.add(45, 154);
        xy.add(46, 155);
        xy.add(47, 155);
        xy.add(48, 156);
        xy.add(49, 157);
        xy.add(50, 158);
        xy.add(51, 158);
        xy.add(52, 159);
        xy.add(53, 160);
        xy.add(54, 161);
        xy.add(55, 162);
        xy.add(56, 163);
        xy.add(57, 163);
        xy.add(58, 164);
        xy.add(59, 165);
        xy.add(60, 165);
        xy.add(61, 166);
        xy.add(62, 167);
        xy.add(63, 168);
        xy.add(64, 169);
        xy.add(65, 169);
        xy.add(66, 170);
        xy.add(67, 171);
        xy.add(68, 172);
        xy.add(69, 172);
        xy.add(70, 173);
        xy.add(71, 174);
        xy.add(72, 174);
        xy.add(73, 175);
        xy.add(74, 176);
        xy.add(75, 177);
        xy.add(76, 178);
        xy.add(77, 178);
        xy.add(78, 179);
        xy.add(79, 180);
        xy.add(80, 180);
        xy.add(81, 181);
        xy.add(82, 182);
       
        xy.add(83, 183);  // <-- 408
        xy.add(84, 184);
        xy.add(85, 184);
        xy.add(86, 184);
        xy.add(87, 184);
        xy.add(88, 184);
        xy.add(89, 184);
        xy.add(90, 184);
        xy.add(91, 184);
        xy.add(92, 184);
        xy.add(93, 184);
        xy.add(94, 184);
        xy.add(95, 184);
        xy.add(96, 184);
        xy.add(97, 184);
        xy.add(98, 184);
        xy.add(99, 184);
        xy.add(100, 184);
        xy.add(101, 184);
        xy.add(102, 184);
        xy.add(103, 184);
        xy.add(104, 184);
        xy.add(105, 184);
        xy.add(106, 184);
        xy.add(107, 184);
        xy.add(108, 184);
        xy.add(109, 184);
        xy.add(110, 184);
        xy.add(111, 183);
        xy.add(112, 183);
        xy.add(113, 183);
        xy.add(114, 183);
        xy.add(115, 183);
        xy.add(116, 183);
        xy.add(117, 183);
        xy.add(118, 183);
        xy.add(119, 183);
        xy.add(120, 183);
        xy.add(121, 183);
        xy.add(122, 183);
        xy.add(123, 183);
        xy.add(124, 183);
        xy.add(125, 183);
        xy.add(126, 183);
        xy.add(127, 183);
        xy.add(128, 183);
        xy.add(129, 183);
        xy.add(130, 183);
        xy.add(131, 183);
        xy.add(132, 183);
        xy.add(133, 183);
        xy.add(134, 183);
        xy.add(135, 184);
        xy.add(136, 185);
        xy.add(135, 186);
        xy.add(135, 187);
        xy.add(135, 188);
        xy.add(135, 189);
        xy.add(135, 190);
        xy.add(135, 191);
        xy.add(135, 192);
        xy.add(135, 193);
        xy.add(135, 194);
        xy.add(135, 195);
        xy.add(135, 196);
        xy.add(135, 197);
        xy.add(135, 198);
        xy.add(135, 199);
        xy.add(135, 200);
        xy.add(135, 201);
        xy.add(135, 202);
        xy.add(136, 203);
        xy.add(137, 204);
        xy.add(138, 205);
        xy.add(139, 206);
        xy.add(140, 206); // <--- 483
        xy.add(141, 207);
        xy.add(142, 207);
        xy.add(143, 208);
        xy.add(144, 208);
        xy.add(145, 208);
        xy.add(146, 209);
        xy.add(147, 209);
        xy.add(148, 210);
        xy.add(149, 210);
        xy.add(150, 211);
        xy.add(151, 211);
        xy.add(152, 211);
        xy.add(153, 212);
        xy.add(154, 212);
        xy.add(155, 212); // <--- 498
         
        xy.add(156, 211);
        xy.add(157, 210);
        xy.add(158, 209);
        xy.add(159, 209);
        xy.add(160, 208);
        xy.add(161, 207);
        xy.add(162, 206);
        xy.add(163, 205);
        xy.add(164, 204);
        xy.add(165, 203);
        xy.add(166, 202);
        xy.add(167, 201);
        xy.add(167, 200);
        xy.add(167, 199);
        xy.add(167, 198);
        xy.add(168, 197);
        xy.add(168, 196);
        xy.add(168, 195);
        xy.add(168, 194);
        xy.add(168, 193);
        xy.add(168, 192);
        xy.add(168, 191);
        xy.add(168, 190);
        xy.add(168, 189);
        xy.add(168, 188);
        xy.add(168, 187);
        xy.add(168, 186);
        xy.add(168, 185);
        xy.add(169, 184);
        xy.add(169, 183);
        xy.add(170, 182);
        xy.add(171, 181);
        xy.add(170, 180);
        xy.add(169, 179);
        xy.add(169, 178);
        xy.add(169, 177);
        xy.add(169, 176);
        xy.add(169, 175);
        xy.add(169, 174);
        xy.add(170, 173);
        xy.add(170, 172);
        xy.add(170, 171);
        xy.add(170, 170);
        xy.add(170, 169);
        xy.add(170, 168);
        xy.add(170, 167);
        xy.add(170, 166);
        xy.add(170, 165);
        xy.add(170, 164);
        xy.add(171, 163);
        xy.add(171, 162);
        xy.add(170, 161);
        xy.add(169, 160);
        xy.add(168, 161);
        xy.add(167, 162);
        xy.add(166, 162);
        xy.add(165, 163);
        xy.add(164, 164);
        xy.add(163, 165);
        xy.add(162, 166);
        xy.add(161, 167);
        xy.add(160, 167);
        xy.add(159, 168);
        xy.add(158, 169);
        xy.add(157, 169);
        xy.add(156, 169);
        xy.add(155, 169);
        xy.add(154, 168); // <-- 566
        
        xy.add(153, 168);
        xy.add(152, 168);
        xy.add(151, 168);
        xy.add(150, 167);
        xy.add(149, 167);
        xy.add(148, 166);
        xy.add(147, 166);
        xy.add(146, 166);
        xy.add(145, 165);
        xy.add(144, 165);
        xy.add(143, 165);
        xy.add(142, 164);
        xy.add(141, 164);
        xy.add(140, 164);
        xy.add(139, 163);
        xy.add(138, 162);
        xy.add(137, 161); // <-- 583
        
        xy.add(138, 160);
        xy.add(139, 159);
        xy.add(140, 158);
        xy.add(141, 157);
        xy.add(142, 156);
        xy.add(143, 155);
        xy.add(144, 154);
        xy.add(145, 154);
        xy.add(146, 153);
        xy.add(147, 152);
        xy.add(148, 152);
        xy.add(149, 151);
        xy.add(150, 152); // <--- 596
       
        xy.add(151, 152);
        xy.add(152, 152);
        xy.add(153, 153);
        xy.add(154, 153);
        xy.add(155, 153);
        xy.add(156, 153);
        xy.add(157, 154);
        xy.add(158, 154);
        xy.add(159, 155);
        xy.add(160, 155);
        xy.add(161, 155);
        xy.add(162, 156);
        xy.add(163, 156);
        xy.add(164, 156);
        xy.add(165, 157);
        xy.add(166, 157);
        xy.add(167, 157);
        xy.add(168, 158); // <-- 614
        xy.add(169, 159);
        return xy;
    }
}
