package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SkylineExtractorTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SkylineExtractorTest() {
    }
    
    public void estOrderByProximity() throws Exception {
                
        int x0 = 468;
        int y0 = 98;
        List<PairInt> expected = new ArrayList<PairInt>();
        expected.add(new PairInt(x0, y0));
        expected.add(new PairInt(469, 97));
        expected.add(new PairInt(469, 96));
        expected.add(new PairInt(469, 95));
        expected.add(new PairInt(470, 95));
        expected.add(new PairInt(471, 94));
        expected.add(new PairInt(472, 93));
        expected.add(new PairInt(473, 92));
        expected.add(new PairInt(474, 92));
        expected.add(new PairInt(475, 91));
        expected.add(new PairInt(476, 92));
        expected.add(new PairInt(477, 91));
        expected.add(new PairInt(478, 92));
        expected.add(new PairInt(477, 93));        
        expected.add(new PairInt(476, 94));
        expected.add(new PairInt(475, 95));
        expected.add(new PairInt(474, 95));
        expected.add(new PairInt(473, 95));
        expected.add(new PairInt(472, 96));
        expected.add(new PairInt(471, 97));
        expected.add(new PairInt(470, 98));
        expected.add(new PairInt(469, 98));
        
        Set<PairInt> points = new HashSet<PairInt>(expected);
        
        SkylineExtractor instance = new SkylineExtractor();
        List<PairInt> result = instance.orderByProximity(points);
        
        assertTrue(result.size() == expected.size());
        
        for (int i = 0; i < result.size(); ++i) {
            PairInt r = result.get(i);
            PairInt e = expected.get(i);
            assertTrue(r.equals(e));
        }
        
    }
    
    public void testOrderByIncreasingXThenY() throws Exception {
        
        int[] x = new int[]{7,  6,   5, 4,  4, 3, 2,  1, 0};
        int[] y = new int[]{10, 10, 10, 0, 10, 4, 10, 10, 10};
        
        int[] expectedX = new int[]{0,   1,  2,  3,  4,  4,  5,  6,  7};
        int[] expectedY = new int[]{10, 10, 10,  4,  0, 10, 10, 10, 10};
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < x.length; ++i) {
            points.add(new PairInt(x[i], y[i]));
        }
        
        SkylineExtractor skylineExtractor = new SkylineExtractor();
        
        DoubleLinkedCircularList outputList = new DoubleLinkedCircularList();
        Map<Integer, HeapNode> outputMap = new
            HashMap<Integer, HeapNode>();
        
        skylineExtractor.orderByIncreasingXThenY(points, outputList, outputMap);
        
        int i = 0;
        HeapNode cn = outputList.getSentinel().getLeft();
        while (cn != outputList.getSentinel()) {
            PairInt p = (PairInt)cn.getData();
            assertTrue(p.getX() == expectedX[i]);
            assertTrue(p.getY() == expectedY[i]);
            ++i;
            cn = cn.getLeft();
        }
        
        assertTrue(outputMap.size() == 8);
    }
    
    public void testPopulatePolygon() throws Exception {
        
        SkylineExtractor skylineExtractor = new SkylineExtractor();
        
        /*
        y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        y = 10 + 0
        */
        float[] polynomialCoeff = new float[] {10.0f, 0.0f, 0.f};
        float[] x = new float[]{0,    1,     2,     3};
        float[] y = new float[]{10f,  10.0f, 10.0f, 10.0f};
        
        float dist = 2;
        
        Set<PairInt> points = new HashSet<PairInt>();
        points.add(new PairInt(1, 11));
        points.add(new PairInt(2, 11));
        points.add(new PairInt(3, 11));
        points.add(new PairInt(1, 9));
        points.add(new PairInt(2, 9));
        points.add(new PairInt(3, 9));

        /*
        y = c0*1 + c1*x[i] + c2*x[i]*x[i]
           
        //    0     1     2     3      4    5     6   7   8
        outX=[0,    1,    2,    3,      3,  2,    1,  0,  0
        outY=[12,   12,   12,   12,     8,  8,    8,  8,  12

        1 2 3 4
        0     5
        9 8 7 6  n=11
        */
        float[] outputXPoly = new float[(2 * x.length) + 1];
        float[] outputYPoly = new float[outputXPoly.length];

        skylineExtractor.populatePolygon(x, y, dist, outputXPoly, outputYPoly, 
            polynomialCoeff, 200, 200);

        float[] expectedOutputX = new float[]{
            0,    1,    2,    3,      3,  2,    1,  0,  0};

        float[] expectedOutputY = new float[]{
            12,   12,   12,   12,     8,  8,    8,  8,  12};

        for (int i = 0; i < outputXPoly.length; i++) {
            assertTrue(Math.abs(outputXPoly[i] - expectedOutputX[i]) < 0.1);
            assertTrue(Math.abs(outputYPoly[i] - expectedOutputY[i]) < 0.1);
        }

        int n = skylineExtractor.nPointsInPolygon(points, outputXPoly,
            outputYPoly);
        
        assertTrue(n == points.size());
        
        points.add(new PairInt(100, 100));
        
        n = skylineExtractor.nPointsInPolygon(points, outputXPoly,
            outputYPoly);
        
        assertTrue(n == 6);
    }
   
    public void testPopulatePolygon2() throws Exception {
        
        SkylineExtractor skylineExtractor = new SkylineExtractor();
        
        /*
        y = x^2 + 3*x
        
        */
        
        float[] polynomialCoeff = new float[] {0.0f, 3.0f, 1.0f};
        
        float[] x = new float[]{-3.5f, -3.f,  -1.5f, 0.f, 0.5f};
        float[] y = new float[]{ 2.f,   0.f, -2.25f, 0.f, 2.0f};
        
        float dist = 1;
        
        float[] outputXPoly = new float[(2 * x.length) + 1];
        float[] outputYPoly = new float[outputXPoly.length];

        skylineExtractor.populatePolygon(x, y, dist, outputXPoly, outputYPoly, 
            polynomialCoeff, 200, 200);

        for (int i = 0; i < outputXPoly.length; i++) {
            System.out.println("i=" + i 
                + " (" + outputXPoly[i] + " , " + outputYPoly[i] + ")");
        }
        
        int z = 1;
    }
    
    public static void main(String[] args) {
        
        try {
            SkylineExtractorTest test = new SkylineExtractorTest();

            //test.estOrderByProximity();
            //test.testOrderByIncreasingXThenY();
            //test.testPopulatePolygon();
            test.testPopulatePolygon2();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
