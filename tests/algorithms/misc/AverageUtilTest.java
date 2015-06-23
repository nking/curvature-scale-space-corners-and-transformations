package algorithms.misc;

import algorithms.util.PairIntArray;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class AverageUtilTest extends TestCase {
    
    public AverageUtilTest() {
    }
        
    public void testCalculate1() throws Exception {
        
        /*
        k=3
        curve=10 4's
        2  3  4  3  2  4  4  4  4  4
           5  9  10 9  9  10 12 12 12
        
        2  2  3  3  3  3  3  4  4  4
        */
        AverageUtil avgUtil = new AverageUtil();
        
        int[] curveY = new int[]{2, 3, 4, 3, 2, 4, 4, 4, 4, 4};
        int kPoints = 3;
        int[] result = avgUtil.calculateBoxCarAverage(curveY, kPoints);
                
        int[] expected = new int[]{2,  3,  3,  3,  3,  3,  3,  4,  4,  4};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
    }
    
    public void testCalculate2() throws Exception {
        
        /*        
        k=5
        x:  0 1 2 3 4 5 6 7
        
        y:  2 3 4 3 2 4 4 4
            2 2 4 3 2 3 3 3
        */
        
        AverageUtil avgUtil = new AverageUtil();
        
        PairIntArray curve = new PairIntArray(8);
        curve.add(0, 2);
        curve.add(1, 3);
        curve.add(2, 4);
        curve.add(3, 3);
        curve.add(4, 2);
        curve.add(5, 4);
        curve.add(6, 4);
        curve.add(7, 4);
        
        int kPoints = 5;
        
        int[] expectedX = new int[]{0, 1, 2, 3, 4, 5, 6, 7};
        int[] expectedY = new int[]{2, 3, 3, 3, 3, 3, 3, 3};
        
        PairIntArray result = avgUtil.calculateBoxCarAverage(curve, kPoints);
        
        assertTrue(result.getN() == expectedX.length);
        
        for (int i = 0; i < result.getN(); ++i) {
            int x = result.getX(i);
            int y = result.getY(i);
            assertTrue(x == expectedX[i]);
            assertTrue(y == expectedY[i]);
        }
    }
    
    public void testCalculate3() throws Exception {
        
        // test exceptions
        
        AverageUtil avgUtil = new AverageUtil();
        
        PairIntArray curve = null;
        
        int kPoints = 5;
        
        boolean caughtException = false;
        try {
            PairIntArray result = avgUtil.calculateBoxCarAverage(curve, kPoints);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        //-----------------------------------
        curve = new PairIntArray();
        curve.add(0, 2);
        curve.add(1, 3);
        curve.add(2, 4);
        caughtException = false;
        try {
            PairIntArray result = avgUtil.calculateBoxCarAverage(curve, kPoints);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        //------------------------------------
        int[] curveY = null;
        caughtException = false;
        try {
            int[] result = avgUtil.calculateBoxCarAverage(curveY, kPoints);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        //------------------------------------
        curveY = new int[]{1, 2, 2};
        caughtException = false;
        try {
            int[] result = avgUtil.calculateBoxCarAverage(curveY, kPoints);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
    }
    
    public void testBin0() throws Exception {
        
        /*
        k=3
        curve=10 4's
        4 4 4 4 4 4 4 4 4 4
         
            4     4     4 1
        */
        AverageUtil avgUtil = new AverageUtil();
        
        int[] curveY = new int[]{4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
        int kPoints = 3;
        int[] result = avgUtil.bin(curveY, kPoints);
                
        int[] expected = new int[]{4, 4, 4, 1};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
    }
    
    public void testBin1() throws Exception {
        
        /*
        k=2
        curve=10 4's
        4 4 4 4 4 4 4 4 4 4
         
          4   4   4   4   4
        */
        AverageUtil avgUtil = new AverageUtil();
        
        int[] curveY = new int[]{4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
        int kPoints = 2;
        int[] result = avgUtil.bin(curveY, kPoints);
                
        int[] expected = new int[]{4, 4, 4, 4, 4};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
    }
    
    public void tsetBin2() throws Exception {
        
        /*
        k=5
        x:  0 1 2 3 4 5 6 7 
            - - - - -|- - - - -
                2       3
            
        y:  2 3 4 3 2 4 4 4
            
            - - - - -|- - - - -
                2        2
        */
        
        int[] expectedX = new int[]{2, 3};
        int[] expectedY = new int[]{2, 2};
        
        AverageUtil avgUtil = new AverageUtil();
        
        PairIntArray curve = new PairIntArray(8);
        curve.add(0, 2);
        curve.add(1, 3);
        curve.add(2, 4);
        curve.add(3, 3);
        curve.add(4, 2);
        curve.add(5, 4);
        curve.add(6, 4);
        curve.add(7, 4);
        
        int kPoints = 5;
        
        PairIntArray result = avgUtil.bin(curve, kPoints);
        
        assertTrue(result.getN() == expectedX.length);
        
        for (int i = 0; i < result.getN(); ++i) {
            int x = result.getX(i);
            int y = result.getY(i);
            assertTrue(x == expectedX[i]);
            assertTrue(y == expectedY[i]);
        }
    }
    
    public void tsetBin3() throws Exception {
        
        /*
        k=2
        x:  0 1 2 3 4 5 6 7
            - -|- -|- -|- -|
             0   2   4   6             
            
        y:  2 3 4 3 2 4 4 4
            - -|- -|- -|- -|
              4  3   3   4          
        */
        
        int[] expectedX = new int[]{0, 2, 4, 6};
        int[] expectedY = new int[]{4, 3, 3, 4};
        
        AverageUtil avgUtil = new AverageUtil();
        
        PairIntArray curve = new PairIntArray(8);
        curve.add(0, 2);
        curve.add(1, 3);
        curve.add(2, 4);
        curve.add(3, 3);
        curve.add(4, 2);
        curve.add(5, 4);
        curve.add(6, 4);
        curve.add(7, 4);
        
        int kPoints = 2;
        
        PairIntArray result = avgUtil.bin(curve, kPoints);
        
        assertTrue(result.getN() == expectedX.length);
        
        for (int i = 0; i < result.getN(); ++i) {
            int x = result.getX(i);
            int y = result.getY(i);
            assertTrue(x == expectedX[i]);
            assertTrue(y == expectedY[i]);
        }
    }
}
