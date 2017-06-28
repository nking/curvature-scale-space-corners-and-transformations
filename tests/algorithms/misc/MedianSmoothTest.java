package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.misc.MedianSmooth.SortedVector;
import algorithms.util.PairIntArray;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MedianSmoothTest extends TestCase {
    
    public MedianSmoothTest() {
    }
    
    public void testSortedVector() throws Exception {
        
        int k = 5;
        
        SortedVector sVec = new SortedVector(k);
        
        assertTrue(sVec.a.length == k);
        
        assertTrue(sVec.n == 0);
        
        sVec.append(5);
        sVec.append(4);
        assertTrue(sVec.n == 2);
        assertTrue(sVec.a[0] == 5);
        assertTrue(sVec.a[1] == 4);
        sVec.append(3);
        sVec.append(2);
        sVec.append(1);
        assertTrue(sVec.n == 5);
        
        boolean caughtException = false;
        try {
            sVec.append(6);
        } catch(Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        // assert the side effect of ascending sort is true
        for (int i = 0; i < 5; i++) {
            assertTrue(sVec.a[i] == (i + 1));
        }
        
        //[1,2,3,4,5]
        sVec.remove(5);
        assertTrue(sVec.a[4] == Integer.MAX_VALUE);
        assertTrue(sVec.n == 4);
        assertTrue(sVec.availSlot == 4);
        
        caughtException = false;
        try {
            sVec.remove(4);
        } catch(Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        caughtException = false;
        try {
            sVec.remove(6);
        } catch(Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
      
        caughtException = false;
        try {
            sVec.getMedian();
        } catch(Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
          
        //[1,2,3,4,inf]
        sVec.insertIntoOpenSlot(5);
        assertTrue(sVec.a[4] == 5);
        assertTrue(sVec.n == 5);
        assertTrue(sVec.getMedian() == 3);
        
        //[1,2,3,4,5]
        sVec.remove(5);
        assertTrue(sVec.a[4] == Integer.MAX_VALUE);
        assertTrue(sVec.n == 4);
        assertTrue(sVec.availSlot == 4);
        
        //[1,2,3,4,inf]
        sVec.insertIntoOpenSlot(0);
        assertTrue(sVec.a[4] == 4);
        assertTrue(sVec.n == 5);
        for (int i = 0; i < 5; i++) {
            assertTrue(sVec.a[i] == i);
        }
        assertTrue(sVec.getMedian() == 2);
        
        //[0,1,2,3,4]
        sVec.remove(0);
        assertTrue(sVec.a[0] != 0);
        assertTrue(sVec.n == 4);
        assertTrue(sVec.availSlot == 0);
        
        //[#,1,2,3,4]
        sVec.insertIntoOpenSlot(2);
        assertTrue(sVec.a[4] == 4);
        assertTrue(sVec.n == 5);
        //[1,2,2,3,4]
        assertTrue(sVec.a[0] == 1);
        assertTrue(sVec.a[1] == 2);
        assertTrue(sVec.a[2] == 2);
        assertTrue(sVec.a[3] == 3);
        assertTrue(sVec.a[4] == 4);
        assertTrue(sVec.getMedian() == 2);
        
        //[1,2,2,3,4]
        sVec.remove(2);
        assertTrue(sVec.n == 4);
        assertTrue(sVec.availSlot == 1 || sVec.availSlot == 2);
        //[1,2,#,3,4]
        sVec.insertIntoOpenSlot(0);
        assertTrue(sVec.a[4] == 4);
        assertTrue(sVec.n == 5);
        //[0,1,2,3,4]
        for (int i = 0; i < 5; i++) {
            assertTrue(sVec.a[i] == i);
        }
        assertTrue(sVec.getMedian() == 2);
        
        
        sVec = new SortedVector(4);
        for (int i = 3; i > -1; i--) {
            sVec.append(i);
        }
        //[0,1,2,3]
        assertTrue(sVec.a[0] == 0);
        assertTrue(sVec.a[1] == 1);
        assertTrue(sVec.a[2] == 2);
        assertTrue(sVec.a[3] == 3);
        assertTrue(sVec.getMedian() == 1);
        
        sVec = new SortedVector(3);
        sVec.append(3);
        sVec.append(5);
        sVec.insertIntoOpenSlot(4);
        assertTrue(sVec.n == 3);
        assertTrue(sVec.a[0] == 3);
        assertTrue(sVec.a[1] == 4);
        assertTrue(sVec.a[2] == 5);
        assertTrue(sVec.getMedian() == 4);
        //[3,4,5]
        sVec.remove(3);
        //[#,4,5]  insert 6
        sVec.insertIntoOpenSlot(6);
        assertTrue(sVec.n == 3);
        assertTrue(sVec.a[0] == 4);
        assertTrue(sVec.a[1] == 5);
        assertTrue(sVec.a[2] == 6);
        assertTrue(sVec.getMedian() == 5);
        
        
        //[1, 2, 3, 404, 459]
        sVec = new SortedVector(5);
        sVec.append(1);
        sVec.append(2);
        sVec.append(3);
        sVec.append(404);
        sVec.append(459);
        //[1, #, 3, 404, 459]  insert 407
        sVec.remove(2);
        sVec.insertIntoOpenSlot(407);
        assertTrue(sVec.a[0] == 1);
        assertTrue(sVec.a[1] == 3);
        assertTrue(sVec.a[2] == 404);
        assertTrue(sVec.a[3] == 407);
        assertTrue(sVec.a[4] == 459);
        
    }
    
    public void testCalculate2() throws Exception {
        
        /*
        x:  0 1 2 3 4 5 6 7
        k=5
        averages for x:        
            0 1 2 3  4  5  6  7    
        sum         10 15 20 25
        avg          2  3  4  5
        
        y:  2 3 4 3  2  4  4  4
        y:           3  3  4  4        
        */
        
        int[] expectedX = new int[]{2, 3, 4, 5};
        int[] expectedY = new int[]{3, 3, 4, 4};
        
        MedianSmooth interp = new MedianSmooth();
        
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
        
        PairIntArray result = interp.calculate(curve, kPoints);
        
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
        
        MedianSmooth interp = new MedianSmooth();
        
        PairIntArray curve = null;
        
        int kPoints = 5;
        
        boolean caughtException = false;
        try {
            PairIntArray result = interp.calculate(curve, kPoints);
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
            PairIntArray result = interp.calculate(curve, kPoints);
        } catch (Throwable t) {
            caughtException = true;
        }
        assertTrue(caughtException);
        
        
    }
    
    public void testMedianFilter2D() throws Exception {
        
        MedianSmooth medianSmooth = new MedianSmooth();
        
        /*
          4
          3
          2 7 8 9       2   8 9 10
          1 4 5 6       1   5 6 7 
          0 1 2 3       0   2 3 4
            0 1 2 3 4     0 1 2 3 
        
          3 10 11 12          3   11 12  13  
          2 7  8  9           2    8  9  10
          1 4  5  6           1    5  6  7 
          0                   0        
            0  1  2  3  4       0  1  2  3 
        */
        
        GreyscaleImage input = new GreyscaleImage(4, 4);
        input.setValue(0, 0, 1);
        input.setValue(1, 0, 5);
        input.setValue(2, 0, 3);
        input.setValue(0, 1, 4);
        input.setValue(1, 1, 2);
        input.setValue(2, 1, 6);
        input.setValue(0, 2, 7);
        input.setValue(1, 2, 8);
        input.setValue(2, 2, 9);
        
        input.setValue(3, 0, 4);
        input.setValue(3, 1, 7);
        input.setValue(3, 2, 10);
        
        input.setValue(0, 3, 10);
        input.setValue(1, 3, 11);
        input.setValue(2, 3, 12);
        
        input.setValue(3, 3, 13);
        
        GreyscaleImage out = medianSmooth.calculate(input, 3, 3);
                
        assertEquals(5, out.getValue(1, 1));
        assertEquals(6, out.getValue(2, 1));
        assertEquals(8, out.getValue(1, 2));
        assertEquals(9, out.getValue(2, 2));
    }
}
