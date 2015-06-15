package algorithms.misc;

import algorithms.misc.MedianInterpolation.SortedVector;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MedianInterpolationTest extends TestCase {
    
    public MedianInterpolationTest() {
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
    
    public void testInterpolation0() throws Exception {
        
        /*
        k=3
        curve=10 4's
        expectedN = 10 - kPoints + 1
        4 4 4 4 4 4 4 4 4 4
            4 4 4 4 4 4 4 4
            0 1 2 3 4 5 6 7
        */
        MedianInterpolation interp = new MedianInterpolation();
        
        int[] curveY = new int[]{4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
        int kPoints = 3;
        int[] result = interp.interpolate(curveY, kPoints);
                
        int[] expected = new int[]{4, 4, 4, 4, 4, 4, 4, 4};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
    }
    
    public void testInterpolation1() throws Exception {
        
        /*
        k=3
        curve=10 4's
        expectedN = 10 - kPoints + 1
        2 3 4 3 2 4 4 4 4 4
            3 3 3 3 4 4 4 4
            0 1 2 3 4 5 6 7
        */
        MedianInterpolation interp = new MedianInterpolation();
        
        int[] curveY = new int[]{2, 3, 4, 3, 2, 4, 4, 4, 4, 4};
        int kPoints = 3;
        int[] result = interp.interpolate(curveY, kPoints);
                
        int[] expected = new int[]{3, 3, 3, 3, 4, 4, 4, 4};
        
        assertTrue(result.length == expected.length);
        
        assertTrue(Arrays.equals(result, expected));
        
    }
}
