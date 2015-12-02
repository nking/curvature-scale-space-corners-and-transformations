package algorithms.imageProcessing;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FixedSizeSortedIntVectorTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testAdd() throws Exception {
        
        FixedSizeSortedIntVector sortedVector = 
            new FixedSizeSortedIntVector(4);
        
        assertTrue(sortedVector.add(7));
        assertTrue(sortedVector.add(6));
        assertTrue(sortedVector.add(5));
        assertTrue(sortedVector.add(4));
        
        int[] values = sortedVector.getArray();
        
        assertNotNull(values);
        
        assertEquals(4, values[0]);
        assertEquals(5, values[1]);
        assertEquals(6, values[2]);
        assertEquals(7, values[3]);
        
        assertFalse(sortedVector.add(10));

        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertEquals(4, values[0]);
        assertEquals(5, values[1]);
        assertEquals(6, values[2]);
        assertEquals(7, values[3]);
        
        assertTrue(sortedVector.add(3));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertEquals(3, values[0]);
        assertEquals(4, values[1]);
        assertEquals(5, values[2]);
        assertEquals(6, values[3]);
        
        //----
        sortedVector = new FixedSizeSortedIntVector(4);
        
        sortedVector.add(6);
        sortedVector.add(4);
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        
        assertEquals(4, values[0]);
        assertEquals(6, values[1]);
        
        //test whether can have an array with nulls in it and still use
        //Arrays.binarySearch for non-null indexes
   
        int[] a = new int[4];
        a[0] = 6;
        a[1] = 7;
        int idx = Arrays.binarySearch(a, 0, 2, 6);
        assertTrue(idx == 0);
        
        sortedVector = new FixedSizeSortedIntVector(4);
        
        sortedVector.add(7);
        sortedVector.add(6);
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        assertEquals(2, sortedVector.getNumberOfItems());
        
        assertEquals(6, values[0]);
        assertEquals(7, values[1]);
        
        sortedVector.add(4);
        sortedVector.add(5);
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        assertEquals(4, sortedVector.getNumberOfItems());
        
        assertEquals(4, values[0]);
        assertEquals(5, values[1]);
        assertEquals(6, values[2]);
        assertEquals(7, values[3]);
        
        sortedVector.add(10);
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        assertEquals(4, sortedVector.getNumberOfItems());
        
        assertEquals(4, values[0]);
        assertEquals(5, values[1]);
        assertEquals(6, values[2]);
        assertEquals(7, values[3]);
    }
    
    public void testRandom() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed=1393807938003554000l;
        sr.setSeed( seed );
        
        int n = 1000;
        int sz = 100;
                
        List<Integer> sorted = new ArrayList<Integer>();
        
        FixedSizeSortedIntVector sVec = new FixedSizeSortedIntVector(sz);
                
        for (int i = 0; i < n; ++i) {
            
            int number = sr.nextInt(sz);
             
            sorted.add(Integer.valueOf(number));
            
            sVec.add(number);
            
            if (i >= sz) {
                
                Collections.sort(sorted);
                sorted.remove(sorted.size() - 1);
                
                // compare contents
                for (int j = 0; j < sorted.size(); ++j) {
                    Integer number2 = sorted.get(j);
                    assertEquals(number2.intValue(), sVec.getValue(j));
                }
            }
        }
    }
}
