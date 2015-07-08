package algorithms.imageProcessing;

import junit.framework.TestCase;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class FixedSizeSortedVectorTest extends TestCase {
    
    public void testAdd() throws Exception {
        
        FixedSizeSortedVector<Integer> sortedVector = 
            new FixedSizeSortedVector(4);
        
        sortedVector.add(Integer.valueOf(7));
        sortedVector.add(Integer.valueOf(6));
        sortedVector.add(Integer.valueOf(5));
        sortedVector.add(Integer.valueOf(4));
        
        Integer[] values = sortedVector.getArray();
        
        assertNotNull(values);
        
        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        sortedVector.add(Integer.valueOf(10));

        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        sortedVector.add(Integer.valueOf(3));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertTrue(values[0].intValue() == 3);
        assertTrue(values[1].intValue() == 4);
        assertTrue(values[2].intValue() == 5);
        assertTrue(values[3].intValue() == 6);
    }
    
}
