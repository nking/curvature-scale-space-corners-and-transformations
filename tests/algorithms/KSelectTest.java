package algorithms;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class KSelectTest {
    
    public KSelectTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testSelect() {
        
        KSelect kSelect = new KSelect();
        
                            //  0  1  2  3  4  5  6  7  8  9
                            //              |  |
        float[] a = new float[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};

        float median = kSelect.findMedianOfMedians(a, 0, a.length - 1);
        assertTrue(median == 4);

        a = new float[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};
        median = kSelect.selectKth(a, 0, a.length - 1, 4);
        assertTrue(median == 4);

        // select kth smallest integer
        float kth = -1;
        for (int i = 0; i < a.length; i++) {
            a = new float[]{4, 5, 9, 0, 3, 1, 6, 2, 7, 8};
            kth = kSelect.selectKth(a, 0, a.length - 1, i);
            assertTrue(kth == i);
        }
    }
}
