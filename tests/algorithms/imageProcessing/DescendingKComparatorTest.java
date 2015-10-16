package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DescendingKComparatorTest extends TestCase {
    
    public DescendingKComparatorTest() {
    }
    
    public void test0() throws Exception {
        
        DescendingKComparator comp = new DescendingKComparator();
        
        CornerRegion cr0 = new CornerRegion(100, 3, 1);
        cr0.set(1, 0.5f, 10, 20);
        
        CornerRegion cr1 = new CornerRegion(50, 3, 1);
        cr1.set(1, 1.5f, 20, 40);
        
        CornerRegion cr2 = null;
        
        assertTrue(comp.compare(cr0, cr1) == 1);
        
        assertTrue(comp.compare(cr2, cr0) == 1);
                
        assertTrue(comp.compare(cr1, cr0) == -1);
        
        assertTrue(comp.compare(cr0, cr2) == -1);
        
        assertTrue(comp.compare(cr0, cr0) == 0);
        
        assertTrue(comp.compare(cr2, cr2) == 0);
        
        List<CornerRegion> list = new ArrayList<CornerRegion>();
        list.add(cr0);
        list.add(cr1);
        
        Collections.sort(list, comp);
        
        assertEquals(cr1, list.get(0));
    }
}
