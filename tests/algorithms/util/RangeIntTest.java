package algorithms.util;

import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class RangeIntTest extends TestCase {
    
    public RangeIntTest() {
    }

    public void test0() {
        RangeInt range = new RangeInt(1, 2);
        assertTrue(range.getStart() == 1);
        assertTrue(range.getStop() == 2);
        
        RangeInt range0 = new RangeInt(range);
        assertTrue(range0.getStart() == 1);
        assertTrue(range0.getStop() == 2);
        assertTrue(range.getStart() == 1);
        assertTrue(range.getStop() == 2);
    
        range.resetBoundsIfNeeded(1, 2);
        assertTrue(range.getStart() == 1);
        assertTrue(range.getStop() == 2);
        
        range.resetBoundsIfNeeded(2, 2);
        assertTrue(range.getStart() == 2);
        assertTrue(range.getStop() == 2);
    }
    
}
