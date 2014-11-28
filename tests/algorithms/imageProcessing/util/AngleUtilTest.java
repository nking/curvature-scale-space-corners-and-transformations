package algorithms.imageProcessing.util;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class AngleUtilTest {
    
    public AngleUtilTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSubtract() {
        
        double diffX1 = 0;
        double diffY1 = 0;
        double diffX2 = 0;
        double diffY2 = 0;
        
        AngleUtil instance = new AngleUtil();
        
        double result = instance.subtract(diffX1, diffY1, diffX2, diffY2);
        
        
    }
    
}
