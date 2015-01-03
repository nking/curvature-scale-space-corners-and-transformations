package algorithms.util;

import static junit.framework.Assert.assertTrue;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PairFloatTest {
    
    public PairFloatTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGetSetX() {
        float xPoint = 3;
        float yPoint = 100;
        
        PairFloat instance = new PairFloat();
        instance.setX(xPoint);
        instance.setY(yPoint);
        
        assertTrue(instance.getX() == xPoint);
        assertTrue(instance.getY() == yPoint);
    }
    
   
}
