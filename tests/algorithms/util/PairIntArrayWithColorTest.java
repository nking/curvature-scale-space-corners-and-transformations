package algorithms.util;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PairIntArrayWithColorTest {
    
    public PairIntArrayWithColorTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGetColor() {
        
        PairIntArray xy = new PairIntArray();
        xy.add(10, 20);
        
        PairIntArrayWithColor xyColor = new PairIntArrayWithColor(xy);
        xyColor.setColor(1);
        
        assertTrue(xyColor.getColor() == 1);
        
        assertTrue(xyColor.getN() == 1);
        
        assertTrue((xyColor.getX(0) == 10) && (xyColor.getY(0) == 20));
    }

}
