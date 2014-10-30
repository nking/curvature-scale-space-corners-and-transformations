package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PairIntTest {
    
    public PairIntTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGetSetX() {
        int xPoint = 3;
        int yPoint = 100;
        PairInt instance = new PairInt();
        instance.setX(xPoint);
        instance.setY(yPoint);
        
        assertTrue(instance.getX() == xPoint);
        assertTrue(instance.getY() == yPoint);
    }
    
}
