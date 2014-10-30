package algorithms.imageProcessing;

import java.security.SecureRandom;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PairFloatArrayTest {
    
    public PairFloatArrayTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void test() {
        
        float xPoint = 1;
        float yPoint = 10;
        
        PairFloatArray instance = new PairFloatArray();
        
        assertTrue(instance.getN() == 0);
        
        instance.add(xPoint, yPoint);
        
        assertTrue(instance.getN() == 1);
        
        assertTrue(instance.getX(0) == xPoint);
        assertTrue(instance.getY(0) == yPoint);
        
        assertTrue(instance.getX()[0] == xPoint);
        assertTrue(instance.getY()[0] == yPoint);
        
        instance.set(0, 2*xPoint, 2*yPoint);
        assertTrue(instance.getX(0) == 2*xPoint);
        assertTrue(instance.getY(0) == 2*yPoint);
    }

    @Test
    public void test2() {
        
        PairFloatArray instance = new PairFloatArray();
        
        SecureRandom sr = new SecureRandom();
        sr.setSeed(System.currentTimeMillis());
        
        int nr = 51;
        float[] x = new float[nr];
        float[] y = new float[nr];
        for (int i = 0; i < nr; i++) {
            x[i] = sr.nextFloat();
            y[i] = sr.nextFloat();
            instance.add(x[i], y[i]);
        }
        
        for (int i = 0; i < x.length; i++) {
            assertTrue(instance.getX(i) == x[i]);
            assertTrue(instance.getY(i) == y[i]);
            assertTrue(instance.getX()[i] == x[i]);
            assertTrue(instance.getY()[i] == y[i]);
        }
        
        instance.set(nr/2, 100, 101);
        assertTrue(instance.getX(nr/2) == 100);
        assertTrue(instance.getY(nr/2) == 101);
        
        assertTrue(instance.toString().contains("x="));
    }

}
