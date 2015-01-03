package algorithms.util;

import java.security.SecureRandom;
import static junit.framework.Assert.assertTrue;
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
        
        PairFloatArray copied = instance.copy();
        assertTrue(copied.getN() == instance.getN());
        for (int i = 0; i < x.length; i++) {
            assertTrue(copied.getX(i) == x[i]);
            assertTrue(copied.getY(i) == y[i]);
            assertTrue(copied.getX()[i] == x[i]);
            assertTrue(copied.getY()[i] == y[i]);
        }
        
        instance.set(nr/2, 100, 101);
        assertTrue(instance.getX(nr/2) == 100);
        assertTrue(instance.getY(nr/2) == 101);
        
        assertTrue(instance.toString().contains("x="));
        
    }

    public void testRemoveRange() throws Exception {
        
        PairFloatArray xy = new PairFloatArray();
        for (int i = 0; i < 7; i++) {
            xy.add(i, i);
        }
        
        xy.removeRange(0, 2);
        
        assertTrue(xy.getN() == 4);
        
        assertTrue(xy.getX(0) == 3);
        assertTrue(xy.getY(0) == 3);
        assertTrue(xy.getX(1) == 4);
        assertTrue(xy.getY(1) == 4);
        assertTrue(xy.getX(2) == 5);
        assertTrue(xy.getY(2) == 5);
        assertTrue(xy.getX(3) == 6);
        assertTrue(xy.getY(3) == 6);
        
        //=======
        xy = new PairFloatArray();
        for (int i = 0; i < 7; i++) {
            xy.add(i, i);
        }
        
        xy.removeRange(4, 6);
        
        assertTrue(xy.getN() == 4);
        
        assertTrue(xy.getX(0) == 0);
        assertTrue(xy.getY(0) == 0);
        assertTrue(xy.getX(1) == 1);
        assertTrue(xy.getY(1) == 1);
        assertTrue(xy.getX(2) == 2);
        assertTrue(xy.getY(2) == 2);
        assertTrue(xy.getX(3) == 3);
        assertTrue(xy.getY(3) == 3);
        
        xy = new PairFloatArray();
        xy.add(0, 10);
        xy.add(12, 23);
        xy.add(23, 31);
        xy.add(34, 84);
        int i = 2;
        int r1 = 31;
        xy.set(i - 1, xy.getX(i - 1), r1);
        xy.removeRange(i, i);
        assertTrue(xy.getN() == 3);
        assertTrue((xy.getX(0) == 0) && (xy.getY(0) == 10));
        assertTrue((xy.getX(1) == 12) && (xy.getY(1) == 31));
        assertTrue((xy.getX(2) == 34) && (xy.getY(2) == 84));

    }
    
    public void testToPairIntArray() throws Exception {
        
        PairFloatArray xy = new PairFloatArray();
        for (int i = 0; i < 7; i++) {
            xy.add(i, i);
        }
        
        PairIntArray xyF = xy.toPairIntArray();
        for (int i = 0; i < 7; i++) {
            assertTrue(Math.abs(xy.getX(i) - xyF.getX(i)) < 0.001f);
            assertTrue(Math.abs(xy.getY(i) - xyF.getY(i)) < 0.001f);
        }
    }
}
