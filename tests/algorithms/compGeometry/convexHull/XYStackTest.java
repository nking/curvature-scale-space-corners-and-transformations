package algorithms.compGeometry.convexHull;

import java.security.SecureRandom;
import java.util.logging.Logger;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 */
public class XYStackTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public XYStackTest(String testName) {
        super(testName);
    }

    public void testLIFO() {

        SecureRandom sr = new SecureRandom();
        sr.setSeed(System.currentTimeMillis());

        int ntries = 1;

        for (int i = 0; i < ntries; i++) {

            int n = 1000;

            float[] x = createRandomNumbers(sr, n);
            float[] y = createRandomNumbers(sr, n);

            XYStack stack = new XYStack();

            for (int j = 0; j < n; j++) {
                stack.push(x[j], y[j]);
            }

            float xCheck0 = x[n-1];
            float yCheck0 = y[n-1];

            float xCheck1 = x[n-2];
            float yCheck1 = y[n-2];

            float xCheck2 = x[n-3];
            float yCheck2 = y[n-3];

            float xCheck3 = x[n-4];
            float yCheck3 = y[n-4];

            float[] result0 = stack.pop();
            float[] result1 = stack.pop();
            float[] result2 = stack.pop();
            float[] result3 = stack.pop();

            assertTrue( Math.abs(xCheck0 - result0[0]) < 0.01);
            assertTrue( Math.abs(yCheck0 - result0[1]) < 0.01);

            assertTrue( Math.abs(xCheck1 - result1[0]) < 0.01);
            assertTrue( Math.abs(yCheck1 - result1[1]) < 0.01);

            assertTrue( Math.abs(xCheck2 - result2[0]) < 0.01);
            assertTrue( Math.abs(yCheck2 - result2[1]) < 0.01);

            assertTrue( Math.abs(xCheck3 - result3[0]) < 0.01);
            assertTrue( Math.abs(yCheck3 - result3[1]) < 0.01);
        }
    }
    
    public void testEmptyStack() {
        XYStack stack = new XYStack();
        assertTrue(Float.isInfinite(stack.peekTopX()));
        assertTrue(Float.isInfinite(stack.peekTopY()));
        assertTrue(Float.isInfinite(stack.peekNextToTopX()));
        assertTrue(Float.isInfinite(stack.peekNextToTopY()));
        assertNull(stack.pop());
        assertTrue(stack.isEmpty());
        
        stack = new XYStack(100);
        for (int i = 0; i < 10; i++) {
            stack.push(i, i);
        }
        stack.compressArrays();
        assertTrue(stack.getNPoints() == 10);
        assertTrue(stack.x.length == 10);
    }

    public void testPush() {
        XYStack stack = new XYStack(10);
        for (int i = 0; i < 11; i++) {
            stack.push(i, i);
        }
        assertTrue(stack.getNPoints() == 11);
        assertTrue(stack.x.length > 10);
    }

    protected float[] createRandomNumbers(SecureRandom sr, int n) {
        float[] x = new float[n];
        for (int i = 0; i < n; i++) {
            x[i] = sr.nextFloat();
        }
        return x;
    }

    /**
     * Test suite
     * @return static Test
    */
    public static Test suite(){
        return new TestSuite(XYStackTest.class);
    }

    /**
     * Set up a Junit test runner
     * @param args Not used.
    */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }

}
