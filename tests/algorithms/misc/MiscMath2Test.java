package algorithms.misc;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MiscMath2Test extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void test() {
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            int nb = MiscMath.numberOfBits(v);
            
            //System.out.println("i=" + i + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            v *= -1;
            int nb = MiscMath.numberOfBits(v);
            
            //System.out.println("i=" + i + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    
    public void test1() {
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            int nb = MiscMath0.numberOfBitsWOB(v);
            
            //System.out.println("i=" + i + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
        for (int i = 0; i < 31; ++i) {
            int v = 1 << i;
            String bs = Integer.toBinaryString(v);
            v *= -1;
            int nb = MiscMath0.numberOfBitsWOB(v);
            
            //System.out.println("bs=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    
    public void test2() {
        for (int i = 1; i < 63; ++i) {
            long v = 1L << i;
            int nb = MiscMath.numberOfBits(v);
            String bs = Long.toBinaryString(v);
            //System.out.println("i=" + i + " bs=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }

    public void test3() {
        for (int i = 0; i < 31; ++i) {
            int nb = MiscMath.numberOfBits(i);
            String bs = Integer.toBinaryString(i);
            //System.out.println("v=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
    public void test4() {
        for (int i = 0; i < 31; ++i) {
            int nb = MiscMath.numberOfBits((long)i);
            String bs = Integer.toBinaryString(i);
            //System.out.println("v=" + bs + " nb=" + nb);
            assertEquals(bs.length(), nb);
        }
    }
}

