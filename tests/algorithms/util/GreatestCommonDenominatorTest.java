package algorithms.util;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GreatestCommonDenominatorTest extends TestCase {

    public GreatestCommonDenominatorTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
    
    public void testEuclid0() {
        System.out.println("testEuclid0");
        int a = 11;
        int nIter = 0;
        for (int b = (a-1); b > 0; b--) {
            GreatestCommonDenominator.count = 0;
            int result = GreatestCommonDenominator.euclid(a, b);
            System.out.println("count=" + GreatestCommonDenominator.count + " result=" + result);
            nIter += GreatestCommonDenominator.count;
            if (result == 1) {
                // a is possibly prime
            } else {
                // a is not prime
            }
        }
        System.out.println("nIter=" + nIter + " to prove that 11 is prime");
    }

    /**
     * Test of euclid method, of class GreatestCommonDenominator.
     */
    public void testEuclid() {
        System.out.println("testEuclid");
        int a = 30;
        int b = 21;
        int expResult = 3;
        int result = GreatestCommonDenominator.euclid(a, b);
        assertTrue(expResult == result);
    }

    public void testExtendedEuclid() {
        System.out.println("testExtendedEuclid");
        int a = 99;
        int b = 78;
        int[] expResult = new int[]{3, -11, 14};
        int[] result = GreatestCommonDenominator.extendedEuclid(a, b);
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(78, 21);
        expResult = new int[]{3, 3, -11};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(21, 15);
        expResult = new int[]{3, -2, 3};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(15, 6);
        expResult = new int[]{3, 1, -2};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(6, 3);
        expResult = new int[]{3, 0, 1};
        assertTrue( Arrays.equals(expResult, result));
        
        result = GreatestCommonDenominator.extendedEuclid(3, 0);
        expResult = new int[]{3, 1, 0};
        assertTrue( Arrays.equals(expResult, result));
    }

}
