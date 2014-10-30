package algorithms.util;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class VeryLargeNumberTest {
    
    public VeryLargeNumberTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testConstructor() {
                
        // this one could require change w/ implementation
        
        VeryLargeNumber instance = new VeryLargeNumber(1294);
        
        int[] a = instance.getInternalArray();
   
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(a[0] == 1);
            assertTrue(a[1] == 2);
            assertTrue(a[2] == 9);
            assertTrue(a[3] == 4);
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(a[0] == 1294);
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        instance = new VeryLargeNumber(Integer.MAX_VALUE);
        a = instance.getInternalArray();
   
        if (VeryLargeNumber.BASE == 10) {
            //2147483647
            assertTrue(a[0] == 2);
            assertTrue(a[1] == 1);
            assertTrue(a[2] == 4);
            assertTrue(a[3] == 7);
            assertTrue(a[4] == 4);
            assertTrue(a[5] == 8);
            assertTrue(a[6] == 3);
            assertTrue(a[7] == 6);
            assertTrue(a[8] == 4);
            assertTrue(a[9] == 7);
            assertTrue(instance.getInternalArraySize() == 10);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(a[0] == 2);
            assertTrue(a[1] == 1);
            assertTrue(instance.getInternalArraySize() == 2);
        }
    }
    
    @Test
    public void testResetTo() {
                
        VeryLargeNumber instance = new VeryLargeNumber(1294);
        
        VeryLargeNumber instance2 = new VeryLargeNumber(-9876);
        
        instance.resetTo(instance2);
        
        assertTrue(instance.compareTo(instance2) == 0);
        
        assertTrue(instance.getInternalArraySize() == 
            instance2.getInternalArraySize());
        
        assertTrue(instance.isPositive() == instance2.isPositive());
        
        assertFalse(instance2.isPositive());
    }
    
    @Test
    public void testCompareTo() {
                
        VeryLargeNumber instance = new VeryLargeNumber(10);
        
        VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();
        
        int comp = instance.compareTo(maxLong);
        
        assertTrue(comp < 0);
        
        comp = maxLong.compareTo(instance);
        
        assertTrue(comp > 0);
        
        comp = maxLong.compareTo(maxLong);
        
        assertTrue(comp == 0);
        
        
        instance = new VeryLargeNumber(-10);
        
        VeryLargeNumber instance1 = new VeryLargeNumber(10);
        
        comp = instance.compareTo(instance1);
        
        assertTrue(comp < 0);
        
        comp = instance1.compareTo(instance);
        
        assertTrue(comp > 0);
        
        comp = instance.compareTo(instance);
        
        assertTrue(comp == 0);
        
        
        instance = new VeryLargeNumber(-100);
        
        instance1 = new VeryLargeNumber(10);
        
        comp = instance.compareTo(instance1);
        
        assertTrue(comp < 0);

        comp = instance1.compareTo(instance);
        
        assertTrue(comp > 0);
        
        comp = instance.compareTo(instance);
        
        assertTrue(comp == 0);
    }

    @Test
    public void testAdd() {
                
        //----------------------------------------------------------------
        // add large positive w/ smaller positive
        VeryLargeNumber instance = new VeryLargeNumber(1294);
        
        assertTrue(instance.isPositive());
        
        VeryLargeNumber addThis = new VeryLargeNumber(10);
        
        instance.add(addThis);
        
        VeryLargeNumber expected = new VeryLargeNumber(1304);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.isPositive());
             
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        //----------------------------------------------------------------
        // same test, but swap order
        instance = new VeryLargeNumber(10);
        
        assertTrue(instance.isPositive());
        
        addThis = new VeryLargeNumber(1294);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(1304);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertTrue(instance.isPositive());
        
        //----------------------------------------------------------------
        // add larger negative w/ smaller positive
        instance = new VeryLargeNumber(-1294);
        
        assertFalse(instance.isPositive());
        
        addThis = new VeryLargeNumber(10);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(-1284);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertFalse(instance.isPositive());
        
        //----------------------------------------------------------------
        // add larger positive w/ smaller negative
        instance = new VeryLargeNumber(1294);
        
        assertTrue(instance.isPositive());
        
        addThis = new VeryLargeNumber(-10);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(1284);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertTrue(instance.isPositive());
        
        //-------------------------------------------------------------
        // add a negative number to 0 with same nLen
        instance = new VeryLargeNumber(0);
        
        addThis = new VeryLargeNumber(-2);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(-2);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
        
        assertFalse(instance.isPositive());
        
        //-------------------------------------------------------------
        // small negative number gets added to a larger positive number, both with same nLen
        instance = new VeryLargeNumber(-2);
        
        addThis = new VeryLargeNumber(4);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(2);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
        
        assertTrue(instance.isPositive());
        
        //-------------------------------------------------------------
        // add to a positive number near carry over
        instance = new VeryLargeNumber(9999);
        
        addThis = new VeryLargeNumber(10);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(10009);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 5);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertTrue(instance.isPositive());
        
        //-------------------------------------------------------------
        // add to a negative number near carry over
        instance = new VeryLargeNumber(-9999);
        
        addThis = new VeryLargeNumber(-10);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(-10009);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 5);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertFalse(instance.isPositive());
        
        //-------------------------------------------------------------
        instance = new VeryLargeNumber(-9999);
        
        assertFalse(instance.isPositive());
        
        addThis = new VeryLargeNumber(10);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(-9989);
        
        assertFalse(instance.isPositive());
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        //-------------------------------------------------------------
        instance = new VeryLargeNumber(-100);
        
        assertFalse(instance.isPositive());
        
        addThis = new VeryLargeNumber(10);
        
        instance.add(addThis);
        
        expected = new VeryLargeNumber(-90);
        
        assertFalse(instance.isPositive());
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 2);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        //-------------------------------------------------------------
        instance = new VeryLargeNumber(-9999);
        
        addThis = new VeryLargeNumber(-9999);
        
        instance.add(addThis);
        
        assertFalse(instance.isPositive());
        
        expected = new VeryLargeNumber(-19998);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 5);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertFalse(instance.isPositive());
        
        //-------------------------------------
        // test that adding to an int until the number is larger than a long
        // will hold expected results
        // max long = (2^63)-1
        VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();
        
        instance = new VeryLargeNumber(Integer.MAX_VALUE);
        
        instance.add(maxLong);
        
        assertTrue(instance.isPositive());
        
        //922337203 9002259454
        expected = new VeryLargeNumber(0);
        if (VeryLargeNumber.BASE == 10) {
            expected.setInternalArray(
                new int[]{
                9, 2, 2, 3, 3, 7, 2, 0, 3,   9, 0, 0, 2, 2, 5, 9, 4, 5, 4
            }, 19, true);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            //9 223372039 002259454
            expected.setInternalArray(
                new int[]{8, 18, 8}, 3, true);
        }
        
        assertTrue(instance.compareTo(expected) == 0);
        
        addThis = new VeryLargeNumber(Integer.MIN_VALUE);
        
        assertFalse(addThis.isPositive());
        
        instance.add(addThis);
        
        int comp = instance.compareTo(expected);
        
        assertTrue(comp == -1);
        
        assertTrue(expected.compareTo(instance) == 1);
        
        //9223372036854775806
        expected = new VeryLargeNumber(0);
        
        if (VeryLargeNumber.BASE == 10) {
            expected.setInternalArray(
                new int[]{
                9, 2, 2, 3, 3, 7, 2, 0, 3,   6, 8, 5, 4, 7, 7, 5, 8, 0, 6
            }, 19, true);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            //9223372036854775806
            expected.setInternalArray(new int[]{8, 16, 6}, 3, true);
        }
        
        assertTrue(instance.compareTo(expected) == 0);
    }

    @Test
    public void testSubtract() {
                
        VeryLargeNumber instance = new VeryLargeNumber(1294);
        
        VeryLargeNumber subtrThis = new VeryLargeNumber(10);
        
        instance.subtract(subtrThis);
        
        VeryLargeNumber expected = new VeryLargeNumber(1284);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        //-------------------------------------------------------------
        VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();
        
        instance = new VeryLargeNumber(Integer.MAX_VALUE);
        
        String str0 = instance.toString();        
        String strA = maxLong.toString();
        
        assertEquals("9223372036854775807", strA);
        assertEquals("2147483647", str0);
        
        //((1<<63)-1) + ((1<<31)-1) = 9223372039002259454L
        instance.add(maxLong);
        
        // overflows, so we get '9223372019674906632 + 9223372034707292178'
        str0 = instance.toString();
        
        // 9223372039002259454L - 10 = 9223372039002259444
        instance.subtract(subtrThis);
        expected = new VeryLargeNumber(0);
        
        //922337203 9002259444
        if (VeryLargeNumber.BASE == 10) {
            expected.setInternalArray(
                new int[]{
                9, 2, 2, 3, 3, 7, 2, 0, 3,   9, 0, 0, 2, 2, 5, 9, 4, 4, 4
            }, 19, true);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            //9 223372039 002259444
            expected.setInternalArray(new int[]{8, 17, 1073741821}, 3, true);
        }
        String str1 = instance.toString();
        String str2 = expected.toString();
        assertTrue(instance.compareTo(expected) == 0);
        
        //-------------------------------------------------------------
        //1004-10
        instance = new VeryLargeNumber(1004);
        
        subtrThis = new VeryLargeNumber(10);
        
        instance.subtract(subtrThis);
        
        expected = new VeryLargeNumber(0);
        
        if (VeryLargeNumber.BASE == 10) {
            expected.setInternalArray(new int[]{9, 9, 4}, 3, true);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            expected.setInternalArray(new int[]{994}, 1, true);
        }
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 3);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        assertTrue(instance.compareTo(expected) == 0);        
        
        //-------------------------------------------------------------
        //4 - 10
        instance = new VeryLargeNumber(4);
        
        subtrThis = new VeryLargeNumber(10);
        
        instance.subtract(subtrThis);
        
        assertFalse(instance.isPositive());
        
        expected = new VeryLargeNumber(0);
        
        // -6
        expected.setInternalArray(new int[]{6}, 1, false);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
        
        //-------------------------------------------------------------
        //10 - 10
        instance = new VeryLargeNumber(10);
        
        subtrThis = new VeryLargeNumber(10);
        
        instance.subtract(subtrThis);
        
        assertTrue(instance.isPositive());
        
        expected = new VeryLargeNumber(0);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
    }
    
    @Test
    public void testIncrement() {
                
        VeryLargeNumber instance = new VeryLargeNumber(1294);
        
        instance.increment();
        
        VeryLargeNumber expected = new VeryLargeNumber(1295);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        //-------------------------------------------------------------
        instance = new VeryLargeNumber(-1294);
        
        instance.increment();
        
        expected = new VeryLargeNumber(-1293);
        
        assertTrue(instance.compareTo(expected) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        //------------------------------------------------------------
        instance = new VeryLargeNumber(-1);
        
        instance.increment();
        
        expected = new VeryLargeNumber(0);
                
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
        
        assertTrue(instance.isPositive());
        
        instance.increment();
        
        expected = new VeryLargeNumber(1);
                
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
        
        //------------------------------------------------------------
        instance = new VeryLargeNumber(0);
        
        instance.increment();
        
        expected = new VeryLargeNumber(1);
                
        assertTrue(instance.compareTo(expected) == 0);
        
        assertTrue(instance.getInternalArraySize() == 1);
        
        //------------------------------------------------------------
        instance = new VeryLargeNumber(0);
        for (int i = 0; i < 110; i++) {
            instance.increment();
        }
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 3);
        } else if (VeryLargeNumber.BASE == 1073741823) {
            //1073741823
            assertTrue(instance.getInternalArraySize() == 1);
        }
        
        expected = new VeryLargeNumber(110);
        assertTrue(instance.compareTo(expected) == 0);
        
        //------------------------------------------------------------
        instance = new VeryLargeNumber(-110);
        for (int i = 0; i < 110; i++) {
            instance.increment();
        }
        assertTrue(instance.isZero());
        expected = new VeryLargeNumber(0);
        assertTrue(instance.compareTo(expected) == 0);
        assertTrue(instance.getInternalArraySize() == 1);
        
        //------------------------------------------------------------
        instance = new VeryLargeNumber(-10);
        for (int i = 0; i < 15; i++) {
            instance.increment();
        }
        expected = new VeryLargeNumber(5);
        assertTrue(instance.compareTo(expected) == 0);
        assertTrue(instance.getInternalArraySize() == 1);
    }
    
    @Test
    public void testClone() throws Exception {

        VeryLargeNumber instance = new VeryLargeNumber(1294);
        
        VeryLargeNumber clone = instance.clone();
                
        assertTrue(instance.compareTo(clone) == 0);
        
        if (VeryLargeNumber.BASE == 10) {
            assertTrue(instance.getInternalArraySize() == 4);
        } else if (VeryLargeNumber.BASE == ((1<<30)-1)) {
            assertTrue(instance.getInternalArraySize() == 1);
        }
        assertTrue(instance.getInternalArraySize() == 
            clone.getInternalArraySize());
        
        //-------------------------------------------------------------
        instance = new VeryLargeNumber(-1294);
        
        clone = instance.clone();
                
        assertTrue(instance.compareTo(clone) == 0);
        
        assertTrue(instance.getInternalArraySize() == 
            clone.getInternalArraySize());
    }
    
    @Test
    public void testDivideByAndPrint() {
        
        VeryLargeNumber divisor = new VeryLargeNumber(10);        
        VeryLargeNumber instance = new VeryLargeNumber(1294);
        String str = instance.divideByAndPrint(divisor);
        //System.out.println("1294/10=" + str);
        assertEquals("129.4", str);
        
        //---------------------------------------------------------
        // test for divisor is negative
        instance = new VeryLargeNumber(1294);
        divisor = new VeryLargeNumber(-10);
        str = instance.divideByAndPrint(divisor);        
        //System.out.println("1294/(-10)=" + str);
        assertEquals("-129.4", str);
        
        //---------------------------------------------------------
        // test for instance is negative
        instance = new VeryLargeNumber(-1294);
        divisor = new VeryLargeNumber(10);
        str = instance.divideByAndPrint(divisor);
        //System.out.println("-1294/10=" + str);
        assertEquals("-129.4", str);
        
        //---------------------------------------------------------
        // test for instance and divisor are negative
        instance = new VeryLargeNumber(-1294);
        divisor = new VeryLargeNumber(-10);
        str = instance.divideByAndPrint(divisor);        
        //System.out.println("-1294/(-10)=" + str);
        assertEquals("129.4", str);
        
        //---------------------------------------------------------
        // test for expected remainder is zero
        instance = new VeryLargeNumber(100);
        divisor = new VeryLargeNumber(10);
        str = instance.divideByAndPrint(divisor);
        //System.out.println("100/10=" + str);
        assertTrue(str.startsWith("10.0"));
        
        //---------------------------------------------------------
        instance = new VeryLargeNumber(1294);
        divisor = new VeryLargeNumber(7);
        str = instance.divideByAndPrint(divisor);
        assertEquals("184.8571428571428571", str);
        //                 123 1234567890123456
        
        //System.out.println("1294/7=" + str);
     
        //---------------------------------------------------------
        // very large number test
        //9223372036854775807
        /*VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();
        
        instance = new VeryLargeNumber(Integer.MAX_VALUE);
        instance.add(maxLong);
        // ===> 9223372039002259454
        
        divisor = new VeryLargeNumber(Integer.MAX_VALUE);
        divisor.add(new VeryLargeNumber(Integer.MAX_VALUE));
        // ==> 4294967294
        
        // 9223372039002259454 / 4294967294
        str = instance.divideByAndPrint(divisor);
        
        System.out.println("9223372039002259454 / 4294967294=" + str);
        
        assertTrue(str.startsWith("2147483649.5"));
        */
    }

}
