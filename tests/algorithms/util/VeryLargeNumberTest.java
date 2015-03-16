package algorithms.util;

import java.util.Arrays;
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
        
        assertFalse(instance.isOdd());
        
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
        
        assertTrue(instance.isOdd());
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
        
        assertFalse(instance.isOdd());
        assertFalse(instance2.isOdd());
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
        
        if (VeryLargeNumber.BASE == 1073741823) {
            
            instance = new VeryLargeNumber(1073741823);
            
            assertTrue(instance.getInternalArraySize() == 2);
            
            addThis = new VeryLargeNumber(1073741823);
            
            instance.add(addThis);
            
            assertTrue(instance.getInternalArraySize() == 2);
            
            expected = new VeryLargeNumber(2147483646);
            
            assertTrue(expected.compareTo(instance) == 0);
        }         
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
        
        assertTrue(instance.isOdd());
        
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
    
    //@Test
    public void testInverseByEuclidean() {
        
        VeryLargeNumber instance = new VeryLargeNumber(4);
        
        VeryLargeNumber divisor = new VeryLargeNumber(4);        
        instance = new VeryLargeNumber(1);
        String str = instance.divideByAndPrint(divisor);
        
        divisor = new VeryLargeNumber(-4);        
        instance = new VeryLargeNumber(1);
        str = instance.divideByAndPrint(divisor);
        
        instance = new VeryLargeNumber(4);
        double inverted = instance.inverseByEuclidean(instance);
        assertTrue(Math.abs(inverted - 0.25) < 0.01);
        
        instance = new VeryLargeNumber(-4);
        inverted = instance.inverseByEuclidean(instance);
        assertTrue(Math.abs(inverted - -0.25) < 0.01);        
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
    }
    
    @Test
    public void testSplitAt() {
        
        VeryLargeNumber number;
        VeryLargeNumber[] split;
        
        int m = 3;
        
        number = new VeryLargeNumber(0);
        number.setInternalArray(new int[]{1, 2, 3, 4, 5}, 5, true);
        split = number.splitAt(m);
        assertTrue(split[0].getInternalArraySize() == 2);
        assertTrue(Arrays.equals(split[0].getInternalArray(), new int[]{1, 2}));
        assertTrue(split[1].getInternalArraySize() == 3);
        assertTrue(Arrays.equals(split[1].getInternalArray(), new int[]{3, 4, 5}));
        
        number = new VeryLargeNumber(0);
        number.setInternalArray(new int[]{6, 7, 8, 9}, 4, true);
        split = number.splitAt(m);
        assertTrue(split[0].getInternalArraySize() == 1);
        assertTrue(Arrays.equals(split[0].getInternalArray(), new int[]{6}));
        assertTrue(split[1].getInternalArraySize() == 3);
        assertTrue(Arrays.equals(split[1].getInternalArray(), new int[]{7, 8, 9}));
        
    }
    
    @Test
    public void testKaratasuba() throws CloneNotSupportedException {
        
        VeryLargeNumber num1 = new VeryLargeNumber(0);
        VeryLargeNumber num2 = new VeryLargeNumber(0);
        VeryLargeNumber result = VeryLargeNumber.karatsuba(num1, num2);
        
        VeryLargeNumber expected = new VeryLargeNumber(0);
        assertTrue(result.compareTo(expected) == 0);
        
        num1 = new VeryLargeNumber(12345);
        num2 = new VeryLargeNumber(6789);
        result = VeryLargeNumber.karatsuba(num1, num2);
        expected = new VeryLargeNumber(0);
        //83810205
        expected.setInternalArray(new int[]{83810205}, 1, true);        
        assertTrue(result.compareTo(expected) == 0); 
        
        //-----------------------------------------------
        num1 = new VeryLargeNumber(0);
        num1.setInternalArray(new int[]{1, 1}, 2, true);
        num2 = new VeryLargeNumber(0);
        num2.setInternalArray(new int[]{1, 1}, 2, true);
        result = VeryLargeNumber.karatsuba(num1, num2);
        expected = new VeryLargeNumber(0);
        //1152921504606846976
        expected.setInternalArray(new int[]{1, 2, 1}, 3, true);        
        assertTrue(result.compareTo(expected) == 0); 
        
        //------------------------------------------------
               
        num1 = new VeryLargeNumber(Integer.MAX_VALUE);
        num2 = new VeryLargeNumber(Integer.MAX_VALUE);
        result = VeryLargeNumber.karatsuba(num1, num2);
        expected = new VeryLargeNumber(0);
        //4611686014132420609
        expected.setInternalArray(new int[]{4, 4, 1}, 3, true);
        assertTrue(result.compareTo(expected) == 0);
        
        num1 = new VeryLargeNumber(Integer.MIN_VALUE);
        num2 = new VeryLargeNumber(Integer.MAX_VALUE);
        result = VeryLargeNumber.karatsuba(num1, num2);
        expected = new VeryLargeNumber(0);
        //-4611686016279904256
        expected.setInternalArray(new int[]{4, 6, 2}, 3, false);
        assertTrue(result.compareTo(expected) == 0);
        
        
        num1 = VeryLargeNumber.createMaxLong();
        num2 = new VeryLargeNumber(0);
        num2.setInternalArray(new int[]{1, 1}, 2, true);
        result = VeryLargeNumber.karatsuba(num1, num2);
        assertNotNull(result);
        expected = new VeryLargeNumber(0);
        //9903520314283042198119251968L
        expected.setInternalArray(new int[]{8, 24, 23, 7}, 4, true);
        assertTrue(result.compareTo(expected) == 0);
          
        num1 = VeryLargeNumber.createMaxLong();
        num2 = new VeryLargeNumber(Integer.MAX_VALUE);
        result = VeryLargeNumber.karatsuba(num1, num2);
        assertNotNull(result);
        expected = new VeryLargeNumber(0);
        //19807040619342712359383728129L
        expected.setInternalArray(new int[]{16, 40, 30, 7}, 4, true);
        assertTrue(result.compareTo(expected) == 0);
        
        
        num1 = VeryLargeNumber.createMaxLong();
        num2 = VeryLargeNumber.createMaxLong();
        result = VeryLargeNumber.karatsuba(num1, num2);
        assertNotNull(result);
        //85070591730234615847396907784232501249L
        expected = new VeryLargeNumber(0);      
        expected.setInternalArray(new int[]{64, 256, 368, 224, 49}, 5, true);
        assertTrue(result.compareTo(expected) == 0);
        
        num1 = result;
        num2 = VeryLargeNumber.createMaxLong();
        result = VeryLargeNumber.karatsuba(num1, num2);
        assertNotNull(result);
        //784637716923335095224261902710254454442933591094742482943L
        expected = new VeryLargeNumber(0);      
        expected.setInternalArray(new int[]{512, 3072, 7488, 9472, 6552, 2352, 343}, 7, true);
        assertTrue(result.compareTo(expected) == 0);
        
        //TODO: not correct yet.
        // this number is completely backwards, but has expected components
        num1 = result;
        num2 = result.clone();
        result = VeryLargeNumber.karatsuba(num1, num2);
        assertNotNull(result);
        //615656346818663737291362432329573325363859439854215515892590030606645220018700810278654093633374614700204645941249L
        expected = new VeryLargeNumber(0);      
        expected.setInternalArray(new int[]{
            262144, 3145728, 17104896, 55705600, 120975360, 184516608, 
            202643456, 161452032, 92621760, 37318400, 10026576, 1613472, 
            117649}, 13,
            true);
        /*
        int[] pr = Arrays.copyOf(result.getInternalArray(), 
            result.getInternalArraySize());
        System.out.println("result=" + Arrays.toString(pr));
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < (pr.length - 1); i++) {
            if (sb.length() > 0) {
                sb.append(" + ");
            }
            sb.append(Integer.toString(pr[i]));
            int nb = pr.length - i - 1;
            sb.append("*math.pow(base,").append(Integer.toString(nb)).append(")");
        }
        sb.append(" + ").append(Integer.toString(pr[(pr.length - 1)]));
        System.out.println(sb.toString());
        */
        assertTrue(result.compareTo(expected) == 0);
        
    }
    
    @Test
    public void testPow() {
        
        VeryLargeNumber instance = new VeryLargeNumber(2);
        VeryLargeNumber result = instance.pow(-1);
        double r = result.getDoubleValueIfExists();
        double exp = 0.5;
        assertTrue(Math.abs(r - 0.5) < 0.01*exp);
        
        instance = new VeryLargeNumber(Integer.MAX_VALUE);
        result = instance.pow(-3);
        r = result.getDoubleValueIfExists();
        exp = 1.0097419600934883e-28;
        assertTrue(Math.abs(r - 1.0097419600934883e-28) < 0.01*exp);
        
        instance = new VeryLargeNumber(2);
        result = instance.pow(0);
        VeryLargeNumber expected = new VeryLargeNumber(1);
        assertTrue(expected.compareTo(result) == 0);
        
        instance = new VeryLargeNumber(2);
        result = instance.pow(1);
        expected = new VeryLargeNumber(2);
        assertTrue(expected.compareTo(result) == 0);
        
        instance = new VeryLargeNumber(2);
        result = instance.pow(4);
        expected = new VeryLargeNumber(16);
        assertTrue(expected.compareTo(result) == 0);
        
        instance = new VeryLargeNumber(2);
        result = instance.pow(5);
        expected = new VeryLargeNumber(32);
        assertTrue(expected.compareTo(result) == 0);
        
        
        //9223372036854775807
        instance = VeryLargeNumber.createMaxLong();
        result = instance.pow(5);
        //6.674959487252844e+94
        //66749594872528440038659400431137192519314960677663778810642210898138927029830146832929608695807L
        VeryLargeNumber expected2 = new VeryLargeNumber(0);
        expected.setInternalArray(new int[]{
            32768, 327680, 1454080, 3768320, 6312960, 7141376, 5523840, 
            2885120, 974120, 192080, 16807}, 11,
            true);
        assertTrue(expected.compareTo(result) == 0);
    }
    
    //@Test
    public void estDivideByAndPrint2() {
        
        //---------------------------------------------------------
        // very large number test
        //9223372036854775807
        VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();
        
        VeryLargeNumber instance = new VeryLargeNumber(Integer.MAX_VALUE);
        instance.add(maxLong);
        // ===> 9223372039002259454
        
        VeryLargeNumber divisor = new VeryLargeNumber(Integer.MAX_VALUE);
        divisor.add(new VeryLargeNumber(Integer.MAX_VALUE));
        // ==> 4294967294
        
        // 9223372039002259454 / 4294967294
        String str = instance.divideByAndPrint(divisor);
        
        System.out.println("9223372039002259454 / 4294967294=" + str);
        
        assertTrue(str.startsWith("2147483649.5"));
    
    }
    
    public void testMoveUpIfStartsWithZeros() {
        
        int[] a = new int[]{0, 0, 0, 1, 2, 3, 4, 5};
                
        int len = VeryLargeNumber.moveUpIfStartsWithZeros(a, a.length);
        
        assertTrue(Arrays.equals(new int[]{1, 2, 3, 4, 5, 0, 0, 0}, a));
        
        assertTrue(len == 5);
        
        a = new int[]{1, 2, 3, 4, 5};
        len = VeryLargeNumber.moveUpIfStartsWithZeros(a, a.length);
        
        assertTrue(Arrays.equals(new int[]{1, 2, 3, 4, 5}, a));
        
        assertTrue(len == 0);
    }
    
    public void testShiftDown() {
        
        int[] a = new int[]{1, 2, 3, 4, 5, 0, 0, 0};
                
        VeryLargeNumber.moveDown(a, 3);
        
        assertTrue(Arrays.equals(new int[]{0, 0, 0, 1, 2, 3, 4, 5}, a));
        
    }
    
    @Test
    public void testMultiply() {
        
        VeryLargeNumber v = new VeryLargeNumber(0);
        VeryLargeNumber result = v.multiplySmall(new VeryLargeNumber(0));
        
        VeryLargeNumber expected = new VeryLargeNumber(0);
        
        assertTrue(result.compareTo(expected) == 0);
        
        // ---------------------
        v = new VeryLargeNumber(1);
        result = v.multiplySmall(new VeryLargeNumber(2));
        
        expected = new VeryLargeNumber(2);
        assertTrue(result.compareTo(expected) == 0);
        
        //-------------------------------------------------
        v = new VeryLargeNumber(1111);
        result = v.multiplySmall(new VeryLargeNumber(234));
        
        expected = new VeryLargeNumber(259974);
        assertTrue(result.compareTo(expected) == 0);
        
        //-------------------------------------------------
        v = new VeryLargeNumber(-1111);
        result = v.multiplySmall(new VeryLargeNumber(-234));
        
        expected = new VeryLargeNumber(259974);
        assertTrue(result.compareTo(expected) == 0);
        
        //-------------------------------------------------
        v = new VeryLargeNumber(-1111);
        result = v.multiplySmall(new VeryLargeNumber(234));
        
        expected = new VeryLargeNumber(-259974);
        assertTrue(result.compareTo(expected) == 0);
        
        //-------------------------------------------------
        v = new VeryLargeNumber(1111);
        result = v.multiplySmall(new VeryLargeNumber(-234));
        
        expected = new VeryLargeNumber(-259974);
        assertTrue(result.compareTo(expected) == 0);
        
        //-------------------------------------------------
        //9223372036854775807L * 2147483647 = 19807040619342712359383728129L
        VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();        
        v = new VeryLargeNumber(Integer.MAX_VALUE);
        result = v.multiplySmall(maxLong);
        
        /*
        19807040619342712359383728129 % base = 7 <===
        19807040619342712359383728129 /= base = 18446744082299486214L
        7, 30, 40, 16, 0
        */
        expected = new VeryLargeNumber(0);
        expected.setInternalArray(
            new int[]{16, 40, 30, 7}, 4, true);
        assertTrue(result.compareTo(expected) == 0);
    }
    
}
