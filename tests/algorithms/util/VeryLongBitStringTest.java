package algorithms.util;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class VeryLongBitStringTest extends TestCase {
    
    protected Logger log = Logger.getLogger(VeryLongBitStringTest.class.getName());
    
    public VeryLongBitStringTest(String testName) {
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

    public void testGetRowNumber() {
        System.out.println("getRowNumber");
        
        VeryLongBitString bitstring = new VeryLongBitString(2);
        assertTrue(bitstring.getRowNumber(1) == 0);
        
        bitstring = new VeryLongBitString(64);
        assertTrue(bitstring.getRowNumber(1) == 0);
        
        bitstring = new VeryLongBitString(Integer.MAX_VALUE);
        assertTrue(bitstring.getRowNumber(65) == 1);
        
        assertTrue(bitstring.getRowNumber((64*5) + 2) == 5);
        
        assertTrue(bitstring.getRowNumber((64*5) -1) == 4);
        
        int itemIdx = bitstring.getRowNumber(1);
        int bitIdx = bitstring.getBitIdx(1, itemIdx);
        assertTrue(itemIdx == 0);
        assertTrue(bitIdx == 1);
        
        // bits are 0 through 63, 64 through 127
        itemIdx = bitstring.getRowNumber(65);
        bitIdx = bitstring.getBitIdx(65, itemIdx);
        assertTrue(itemIdx == 1);
        assertTrue(bitIdx == 1);
        
        int nBit = 64 + (64/2);
        itemIdx = bitstring.getRowNumber(nBit);
        bitIdx = bitstring.getBitIdx(nBit, itemIdx);
        assertTrue(itemIdx == 1);
        assertTrue(bitIdx == (64/2));
        
    }
    
    public void testSetBit() throws Exception {
        
        System.out.println("testSetBit");
        
        SecureRandom secureRand = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed=1393807938003554000l;
        secureRand.setSeed( seed );
        log.info("seed=" + seed);
        
        int nTests = 100;
        
        // make nTests random bitstrings of size up to value 10^3
        
        // similar to the programming pearls example 1 bitvector, use a very 
        // large bitstring to hold each number as a bitstring itself.  so the 
        // number 1 million would be all 0's and only the 1millionth-1 bit set.
        
        int max = (int)Math.pow(10, 3);
        
        LinkedHashSet<Integer> numbers = new LinkedHashSet<Integer>();
       
        VeryLongBitString[] bitstrings = new VeryLongBitString[nTests];
        for (int i = 0; i < bitstrings.length; i++) {
            
            int r = secureRand.nextInt(max);
            
            while (numbers.contains(Integer.valueOf(r))) {
                r = secureRand.nextInt(max);
            }
         
            numbers.add(Integer.valueOf(r));
            
            bitstrings[i] = new VeryLongBitString(max);
            
            bitstrings[i].setBit(r);
        }
        
        // ===== test isSet bit and isNotSet bit ========
        
        // check the test bits... numbers order is preserved so iterate over it
        Iterator<Integer> iter = numbers.iterator();
        int count = 0;
        while (iter.hasNext()) {
            
            Integer r = iter.next();
            
            // check that bit r is set and the rest are not in each bitstring
            VeryLongBitString bitstring = bitstrings[count];
            
            assertTrue(bitstring.isSet(r.intValue()));
            
            long capacity = bitstring.getCapacity();
            
            assertTrue(capacity >= r.intValue());
            
            for (int i = 0; i < (int)capacity; i++) {
                if (i != r.intValue()) {
                    assertTrue(bitstring.isNotSet(i));
                }
            }
            
            count++;
        }
        
        // ===== test clear bit ========
        iter = numbers.iterator();
        count = 0;
        while (iter.hasNext()) {
            
            Integer r = iter.next();
            
            // check that bit r is set and the rest are not in each bitstring
            VeryLongBitString bitstring = bitstrings[count];
                        
            bitstring.clearBit(r.intValue());
            
            long capacity = bitstring.getCapacity();
            
            assertTrue(capacity >= r.intValue());
            
            for (int i = 0; i < (int)capacity; i++) {
                assertTrue(bitstring.isNotSet(i));
            }
            
            count++;
        }
        
    }

    public void testSetBit2() throws Exception {
        
        System.out.println("testSetBit2");
        long seed = System.currentTimeMillis();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);
        
        // test for the more common use of using the bitstring form
        // of a number
        
        int nTests = 100;
        
        int maxNBits = 3*64;
        
        // testing w/ BigInteger as a bootstrap
        
        List<BigInteger> numbers = new ArrayList<BigInteger>();
        
        VeryLongBitString[] bitstrings = new VeryLongBitString[nTests];
        for (int i = 0; i < bitstrings.length; i++) {
            
            bitstrings[i] = new VeryLongBitString(maxNBits);
            
            BigInteger r = BigInteger.probablePrime(maxNBits, rand);
            
            assertTrue(numbers.add(r));
            
            for (int ii = 0; ii < maxNBits; ii++) {
                if (r.testBit(ii)) {
                    bitstrings[i].setBit(ii);
                }
            }
        }
        
        // test bits were set and not set
        for (int i = 0; i < numbers.size(); i++) {
            
            BigInteger r = numbers.get(i);
            
            VeryLongBitString bitstring = bitstrings[i];
            
            for (int ii = 0; ii < maxNBits; ii++) {
                if (r.testBit(ii)) {
                    assertTrue(bitstring.isSet(ii));
                    assertFalse(bitstring.isNotSet(ii));
                } else {
                    assertTrue(bitstring.isNotSet(ii));
                    assertFalse(bitstring.isSet(ii));
                }
            }
        }        
    }

    public void testToggleBit() throws Exception {
        
        System.out.println("testToggleBit");
        
        SecureRandom secureRand = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed=1393807938003554000l;
        secureRand.setSeed(seed);
        Random rand = new Random(seed);
        log.info("seed=" + seed);
        
        int nTests = 100;
        int maxNBits = (2 * 64) + (64/2);
        VeryLongBitString bitstring = null;
               
        for (int i = 0; i < nTests; i++) {
            
            bitstring = new VeryLongBitString(maxNBits);
            
            int chkBit = secureRand.nextInt(maxNBits);
           
            BigInteger r = BigInteger.probablePrime(maxNBits, rand);
                        
            for (int ii = 0; ii < maxNBits; ii++) {
                if (r.testBit(ii)) {
                    bitstring.setBit(ii);
                }
            }
            
            boolean bitIsSet = bitstring.isSet(chkBit);
            
            if (!bitIsSet) {
                bitstring.toggleBit(chkBit);
                assertTrue(bitstring.isSet(chkBit));
            } else {
                bitstring.toggleBit(chkBit);
                assertTrue(bitstring.isNotSet(chkBit));
            }            
        }
        
        // test clear all bits
        bitstring.clearAllBits();
        for (int i = 0; i < maxNBits; i++) {
            assertTrue(bitstring.isNotSet(i));
        }
        
    }

    public void testCopy() throws Exception {
        
        System.out.println("testCopy");
        
        long seed = System.currentTimeMillis();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);
        
        int nTests = 100;
        
        int maxNBits = 3*64;
        
        // testing w/ BigInteger as a bootstrap
                
        for (int i = 0; i < nTests; i++) {
            
            VeryLongBitString b = new VeryLongBitString(maxNBits);
            
            BigInteger r = BigInteger.probablePrime(maxNBits, rand);
                        
            for (int ii = 0; ii < maxNBits; ii++) {
                if (r.testBit(ii)) {
                    b.setBit(ii);
                    assertTrue(b.isSet(ii));
                } else {
                    assertTrue(b.isNotSet(ii));
                }
            }
            
            VeryLongBitString c = b.copy();
            
            for (int ii = 0; ii < maxNBits; ii++) {
                if (r.testBit(ii)) {
                    assertTrue(c.isSet(ii));
                } else {
                    assertTrue(c.isNotSet(ii));
                }
            }
        }
        
    }
    
    public void testResetTo() throws Exception {
        
        System.out.println("testResetTo");
        
        long seed = System.currentTimeMillis();
        System.out.println("seed=" + seed);
        Random rand = new Random(seed);
        
        int nTests = 100;
        
        int maxNBits = 3*64;
        
        VeryLongBitString a = new VeryLongBitString(maxNBits);
                        
        for (int i = 0; i < nTests; i++) {
            
            VeryLongBitString b = new VeryLongBitString(maxNBits);
            
            BigInteger r = BigInteger.probablePrime(maxNBits, rand);
                        
            for (int ii = 0; ii < maxNBits; ii++) {
                if (r.testBit(ii)) {
                    b.setBit(ii);
                    assertTrue(b.isSet(ii));
                } else {
                    assertTrue(b.isNotSet(ii));
                }
            }
            
            b.resetAllTo(a);
            
            for (int ii = 0; ii < maxNBits; ii++) {
                if (a.isSet(ii)) {
                    assertTrue(b.isSet(ii));
                } else {
                    assertTrue(b.isNotSet(ii));
                }
            }
        }
    }
    
    public void testEquals() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        sr.setSeed(seed);
                
        for (int i = 0; i < 1000; ++i) {
            int sz = 10 + sr.nextInt(3000);
            
            VeryLongBitString vlbs1 = new VeryLongBitString(sz);
            VeryLongBitString vlbs2 = new VeryLongBitString(sz);
            
            for (int j = 0; j < 10; ++j) {
                int v = sr.nextInt(sz - 1);
                vlbs1.setBit(v);
                vlbs2.setBit(v);
            }
            
            assertTrue(vlbs1.equals(vlbs2));
        }
    }
    
    public void testGetSetBits() {
        
        VeryLongBitString bitstring = new VeryLongBitString(16);
        
        int[] a = new int[]{0, 2, 3, 10, 15};
        for (int i = 0; i < a.length; ++i) {
            bitstring.setBit(a[i]);
        }
        
        int[] s = bitstring.getSetBits();
        
        assertEquals(a.length, s.length);
        
        for (int i = 0; i < a.length; ++i) {
            assertEquals(a[i], s[i]);
        }
        
        // -------
        bitstring = new VeryLongBitString(128);
        
        a = new int[]{0, 2, 3, 10, 15, 65, 88, 99, 121, 127};
        for (int i = 0; i < a.length; ++i) {
            bitstring.setBit(a[i]);
        }
        
        s = bitstring.getSetBits();
        
        assertEquals(a.length, s.length);
        
        for (int i = 0; i < a.length; ++i) {
            assertEquals(a[i], s[i]);
        }
    }
    
    public void testOr() {
        
        //0 1 1 0   
        VeryLongBitString bs6 = new VeryLongBitString(4);
        bs6.setBit(1);
        bs6.setBit(2);
        
        //1 0 1 0
        VeryLongBitString bs10 = new VeryLongBitString(4);
        bs10.setBit(1);
        bs10.setBit(3);
        
        VeryLongBitString bs14 = bs6.or(bs10);
        
        int[] setBits = bs14.getSetBits();
        
        assertTrue(Arrays.equals(new int[]{1, 2, 3}, setBits));
    }
    
    public void testAnd() {
        
        // 86                    74                  66
        // 0 1 0 1 0 1 1 0   &   0 1 0 0 1 0 1 0  =  0 1 0 0 0 0 1 0
        
            
        VeryLongBitString bs1 = new VeryLongBitString(8);
        bs1.setBit(1); bs1.setBit(2); bs1.setBit(4); bs1.setBit(6);
        
        VeryLongBitString bs2 = new VeryLongBitString(8);
        bs2.setBit(1); bs2.setBit(3); bs2.setBit(6);
        
        VeryLongBitString bs3 = new VeryLongBitString(8);
        bs3.setBit(1); bs3.setBit(6);
        
        VeryLongBitString r = bs1.and(bs2);
        assertEquals(bs3, r);
        
        r = bs2.and(bs1);
        assertEquals(bs3, r);
    }
    
    public void testDifference() {
        // 86                    74                  66
        // 0 1 0 1 0 1 1 0   &   0 1 0 0 1 0 1 0  =  0 1 0 0 0 0 1 0
        
        // 86                   66                   20
        // 0 1 0 1 0 1 1 0   &  ~0 1 0 0 0 0 1 0  =  0 0 0 1 0 1 0 0
        
        // 74                   66                   8
        // 0 1 0 0 1 0 1 0   &  ~0 1 0 0 0 0 1 0  =  0 0 0 0 1 0 0 0
        
        VeryLongBitString bs1 = new VeryLongBitString(8);
        bs1.setBit(1); bs1.setBit(2); bs1.setBit(4); bs1.setBit(6);
        
        VeryLongBitString bs2 = new VeryLongBitString(8);
        bs2.setBit(1); bs2.setBit(3); bs2.setBit(6);
        
        // 20
        VeryLongBitString bs3 = new VeryLongBitString(8);
        bs3.setBit(2); bs3.setBit(4);
        
        // 8
        VeryLongBitString bs4 = new VeryLongBitString(8);
        bs4.setBit(3); 
        
        VeryLongBitString r = bs1.difference(bs2);
        assertEquals(bs3, r);
        
        r = bs2.difference(bs1);
        assertEquals(bs4, r);
    }
    
    public void testNextHighestBitSet() {
        
        // 0 1 0 1 1 0
        VeryLongBitString bs = new VeryLongBitString(256);
        bs.setBit(1);
        bs.setBit(2);
        bs.setBit(4);
        bs.setBit(34);
        bs.setBit(37);
        
        int b = bs.nextHighestBitSet(1);
        //System.out.println("b=" + b);
        assertEquals(2, b);
        
        b = bs.nextHighestBitSet(34);
        //System.out.println("b=" + b);
        assertEquals(37, b);
    }
    
    public void testLeastBitSet() {
        
        VeryLongBitString bs = new VeryLongBitString(256);
        bs.setBit(2);
        bs.setBit(3);
        long lb = bs.leastBitSet();
        //System.out.println("leastBitSet=" + lb);
        assertEquals(2, lb);
        
        bs = new VeryLongBitString(256);
        bs.setBit(5);
        bs.setBit(7);
        lb = bs.leastBitSet();
        //System.out.println("leastBitSet=" + lb);
        assertEquals(5, lb);
    }
    
    public void testHighestBitSet() {
        
        VeryLongBitString bs = new VeryLongBitString(256);
        bs.setBit(2);
        bs.setBit(3);
        long lb = bs.highestBitSet();
        //System.out.println("leastBitSet=" + lb);
        assertEquals(3, lb);
        
        bs = new VeryLongBitString(256);
        bs.setBit(5);
        bs.setBit(7);
        lb = bs.highestBitSet();
        //System.out.println("leastBitSet=" + lb);
        assertEquals(7, lb);
    }
    
    public void testNextHighestBitSet_2() {
        
        // 0 1 0 1 1 0
        VeryLongBitString bs = new VeryLongBitString(256);
        bs.setBit(201);
        bs.setBit(82);
        bs.setBit(67);
        bs.setBit(73);
        bs.setBit(62);
        bs.setBit(59);
        bs.clearBit(59);
        
        int b = bs.nextHighestBitSet(59);
        //System.out.println("b=" + b);
        assertEquals(62, b);
        
        bs.setBit(63);
        bs.setBit(59);
        bs.clearBit(59);
        
        b = bs.nextHighestBitSet(59);
        //System.out.println("b=" + b);
        assertEquals(62, b);
    
        bs = new VeryLongBitString(256);
        bs.setBit(61);
        bs.setBit(63);
        bs.setBit(64);
        bs.setBit(65);
        b = bs.nextHighestBitSet(62);
        assertEquals(63, b);
        b = bs.nextHighestBitSet(64);
        assertEquals(65, b);
        b = bs.nextHighestBitSet(63);
        System.out.println("b=" + b);
        assertEquals(64, b);
        
        bs.setBit(61 + 64);
        bs.setBit(63 + 64);
        bs.setBit(64 + 64);
        bs.setBit(65 + 64);
        b = bs.nextHighestBitSet(62 + 64);
        assertEquals(63 + 64, b);
        b = bs.nextHighestBitSet(64 + 64);
        assertEquals(65 + 64, b);
        b = bs.nextHighestBitSet(63 + 64);
        System.out.println("b=" + b);
        assertEquals(64 + 64, b);
    }
}
