package algorithms;

import algorithms.misc.MiscMath;
import java.math.BigInteger;

/**
Class to iterate over every combination of subsets within n objects in an 
ordered manner.

The class uses Gosper's hack from
  http://read.seas.harvard.edu/cs207/2012/

 * @author nichole
 */
public class SubsetChooser {
    
    private final int n;

    private final int k;

    private long x64;
    
    private final long highBit64;
    
    private BigInteger x;
    
    private final BigInteger highBit;

    private long count = 0;

    private final long np;

    /**
     * constructor with the number of indexes to choose from, n, and the size of
     * the subset, k.
     * @param n the number of indexes that the selector will choose from
     * @param k the subset size of selected indexes.  the maximum value currently
     * accepted is 12.
     */
    public SubsetChooser(int n, int k) {
        
        if (k > 13) {
            throw new IllegalArgumentException(
                "currently, class can only handle k < 13, but changes to accomodate larger could be made");
        }
        if (n < 1) {
            throw new IllegalArgumentException("n must be larger than 0");
        }
        if (k < 1) {
            throw new IllegalArgumentException("k must be larger than 0");
        }
        if (k > n) {
            throw new IllegalArgumentException("k must be less than or equal to n");
        }
        
        this.n = n;

        this.k = k;

        count = 1;

        np = MiscMath.computeNDivNMinusK(n, k)/MiscMath.factorial(k);

        if (n < 64) {
            
            highBit64 = 1L << n;

            x64 = (1L << k) - 1;

            x = null;
            
            highBit = null;
            
        } else {
            
            byte[] val = MiscMath.writeToBigEndianBytes((1L << k) - 1);
            x = new BigInteger(val);
            
            BigInteger hb = new BigInteger(MiscMath.writeToBigEndianBytes(1));
            highBit = hb.shiftLeft(n);  
            
            highBit64 = Long.MAX_VALUE;
            x64 = Long.MAX_VALUE;
        }
    }

    /**
     * given a constructed array, populates it with the next selected subset
     * of indexes and returns the number of values placed in the subset.
     * Returns a -1 when there are no more subsets to return;
     * @param outputIndexes
     * @return
     */
    public int getNextSubset(int[] outputIndexes) {

        if (outputIndexes == null || outputIndexes.length != k) {
            throw new IllegalArgumentException(
                "outputIndexes cannot be null and has to be size k");
        }
        
        if (n < 64) {
            return getNextSubset64(outputIndexes);
        } else {
            return getNextSubsetBigInteger(outputIndexes);
        }
        
    }
    
    protected int getNextSubsetBigInteger(int[] outputIndexes) {

        if (outputIndexes == null || outputIndexes.length != k) {
            throw new IllegalArgumentException(
                "outputIndexes cannot be null and has to be size k");
        }
        
        if (count > np) {
            return -1;
        }

        int nValues = selectBigInteger(outputIndexes);

        x = nextSubsetBigInteger(x);

        return nValues;
    }
    
    /**
     * given a constructed array, populates it with the next selected subset
     * of indexes and returns the number of values placed in the subset.
     * Returns a -1 when there are no more subsets to return;
     * @param outputIndexes
     * @return
     */
    protected int getNextSubset64(int[] outputIndexes) {

        if (outputIndexes == null || outputIndexes.length != k) {
            throw new IllegalArgumentException(
                "outputIndexes cannot be null and has to be size k");
        }

        if (count > np) {
            return -1;
        }

        int nValues = select64(outputIndexes);

        x64 = nextSubset64(x64);

        return nValues;
    }

    /**
     * highbit is 1 << entire_set_size (not subset size, k).
     *
     * @param x
     * @param highbit
     * @return
     */
    private long nextSubset64(long x0) {

        long y = x0 & -x0;  // = the least significant one bit of x
        long c = x0 + y;

        x0 = c + (((c ^ x0) / y) >> 2);

        count++;

        return x0;
    }

    /**
     * highbit is 1 << entire_set_size (not subset size, k).
     *
     * @param x
     * @param highbit
     * @return
     */
    private BigInteger nextSubsetBigInteger(BigInteger x0) {

        BigInteger y = x0.and(x0.negate()); // = the least significant one bit of x
        BigInteger c = x0.add(y);
        
        //x0 = c + (((c ^ x0) / y) >> 2);
        BigInteger tmp = c.xor(x0).divide(y).shiftRight(2);
        
        x0 = c.add(tmp);

        count++;

        return x0;
    }

    protected int select64(int[] selected) {

        // interpret the bit string x:  1 is 'selected' and 0 is not

        /*
        String str = Long.toBinaryString(x);
        while (str.length() < n) {
            str = "0" + str;
        }
        System.out.format("%d\t%10s\n", x, str);
        */

        int nBits = 0;
        int nOneBits = 0;
        long xp = x64;
        while (xp > 0) {
            if ((xp & 1) == 1) {
                int idx2 = n - 1 - nBits;
                selected[nOneBits] = idx2;
                nOneBits++;
            }
            xp = xp >> 1;
            nBits++;
        }

        return nOneBits;
    }

    protected int selectBigInteger(int[] selected) {

        // interpret the bit string x:  1 is 'selected' and 0 is not

        /*
        String str = Long.toBinaryString(x);
        while (str.length() < n) {
            str = "0" + str;
        }
        System.out.format("%d\t%10s\n", x, str);
        */

        int nBits = 0;
        int nOneBits = 0;
        BigInteger xp = new BigInteger(x.toByteArray());
        while (xp.compareTo(BigInteger.ZERO) > 0) {
            if (xp.and(BigInteger.ONE).equals(BigInteger.ONE)) {
                int idx2 = n - 1 - nBits;
                selected[nOneBits] = idx2;
                nOneBits++;
            }
            xp = xp.shiftRight(1);
            nBits++;
        }

        return nOneBits;
    }
}
