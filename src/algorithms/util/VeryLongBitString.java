package algorithms.util;

import algorithms.misc.MiscMath;
import java.util.Arrays;

/**
 * 
 * creating a modifiable holder for a very long bit string of fixed size nbits. 
 * the internal array type is long which means each item is a 64 bit holder
 * (128 bits on the stack of a 64 bit platform).  
 * the number of items in an array is limited to 32 bits, so that's a capacity
 * of bits in the billions, if the JVM has enough memory for the maximum size
 * of the long array.  
 * 
 * @author nichole
 */
public final class VeryLongBitString {
    
    private final long[] bitstring;
    
    protected final long nBits;
    
    protected static final long itemBitLength = 64L;
    
    protected static final long nMaxBits = Integer.MAX_VALUE * itemBitLength;
    
    protected final long capacityBits;
   
    public VeryLongBitString(long nBits) {
        
        if (nBits > nMaxBits) {
            throw new IllegalArgumentException("cannot hold more than " + nMaxBits + " bits");
        }
        
        int nElements = getRowNumber(nBits) + 1;
        
        bitstring = new long[nElements];
        
        capacityBits = nElements * itemBitLength;
        
        this.nBits = nBits;
    }
    
    protected VeryLongBitString(long[] bitstrings, long nBits) {
        
        if (nBits > nMaxBits) {
            throw new IllegalArgumentException("cannot hold more than " + nMaxBits + " bits");
        }
        
        bitstring = Arrays.copyOf(bitstrings, bitstrings.length);
        
        capacityBits = bitstring.length * itemBitLength;
        
        this.nBits = nBits;
    }
    
    public long getCapacity() {
        return capacityBits;
    }
    
    public void setBit(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        bitstring[idx] |= (1L << bitIdx);
        
    }
    
    public void clearBit(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        bitstring[idx] &= ~(1L << bitIdx);
    }
    
    public void toggleBit(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        bitstring[idx] ^= (1L << bitIdx);
    }
    
    public boolean isSet(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        return ((bitstring[idx] & (1L << bitIdx)) != 0);
    }
    
    public boolean isNotSet(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        return ((bitstring[idx] & (1L << bitIdx)) == 0);
    }
         
    protected int getRowNumber(long n) {
        
        int nthElement = (int)(n/itemBitLength);
        
        return nthElement;
    }
    
    /**
     * given the arrayIdx, calculate the bit position within that item
     * for the very large bitstring position nthBit.  Note that if the
     * result is out of bounds, -1 is returned and the invoker should
     * handle that.
     * 
     * @param nthBit
     * @param arrayIdx
     * @return 
     */
    int getBitIdx(long nthBit, int arrayIdx) {
        
        int bitIdx = (int)(nthBit - (arrayIdx * itemBitLength));
        
        if ((bitIdx < 0) || (bitIdx > (itemBitLength - 1))) {
            return -1;
        }
        
        return bitIdx;
    }
    
    public void clearAllBits() {
        Arrays.fill(bitstring, 0);
    }
    
    public VeryLongBitString copy() {
        
        VeryLongBitString c = new VeryLongBitString(this.bitstring, this.nBits);
        
        return c;
    }
    
    public void resetAllTo(VeryLongBitString other) {
        if (other.nBits != nBits) {
            throw new IllegalArgumentException("nBits must be the same in both to use this method");
        }
        System.arraycopy(other.bitstring, 0, bitstring, 0, other.bitstring.length);
    }
    
    /**
     * get the '01...' bitstring representation of this object
     * @return 
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = (bitstring.length - 1); i > -1 ; i--) {
            sb.append(Long.toBinaryString(bitstring[i]));
        }
        return sb.toString();
    }
    
    /**
     * get a list of the bit numbers that are set.
     * @return 
     */
    public int[] getSetBits() {
        
        int n = 0;
        for (int i = 0; i < bitstring.length; ++i) {
            n += Long.bitCount(bitstring[i]);
        }
        
        int[] setBits = new int[n];
        for (int i = 0; i < bitstring.length; ++i) {
            long b = bitstring[i];
            int count = 0;
            n = 0;
            while (b > 0) {
                long idx = (b & 1L);
                if (idx == 1) {
                    setBits[n] = count;
                    n++;
                }
                b >>= 1L;
                count++;
            }
        }
        
        return setBits;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof VeryLongBitString)) {
            return false;
        }
        
        VeryLongBitString other = (VeryLongBitString)obj;
        
        if (nBits != other.nBits) {
            return false;
        }
        
        return Arrays.equals(this.bitstring, other.bitstring);
    }    

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 11 * hash + Arrays.hashCode(this.bitstring);
        hash = 11 * hash + (int) (this.nBits ^ (this.nBits >>> 32));
        return hash;
    }
    
    /**
     * approximate the memory (in Bytes) used by this instance and its member variables, and 
     * return that portion on the stack.
     * @return 
     */
    protected long approximateMemoryUsed_Stack() {
       
        String arch = System.getProperty("sun.arch.data.model");
        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
        int nbits = (is32Bit) ? 32 : 64;
        
        /*
        instance contains members:
            long[], long, long, long
        
        mem Stack:
           3 long primitives (a long is 2 times word size)
           1 array ref
        mem Heap:
           object overhead
           sum of each item in the long array
        */
        
        long sumBits = 3*(nbits*2) + 32;
        
        long sumBytes = (sumBits/8);
        if (sumBytes*8 > sumBits) {
            sumBytes++;
        }
        
        return sumBytes;
    }
    
    /**
     * approximate the memory (in Bytes) used by this instance and its member variables, and 
     * return that portion on the heap.
     * @return 
     */
    protected long approximateMemoryUsed_Heap() {
        /*
        see comments within approximateMemoryUsed_Stack
        
        mem Heap:
           object overhead
           sum of each item in the long array  (a long is 2*word size)
        */
        
        String arch = System.getProperty("sun.arch.data.model");
        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
        int nbits = (is32Bit) ? 32 : 64;
        
        long sumBytes = 16;
        if (bitstring != null) {
            sumBytes += bitstring.length * (2 * nbits);
        }
        long padding = (sumBytes % 8);
        sumBytes += padding;
        
        return sumBytes;
    }

}
