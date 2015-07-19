package algorithms;

/**
Uses Gosper's hack from
  http://read.seas.harvard.edu/cs207/2012/
     
 * @author nichole
 */
public class SubsetChooser {
    
    private final int n;
    
    private final int k;
    
    private final long highBit;
    
    private final long check;
    
    private long x;
        
    /**
     * constructor with the number of indexes to choose from, n, and the size of
     * the subset, k
     * @param n
     * @param k 
     */
    public SubsetChooser(int n, int k) {
        
        this.n = n;
        
        this.k = k;
        
        highBit = 1L << n;
        
        x = (1L << k) - 1;
        
        check = 1L << n;
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
        
        long nextX = nextSubset64(x);
        
        if ((nextX & check) != 0)  {
            return -1;
        }
        
        int nValues = select(outputIndexes);
        
        x = nextX;
                
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
        x0 = (((x0 ^ c) >> 2) / y) | c;
        if ((x0 & highBit) > 0) {
            x0 = ((x0 & (highBit - 1)) << 2) | 3;
        }
        return x0;
    }
    
    protected int select(int[] selected) {
        
        // interpret the bit string:  1 is 'selected' and 0 is not
        
        /*String str = Long.toBinaryString(x);
        while (str.length() < n) {
            str = "0" + str;
        }
        System.out.format("%d\t%10s\n", x, str);
        */
            
        int nBits = 0;
        int nOneBits = 0;
        long xp = x;
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

}
