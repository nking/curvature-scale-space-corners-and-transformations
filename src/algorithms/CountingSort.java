package algorithms;

import algorithms.util.PairInt;
import java.util.Arrays;

/**
 * a sort for integers in range 0 to k that has an O(N) runtime at the expense 
 * of space where N is the maximum of (size of array to be sorted, maximum value in array).
 * <pre>
 * To use this algorithm:
 *    (1) numbers must be positive.
 *    (2) the maximum number in the array should probably not be much greater
 *        than 10,000,000 unless the jvm settings for maximum stack size
 *        are increased.  An internal long array of size maximum of array values
 *        is constructed and that consumes memory which also affects
 *        performance for max &gt; 10,000,000.
 * 
 * implemented from Cormen et al. "Introduction to Algorithms"
 * </pre>
 * @author nichole
 */
public class CountingSort {
    
    /**
     * sort the members of a with values less than max and return as an
     * ordered array.
     * <pre>
     * Note that the numbers have to be positive, so if negative numbers are
     * in the array, the invoker needs to add a number to bring the values
     * to &gt; == 0 and then subtract that after the sort.
     * </pre>
     * @param a
     * @param max 
     * @return  
     */
    public static int[] sort(int[] a, int max) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (max <= 0) {
            throw new IllegalArgumentException("max must be > 0");
        }
        
        if (a.length < 2) {
            return Arrays.copyOf(a, a.length);
        }
        
        long[] c = new long[max + 1];

        // c holds frequency of each number by index, e.g. c[0] holds the number of 0's in a
        for (int i = 0; i < a.length; i++) {
            
            int idx = a[i];
            
            c[idx]++;
        }
                
        // cumulative sum to end of array c.  the last item in c holds the 
        // total number of items in 'a' less than or equal to max
        for (int i = 1; i < c.length; i++) {
            c[i] = c[i] + c[i - 1];
        }
        
        int[] b = new int[a.length];
                
        // use the order imposed by c to write the values of a into b.  c holds
        // frequency, that is place markers too so that is updated as b is written
        for (int i = (a.length - 1); i > -1; i--) {
            
            int aIdx = a[i];
                        
            int cfa = (int)c[aIdx];
            b[cfa - 1] = aIdx;
            c[aIdx]--;
        }

        return b;
    }
    
    /**
     * sort the members of a with values less than max and apply the same
     * changes of item position to b.
     * <pre>
     * Note that the numbers have to be positive, so if negative numbers are
     * in the array, the invoker needs to add a number to bring the values
     * to &gte; == 0 and then subtract that after the sort.
     * </pre>
     * @param a
     * @param b
     * @param max 
     */
    public static void sort(int[] a, int[] b, int max) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "the lengths of a and b must be the same");
        }
        if (max <= 0) {
            throw new IllegalArgumentException("max must be > 0");
        }
        if (a.length < 2) {
            return;
        }
        
        long[] c = new long[max + 1];

        // c holds frequency of each number by index, e.g. c[0] holds the number of 0's in a
        for (int i = 0; i < a.length; i++) {
            
            int idx = a[i];
            
            c[idx]++;
        }
                
        // cumulative sum to end of array c.  the last item in c holds the 
        // total number of items in 'a' less than or equal to max
        for (int i = 1; i < c.length; i++) {
            c[i] = c[i] + c[i - 1];
        }
        
        int[] aa = new int[a.length];
        int[] bb = new int[a.length];
                
        // use the order imposed by c to write the values of a into aa.  c holds
        // frequency, that is place markers too so that is updated as aa is written
        for (int i = (a.length - 1); i > -1; i--) {
            
            int aIdx = a[i];
                        
            int cfa = (int)c[aIdx];
            aa[cfa - 1] = aIdx;
            bb[cfa - 1] = b[i];
            c[aIdx]--;
        }

        System.arraycopy(aa, 0, a, 0, a.length);
        System.arraycopy(bb, 0, b, 0, b.length);
    }
    
    /**
     * apply a descending sort to the members of a with values less than max 
     * and apply the same changes of item position to b.
     * <pre>
     * Note that the numbers have to be positive, so if negative numbers are
     * in the array, the invoker needs to add a number to bring the values
     * to &gte; == 0 and then subtract that after the sort.
     * </pre>
     * @param a
     * @param b
     * @param maxValueInA max value in a 
     */
    public static void sortByDecr(int[] a, int[] b, int maxValueInA) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "the lengths of a and b must be the same");
        }
        if (maxValueInA <= 0) {
            throw new IllegalArgumentException("max must be > 0");
        }
        if (a.length < 2) {
            return;
        }
        
        long[] c = new long[maxValueInA + 1];

        // c holds frequency of each number by index, e.g. c[0] holds the number of 0's in a
        for (int i = 0; i < a.length; i++) {
            
            int idx = a[i];
            
            c[idx]++;
        }
                
        // cumulative sum to end of array c.  the last item in c holds the 
        // total number of items in 'a' less than or equal to max
        for (int i = 1; i < c.length; i++) {
            c[i] = c[i] + c[i - 1];
        }
        
        int[] aa = new int[a.length];
        int[] bb = new int[a.length];
                
        int n = a.length;
        
        // use the order imposed by c to write the values of a into aa.  c holds
        // frequency, that is place markers too so that is updated as aa is written
        for (int i = (n - 1); i > -1; i--) {
            
            int aIdx = a[i];
                        
            int cfa = (int)c[aIdx];
            
            int ii = n - (cfa - 1) - 1;
            
            aa[ii] = aIdx;
            bb[ii] = b[i];
            
            c[aIdx]--;
        }

        System.arraycopy(aa, 0, a, 0, a.length);
        System.arraycopy(bb, 0, b, 0, b.length);
    }

    public static void sort(int[] a, Object[] b, int maxA) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "the lengths of a and b must be the same");
        }
        if (maxA <= 0) {
            throw new IllegalArgumentException("maxA must be > 0");
        }
        if (a.length < 2) {
            return;
        }
        
        long[] c = new long[maxA + 1];

        // c holds frequency of each number by index, e.g. c[0] holds the number of 0's in a
        for (int i = 0; i < a.length; i++) {
            
            int idx = a[i];
            
            c[idx]++;
        }
                
        // cumulative sum to end of array c.  the last item in c holds the 
        // total number of items in 'a' less than or equal to max
        for (int i = 1; i < c.length; i++) {
            c[i] = c[i] + c[i - 1];
        }
        
        int[] aa = new int[a.length];
        Object[] bb = new Object[a.length];
                
        // use the order imposed by c to write the values of a into aa.  c holds
        // frequency, that is place markers too so that is updated as aa is written
        for (int i = (a.length - 1); i > -1; i--) {
            
            int aIdx = a[i];
                        
            int cfa = (int)c[aIdx];
            aa[cfa - 1] = aIdx;
            bb[cfa - 1] = b[i];
            c[aIdx]--;
        }

        System.arraycopy(aa, 0, a, 0, a.length);
        System.arraycopy(bb, 0, b, 0, b.length);
    }
}
