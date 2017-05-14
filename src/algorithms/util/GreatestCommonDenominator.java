package algorithms.util;

/**
 * implemented from pseudocode from Cormen et al. 
 * "Introduction to Algorithms", Chap 31
 *
 * @author nichole
 */
public class GreatestCommonDenominator {

    /**
     * return the greatest common denominator of the 2 integers.
     *
     * worse case runtime can never be greater than five times the number of its digits (base 10).
     * avg case runtime is
     * 
     * @param a
     * @param b
     * @return
     */
    public static int euclid(int a, int b) {
        if (b == 0) {
            return a;
        }
        count++;
        return euclid(b, a % b);
    }
    public static int count = 0;

    /**
     * return the greatest common denominator of the 2 integers.
     *
     * worse case runtime can never be greater than five times the 
     * O(log_2(min(m, n)).
     * 
     * @param a
     * @param b
     * @return
     */
    public static long euclid(long a, long b) {
        if (b == 0) {
            return a;
        }
        count++;
        return euclid(b, a % b);
    }
    
    /*
    r.t. complexity of multiplying 2 n bit numbers is
         O(n log n log log n).
    */

    /**
     * 
     * extended euclid r.t. O(size(a) · size(b))?
     * 
     * @param a
     * @param b
     * @return 
     */
    public static int[] extendedEuclid(int a, int b) {
        if (b == 0) {
            return new int[]{a, 1, 0};
        }
        int[] dxy_p = extendedEuclid(b, a % b);
        int[] dxy = new int[] {
            dxy_p[0], dxy_p[2], (dxy_p[1] - (a/b)*dxy_p[2])
        };
        return dxy;
    }
    
    public static long[] extendedEuclid(long a, long b) {
        if (b == 0) {
            return new long[]{a, 1, 0};
        }
        // r.t. complexity of '%' is O(a/b)
        long[] dxy_p = extendedEuclid(b, a % b);
        long[] dxy = new long[] {
            dxy_p[0], dxy_p[2], (dxy_p[1] - (a/b)*dxy_p[2])
        };
        return dxy;
    }
   
}
