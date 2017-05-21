package algorithms.util;

import algorithms.misc.MiscMath;

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
    
    /**
     * solves the equation a * x = b mod n to
     * find the smallest gcd for which a*x + b*y = d where d is a
     * gcd of number n.
     * @param a
     * @param b
     * @param n
     * @return 
     */
    public static int gcdModularLinearEqnSolver(int a, int b, int n) {
        
        /*
        https://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity
        Bézout's identity (also called Bézout's lemma) is a theorem in 
        elementary number theory: let a and b be nonzero integers and let d be 
        their greatest common divisor. Then there exist integers x and y such 
        that 
            ax+by=d.  
        In addition, the greatest common divisor d is the 
        smallest positive integer that can be written as ax + by every integer 
        of the form ax + by is a multiple of the greatest common divisor d.
        The integers x and y are called Bézout coefficients for (a, b); they 
        are not unique. A pair of Bézout coefficients can be computed by the 
        extended Euclidean algorithm.
        */
        
        int min = Integer.MAX_VALUE;
        
        int[] d_xp_yp = extendedEuclid(a, n);
        if ((d_xp_yp[0] != 0) || d_xp_yp[2] != 0) {
            int d = d_xp_yp[0];
            int x0 = d_xp_yp[1] * (b/d) % n;
            for (int i = 0; i < d; ++i) {
                
                int x1 = (x0 + i*(n/d)) % n;
                
                //System.out.println(" " + d);
                
                if (d > 0 && d < min) {
                    min = d;
                }
            }
        }
        return min;
    }
    
    public static int modularExponentiation(int a, int b, int n) {
        int c = 0;
        int d = 1;
        int nbits = MiscMath.numberOfBits(b);
        for (int i = nbits - 1; i >= 0; --i) {
            c *= 2;
            d = (d*d) % n;
            if ((b & (1 << i)) != 0) {
                c++;
                d = (d*a) % n;
            }
        }
        return d;
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
        
        System.out.format("a=%d b=%d (a/b)=%d d=%d x=%d y=%d\n", 
            a, b, (a/b), dxy_p[0], dxy_p[2], (dxy_p[1] - (a/b)*dxy_p[2]));
        
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
