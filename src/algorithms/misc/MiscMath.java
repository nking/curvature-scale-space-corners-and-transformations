package algorithms.misc;

import algorithms.util.PairInt;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MiscMath {

    public static int findPowerOf10_2(float a) {

        if (a == 0) {
            return 0;
        }

        int power = 0;

        if (a <= 1.0f) {
            while (a < 1.0) {
                a *=  10.0;
                power--;
            }
        } else {
            // precision errors in multiplication here are trouble for non base2 numbers such as powers of 10
            while (a >= 1.0) {
                a /= 10.0;
                power++;
            }
            power--;
        }

        return power;
    }

    public static int findPowerOf10(float a) {

        if (a == 0) {
            return 0;
        }
        if (a < 0.f) {
            a *= -1.0f;
        }
        double b = Math.log10(a);
        if (b > 0) {
            return (int)b;
        } else {
            if (b >= 1) {
                return (int)Math.round(b);
            } else if (b > -1) {
                // fractions between -1 and +1
                 return findPowerOf10_2(a);
            } else {
                return (int)Math.round(b);
            }
        }
    }

    public static int findYMinIndex(float[] a) {
        float min = Float.MAX_VALUE;
        int index = -1;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i])) {
                min = a[i];
                index = i;
            }
        }
        return index;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static float findMax(float[] a) {
        float max = Float.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) 
                && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
            }
        }
        return max;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findMax(int[] a) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) 
                && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
            }
        }
        return max;
    }

    public static float findMin(float[] a) {
        float min = Float.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i])) {
                min = a[i];
            }
        }
        return min;
    }

    public static int findMin(int[] a) {
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] < min) {
                min = a[i];
            }
        }
        return min;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(float[] a) {
        if (a == null) {
            return -1;
        }
        float max = Float.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) 
                && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(int[] a) {
        if (a == null) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && (a[i] < Integer.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static float[] calculateOuterRoundedMinAndMax(float[] a) {

        // find the powers of 10 for the data min and max
        float xmin = Float.MAX_VALUE;
        float xmax = Float.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > xmax) {
                xmax = a[i];
            }
            if (a[i] < xmin) {
                xmin = a[i];
            }
        }

        xmax = MiscMath.roundUpByLargestPower(xmax);

        xmin = MiscMath.roundDownByLargestPower(xmin);

        // xmax > 1 and xmin is between 0 and 1, round xmin down
        if ((xmax > 1) && (xmin > 0) && (xmin < 1.0)) {
            xmin = 0;
        }

        return new float[]{xmin, xmax};
    }

    /**
     * Round number down to in the largest power.
     * For example,
     *     roundDownByLargestPower(3.1) returns 3.0
     *     roundDownByLargestPower(-3.1) returns -4.0
     *     roundDownByLargestPower(31.1) returns 30.0
     *
     * @param f
     * @return
     */
    public static float roundDownByLargestPower(float f) {

        if (f == 0) {
            return 0;
        }

        int pow = findPowerOf10(f);
        double pow10 = Math.pow(10, pow);

        double r;
        if (f > 0) {
            double m = f % pow10;
            r = f - m;
        } else {
            int m = (int)(f/pow10);
            r = pow10 * (m - 1);
        }

        return (float)r;
    }

    /**
     Round number up to the next digit in the largest power.

        MiscMath.roundUpByLargestPower(31.1f) == 40.0f;

        MiscMath.roundUpByLargestPower(0.11f) == 1.0f;

        MiscMath.roundUpByLargestPower(-0.011f) == -0.02f;

        MiscMath.roundUpByLargestPower(310.1f) == 400.0f;

        MiscMath.roundUpByLargestPower(-3.1) == -4.0f;

        MiscMath.roundUpByLargestPower(3.1) == 4.0f;

    @param f
    @return
    */
    public static float roundUpByLargestPower(float f) {

        int pow = findPowerOf10(f);

        double pow10 = Math.pow(10, pow);

        int d = (int)(f/pow10);
        int m;
        if (f == pow10) {
            m = d;
        } else if (f > 0) {
            // residual ?
            if (f > (d*pow10)) {
                m = d + 1;
            } else {
                m = d;
            }
        } else {
            // decimals
            m = d - 1;
        }

        float r = (float) (m * pow10);

        return r;
    }
    
    /*
     http://en.wikipedia.org/wiki/Poisson_distribution
         algorithm poisson random number (Knuth):
     
     @param sr an instance of secure random which is a strong random number generator
     @param lambda a value less than 4 roughly looks like the tested dataset histograms
     */
    public static int poissonRandom(SecureRandom sr, int lambda) {

        double L = Math.exp(-lambda);

        int k = 0;

        double p = 1;

        do {
            k = k + 1;

            double u = sr.nextDouble();

            p = p * u;

        } while (p > L);

        return k - 1;
    }
  
    /**
     * use the taylor series to approximate the natural log of a number.
     * Note that internally, the variable x which is number - 1.
     * 
     * ln(1 + x) is computed, so number = 1 + x and therefore x = number - 1;
     * 
     * ln(1+x) = x - (x^2)/2 + (x^3)/3 - ...
     * 
     * @param number
     * @param n
     * @return
     */
    protected static double taylor(double number, int n) {
        double x = number - 1.0;
        
        double sum = 0;
        for (int i = 1; i <= n; i++) {
            double f = Math.pow(x, i)/(float)i;
            if (i % 2 == 0) {
                sum -= f;
                //System.out.print("     MINUS i=" + i);
            } else {
                sum += f;
                //System.out.print("     PLUS i=" + i);
            }            
        }
        
        return sum;
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(int[] x) {
        
        int n = x.length;
        long sumX = 0;
        for (int i = 0; i < n; i++) {
            sumX += x[i];
        }
        
        double avgX = (double)sumX/(double)n;
        
        sumX = 0;
        for (int i = 0; i < n; i++) {
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = (Math.sqrt(sumX/(n - 1.0f)));
        
        return new double[]{avgX, stdDevX};
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static float[] getAvgAndStDev(float[] x) {
        
        int n = x.length;
        double sumX = 0;
        for (int i = 0; i < n; i++) {
            sumX += x[i];
        }
        
        float avgX = (float)(sumX/(float)n);
        
        sumX = 0;
        for (int i = 0; i < n; i++) {
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        float stdDevX = (float)(Math.sqrt(sumX/(n - 1.0f)));
        
        return new float[]{avgX, stdDevX};
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(double[] x) {
        
        int n = x.length;
        double sumX = 0;
        for (int i = 0; i < n; i++) {
            sumX += x[i];
        }
        
        double avgX = sumX/(double)n;
        
        sumX = 0;
        for (int i = 0; i < n; i++) {
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = Math.sqrt(sumX/(n - 1.0f));
        
        return new double[]{avgX, stdDevX};
    }
    
    public static int[] add(int[] x, int amount) {
        
        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        
        int[] out = new int[x.length];
        
        int n = x.length;
        
        for (int i = 0; i < n; i++) {
            out[i] = x[i] + amount;
        }
        
        return out;
    }
    
    /**
     * compute n!/(n-k)!... needed for large numbers
     *
     * @param n
     * @param k
     * @return
     */
    public static long computeNDivNMinusK(int n, int k) {

        if (n == k) {
            return 1;
        }

        long result = 1;
        for (int i = n; i > (n-k); i--) {
            result *= i;
        }
        return result;
    }
    
    /**
     * compute n!/(n-k)!... needed for large numbers
     *
     * @param n
     * @return
     */
    public static long factorial(int n) {

        if (n < 3) {
            return n;
        }

        long result = 1;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
    
    /**
     * return the next subset after x, represented as a bitstring.
     * if x is null, the first subset is returned, else the subset
     * proceeding it.
     * The set bits are the selected bits in the returned string, that is
     * the set bits should represent the k subset bits in from nBits total
     * to select from.
     * 
     * One can quickly see the results using Long.toBinaryString(x).
     * 
     * @param nBits
     * @param kOnes
     * @param x
     * @return 
     */
    public static Long getNextSubsetBitstring(long nBits, long kOnes, Long x) {
        
        if (x == null) {
            
            long t = (1L << kOnes) - 1;
            
            return Long.valueOf(t);
        }
        
        long t = x.longValue();
        
        boolean proceed = ((t & (1L << nBits)) == 0);
        
        if (!proceed) {
            return null;
        }
        
        // 3 different ways to the same result:

        /*
        long lo = x & ~(x - 1);     // lowest one bit

        long lz = (x + lo) & ~x;    // lowest zero bit above lo

        x |= lz;                    // add lz to the set

        x &= ~(lz - 1);             // reset bits below lz

        x |= (lz / lo / 2) - 1;     // put back right number of bits at end
        */

        long a = t & -t;
        long b = t + a;
        t = (((t^b) >>> 2)/a) | b;

        /*
        t = b + (((b ^ t) / a) >> 2);
        */            
        
        return Long.valueOf(t);
    }
    
    /**
     * read the set bits from bitstring and return them as indexes in 
     * outputIndexes.
     * @param bitstring
     * @param outputIndexes 
     */
    public static void readSetBits(Long bitstring, List<Integer> outputIndexes) {
        
        if (bitstring == null) {
            throw new IllegalArgumentException("bitstring cannot be null");
        }
        if (outputIndexes == null) {
            throw new IllegalArgumentException("outputIndexes cannot be null");
        }
        
        long b = bitstring.longValue();
        
        int count = 0;
        
        while (b > 0) {
            
            long idx = (b & 1L);
            
            if (idx == 1) {
                outputIndexes.add(Integer.valueOf(count));
            }
            
            b >>= 1L;
            
            count++;
        }
    }
    
     /**
      * solve for the roots of equation a0 * x^3 + a1 * x^2 + a2 * x + a4 = 0;
      * 
      * most of the method is adapted from 
      * http://www.csse.uwa.edu.au/~pk/research/matlabfns/Misc/cubicroots.m

        Copyright (c) 2008 Peter Kovesi
        School of Computer Science & Software Engineering
        The University of Western Australia
        pk at csse uwa edu au
        http://www.csse.uwa.edu.au/

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in 
        all copies or substantial portions of the Software.

        The Software is provided "as is", without warranty of any kind.

        Nov 2008
         
      * @param a
      * @param b
      * @param c
      * @param d
      * @return 
      */
    public static double[] solveCubicRoots(double a, double b, double c, 
        double d) {
        
        b /= a; 
        c /= a; 
        d /= a;
        double bOn3 = b/3.;
    
        double q = (3.*c - b*b)/9.;
        double r = (9.*b*c - 27.*d - 2.*b*b*b)/54.;
        double discriminant = q*q*q + r*r;
            
        if (discriminant >= 0) {
            
            // discriminant > 0: one root is real and two complex conjugate
            // discriminant == 0 : all roots are real and at least two are equal
            
            double s = realcuberoot(r + Math.sqrt(discriminant));   
            double t = realcuberoot(r - Math.sqrt(discriminant));
    
            double root = s + t - bOn3;
        
            return new double[]{root};
            
        } else {

            //  all roots are real and unequal
            
            double rho = Math.sqrt(r*r - discriminant); 
            double cubeRootrho = realcuberoot(rho);
            double thetaOn3 = Math.acos(r/rho)/3.;
        
            double crRhoCosThetaOn3 = cubeRootrho * Math.cos(thetaOn3);
            double crRhoSinThetaOn3 = cubeRootrho * Math.sin(thetaOn3);   

            double[] root = new double[3];
            root[0] = 2. * crRhoCosThetaOn3 - bOn3;
            root[1] = -crRhoCosThetaOn3 - bOn3 - Math.sqrt(3.)*crRhoSinThetaOn3;
            root[2] = -crRhoCosThetaOn3 - bOn3 + Math.sqrt(3.)*crRhoSinThetaOn3;
            
            return root;
        }   
    }
   
    public static double realcuberoot(double x) {
        double sign = (x < 0.) ? -1 : 1;
        double y = sign * Math.pow(Math.abs(x), 1./3.);
        return y;
    }
    
    /**
     * choose selected.length unique random numbers between 0 and nMax and
     * populate the array 'selected' with them.
     * 
     * Note that the number of ways to make a subset of 7 distinct numbers
     * out of nMax numbers is nMax!/(nMax-7)!
     * So the probability of selecting the same 7 randomly from nMax is very
     * low when nMax is much larger than selected.length, so this method does
     * not check that the same 7 were not drawn before, but the invoker 
     * should for small nMax.
     * 
     * @param sr class to generate semi random numbers
     * @param selected the output array to populate with unique random numbers 
     * valued 0 to nMax, exclusive.
     * @param nMax the upper limit from which to choose numbers from for 
     * 'selected', where the possible numbers are 0 through nMax-1.
     */
    public static void chooseRandomly(SecureRandom sr, int[] selected, int nMax) {
        
        if (nMax < selected.length) {
            throw new IllegalArgumentException("cannot draw " + 
                Integer.toString(selected.length) + " distinct random numbers "
                + " from only " + Integer.toString(nMax) + " numbers");
        }
        
        // TODO: ideally, would like to be able to predict the ith iteration of 
        // a subset chooser (@see getNextSubsetBitstring)
      
        if (selected.length == nMax) {
            for (int i = 0; i < selected.length; i++) {
                selected[i] = i;
            }
            return;
        }
        
        int nLimit = selected.length * 3;
        // when nMax is smaller than some limit, choose randomly from 
        // unchosen numbers
        if (nMax < nLimit) {
            
            // populate numbers to choose from:
            List<Integer> numbers = new ArrayList<Integer>();
            for (int i = 0; i < nMax; i++) {
                numbers.add(Integer.valueOf(i));
            }
            
            for (int i = 0; i < selected.length; i++) {
                int selIdx = sr.nextInt(numbers.size());
                Integer sel = numbers.get(selIdx);
                selected[i] = sel.intValue();
                numbers.remove(sel);
            }
            
            return;
        }
        
        for (int i = 0; i < selected.length; i++) {
            int sel = sr.nextInt(nMax);
            while (contains(selected, i, sel)) {
                sel = sr.nextInt(nMax);
            }
            selected[i] = sel;
        }
    }
    
    private static boolean contains(int[] values, int lastIdx, int valueToCheck) {
        for (int i = 0; i < lastIdx; i++) {
            if (values[i] == valueToCheck) {
                return true;
            }
        }
        return false;
    }
    
     /**
     * for the given histogram, returns the indexes of the primary peak and
     * any peaks which are larger than frac*maxPeak above their neighboring
     * values.
     * @param h
     * @param frac the fraction of the y value of the maximum peak that is used
     * as a critical limit that any other peaks must have in excess of their
     * neighboring points.  For example, 0.1.
     * @return 
     */
    public static List<Integer> findStrongPeakIndexes(HistogramHolder h, float frac) {
        
        float[] x = h.getXHist();
        int[] y = h.getYHist();
        
        if (x == null || y == null || x.length == 0 || y.length == 0) {
            return null;
        }
        
        int yPeakIdx = MiscMath.findYMaxIndex(y);
        
        if (yPeakIdx == -1) {
            return null;
        }
        
        int yMaxPeak = y[yPeakIdx];
        
        /*
        storing the minima and maxima in the same array list.
        the minima have -1*index within k
        and the maxima keep their positive values of the index within k.
        */
        List<Integer> minMaxIndexes = new ArrayList<Integer>();
        
        float lastY = y[0];
        boolean incr = true;
        for (int ii = 1; ii < y.length; ii++) {
            if ((y[ii] < lastY) && incr) {
                minMaxIndexes.add(Integer.valueOf(ii - 1));
                incr = false;
            } else if ((y[ii] > lastY) && !incr) {
                minMaxIndexes.add(Integer.valueOf(-1*(ii - 1)));
                incr = true;
            }
            lastY = y[ii];
        }
        
        // for the histograms of euclidean point combination differences,
        // this should usually be singly peaked
        if (minMaxIndexes.size() == 1) {
            return minMaxIndexes;
        }
        
        if (incr) {
            // add the last point
             minMaxIndexes.add(Integer.valueOf(y.length - 1));
        }
        
        float limit = frac * yMaxPeak;
        
        // find peaks where y[ii] is > limit above adjacent local minima
        
        List<Integer> peaks = new ArrayList<Integer>();

        for (int ii = 0; ii < minMaxIndexes.size(); ii++) {

            int idx = minMaxIndexes.get(ii).intValue();

            if (idx > -1) {
                // this is maximum
                
                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = y[-1*idx2];
                        float diff = y[idx] - compare;
                        if (diff >= limit) {
                            peaks.add(Integer.valueOf(idx));
                            found = true;
                        }
                        break;
                    }
                }
                if (found) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIndexes.size(); iii++) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = y[-1*idx2];
                        float diff = y[idx] - compare;
                        if (diff >= limit) {
                            peaks.add(Integer.valueOf(idx));
                        }
                        break;
                    }
                }
            }
        }
       
        return peaks;
    }

    public static float[] extractAllXForYAboveHalfMax(HistogramHolder h) {
        
        if (h == null) {
            return null;
        }
        
        int peakIdx = MiscMath.findYMaxIndex(h.getYHist());
        int peakY = h.getYHist()[peakIdx];
        int crit = peakY/2;
        
        int i0 = -1;
        int i1 = -1;
        
        for (int i = 0; i < h.getYHist().length; i++) {
            int y = h.getYHist()[i];
            if (y > crit) {
                if (i0 == -1) {
                    i0 = i;
                }
                i1 = i;
            }
        }
        
        int n = i1 - i0 + 1;
        
        float[] xs = new float[n];
        int count = 0;
        for (int i = i0; i <= i1; i++) {
            xs[count] = h.getXHist()[i];
            count++;
        }
        
        return xs;
    }

    public static int findLastZeroIndex(HistogramHolder h) {
        int n = h.getXHist().length;
        int lastZeroIdx = n - 1;
        for (int i = (n - 1); i > -1; i--) {
            int y = h.getYHist()[i];
            if (y > 0) {
                break;
            }
            lastZeroIdx = i;
        }
        return lastZeroIdx;
    }

    public static int[] findMinMaxXY(Set<PairInt> points) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }

    public static int numberOfBits(long i) {
        
        if (i == 0) {
            return 1;
        }
        
        if (i < 0) {
            i *= -1;
        }
        
        int count = 0;
        while (i > 0) {
            i >>= 1L;
            count++;
        }
        return count;
        
    }
    
    /**
     * 
     * @param v
     * @return 
     */
    public static int bitReverse(int v, int nBits) {

        int r = v;
                
        int rev = 0;
        for (int i = 0; i < nBits; i++) {
            rev = (rev << 1) | (r & 1);
            r >>= 1;
        }
        
        return rev;
        
    }
    
    public static boolean isAPowerOf2(int n) {
        // n XOR n-1
        return ((n == 0) || ((n & (n - 1)) == 0));
    }
    
    /**
     * calculate the angle of theta with respect to an origin of (0, 0) for the
     * point (x, y).
     * The angles are
     * <pre>    90(=pi/2)
     *            |
     *            |
     *   180 ----------- 0
     *   (=pi)    |
     *            | 
     *          270(=3pi/2)
     * </pre>
     * 
     * @param x
     * @param y
     * @return degrees in radians
     */
    public static double calculatePolarTheta(float x, float y) {
        
        if (x == 0) {
            if (y < 0) {
                return 1.5 * Math.PI;
            }
            return Math.PI/2.;
        }
        
        double div = y/x;
        
        double theta = Math.atan(div);
        
        if (x > 0) {
            // if y > 0, no change needed
            if (y == 0) {
                theta = 0;
            } else if (y < 0) {
                theta = 2*Math.PI + theta;
            }
        } else {
            // x < 0
            if (y == 0) {
                theta = Math.PI;
            } else {
                theta = Math.PI + theta;
            }
        }
        
        return theta;
    }

    /**
     * write an array of the rotation angles that start with rotStart and
     * increment by rotDelta until rotStop, inclusive.
     * 
     * @param rotStart start of rotation angles in units of degrees
     * @param rotStop stop of rotation angles in units of degrees
     * @param rotDelta interval between rotation angles in units of degrees
     * @return 
     */
    public static int[] writeDegreeIntervals(int rotStart, int rotStop, 
        int rotDelta) {
        
        int rotStopMinusStart = (rotStop < rotStart) ? 
            ((360 - rotStart) + rotStop) : (rotStop - rotStart);
        
        int nRot = (rotStopMinusStart/rotDelta) + 1;
        int[] rotation = new int[nRot];
        for (int i = 0; i < nRot; ++i) {
            float r = rotStart + (i * rotDelta);
            if (r > 359) {
                r -= 360;
            }
            rotation[i] = (int)r;
        }
        
        return rotation;
    }
    
}
