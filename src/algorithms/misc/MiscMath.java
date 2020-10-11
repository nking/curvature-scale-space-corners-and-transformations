package algorithms.misc;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.FurthestPair;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceContour;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.matching.ORBMatcher;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.set.TIntSet;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MiscMath {

    public static int findPowerOf10_2(float a) {
        return MiscMath0.findPowerOf10_2(a);
    }

    public static int findPowerOf10(float a) {
        return MiscMath0.findPowerOf10(a);
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
        return MiscMath0.findMax(a);
    }
    
    public static double findMax(double[] a) {
        double max = Double.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }
    
    public static float findMax(float[] a, int nIndexes) {
        float max = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < nIndexes; i++) {
            if (a[i] > max) {
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
        return MiscMath0.findMax(a);
    }
    
    public static float findMax(float[][] img) {
        
        return MiscMath0.findMax(img);
    }
    
    public static int findMax(int[][] img) {
        return MiscMath0.findMax(img);
    }

    public static int findMax(int[] a, int nIndexes) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < nIndexes; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) 
                && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
            }
        }
        return max;
    }

    public static float findMin(float[][] img) {
        
        return MiscMath0.findMin(img);
    }
    
    public static double findMin(double[][] img) {
        
        return MiscMath0.findMin(img);
    }
    
    public static double findMax(double[][] img) {
        
        return MiscMath0.findMax(img);
    }
    
    public static int findMin(int[][] img) {
        return MiscMath0.findMin(img);
    }

    public static float findMin(float[] a) {
        return MiscMath0.findMin(a);
    }
    
    public static float findMin(float[] a, int nIndexes) {
        float min = Float.POSITIVE_INFINITY;
        for (int i = 0; i < nIndexes; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i])) {
                min = a[i];
            }
        }
        return min;
    }

    public static int findMin(int[] a) {
        return MiscMath0.findMin(a);
    }
    
    public static int findMin(int[] a, int nIndexes) {
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < nIndexes; i++) {
            if (a[i] < min) {
                min = a[i];
            }
        }
        return min;
    }
    
    public static int findMin(GreyscaleImage img, Set<PairInt> points) {
        int min = Integer.MAX_VALUE;
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            int v = img.getValue(x, y);
            if (v < min) {
                min = v;
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
        return MiscMath0.findYMaxIndex(a);
    }

    public static int findYMaxIndex(int[] a) {
        return MiscMath0.findYMaxIndex(a);
    }
    
    public static int findYMaxIndex(long[] a) {
        return MiscMath0.findYMaxIndex(a);
    }
    
    public static int[] findYMaxIndexes(int[] a) {
        if (a == null || a.length == 0) {
            return null;
        }
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && (a[i] < Integer.MAX_VALUE)) {
                max = a[i];
            }
        }
        if (max == Integer.MIN_VALUE) {
            return null;
        }
        int count = 0;
        for (int i = 0; i < a.length; i++) {
            if (a[i] == max) {
                count++;
            }
        }
        int[] indexes = new int[count];
        count = 0;
        for (int i = 0; i < a.length; i++) {
            if (a[i] == max) {
                indexes[count] = i;
                count++;
            }
        }
        return indexes;
    }
    
    public static float[] calculateOuterRoundedMinAndMax(float[] a) {
        return MiscMath0.calculateOuterRoundedMinAndMax(a);
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
        return MiscMath0.roundDownByLargestPower(f);
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
        return MiscMath0.roundUpByLargestPower(f);
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
        
        return MiscMath0.getAvgAndStDev(x, x.length);
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(int[] x, int length) {
        return MiscMath0.getAvgAndStDev(x, length);
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(List<Double> x) {
        
        int length = x.size();
        
        double sumX = 0;
        for (int i = 0; i < length; i++) {
            sumX += x.get(i);
        }
        
        double avgX = sumX/(double)length;
        
        sumX = 0;
        for (int i = 0; i < length; i++) {
            double diffX = x.get(i) - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = (Math.sqrt(sumX/((double)length - 1.)));
        
        return new double[]{avgX, stdDevX};
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return float[]{avg, stDev}
     */
    public static float[] getAvgAndStDev(float[] x) {        
        return getAvgAndStDev(x, x.length);
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @param length the number of indexes to use in x.  x can be longer than
     * length
     * @return float[]{avg, stDev}
     */
    public static float[] getAvgAndStDev(float[] x, int length) {
        
        int n = length;
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
        float stdDevX = (float)(Math.sqrt(sumX/((float)n - 1.0f)));
        
        return new float[]{avgX, stdDevX};
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @param length the number of indexes to use in x.  x can be longer than
     * length
     * @param sentinel value to flag that an item should not be included in
     * calculations.
     * @return float[]{avg, stDev}
     */
    public static float[] getAvgAndStDevIgnoreForSentinel(int[] x, int length, 
        final int sentinel) {
        
        int n = length;
        double sumX = 0;
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (x[i] == sentinel) {
                continue;
            }
            sumX += x[i];
            count++;
        }
        
        float avgX = (float)(sumX/(float)count);
        
        sumX = 0;
        for (int i = 0; i < n; i++) {
            if (x[i] == sentinel) {
                continue;
            }
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        float stdDevX = (float)(Math.sqrt(sumX/(count - 1.0f)));
        
        return new float[]{avgX, stdDevX};
    }
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @param length the number of indexes to use in x.  x can be longer than
     * length
     * @param sentinel value to flag that an item should not be included in
     * calculations.
     * @return float[]{avg, stDev}
     */
    public static float[] getAvgAndStDevIgnoreForSentinel(float[] x, int length, 
        final float sentinel) {
        
        int n = length;
        double sumX = 0;
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (x[i] == sentinel) {
                continue;
            }
            sumX += x[i];
            count++;
        }
        
        float avgX = (float)(sumX/(float)count);
        
        sumX = 0;
        for (int i = 0; i < n; i++) {
            if (x[i] == sentinel) {
                continue;
            }
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        float stdDevX = (float)(Math.sqrt(sumX/(count - 1.0f)));
        
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
    
    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @return 
     */
    public static double[] getAvgAndStDev(GreyscaleImage img, Set<PairInt> points) {
        
        int n = points.size();
        double sum = 0;
        for (PairInt p : points) {
            sum += img.getValue(p.getX(), p.getY());
        }
        
        double avg = sum/(double)n;
        
        sum = 0;
        for (PairInt p : points) {
            int v = img.getValue(p.getX(), p.getY());
            double diffX = v - avg;
            sum += (diffX * diffX);
        }
        double stdDevX = Math.sqrt(sum/((double)n - 1.0));
        
        return new double[]{avg, stdDevX};
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

        return MiscMath0.computeNDivNMinusK(n, k);
    }
    
    /**
     * compute n!/k!(n-k)!.  Note that if n or k are larger than 12,
     * computeNDivKTimesNMinusKBigInteger is used and in that case,
     * if the result is larger than Long.MAX_VALUE, LONG.MAX_VALUE is
     * returned.
     *
     * @param n
     * @param k
     * @return
     */
    public static long computeNDivKTimesNMinusK(int n, int k) {

        if (n == k) {
            return 1;
        }
        
        if (k > 12 || n > 12) {
            return computeNDivKTimesNMinusKBigInteger(n, k);
        }

        double result = 1;
        for (int i = n; i > (n-k); i--) {
            result *= i;
        }
        double divisor = factorial(k);
        
        result = result/divisor;
        
        return Math.round(result);
    }
    
    /**
     * compute n!/k!(n-k)!.  Note that if n or k are larger than 12,
     * computeNDivKTimesNMinusKBigIntegerExact is used and in that case,
     * if the result is larger than Long.MAX_VALUE an exception is thrown.
     *
     * @param n
     * @param k
     * @return
     * @throws ArithmeticException thrown when result is out of range of type long
     */
    public static long computeNDivKTimesNMinusKExact(int n, int k) {

        return MiscMath0.computeNDivKTimesNMinusKExact(n, k);
    }
    
    /**
     * compute n!
     *
     * @param n
     * @return
     */
    public static long factorial(int n) {

        return MiscMath0.factorial(n);
    }
    
    /**
     * compute n!
     *
     * @param n
     * @return
     */
    public static BigInteger factorialBigInteger(int n) {

        BigInteger result = BigInteger.ONE;

        for (int i = 2; i <= n; i++) {
            
            byte[] bytes = MiscMath.writeToBigEndianBytes(i);
            
            BigInteger v = new BigInteger(bytes);
            
            result = result.multiply(v);
        }
        
        return result;
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
        return MiscMath0.findStrongPeakIndexes(h, frac);
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
    
    public static int findFirstNonZeroIndex(HistogramHolder h) {
        int n = h.getXHist().length;
        int firstNonZero = -1;
        for (int i = 0; i < n; i++) {
            int y = h.getYHist()[i];
            if (y > 0) {
                firstNonZero = i;
                break;
            }
        }
        return firstNonZero;
    }
    
    public static int findFirstNonZeroIndex(int[] h) {
        int n = h.length;
        int firstNonZero = -1;
        for (int i = 0; i < n; i++) {
            int y = h[i];
            if (y > 0) {
                firstNonZero = i;
                break;
            }
        }
        return firstNonZero;
    }
    
    public static int findLastNonZeroIndex(HistogramHolder h) {
        int n = h.getXHist().length;
        int idx = -1;
        for (int i = (n - 1); i > -1; --i) {
            int y = h.getYHist()[i];
            if (y > 0) {
                idx = i;
                break;
            }
        }
        return idx;
    }

    public static int findLastNonZeroIndex(int[] h) {
        int n = h.length;
        int idx = -1;
        for (int i = (n - 1); i > -1; --i) {
            int y = h[i];
            if (y > 0) {
                idx = i;
                break;
            }
        }
        return idx;
    }
    
    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(Collection<PairInt> points) {
        return MiscMath0.findMinMaxXY(points);
    }
    
    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(TIntSet pixelIdxs, int imgWidth) {
        return MiscMath0.findMinMaxXY(pixelIdxs, imgWidth);
    }
    
    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(PairIntArray points) {
        return MiscMath0.findMinMaxXY(points);
    }
    
    /**
     * determine the number of bits, that is, the msb position + 1.
     * Note that a value of 0 returns a bit length of 1.
     * @param v
     * @return 
     */
    public static int numberOfBits(int v) {
        
        return MiscMath0.numberOfBits(v);
    }
  
    /**
     * determine the number of bits, that is the msb position + 1.
     * Note that a value of 0 returns a bit length of 1.
     * @param v
     * @return 
     */
    public static int numberOfBits(long v) {
        
        return MiscMath0.numberOfBits(v);
    }
        
    /**
     * 
     * @param v
     * @return 
     */
    public static int bitReverse(int v, int nBits) {
        return MiscMath0.bitReverse(v, nBits);
    }
    
    public static boolean isAPowerOf2(int n) {
        return MiscMath0.isAPowerOf2(n);
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
     * @return theta in radians
     */
    public static double calculatePolarTheta(double x, double y) {
       
        double theta = Math.atan2(y, x);
        
        if (Double.isFinite(theta) && (theta < 0)) {
            theta += 2. * Math.PI;
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
    public static float[] writeDegreeIntervals(float rotStart, float rotStop, 
        float rotDelta) {
        
        float rotStopMinusStart = (rotStop < rotStart) ? 
            ((360 - rotStart) + rotStop) : (rotStop - rotStart);
        
        int nRot = (int)((rotStopMinusStart/rotDelta) + 1);
        int np = 0;
        if (((rotStopMinusStart + 1) % rotDelta) != 0) {
            np++;
        }
        float[] rotation = new float[nRot + np];
        for (int i = 0; i < nRot; ++i) {
            float r = rotStart + (i * rotDelta);
            if (r > 359) {
                r -= 360;
            }
            rotation[i] = (int)r;
        }
        if (np > 0) {
            rotation[rotation.length - 1] = rotStop;
        }
        
        return rotation;
    }

    public static List<Integer> findStrongPeakIndexesDescSort(
        HistogramHolder hist, float fracMax) {

        return MiscMath0.findStrongPeakIndexesDescSort(hist, fracMax);
    }
    
    /**
     * write value to a byte array in big endian, that is LSB in highest order bit
     * (MSB is in lowest memory address).
     * these are signed values stored as twos complement and can be input
     * to BigInteger's constructor.
     * @param value
     * @return 
     */
    public static byte[] writeToBigEndianBytes(long value) {
    
        return MiscMath0.writeToBigEndianBytes(value);
    }
   
    /**
     * create a random shuffle of the source array.
     * It uses wikipedia pseudo code for the In/Out method.
     * (add reference here).
     * @param source
     * @return 
     */
    public static PairIntArray shuffle(PairIntArray source) throws NoSuchAlgorithmException {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        PairIntArray out = new PairIntArray();
        
        int n = source.getN();
        
        for (int i = 0; i < n; ++i) {
             
            int j = (i == 0) ? 0 : sr.nextInt(out.getN());
            
            if (j == out.getN()) {
                out.add(source.getX(i), source.getY(i));
            } else {
                out.add(out.getX(j), out.getY(j));
                out.set(j, source.getX(i), source.getY(i));
            }
        }
        
        return out;
    }
    
    public static PairIntArray get20NeighborOffsets() {
        
        return MiscMath0.get20NeighborOffsets();
    }
    
    public static int calculateSumXOR(int[] d1, int[] d2) {
        
        assert(d1.length == d2.length);
        
        int sum = 0;
        
        for (int i = 0; i < d1.length; ++i) {
            sum += (d1[i] ^ d2[i]);
        }
                
        return sum;
    }
    
    public static float calculateSSD(int[] a, int[] b, int sentinel) {
        
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "a and b must have the same lengths");
        }
                
        int count = 0;
        double sum = 0;
        
        for (int i = 0; i < a.length; ++i) {
            
            int v1 = a[i];
            if (v1 == sentinel) {
                continue;
            }
            int v2 = b[i];
            if (v2 == sentinel) {
                continue;
            }
            
            int diff = v1 - v2;
            
            sum += (diff * diff);
            
            count++;
        }
        sum /= (double)count;
                
        return (float)sum;
    }
    
    public static float calculateAngular360SSD(int[] a, int[] b, int sentinel) {
        
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "a and b must have the same lengths");
        }
                
        int count = 0;
        
        double sum = 0;
        
        for (int i = 0; i < a.length; ++i) {
            
            int v1 = a[i];
            if (v1 == sentinel) {
                continue;
            }
            int v2 = b[i];
            if (v2 == sentinel) {
                continue;
            }
            
            float diff = AngleUtil.getAngleDifference(v1, v2);
            
            sum += (diff * diff);
            
            count++;
        }
        
        sum /= (double)count;
                
        return (float)sum;
    }
   
    public static float calculateSSD(float[] a, float[] b, float sentinel) {
        
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "a and b must have the same lengths");
        }
                
        int count = 0;
        double sum = 0;
        
        for (int i = 0; i < a.length; ++i) {
            
            float v1 = a[i];
            if (v1 == sentinel) {
                continue;
            }
            float v2 = b[i];
            if (v2 == sentinel) {
                continue;
            }
            
            float diff = v1 - v2;
            
            sum += (diff * diff);
            
            count++;
        }
        sum /= (double)count;
                
        return (float)sum;
    }
    
    /**
     Calculate the cosine similarity between 2 descriptors.
     Calculates a dot b/(magnitude(a) * magnitude(b).
     The equation and use are from the paper
     "Generalized RANSAC framework for relaxed correspondence problems"
     by Zhang and Kosecka.
     The authors use the similarity with a threshold of 0.95 to filter out 
     dissimilar matches.
     * @param a
     * @param b
     * @param sentinel
     * @return 
     */
    public static float calculateCosineSimilarity(float[] a, float[] b, float sentinel) {
        
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "a and b must have the same lengths");
        }
        
        /*
                       a[0]*b[0] + a[1]*b[1] + a[2]*b[2] ...
        ------------------------------------------------------------------------------
        sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]...) * sqrt(b[0]*b[0] + b[1]*b[1] + ...)
        */     
        int count = 0;
        double aDotB = 0;
        double aMagn = 0;
        double bMagn = 0;
        for (int i = 0; i < a.length; ++i) {
            float a1 = a[i];
            float b1 = a[i];
            if (a1 == sentinel || b1 == sentinel) {
                continue;
            }
            aDotB += (a1 * b1); 
            aMagn += (a1 * a1);
            bMagn += (b1 * b1);
            count++;
        }
        aMagn = Math.sqrt(aMagn);
        bMagn = Math.sqrt(bMagn);
        double cosSim = aDotB/(aMagn*bMagn);
        
        return (float)cosSim;
    }
    
    /**
     Calculate the cosine similarity between 2 descriptors.
     Calculates a dot b/(magnitude(a) * magnitude(b).
     The equation and use are from the paper
     "Generalized RANSAC framework for relaxed correspondence problems"
     by Zhang and Kosecka.
     The authors use the similarity with a threshold of 0.95 to filter out 
     dissimilar matches.
     * @param a
     * @param b
     * @param sentinel
     * @return 
     */
    public static float calculateCosineSimilarity(int[] a, int[] b, int sentinel) {
        
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "a and b must have the same lengths");
        }
        
        /*
                       a[0]*b[0] + a[1]*b[1] + a[2]*b[2] ...
        ------------------------------------------------------------------------------
        sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]...) * sqrt(b[0]*b[0] + b[1]*b[1] + ...)
        */     
        int count = 0;
        double aDotB = 0;
        double aMagn = 0;
        double bMagn = 0;
        for (int i = 0; i < a.length; ++i) {
            float a1 = a[i];
            float b1 = a[i];
            if (a1 == sentinel || b1 == sentinel) {
                continue;
            }
            aDotB += (a1 * b1); 
            aMagn += (a1 * a1);
            bMagn += (b1 * b1);
            count++;
        }
        aMagn = Math.sqrt(aMagn);
        bMagn = Math.sqrt(bMagn);
        double cosSim = aDotB/(aMagn*bMagn);
        
        return (float)cosSim;
    }
    
    /**
     * Determine the sum squared error within this array using 
     * auto-correlation and the assumption that the value at the middle index 
     * is the value from the original central pixel.  Values the same as the
     * sentinel are ignored and not included in the calculation.
     * @return 
     */
    public static float sumSquaredError(int[] a, int sentinel, int centralIdx) {
        
        int n = a.length;
        
        int vc = a[centralIdx];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        int count = 0;
        
        double sum = 0;
        
        for (int i = 0; i < a.length; ++i) {
            
            int v1 = a[i];
            if (v1 == sentinel) {
                continue;
            }
            
            int diff = a[i] - vc;
        
            sum += (diff * diff);
            count++;
        }
        sum /= (double)count;
                       
        return (float)sum;
    }

    /**
     * Determine the sum squared error within this array using 
     * auto-correlation and the assumption that the value at the middle index 
     * is the value from the original central pixel.  Values the same as the
     * sentinel are ignored and not included in the calculation.
     * @return 
     */
    public static float sumSquaredAngular360Error(int[] a, int sentinel, 
        int centralIdx) {
        
        int n = a.length;
        
        int vc = a[centralIdx];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        int count = 0;
        
        double sum = 0;
        
        for (int i = 0; i < a.length; ++i) {
            
            int v1 = a[i];
            if (v1 == sentinel) {
                continue;
            }
            
            float diff = AngleUtil.getAngleDifference(v1, vc);
        
            sum += (diff * diff);
            
            count++;
        }
        
        sum /= (double)count;
                       
        return (float)sum;
    }

    /**
     * Determine the sum squared error within this array using 
     * auto-correlation and a self reference value of the item at
     * centralPixIdx in a.  Any item in a having value sentinel is ignored
     * and not included in the calculation.
     * 
     * @return 
     */
    public static float sumSquaredError(float[] a, float sentinel, int centralPixIdx) {
        
        int n = a.length;
        
        float vc = a[centralPixIdx];
        
        if (vc == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        
        int count = 0;
        
        double sum = 0;
        
        for (int i = 0; i < n; ++i) {
            
            float v1 = a[i];
            if (v1 == sentinel) {
                continue;
            }
            
            float diff = a[i] - vc;
        
            sum += (diff * diff);
            
            count++;
        }
        
        sum /= (double)count;
                       
        return (float)sum;
    }

    /**
     <pre>
     order points such that:
        q ≺ p : (i_q < i_p) ∨ ((i_q == i_p) ∧(j_q < j_p))

          (-1, 1)
          (-1, 0)   p=(0,  0)
          (-1,-1)     (0, -1)
     </pre>
     * runtime complexity is 2*O(N_points) + O(N_points * lg2(O(N_points)))
     * @param points
     * @return 
     */
    public static LinkedHashSet<PairInt> lexicographicallyOrderPointsBySort(
        Set<PairInt> points) {
        
        int n = points.size();
        
        LinkedHashSet<PairInt> ordered = new LinkedHashSet<PairInt>(n);
        
        int[] x = new int[n];
        int[] y = new int[n];
        
        // O(N_points)
        int count = 0;
        for (PairInt p : points) {
            x[count] = p.getX();
            y[count] = p.getY();
            count++;
        }
        
        // O(N_points * lg2(O(N_points)))
        MultiArrayMergeSort.sortBy1stArgThen2nd(x, y);
        
        // O(N_points)
        for (int i = 0; i < n; ++i) {
            
            PairInt p = new PairInt(x[i], y[i]);
            
            ordered.add(p);
        }
        
        return ordered;
    }
    
    /**
     <pre>
     order points such that:
        q ≺ p : (i_q < i_p) ∨ ((i_q == i_p) ∧(j_q < j_p))

          (-1, 1)
          (-1, 0)   p=(0,  0)
          (-1,-1)     (0, -1)
     </pre>
     * runtime complexity is O(w*h)
     * @param points
     * @param w
     * @param h
     * @return 
     */
    public static LinkedHashSet<PairInt> lexicographicallyOrderPointsByScan(
        Set<PairInt> points, int w, int h) {
        
        LinkedHashSet<PairInt> ordered = new LinkedHashSet<PairInt>(w * h);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                
                PairInt p = new PairInt(i, j);
                
                if (points.contains(p)) {
                    ordered.add(p);
                }
            }
        }
        
        return ordered;
    }

    public static long findMaxForByteCompressed(long[] aL, int len, 
        int itemNDatum, int datumNBits, int minDatum) {
        
        long mask = (1L << datumNBits) - 1L;
        
        long max = Long.MIN_VALUE;
        
        for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
            long total = aL[elementIdx];
            for (int i = 0; i < itemNDatum; i++) {
                long v = (total >> (long)(i * datumNBits)) & mask;
                v += minDatum;
                if (v > max) {
                    max = v;
                }
            }
        }
        
        return max;
    }
    
    public static long findMinForByteCompressed(long[] aL, int len, 
        int itemNDatum, int datumNBits, int minDatum) {
                
        long mask = (1L << datumNBits) - 1L;
        
        long min = Long.MAX_VALUE;
        
        for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
            long total = aL[elementIdx];
            for (int i = 0; i < itemNDatum; i++) {
                long v = (total >> (long)(i * datumNBits)) & mask;
                v += minDatum;
                if (v < min) {
                    min = v;
                }
            }
        }
        
        return min;
    }

    public static int findMaxForByteCompressed(int[] a, int len, int itemNDatum, 
        int datumNBits, int minDatum) {
        
        int mask = (1 << datumNBits) - 1;
        
        int max = Integer.MIN_VALUE;
        
        for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
            int total = a[elementIdx];
            for (int i = 0; i < itemNDatum; i++) {
                int v = (total >> (i * datumNBits)) & mask;
                v += minDatum;
                if (v > max) {
                    max = v;
                }
            }
        }
        
        return max;
    }
    
    public static int findMinForByteCompressed(int[] a, int len, int itemNDatum, 
        int datumNBits, int minDatum) {
        
        int mask = (1 << datumNBits) - 1;
        
        int min = Integer.MAX_VALUE;
        
        for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
            int total = a[elementIdx];
            for (int i = 0; i < itemNDatum; i++) {
                int v = (total >> (i * datumNBits)) & mask;
                v += minDatum;
                if (v < min) {
                    min = v;
                }
            }
        }
        
        return min;
    }

    public static float[] findMinXY2(List<List<CurvatureScaleSpaceContour>> crLists) {
        
        float minX = Float.MAX_VALUE;
        float minY = Float.MAX_VALUE;
        
        for (List<CurvatureScaleSpaceContour> list : crLists) {
            for (CurvatureScaleSpaceContour cr : list) {
                float x = cr.getPeakDetails()[0].getXCoord();
                float y = cr.getPeakDetails()[0].getYCoord();
                if (x < minX) {
                    minX = x;
                }
                if (y < minY) {
                    minY = y;
                }
            }
        }
        
        return new float[]{minX, minY};
    }

    public static float[] findMaxXY2(List<List<CurvatureScaleSpaceContour>> crLists) {
        
        float maxX = Float.NEGATIVE_INFINITY;
        float maxY = Float.NEGATIVE_INFINITY;
        
        for (List<CurvatureScaleSpaceContour> list : crLists) {
            for (CurvatureScaleSpaceContour cr : list) {
                float x = cr.getPeakDetails()[0].getXCoord();
                float y = cr.getPeakDetails()[0].getYCoord();
                if (x > maxX) {
                    maxX = x;
                }
                if (y > maxY) {
                    maxY = y;
                }
            }
        }
        
        return new float[]{maxX, maxY};
    }
    
    /**
     * computes n!/(k!(n-k)!) and if result overflows a long, returns Long.MAX_VALUE.
     * 
     * @param n
     * @param k
     * @return 
     */
    protected static long computeNDivKTimesNMinusKBigInteger(int n, int k) {
        
        BigInteger result = computeNDivKTimesNMinusKBigIntegerExact(n, k);
        
        // work around for failure of result.compareTo(BigInteger.valueOf(Long.MAX_VALUE))
        if (result.bitLength() > 62) {
            return Long.MAX_VALUE;
        }
        
        if (result.bitLength() > 63) {
            throw new ArithmeticException("the result will not fit in a long");
        }
        return result.longValue();
    }
    
    /**
     * compute n!/k!(n-k)!
     * @param n
     * @param k
     * @return 
     * @throws ArithmeticException thrown when result is out of range of type long
     */
    protected static BigInteger computeNDivKTimesNMinusKBigIntegerExact(int n, int k) {
        
        return MiscMath0.computeNDivKTimesNMinusKBigIntegerExact(n, k);
    }

    /**
     * rescale a to values between vi and vf, inclusive
     * @param a
     * @param vi
     * @param vf
     * @return 
     */
    public static int[] rescale(double[] a, int vi, int vf) {
        
        double minV = Double.MAX_VALUE;
        double maxV = Double.MIN_VALUE;
        
        for (int i = 0; i < a.length; ++i) {
            double v = a[i];
            if (v < minV) {
                minV = v;
            }
            if (v > maxV) {
                maxV = v;
            }
        }
        double range = maxV - minV;
        
        double scale = ((double)(vf - vi))/range;
 
        int[] scaled = new int[a.length];
        
        for (int i = 0; i < a.length; ++i) {
            
            double v = (a[i] - minV) * scale;
            
            scaled[i] = (int)Math.round(v);
        }
        
        return scaled;
    }
    
    /**
     * rescale a to values between vi and vf, inclusive
     * @param a
     * @param vi
     * @param vf
     * @return 
     */
    public static int[][] rescale(double[][] a, int vi, int vf) {
        
        double minV = Double.MAX_VALUE;
        double maxV = Double.MIN_VALUE;
        
        int[][] out = new int[a.length][];
        
        for (int i = 0; i < a.length; ++i) {
            out[i] = new int[a[i].length];
            for (int j = 0; j < a[i].length; ++j) {
                double v = a[i][j];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
            }
        }
        double range = maxV - minV;
        
        double scale = ((double)(vf - vi))/range;
         
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {            
                double v = (a[i][j] - minV) * scale;
                out[i][j] = (int)Math.round(v);
            }
        }
        
        return out;
    }
    
    /**
     * rescale a to values between vi and vf, inclusive
     * @param a
     * @param vi
     * @param vf
     */
    public static void rescale(GreyscaleImage a, int vi, int vf) {
        
        int minV = Integer.MAX_VALUE;
        int maxV = Integer.MIN_VALUE;
        
        for (int i = 0; i < a.getNPixels(); ++i) {
            int v = a.getValue(i);
            if (v < minV) {
                minV = v;
            }
            if (v > maxV) {
                maxV = v;
            }
        }
        int range = maxV - minV;
        
        double scale = ((double)(vf - vi))/range;
         
        for (int i = 0; i < a.getNPixels(); ++i) {
            
            double v = (a.getValue(i) - minV) * scale;
            
            a.setValue(i, (int)Math.round(v));
        } 
    }
    
    /**
     * rescale a to values between vi and vf, inclusive.  the input array a is
     * not modified.
     * @param a
     * @param vi start of range to scale the values to
     * @param vf stop, inclusive, of range to scale the values to
     * @return 
     */
    public static float[] rescale(float[] a, int vi, int vf) {
        
        float minV = Float.MAX_VALUE;
        float maxV = Float.NEGATIVE_INFINITY;
        
        for (int i = 0; i < a.length; ++i) {
            float v = a[i];
            if (v < minV) {
                minV = v;
            }
            if (v > maxV) {
                maxV = v;
            }
        }
        float range = maxV - minV;
        
        float scale = ((float)(vf - vi))/range;
         
        float[] scaled = new float[a.length];
        
        for (int i = 0; i < a.length; ++i) {
                        
            scaled[i] = (a[i] - minV) * scale;
        }
        
        return scaled;
    }
    
    /**
     * rescale a to values between vi and vf, inclusive.  the input array a is
     * not modified.
     * @param a
     * @param vi start of range to scale the values to
     * @param vf stop, inclusive, of range to scale the values to
     * @return 
     */
    public static float[] rescale(float[] a, float minPossibleA, float maxPossibleA,
        int vi, int vf) {
        
        float range = maxPossibleA - minPossibleA;
        
        float scale = ((float)(vf - vi))/range;
         
        float[] scaled = new float[a.length];
        
        for (int i = 0; i < a.length; ++i) {
                        
            scaled[i] = (a[i] - minPossibleA) * scale;
        }
        
        return scaled;
    }

    public static int findYforX(HistogramHolder hist, float xSrch) {
        
        // looking for closest x to xSrch, then returning the y
        int idx = Arrays.binarySearch(hist.getXHist(), xSrch);
                    
        // if it's negative, (-(insertion point) - 1)
        if (idx < 0) {
            // idx = -*idx2 - 1
            idx = -1*(idx + 1);
        }
        if (idx > (hist.getXHist().length - 1)) {
            idx = hist.getXHist().length - 1;
        }
        
        if (idx < 0) {
            return -1;
        }
        
        return hist.getYHist()[idx];
    }

    /**
     * find the median in the double array of values.  
     * runtime complexity is O(N*lg2(N)).
     * @param a
     * @return 
     */
    public static double findMedian(double[][] a) {
        
        double[] values = new double[a.length*a[0].length];
        int count = 0;
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                values[count] = a[i][j];
                count++;
            }
        }
        
        Arrays.sort(values);
        
        return values[count/2];
    }
    
    /**
     * find the median in the int array of values.  
     * runtime complexity is either O(maxValue( or N*O(lg2(N)), whichever is
     * smaller.
     * @param a
     * @return 
     */
    public static int findMedian(int[][] a) {
        
        // put in an array to sort.
        // if max value < math.log(n)*n, will use counting sort
        int maxV = Integer.MIN_VALUE;
        int[] values = new int[a.length*a[0].length];
        int count = 0;
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                values[count] = a[i][j];
                if (values[count] > maxV) {
                    maxV = values[count];
                }
                count++;
            }
        }
        
        if (maxV < (count * Math.log(count)/Math.log(2))) {
            values = CountingSort.sort(values, maxV);
        } else {
            Arrays.sort(values);
        }
        
        return values[count/2];
    }
    
    public static void applyRescale(double[][] a, double minScaled, double maxScaled) {
        
        double minV = Double.MAX_VALUE;
        double maxV = Double.MIN_VALUE;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                double v = a[i][j];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
            }
        }
        double range = maxV - minV;
        
        double scale = (maxScaled - minScaled)/range;
        
//System.out.println("value 0 is rescaled to value=" + (-minV*scale)
//+ " minV=" + minV + " scale=" + scale);
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                a[i][j] = (a[i][j] - minV) * scale;
            }
        }
    }

    public static void applyRescale(int[][] a, int minScaled, int maxScaled) {
    
        int minV = Integer.MAX_VALUE;
        int maxV = Integer.MIN_VALUE;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                int v = a[i][j];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
            }
        }
        float range = maxV - minV;
        
        float scale = (maxScaled - minScaled)/range;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                a[i][j] = Math.round((a[i][j] - minV) * scale);
            }
        }
    }
    
    public static void applyRescale(float[][] a, float minScaled, 
        float maxScaled) {
        
        float minV = Float.MAX_VALUE;
        float maxV = Float.NEGATIVE_INFINITY;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                float v = a[i][j];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
            }
        }
        float range = maxV - minV;
        
        float scale = (maxScaled - minScaled)/range;
        
//System.out.println("value 0 is rescaled to value=" + (-minV*scale)
//+ " minV=" + minV + " scale=" + scale);
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                a[i][j] = (a[i][j] - minV) * scale;
            }
        }
    }
   
    /**
     *rescale a to between 0 and 255 and write to output
     * @param a
     * @return 
     */
    public static GreyscaleImage rescaleAndCreateImage(float[] a,
        int imageWidth, int imageHeight) {
        
        int nPix = a.length;
        
        if (nPix != (imageWidth * imageHeight)) {
            throw new IllegalArgumentException("a.length must be == to "
                + "imageWidth X imageHeight");
        }
        
        float minV = Float.MAX_VALUE;
        float maxV = Float.NEGATIVE_INFINITY;
        
        for (int i = 0; i < a.length; ++i) {
            float v = a[i];
            if (v < minV) {
                minV = v;
            }
            if (v > maxV) {
                maxV = v;
            }
        }
        float range = maxV - minV;
        
        float scale = 255.f/range;
        
        GreyscaleImage out = new GreyscaleImage(imageWidth, imageHeight);
        
        for (int i = 0; i < imageWidth; ++i) {
            for (int j = 0; j < imageHeight; ++j) {
                int pixIdx = (j * imageWidth) + i;
                int v = Math.round((a[pixIdx] - minV) * scale);
                if (v > 255) {
                    v = 255;
                }
                out.setValue(pixIdx, v);
            }
        }
        
        return out;
    }
    
    /**
     * ignores values of MAX_VALUE and MIN_VALUE during scale calculation
     * and then reassigns those to the calculated min and max of other
     * values.
     * @param a
     * @param minScaled
     * @param maxScaled 
     */
    public static void applyRescale2(float[][] a, float minScaled, 
        float maxScaled) {
        
        float minV = Float.MAX_VALUE;
        float maxV = Float.NEGATIVE_INFINITY;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                float v = a[i][j];
                if (v == Float.MAX_VALUE || v == Float.NEGATIVE_INFINITY) {
                    continue;
                }
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
            }
        }
        float range = maxV - minV;
        
        float scale = (maxScaled - minScaled)/range;
        
//System.out.println("value 0 is rescaled to value=" + (-minV*scale)
//+ " minV=" + minV + " scale=" + scale);
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                float v = a[i][j];
                if (v == Float.MAX_VALUE) {
                    a[i][j] = maxV;
                } else if (v == Float.NEGATIVE_INFINITY) {
                    a[i][j] = minV;
                } else {
                    a[i][j] = (v - minV) * scale;
                }
            }
        }
    }

    public static void applyAbsoluteValue(float[][] a) {

        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                if (a[i][j] < 0) {
                    a[i][j] *= -1;
                }
            }
        }
    }

    public static int calculateThetaAverage(GreyscaleImage img,
        int wrapAround, Set<PairInt> set) {
        
        /*
        -- calc peak frequency of values
        -- calculate average of values by assigning each
           value as the image value or image value
              + wrap around, depending on which is closer
              to the peak frequency
        -- correct the value to fit between 0 and wrap around
        */
        
        float[] values = new float[set.size()];
        int count = 0;
        for (PairInt p : set) {
            values[count] = img.getValue(p);
            count++;
        }
        
        HistogramHolder hist = Histogram.createSimpleHistogram(values, 
            Errors.populateYErrorsBySqrt(values));
        
        int peakValueIdx = MiscMath.findYMaxIndex(hist.getYHist());
        
        final float peakValue = hist.getXHist()[peakValueIdx];
        
        double sum = 0;
        for (int i = 0; i < values.length; ++i) {
            float v = values[i];
            float v2 = v + wrapAround;
            if (Math.abs(v - peakValue) <= Math.abs(v2 - peakValue)) {
                sum += v;
            } else {
                sum += v2;
            }
        }
        sum /= (double)values.length;
        
        int avg = (int)Math.round(sum);
        if (avg > wrapAround) {
            avg -= wrapAround;
        }
        
        return avg;
    }

    /**
     * calculate the mean of the averages, the standard deviation of those
     * averages and the addition of the standard deviations added in quadrature.
     * @param avgs
     * @param stdvs
     * @return [mean][std dev of mean][stdvs added in quadrature]
     */
    public static float[] calcMeanAndStDev(float[] avgs, float[] stdvs) {
    
        float[] avgAndStdv = MiscMath.getAvgAndStDev(avgs);
        
        double sum = 0;
        for (float stdv : stdvs) {
            sum += (stdv * stdv);
        }
        float sq = (float)Math.sqrt(sum/((float)stdvs.length - 1.0f));
        
        return new float[]{avgAndStdv[0], avgAndStdv[1], sq};
    }

    public static float[] calcMeanAndStDevWithWrapAround(float[] a, 
        int wrapAroundValue, float[] stdvs) {
        
        HistogramHolder hist = Histogram.createSimpleHistogram(a, 
            Errors.populateYErrorsBySqrt(a));
        
        int peakValueIdx = MiscMath.findYMaxIndex(hist.getYHist());
        
        final float peakValue = hist.getXHist()[peakValueIdx];
        
        float[] a2 = new float[a.length];
        for (int i = 0; i < a.length; ++i) {
            float v = a[i];
            float v2 = v + wrapAroundValue;
            if (Math.abs(v - peakValue) <= Math.abs(v2 - peakValue)) {
                a2[i] = v;
            } else {
                a2[i] = v2;
            }
        }

        float[] avgAndStdv = MiscMath.getAvgAndStDev(a2);
        if (avgAndStdv[0] > wrapAroundValue) {
            avgAndStdv[0] -= wrapAroundValue;
        }
        
        double sum = 0;
        for (float stdv : stdvs) {
            sum += (stdv * stdv);
        }
        float sq = (float)Math.sqrt(sum/((float)stdvs.length - 1.0f));
        
        return new float[]{avgAndStdv[0], avgAndStdv[1], sq};
    }
    
    /**
     * calculate the distance between furthest points as object
     * size.
     * @param points
     * @return 
     */
    public static int calculateObjectSize(Set<PairInt> points) {
        // O(N*lg_2(N))
        FurthestPair furthestPair = new FurthestPair();
        PairInt[] fp = furthestPair.find(points);
        if (fp == null || fp.length < 2) {            
            throw new IllegalArgumentException("did not find a furthest pair" + " in points");
        }
        double dist = ORBMatcher.distance(fp[0], fp[1]);
        return (int) Math.round(dist);
    }
    
    /**
     * calculate the distance between furthest points as object
     * size.
     * @param points
     * @return 
     */
    public static int calculateObjectSize(PairIntArray points) {
        return calculateObjectSize(Misc.convert(points));
    }
   
    // ---- 32 bit hash methods ---
    /**
     * get index of list from value using the hash.
     * 
     * This method is based upon a method from the libquantum library 
     * file qureg.c which has copyright:
     
       Copyright 2003, 2004 Bjoern Butscher, Hendrik Weimer

       This file is part of libquantum

       libquantum is free software; you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published
       by the Free Software Foundation; either version 3 of the License,
       or (at your option) any later version.

       libquantum is distributed in the hope that it will be useful, but
       WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
       General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with libquantum; if not, write to the Free Software
       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
       MA 02110-1301, USA
     * 
     * =============
     * 
     * the hash is based on public domain FNV31
     * http://www.isthe.com/chongo/src/fnv/hash_32.c
     * 
     * @param list
     * @param value
     * @param hash 
     * @return index of list where value is stored
     */
    public static int getHashIndex(final int[] list, final int value, int[] hash) {
        
        int i, mark = 0;

        i = hashToIndex(value, list.length);
     
        while (hash[i] != -1) {
            if (list[hash[i]] == value) {
                return hash[i];
            }
            i++;
            if (i == list.length) {
                i = 0;
            }
        }

        return -1;
    }
    
    /**
     * get index of list from value using the hash.
     * 
     * This method is based upon a method from the libquantum library 
     * file qureg.c which has copyright:
     
       Copyright 2003, 2004 Bjoern Butscher, Hendrik Weimer

       This file is part of libquantum

       libquantum is free software; you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published
       by the Free Software Foundation; either version 3 of the License,
       or (at your option) any later version.

       libquantum is distributed in the hope that it will be useful, but
       WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
       General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with libquantum; if not, write to the Free Software
       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
       MA 02110-1301, USA
     * 
     * =============
     * 
     * the hash is based on public domain FNV31
     * http://www.isthe.com/chongo/src/fnv/hash_32.c
     * 
     * =======
     * NOTE: this method was made to only use primitives, for use with the
     * quantum methods.
     * 
     * @param list
     * @param value
     * @param hash 
     * @return index of list where value is stored
     */
    public static boolean hashContains(final int[] list, final int value, int[] hash) {
        
        int i, mark = 0;

        i = hashToIndex(value, list.length);
        
        //hash[i] holds list idx;

        int nStart = i;
        
        //TODO: this could be improved. worse case could be O(hash.length)
        while (hash[i] != -1) {
            if (list[hash[i]] == value) {
                return true;
            }
            i++;
            if (i == list.length) {
                i = 0;
            }
            if (i == nStart) {
                return false;
            }
        }

        return false;
    }
    
    /**
     * set list[idx] into hash.
     * 
     * This method is based upon a method from the libquantum library 
     * file qureg.c which has copyright:
     
       Copyright 2003, 2004 Bjoern Butscher, Hendrik Weimer

       This file is part of libquantum

       libquantum is free software; you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published
       by the Free Software Foundation; either version 3 of the License,
       or (at your option) any later version.

       libquantum is distributed in the hope that it will be useful, but
       WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
       General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with libquantum; if not, write to the Free Software
       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
       MA 02110-1301, USA
     * 
     * =============
     * 
     * the hash is based on public domain FNV31
     * http://www.isthe.com/chongo/src/fnv/hash_32.c
     * 
     * @param value
     * @param listLength length of array with values to hash (should be same size
     * as hash)
     * @param hash
     * @return index of hash where list[idx] was stored
     */
    public static int getNextHashIndex(final int value, final int listLength, 
        int[] hash) {
        
        int i, mark = 0;

        i = hashToIndex(value, listLength);
        
        while (hash[i] != -1) {
            i++;
            // if i is > last index
            if (i == listLength) {
                if (mark == 0) {
                    i = 0;
                    mark = 1;
                } else {
                    StackTraceElement[] st = Thread.currentThread().getStackTrace();
                    for (StackTraceElement s : st) {
                        System.out.println(s);
                    }
                    throw new IllegalStateException("hash is full.  i=" + i + 
                         " hash.lrngth=" + listLength);
                }
            }
        }
        
        return i;
    }
    
    /**
     * set list[idx] into hash.
     * 
     * This method is based upon a method from the libquantum library 
     * file qureg.c which has copyright:
     
       Copyright 2003, 2004 Bjoern Butscher, Hendrik Weimer

       This file is part of libquantum

       libquantum is free software; you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published
       by the Free Software Foundation; either version 3 of the License,
       or (at your option) any later version.

       libquantum is distributed in the hope that it will be useful, but
       WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
       General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with libquantum; if not, write to the Free Software
       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
       MA 02110-1301, USA
     * 
     * =============
     * 
     * the hash is based on public domain FNV31
     * http://www.isthe.com/chongo/src/fnv/hash_32.c
     * 
     * @param list
     * @param idx
     * @param hash
     * @return index of hash where list[idx] was stored
     */
    public static int setHashValue(final int[] list, final int idx, int[] hash) {
        
        int i = getNextHashIndex(list[idx], list.length, hash);
        
        hash[i] = idx;     

        return i;
    }
    
    /**
     * return an index for the given value
     * 
     * @param value
     * @param length size of the container in which the result is an index in
     * @return 
     */
    protected static int hashToIndex(int value, int length) {
       
        return fnvRetry(value, length);
    }
    
    /**
     * using a hash based on public domain FNV31
     * http://www.isthe.com/chongo/src/fnv/hash_32.c
     * @param value
     * @param length length of the array value will be
     * hashed to.  no index .geq. length will be returned.
     * @return index from hashing value
     */
    protected static int fnvRetry(int value, int length) {

        if (length == 0) {
            throw new IllegalArgumentException("length must be >0");
        }        
        
        //#define TRUE_HASH_SIZE ((u_int32_t)50000) /* range top plus 1
        int FNV_32_PRIME = 16777619;
        //int FNV1_32_INIT = 2166136261;
        int FNV1_32_INIT = 1083068130; // which is 2166136261 >> 1
        //#define MAX_32BIT ((u_int32_t)0xffffffff) /* largest 32 bit unsigned value
        //#define RETRY_LEVEL ((MAX_32BIT / TRUE_HASH_SIZE) * TRUE_HASH_SIZE)
        int RETRY_LEVEL = (Integer.MAX_VALUE/length) * length;
        
        int hash = fnv_31(value, FNV1_32_INIT);
        //System.out.println("hash=" + hash);
        while (hash >= RETRY_LEVEL) {
            hash = (hash * FNV_32_PRIME) + FNV1_32_INIT;
            //System.out.println("  hash=" + hash);
        }
        
        //NOTE: forcing positive for use case
        hash = (hash ^ (hash >> 31)) + (hash >>> 31);
        
        hash %= length;
        
        return hash;
    }

    /**
     * using a hash based on public domain FNV31
     * http://www.isthe.com/chongo/src/fnv/hash_32.c
     * @param value
     * @param FNV1_31_INIT
     * @return 
     */
    private static int fnv_31(int value, int FNV1_31_INIT) {
       
        int hval = FNV1_31_INIT;
        
        //multiply by the 32 bit FNV magic prime mod 2^32
        //if no opt: hval *= 0x01000193;//FNV_32_PRIME;
        //else:
        hval += (hval<<1) + (hval<<4) + (hval<<7) 
            + (hval<<8) + (hval<<24);

	    // xor the bottom with the current octet
	    hval ^= value;
  
        return hval;
    }
    
    /**
     * for the setBits only, create permutations of the bits that are set, that
     * is creating all combinations of number with just those bits being
     * 0 or 1.
     * @param setBits
     * @return 
     */
    public static int[] permuteTheSetBits(int setBits) {
        
        int i;
        int nBits = MiscMath.numberOfBits(setBits);
        int nSetBits = 0;
        int unsetBits = 0;        
        
        for (i = 0; i < nBits; ++i) {
            if ((setBits & (1 << i)) != 0) {
                nSetBits++;
            } else {
                unsetBits |= (1 << i);
            }
        }
        
        int listLen = (int)Math.pow(2, nSetBits);
        int[] outputList = new int[listLen];
       
        int[] hash = new int[outputList.length];
        for (i = 0; i < hash.length; ++i) {
            hash[i] = -1;
        }
                 
        int[] lastIdx = new int[]{-1};
        
        int number = setBits;
                
        permuteTheSetBits(outputList, hash, lastIdx, number, nSetBits, unsetBits);
            
        return outputList;
    }

    private static void permuteTheSetBits(int[] list, int[] hash, int[] lastIdx, int number, 
        int hiIdx, int unsetBits) {
    
        storeIfUnique(list, hash, lastIdx, number);
        
        while ((hiIdx > -1) && ((unsetBits & (1 << hiIdx)) != 0)) {
            hiIdx--;
        }
        if (hiIdx < 0) {
            return;
        }
        
        // recurse the 2 permutations of bit hiIdx
        
        permuteTheSetBits(list, hash, lastIdx, number, hiIdx - 1, unsetBits);
        
        // unset bit hiIdx
        number &= ~(1 << hiIdx);
        permuteTheSetBits(list, hash, lastIdx, number, hiIdx - 1, unsetBits);
        
    }

    private static void storeIfUnique(int[] list, int[] hash, int[] lastIdx, int value) {
        
        if (!MiscMath.hashContains(list, value, hash)) {
            lastIdx[0]++;
            if (lastIdx[0] > (list.length - 1)) {
                throw new IllegalStateException("index out of bounds");
            }
            list[lastIdx[0]] = value;
            MiscMath.setHashValue(list, lastIdx[0], hash);
        }
    }
    
}
