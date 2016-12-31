package algorithms;

import algorithms.util.PairIntArray;
import java.util.Arrays;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/sorting/MultiArrayMergeSort.java
 * under MIT License (MIT), Nichole King 2013

 * @author nichole
 */
public class MultiArrayMergeSort {
    
    /**
     * sort by increasing value a1 and apply same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArgThen2nd(float[] a1, float[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, 0, a1.length - 1);
    }
    
     /**
     * sort by increasing value a1 and apply same changes to a2.
     * 
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArg(float[] a1, float[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArg(a1, a2, 0, a1.length - 1);
    }
    
    /**
     * sort by increasing value a1 for ties sort by a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArgThen2nd(int[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, 0, a1.length - 1);
    }
    
    /**
     * sort by increasing value a1 for ties sort by a2 and apply same changes 
     * to a3.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also, and its used
     * for sorting when a1's are equal.
     * @param a3 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArgThen2nd(int[] a1, int[] a2, int[] a3) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a3 == null) {
            throw new IllegalArgumentException("a3 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        if (a1.length != a3.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, a3, 0, a1.length - 1);
    }
        
    public static void sortBy1stArgThen2nd(float[] a1, float[] a2, int[] a3) {
        
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a3 == null) {
            throw new IllegalArgumentException("a3 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        if (a1.length != a3.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, a3, 0, a1.length - 1);
    }
    
    /**
     * sort by increasing value y.
     * Ties are further sorted by increasing values of x.
     * runtime is O(N * log_2(N))
     *
     * @param xy array of points to be sorted
     */
    public static void sortByYThenX(PairIntArray xy) {
        
        if (xy == null) {
            throw new IllegalArgumentException("xy cannot be null");
        }
        
        sortByYThenX(xy, 0, xy.getN() - 1);
    }
    
    /**
     * sort so that a is decreasing in value for higher indexes and b is
     * is increasing in value for higher indexes for the same a.  swap
     * operations for a and b logic are performed on c and d too.
     * 
     * @param a
     * @param b
     * @param c
     * @param d 
     */
    public static void sortBy1stDescThen2ndAsc(int[] a, double[] b, 
        Integer[][] c, int[] d) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (d == null) {
            throw new IllegalArgumentException("d cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "number of items in a must be the same as in b");
        }
        if (a.length != c.length) {
            throw new IllegalArgumentException(
            "number of items in a must be the same as in c");
        }
        if (a.length != d.length) {
            throw new IllegalArgumentException(
            "number of items in a must be the same as in d");
        }
        
        sortBy1stDescThen2ndAsc(a, b, c, d, 0, a.length - 1);
    }
    
    /**
     * sort so that a is decreasing in value for higher indexes and b is
     * is increasing in value for higher indexes for the same a.  swap
     * operations for a and b logic are performed on c and d too.
     * 
     * @param a
     * @param b
     * @param c
     */
    public static void sortBy1stDescThen2ndAsc(int[] a, double[] b, int[] c) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "number of items in a must be the same as in b");
        }
        if (a.length != c.length) {
            throw new IllegalArgumentException(
            "number of items in a must be the same as in c");
        }
        
        sortBy1stDescThen2ndAsc(a, b, c, 0, a.length - 1);
    }
    
    /**
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    private static void sortBy1stArgThen2nd(float[] a1, float[] a2, int idxLo, 
        int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArg(float[] a1, float[] a2, int idxLo, 
        int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArg(a1, a2, idxLo, indexMid);
            sortBy1stArg(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArg(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    private static void sortBy1stArgThen2nd(int[] a1, int[] a2, int idxLo, 
        int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param a3 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    private static void sortBy1stArgThen2nd(int[] a1, int[] a2, int[] a3, 
        int idxLo, int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, a3, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, a3, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, a3, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param a3 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    private static void sortBy1stArgThen2nd(float[] a1, float[] a2, int[] a3, 
        int idxLo, int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, a3, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, a3, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, a3, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     * @param xy array of points to be sorted
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    private static void sortByYThenX(PairIntArray xy, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            sortByYThenX(xy, idxLo, indexMid);
            sortByYThenX(xy, indexMid + 1, idxHi);
            mergeByYThenX(xy, idxLo, indexMid, idxHi);
        }
    }

    private static void mergeBy1stArgThen2nd(float[] a1, float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        a2Left[nLeft] = Float.MAX_VALUE;
        a1Left[nLeft] = Float.MAX_VALUE;
        a2Right[nRight] = Float.MAX_VALUE;
        a1Right[nRight] = Float.MAX_VALUE;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l == r) {
                float lx = a2Left[leftPos];
                float rx = a2Right[rightPos];

                if (lx <= rx) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArg(float[] a1, float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        a2Left[nLeft] = Float.MAX_VALUE;
        a1Left[nLeft] = Float.MAX_VALUE;
        a2Right[nRight] = Float.MAX_VALUE;
        a1Right[nRight] = Float.MAX_VALUE;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l <= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArgThen2nd(int[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MAX_VALUE;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l == r) {
                float lx = a2Left[leftPos];
                float rx = a2Right[rightPos];

                if (lx <= rx) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArgThen2nd(int[] a1, int[] a2, int[] a3, 
        int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];
        int[] a3Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];
        int[] a3Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        System.arraycopy(a3, idxLo, a3Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        System.arraycopy(a3, idxMid + 1, a3Right, 0, nRight);
        
        int sentinel = Integer.MAX_VALUE;
        a2Left[nLeft] = sentinel;
        a3Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a3Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l == r) {
                int lx = a2Left[leftPos];
                int rx = a2Right[rightPos];

                if (lx <= rx) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    a3[k] = a3Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    a3[k] = a3Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                a3[k] = a3Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                a3[k] = a3Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArgThen2nd(float[] a1, float[] a2, int[] a3, 
        int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];
        int[] a3Left = new int[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];
        int[] a3Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        System.arraycopy(a3, idxLo, a3Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        System.arraycopy(a3, idxMid + 1, a3Right, 0, nRight);
        
        float sentinel = Float.MAX_VALUE;
        int sentinel2 = Integer.MAX_VALUE;
        a2Left[nLeft] = sentinel;
        a3Left[nLeft] = sentinel2;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a3Right[nRight] = sentinel2;
        a1Right[nRight] = sentinel;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l == r) {
                float lx = a2Left[leftPos];
                float rx = a2Right[rightPos];

                if (lx <= rx) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    a3[k] = a3Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    a3[k] = a3Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                a3[k] = a3Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                a3[k] = a3Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeByYThenX(PairIntArray xy, int idxLo, 
        int idxMid, int idxHi) {
 
        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] yLeft = new int[nLeft + 1];
        int[] xLeft = new int[nLeft + 1];

        int[] yRight = new int[nRight + 1];
        int[] xRight = new int[nRight + 1];

        System.arraycopy(xy.getX(), idxLo, xLeft, 0, nLeft);
        System.arraycopy(xy.getY(), idxLo, yLeft, 0, nLeft);
        
        System.arraycopy(xy.getX(), idxMid + 1, xRight, 0, nRight);
        System.arraycopy(xy.getY(), idxMid + 1, yRight, 0, nRight);
        
        int sentinel = Integer.MAX_VALUE;
        yLeft[nLeft] = sentinel;
        xLeft[nLeft] = sentinel;
        yRight[nRight] = sentinel;
        xRight[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = yLeft[leftPos];
            int r = yRight[rightPos];

            if (l == r) {
                int lx = xLeft[leftPos];
                int rx = xRight[rightPos];

                if (lx <= rx) {
                    xy.getX()[k] = xLeft[leftPos];
                    xy.getY()[k] = yLeft[leftPos];
                    leftPos++;
                } else {
                    xy.getX()[k] = xRight[rightPos];
                    xy.getY()[k] = yRight[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                xy.getX()[k] = xLeft[leftPos];
                xy.getY()[k] = yLeft[leftPos];
                leftPos++;
            } else {
                xy.getX()[k] = xRight[rightPos];
                xy.getY()[k] = yRight[rightPos];
                rightPos++;
            }
        }
    }

    /**
     * sort by decreasing value a1 and apply same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     */
    public static void sortByDecr(int[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortByDecr(a1, a2, 0, a1.length - 1);
              
    }
    
    /**
     * sort by decreasing value a1 and apply same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     */
    public static void sortByDecr(Double[] a1, Float[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortByDecr(a1, a2, 0, a1.length - 1);
              
    }
    
    public static void sortByDecr(float[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortByDecr(a1, a2, 0, a1.length - 1);
    }

    /**
     * sort by decreasing value a1 and apply same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     */
    public static void sortByDecr(double[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortByDecr(a1, a2, 0, a1.length - 1);
              
    }
    
    public static void sortByDecr(int[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    public static void sortByDecr(Double[] a1, Float[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    public static void sortByDecr(float[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    public static void sortByDecr(double[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeByDecr(int[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MIN_VALUE;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeByDecr(Double[] a1, Float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        Float[] a2Left = new Float[nLeft + 1];
        Double[] a1Left = new Double[nLeft + 1];

        Float[] a2Right = new Float[nRight + 1];
        Double[] a1Right = new Double[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        float sentinel = Float.NEGATIVE_INFINITY;
        double sentinel2 = Double.NEGATIVE_INFINITY;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel2;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel2;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            Double l = a1Left[leftPos];
            Double r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeByDecr(float[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        float sentinel = Float.NEGATIVE_INFINITY;
        int sentinel2 = Integer.MIN_VALUE;
        a2Left[nLeft] = sentinel2;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel2;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            double l = a1Left[leftPos];
            double r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeByDecr(double[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        double[] a1Left = new double[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        double[] a1Right = new double[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        double sentinel = Double.NEGATIVE_INFINITY;
        int sentinel2 = Integer.MIN_VALUE;
        a2Left[nLeft] = sentinel2;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel2;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            double l = a1Left[leftPos];
            double r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    protected static void sortBy1stDescThen2ndAsc(int[] a, double[] b, 
        int[] c, int idxLo, int idxHi) {
        
        if (idxLo < idxHi) {

            int idxMid = (idxLo + idxHi) >> 1;
            
            sortBy1stDescThen2ndAsc(a, b, c, idxLo, idxMid);
            
            sortBy1stDescThen2ndAsc(a, b, c, idxMid + 1, idxHi);
            
            mergeBy1stDescThen2ndAsc(a, b, c, idxLo, idxMid, idxHi);
        }
    }
    
    protected static void sortBy1stDescThen2ndAsc(int[] a, double[] b, 
        Integer[][] c, int[] d, int idxLo, int idxHi) {
        
        if (idxLo < idxHi) {

            int idxMid = (idxLo + idxHi) >> 1;
            
            sortBy1stDescThen2ndAsc(a, b, c, d, idxLo, idxMid);
            
            sortBy1stDescThen2ndAsc(a, b, c, d, idxMid + 1, idxHi);
            
            mergeBy1stDescThen2ndAsc(a, b, c, d, idxLo, idxMid, idxHi);
        }
    }
    
    public static void sortBy1stAscThen2ndDesc(double[] a1, int[] a2, 
        Integer[][] a3, int[] a4, int idxLo, int idxHi) {
        
        if (idxLo < idxHi) {

            int idxMid = (idxLo + idxHi) >> 1;
            
            sortBy1stAscThen2ndDesc(a1, a2, a3, a4, idxLo, idxMid);
            
            sortBy1stAscThen2ndDesc(a1, a2, a3, a4, idxMid + 1, idxHi);
            
            mergeBy1stAscThen2ndDesc(a1, a2, a3, a4, idxLo, idxMid, idxHi);
        }
    }
    
    private static void mergeBy1stDescThen2ndAsc(int[] a, double[] b, 
        Integer[][] c, int[] d, int idxLo, int idxMid, int idxHi) {
        
        int[] aLeft = Arrays.copyOfRange(a, idxLo, idxMid + 2);
        double[] bLeft = Arrays.copyOfRange(b, idxLo, idxMid + 2);
        int[] dLeft = Arrays.copyOfRange(d, idxLo, idxMid + 2);
        Integer[][] cLeft = new Integer[(idxMid + 2 - idxLo)][];
        for (int i = 0; i < (cLeft.length - 1); ++i) {
            int idx = i + idxLo;
            cLeft[i] = Arrays.copyOf(c[idx], c[idx].length);
        }        
        
        int[] aRight = Arrays.copyOfRange(a, idxMid + 1, idxHi + 2);
        double[] bRight = Arrays.copyOfRange(b, idxMid + 1, idxHi + 2);
        int[] dRight = Arrays.copyOfRange(d, idxMid + 1, idxHi + 2);
        Integer[][] cRight = new Integer[(idxHi + 2) - (idxMid + 1)][];
        for (int i = 0; i < (cRight.length - 1); ++i) {
            int idx = i + idxMid + 1;            
            cRight[i] = Arrays.copyOf(c[idx], c[idx].length);
        }
        
        aLeft[aLeft.length - 1] = Integer.MIN_VALUE;
        bLeft[bLeft.length - 1] = Double.MIN_VALUE;
        dLeft[dLeft.length - 1] = Integer.MIN_VALUE;
        cLeft[cLeft.length - 1] = new Integer[c[idxLo].length];// not compared, so can be 0's
        
        aRight[aRight.length - 1] = Integer.MIN_VALUE;
        bRight[bRight.length - 1] = Double.MIN_VALUE;
        dRight[dRight.length - 1] = Integer.MIN_VALUE;
        cRight[cRight.length - 1] = new Integer[c[idxMid + 1].length];// not compared, so can be 0's
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = aLeft[leftPos];
            int r = aRight[rightPos];
            if (l > r) {
                a[k] = aLeft[leftPos];
                b[k] = bLeft[leftPos];
                c[k] = cLeft[leftPos];
                d[k] = dLeft[leftPos];
                leftPos++;
            } else if (l == r) {
                // sort for ascending values of b
                double l2 = bLeft[leftPos];
                double r2 = bRight[rightPos];
                if (l2 <= r2) {
                    a[k] = aLeft[leftPos];
                    b[k] = bLeft[leftPos];
                    c[k] = cLeft[leftPos];
                    d[k] = dLeft[leftPos];
                    leftPos++;
                } else {
                    a[k] = aRight[rightPos];
                    b[k] = bRight[rightPos];
                    c[k] = cRight[rightPos];
                    d[k] = dRight[rightPos];
                    rightPos++;
                }
            } else {
                a[k] = aRight[rightPos];
                b[k] = bRight[rightPos];
                c[k] = cRight[rightPos];
                d[k] = dRight[rightPos];
                rightPos++;
            }
        }
        
    }

    private static void mergeBy1stDescThen2ndAsc(int[] a, double[] b, 
        int[] c, int idxLo, int idxMid, int idxHi) {
        
        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] aLeft = new int[nLeft + 1];
        double[] bLeft = new double[nLeft + 1];
        int[] cLeft = new int[nLeft + 1];
        
        int[] aRight = new int[nRight + 1];
        double[] bRight = new double[nRight + 1];
        int[] cRight = new int[nRight + 1];
        
        System.arraycopy(a, idxLo, aLeft, 0, nLeft);
        System.arraycopy(b, idxLo, bLeft, 0, nLeft);
        System.arraycopy(c, idxLo, cLeft, 0, nLeft);
        
        System.arraycopy(a, idxMid + 1, aRight, 0, nRight);
        System.arraycopy(b, idxMid + 1, bRight, 0, nRight);
        System.arraycopy(c, idxMid + 1, cRight, 0, nRight);
        
        int sentinel = Integer.MIN_VALUE;
        double sentinel2 = Double.NEGATIVE_INFINITY;
        bLeft[nLeft] = sentinel2;
        aLeft[nLeft] = sentinel;
        cLeft[nLeft] = sentinel;
        bRight[nRight] = sentinel2;
        aRight[nRight] = sentinel;
        cRight[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = aLeft[leftPos];
            int r = aRight[rightPos];
            if (l > r) {
                a[k] = aLeft[leftPos];
                b[k] = bLeft[leftPos];
                c[k] = cLeft[leftPos];
                leftPos++;
            } else if (l == r) {
                // sort for ascending values of b
                double l2 = bLeft[leftPos];
                double r2 = bRight[rightPos];
                if (l2 <= r2) {
                    a[k] = aLeft[leftPos];
                    b[k] = bLeft[leftPos];
                    c[k] = cLeft[leftPos];
                    leftPos++;
                } else {
                    a[k] = aRight[rightPos];
                    b[k] = bRight[rightPos];
                    c[k] = cRight[rightPos];
                    rightPos++;
                }
            } else {
                a[k] = aRight[rightPos];
                b[k] = bRight[rightPos];
                c[k] = cRight[rightPos];
                rightPos++;
            }
        }
        
    }
    
    private static void mergeBy1stAscThen2ndDesc(double[] a1, int[] a2, 
        Integer[][] a3, int[] a4, int idxLo, int idxMid, int idxHi) {

        double[] aLeft = Arrays.copyOfRange(a1, idxLo, idxMid + 2);
        int[] bLeft = Arrays.copyOfRange(a2, idxLo, idxMid + 2);
        int[] dLeft = Arrays.copyOfRange(a4, idxLo, idxMid + 2);
        Integer[][] cLeft = new Integer[(idxMid + 2 - idxLo)][];
        for (int i = 0; i < (cLeft.length - 1); ++i) {
            int idx = i + idxLo;
            cLeft[i] = Arrays.copyOf(a3[idx], a3[idx].length);
        }        
        
        double[] aRight = Arrays.copyOfRange(a1, idxMid + 1, idxHi + 2);
        int[] bRight = Arrays.copyOfRange(a2, idxMid + 1, idxHi + 2);
        int[] dRight = Arrays.copyOfRange(a4, idxMid + 1, idxHi + 2);
        Integer[][] cRight = new Integer[(idxHi + 2) - (idxMid + 1)][];
        for (int i = 0; i < (cRight.length - 1); ++i) {
            int idx = i + idxMid + 1;            
            cRight[i] = Arrays.copyOf(a3[idx], a3[idx].length);
        }
        
        aLeft[aLeft.length - 1] = Double.MAX_VALUE;
        bLeft[bLeft.length - 1] = Integer.MAX_VALUE;
        dLeft[dLeft.length - 1] = Integer.MAX_VALUE;
        cLeft[cLeft.length - 1] = new Integer[a3[idxLo].length];// not compared, so can be 0's
        
        aRight[aRight.length - 1] = Double.MAX_VALUE;
        bRight[bRight.length - 1] = Integer.MAX_VALUE;
        dRight[dRight.length - 1] = Integer.MAX_VALUE;
        cRight[cRight.length - 1] = new Integer[a3[idxMid + 1].length];// not compared, so can be 0's
        
        int leftPos = 0;
        int rightPos = 0;
        
        for (int k = idxLo; k <= idxHi; k++) {
            double l = aLeft[leftPos];
            double r = aRight[rightPos];
            if (l < r) {
                a1[k] = aLeft[leftPos];
                a2[k] = bLeft[leftPos];
                a3[k] = cLeft[leftPos];
                a4[k] = dLeft[leftPos];
                leftPos++;
            } else if (l == r) {
                // sort for descending values of b
                double l2 = bLeft[leftPos];
                double r2 = bRight[rightPos];
                if (l2 >= r2) {
                    a1[k] = aLeft[leftPos];
                    a2[k] = bLeft[leftPos];
                    a3[k] = cLeft[leftPos];
                    a4[k] = dLeft[leftPos];
                    leftPos++;
                } else {
                    a1[k] = aRight[rightPos];
                    a2[k] = bRight[rightPos];
                    a3[k] = cRight[rightPos];
                    a4[k] = dRight[rightPos];
                    rightPos++;
                }
            } else {
                a1[k] = aRight[rightPos];
                a2[k] = bRight[rightPos];
                a3[k] = cRight[rightPos];
                a4[k] = dRight[rightPos];
                rightPos++;
            }
        }
    }

}
