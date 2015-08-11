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

    private static void mergeBy1stArgThen2nd( float[] a1, float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        float[] a1Left = Arrays.copyOfRange(a1, idxLo, idxMid + 2);
        float[] a2Left = Arrays.copyOfRange(a2, idxLo, idxMid + 2);
        
        float[] a1Right = Arrays.copyOfRange(a1, idxMid + 1, idxHi + 2);
        float[] a2Right = Arrays.copyOfRange(a2, idxMid + 1, idxHi + 2);
        
        a1Left[a1Left.length - 1] = Float.MAX_VALUE;
        a2Left[a2Left.length - 1] = Float.MAX_VALUE;
        a1Right[a1Right.length - 1] = Float.MAX_VALUE;
        a2Right[a2Right.length - 1] = Float.MAX_VALUE;

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
    
    private static void mergeBy1stArgThen2nd(int[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int[] a1Left = Arrays.copyOfRange(a1, idxLo, idxMid + 2);
        int[] a2Left = Arrays.copyOfRange(a2, idxLo, idxMid + 2);
        
        int[] a1Right = Arrays.copyOfRange(a1, idxMid + 1, idxHi + 2);
        int[] a2Right = Arrays.copyOfRange(a2, idxMid + 1, idxHi + 2);
        
        a1Left[a1Left.length - 1] = Integer.MAX_VALUE;
        a2Left[a2Left.length - 1] = Integer.MAX_VALUE;
        a1Right[a1Right.length - 1] = Integer.MAX_VALUE;
        a2Right[a2Right.length - 1] = Integer.MAX_VALUE;

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
    
    private static void mergeByYThenX(PairIntArray xy, int idxLo, 
        int idxMid, int idxHi) {

        int[] yLeft = Arrays.copyOfRange(xy.getY(), idxLo, idxMid + 2);
        int[] xLeft = Arrays.copyOfRange(xy.getX(), idxLo, idxMid + 2);
        
        int[] yRight = Arrays.copyOfRange(xy.getY(), idxMid + 1, idxHi + 2);
        int[] xRight = Arrays.copyOfRange(xy.getX(), idxMid + 1, idxHi + 2);
        
        yLeft[yLeft.length - 1] = Integer.MAX_VALUE;
        xLeft[xLeft.length - 1] = Integer.MAX_VALUE;
        yRight[yRight.length - 1] = Integer.MAX_VALUE;
        xRight[xRight.length - 1] = Integer.MAX_VALUE;

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
    
    public static void sortByDecr(int[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeByDecr(int[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int[] a1Left = Arrays.copyOfRange(a1, idxLo, idxMid + 2);
        int[] a2Left = Arrays.copyOfRange(a2, idxLo, idxMid + 2);
        
        int[] a1Right = Arrays.copyOfRange(a1, idxMid + 1, idxHi + 2);
        int[] a2Right = Arrays.copyOfRange(a2, idxMid + 1, idxHi + 2);
        
        a1Left[a1Left.length - 1] = Integer.MIN_VALUE;
        a2Left[a2Left.length - 1] = Integer.MIN_VALUE;
        a1Right[a1Right.length - 1] = Integer.MIN_VALUE;
        a2Right[a2Right.length - 1] = Integer.MIN_VALUE;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l > r) {
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
    
}
