package algorithms.sorting;

import java.util.Arrays;

/**
 * <pre>
 * merge sort worse case runtime is O(N * log_2(N))
 *
 * The method to sort 2 arrays below uses divide and conquer paradigm and sorts one additional array
 * (with no additional logic) along with the original array which is sorted.
 * the result is a constant factor added to an internal constant for the O(n) for merge of single
 * sort.
 *
 * constructed from pseudo-code in Cormen et al.
 * "Introduction to Algorithms"
 * <pre>
 * 
 * @author Nichole King
 */
public class MultiArrayMergeSort {
    
    /**
     * sort by increasing value a1 and apply same changes to a2.
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
        	throw new IllegalArgumentException("number of items in a1 must be the same as in a2");
        }
        sortBy1stArg(a1, a2, 0, a1.length - 1);
    }

    /**
     * sort by increasing value a1 and apply same changes to a2 and a3.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param a3 array of points to apply a1 sorting to also
     * @param nLength the number of points to sort in a1 starting with item 0
     */
    public static void sortBy1stArg(float[] a1, float[] a2, int[] a3, int nLength) {
        if (a1 == null) {
        	throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
        	throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a3 == null) {
            throw new IllegalArgumentException("a3 cannot be null");
        }
        if ((a1.length != nLength) || (a2.length != nLength) || (a3.length != nLength)) {
            throw new IllegalArgumentException("number of items in a1 and a2 and a3 must be equal to nLength");
        }
        sortBy1stArg(a1, a2, a3, 0, nLength - 1);
    }

    /**
     * sort by increasing value a1 and apply same changes to a2 and a3.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArg(float[] a1, float[] a2, int idxLo, int idxHi) {

        int idxMid = -1;

        if (idxLo < idxHi) {

            idxMid = (idxLo + idxHi)/2;

            sortBy1stArg(a1, a2, idxLo, idxMid);
            sortBy1stArg(a1, a2, idxMid + 1, idxHi);
            merge(a1, a2, idxLo, idxMid, idxHi);
        }
    }

    /**
     * sort by increasing value a1 and apply same changes to a2 and a3.
     * runtime is O(N * log_2(N))
     *
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param a3 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArg(float[] a1, float[] a2, int[] a3, int idxLo, int idxHi) {

        int idxMid = -1;

        if (idxLo < idxHi) {

            idxMid = (idxLo + idxHi)/2;

            sortBy1stArg(a1, a2, a3, idxLo, idxMid);
            sortBy1stArg(a1, a2, a3, idxMid + 1, idxHi);
            merge(a1, a2, a3, idxLo, idxMid, idxHi);
        }
    }

    /**
     * sort by increasing value a1 and apply same changes to a2.
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
        	throw new IllegalArgumentException("number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, 0, a1.length - 1);
    }

    /**
     * @param a1 array of points to be sorted
     * @param a2 array of points to apply a1 sorting to also
     * @param a3 array of points to apply a1 sorting to also
     * @param idxLo starting index of sorting of a1, inclusive
     * @param idxHi stopping index of sorting of a1, inclusive
     * @param incrA2 a2 should be sorted by increasing value also
     */
    private static void sortBy1stArgThen2nd(float[] a1, float[] a2, int idxLo, int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi)/2;

            sortBy1stArgThen2nd(a1, a2, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, idxLo, indexMid, idxHi);
        }
    }

    /**
     * merge array
     *
     * @param a1 array of points to be sorted and merged
     * @param a2 array of points to apply a1 changes to also
     * @param idxLo starting index of merging of a1, inclusive
     * @param idxMid mid point index of merging of a1, inclusive
     * @param idxHi stopping index of merging of a1, inclusive
     */
    private static void merge(float[] a1, float[] a2, int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = idxLo + i;
            a2Left[i] = a2[index];
            a1Left[i] = a1[index];
        }
        for (j = 0; j < nRight; j++) {
            index = idxMid + j + 1;
            a2Right[j] = a2[index];
            a1Right[j] = a1[index];
        }

        a2Left[nLeft] = Float.MAX_VALUE;
        a1Left[nLeft] = Float.MAX_VALUE;
        a2Right[nRight] = Float.MAX_VALUE;
        a1Right[nRight] = Float.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            if (a1Left[i] <= a1Right[j]) {
                a2[k] = a2Left[i];
                a1[k] = a1Left[i];
                i += 1;
            } else {
                a2[k] = a2Right[j];
                a1[k] = a1Right[j];
                j += 1;
            }
        }
    }

    /**
     * merge array
     *
     * @param a1 array of points to be sorted and merged
     * @param a2 array of points to apply a1 changes to also
     * @param a3 array of points to apply a1 changes to also
     * @param idxLo starting index of merging of a1, inclusive
     * @param idxMid mid point index of merging of a1, inclusive
     * @param idxHi stopping index of merging of a1, inclusive
     */
    private static void merge(float[] a1, float[] a2, int[] a3, int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];
        int[] a3Left = new int[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];
        int[] a3Right = new int[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = idxLo + i;
            a2Left[i] = a2[index];
            a1Left[i] = a1[index];
            a3Left[i] = a3[index];
        }
        for (j = 0; j < nRight; j++) {
            index = idxMid + j + 1;
            a2Right[j] = a2[index];
            a1Right[j] = a1[index];
            a3Right[j] = a3[index];
        }

        a2Left[nLeft] = Float.MAX_VALUE;
        a1Left[nLeft] = Float.MAX_VALUE;
        a3Left[nLeft] = Integer.MAX_VALUE;
        a2Right[nRight] = Float.MAX_VALUE;
        a1Right[nRight] = Float.MAX_VALUE;
        a3Right[nRight] = Integer.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            if (a1Left[i] <= a1Right[j]) {
                a2[k] = a2Left[i];
                a1[k] = a1Left[i];
                a3[k] = a3Left[i];
                i += 1;
            } else {
                a2[k] = a2Right[j];
                a1[k] = a1Right[j];
                a3[k] = a3Right[j];
                j += 1;
            }
        }
    }

    /**
     * @param a1 array of points to be sorted and merged
     * @param a2 array of points to apply a1 changes to also
     * @param idxLo first index of first subarray for sorting
     * @param idxMid first index of second subarray for sorting
     * @param idxHi last index of second subarray for sorting
     */
    private static void mergeBy1stArgThen2nd( float[] a1, float[] a2, int idxLo, int idxMid, int idxHi) {

                                                            // cost  times          where n'=(indexHi - indexLo)
                                                            // c01   1
                                                            // c02   1

                                                            // c03   4 * n'/2       for assignment to the subsequent 4 arrays

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

                                                                                     // where M' + P' <= N'; P' is when L==R; M'<N'; P'<N'
        for (int k = idxLo; k <= idxHi; k++) {                                       // c04   n'
            float l = a1Left[leftPos];                                               // c04   n' * 1
            float r = a1Right[rightPos];                                             // c04   n' * 1

            if (l == r) {                                                            // c04   n' * 1
                // compare the a2 values.  lowest a2 value should be moved to the left
                float lx = a2Left[leftPos];                                                // c04   p'
                float rx = a2Right[rightPos];                                               // c04   p'

                if (lx <= rx) {                                                      // 4*c04   p' for total of either block
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    rightPos++;
                }
            
            } else if (l < r) {                                                      // 4*c04   m'-n' for total of either block
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
