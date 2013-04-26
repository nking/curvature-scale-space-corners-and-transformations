package algorithms.sorting;

/**
 * merge sort worse case runtime is O(N * log_2(N))
 *
 * The method to sort 2 arrays below uses divide and conquer paradigm and sorts one additional array
 * (with no additional logic) along with the original array which is sorted.
 * the result is a constant factor added to an internal constant for the O(n) for merge of single
 * sort.
 *
 * The method to sort 2 arrays by one and then the other, sortByYThenX, is approximated as:
 *
 * the algorithm does not sort in place, so does temporarily take up more memory as it sorts.
 *
 * constructed from pseudo-code in Cormen et al.
 * "Introduction to Algorithms"
 *
 * @author Nichole King
 */
public class MultiArrayMergeSort {

	/**
     * sort array using mergesort (divide, conquer, and combine) sorting both by increasing y
     * (the x array is kept parallel to y at all times).
     * runtime is O(N * log_2(N))
     *
     * @param x of x points
     * @param y array of y points
     */
    public static void sortByY(float[] x, float[] y) {
        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
        	throw new IllegalArgumentException("number of items in x must be the same as in y");
        }
        sortByY(x, y, 0, x.length - 1);
    }

    /**
     * sort array using mergesort (divide, conquer, and combine) sorting both by increasing y
     * (the x array is kept parallel to y at all times).
     * runtime is O(N * log_2(N))
     *
     * @param x of x points
     * @param y array of y points
     * @param nLength number of items in x and y to sort, starting with item 0
     */
    public static void sortByY(float[] x, float[] y, int nLength) {
        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length < nLength) {
        	throw new IllegalArgumentException("number of items in x must be at least nLength");
        }
        if (y.length < nLength) {
        	throw new IllegalArgumentException("number of items in y must be at least nLength");
        }
        sortByY(x, y, 0, nLength - 1);
    }

    public static void sortByY(float[] x, float[] y, int[] a, int[] b, int nLength) {
        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length < nLength) {
        	throw new IllegalArgumentException("number of items in x must be at least nLength");
        }
        if (y.length < nLength) {
        	throw new IllegalArgumentException("number of items in y must be at least nLength");
        }
        sortByY(x, y, a, b, 0, nLength - 1);
    }

    public static void sortByY(float[] x, float[] y, int[] a, int nLength) {
        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length < nLength) {
        	throw new IllegalArgumentException("number of items in x must be at least nLength");
        }
        if (y.length < nLength) {
        	throw new IllegalArgumentException("number of items in y must be at least nLength");
        }
        sortByY(x, y, a, 0, nLength - 1);
    }

    public static void sortByY(float[] x, float[] y, int indexLo, int indexHi) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi)/2;

            sortByY(x, y, indexLo, indexMid);
            sortByY(x, y, indexMid + 1, indexHi);
            merge(x, y, indexLo, indexMid, indexHi);
        }
    }

    public static void sortByY(float[] x, float[] y, int[] a, int[] b, int indexLo, int indexHi) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi)/2;

            sortByY(x, y, a, b, indexLo, indexMid);
            sortByY(x, y, a, b, indexMid + 1, indexHi);
            merge(x, y, a, b, indexLo, indexMid, indexHi);
        }
    }

    public static void sortByY(float[] x, float[] y, int[] a, int indexLo, int indexHi) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi)/2;

            sortByY(x, y, a, indexLo, indexMid);
            sortByY(x, y, a, indexMid + 1, indexHi);
            merge(x, y, a, indexLo, indexMid, indexHi);
        }
    }

    public static void sortByYThenX(float[] x, float[] y) {
        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
        	throw new IllegalArgumentException("number of items in x must be the same as in y");
        }
        sortByYThenX(x, y, 0, x.length - 1, true);
    }

    /**
     * @param x array
     * @param y array
     * @param indexLo first index of first subarray for sorting
     * @param indexMid first index of second subarray for sorting
     * @param indexHi last index of second subarray for sorting
     * @param incrX sort by increasing x.  if false, will use decreasing x.
     */
    private static void sortByYThenX(float[] x, float[] y, int indexLo, int indexHi, boolean incrX) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi)/2;

            sortByYThenX(x, y, indexLo, indexMid, incrX);
            sortByYThenX(x, y, indexMid + 1, indexHi, incrX);
            mergeYThenX(x, y, indexLo, indexMid, indexHi, incrX);
        }
    }

    private static void merge( float[] x, float[] y, int indexLo, int indexMid, int indexHi) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        float[] yLeft = new float[nLeft + 1];
        float[] xLeft = new float[nLeft + 1];

        float[] yRight = new float[nRight + 1];
        float[] xRight = new float[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            yLeft[i] = y[index];
            xLeft[i] = x[index];
        }
        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            yRight[j] = y[index];
            xRight[j] = x[index];
        }

        yLeft[nLeft] = Float.MAX_VALUE;
        xLeft[nLeft] = Float.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = indexLo; k <= indexHi; k++) {
            if (yLeft[i] <= yRight[j]) {
                y[k] = yLeft[i];
                x[k] = xLeft[i];
                i += 1;
            } else {
                y[k] = yRight[j];
                x[k] = xRight[j];
                j += 1;
            }
        }
    }

    private static void merge( float[] x, float[] y, int[] a, int[] b, int indexLo, int indexMid, int indexHi) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        float[] yLeft = new float[nLeft + 1];
        float[] xLeft = new float[nLeft + 1];
        int[] aLeft = new int[nLeft + 1];
        int[] bLeft = new int[nLeft + 1];

        float[] yRight = new float[nRight + 1];
        float[] xRight = new float[nRight + 1];
        int[] aRight = new int[nRight + 1];
        int[] bRight = new int[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            yLeft[i] = y[index];
            xLeft[i] = x[index];
            aLeft[i] = a[index];
            bLeft[i] = b[index];
        }
        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            yRight[j] = y[index];
            xRight[j] = x[index];
            aRight[j] = a[index];
            bRight[j] = b[index];
        }

        yLeft[nLeft] = Float.MAX_VALUE;
        xLeft[nLeft] = Float.MAX_VALUE;
        aLeft[nLeft] = Integer.MAX_VALUE;
        bLeft[nLeft] = Integer.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;
        aRight[nRight] = Integer.MAX_VALUE;
        bRight[nRight] = Integer.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = indexLo; k <= indexHi; k++) {
            if (yLeft[i] <= yRight[j]) {
                y[k] = yLeft[i];
                x[k] = xLeft[i];
                a[k] = aLeft[i];
                b[k] = bLeft[i];
                i += 1;
            } else {
                y[k] = yRight[j];
                x[k] = xRight[j];
                a[k] = aRight[j];
                b[k] = bRight[j];
                j += 1;
            }
        }
    }

    private static void merge( float[] x, float[] y, int[] a, int indexLo, int indexMid, int indexHi) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        float[] yLeft = new float[nLeft + 1];
        float[] xLeft = new float[nLeft + 1];
        int[] aLeft = new int[nLeft + 1];

        float[] yRight = new float[nRight + 1];
        float[] xRight = new float[nRight + 1];
        int[] aRight = new int[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            yLeft[i] = y[index];
            xLeft[i] = x[index];
            aLeft[i] = a[index];
        }
        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            yRight[j] = y[index];
            xRight[j] = x[index];
            aRight[j] = a[index];
        }

        yLeft[nLeft] = Float.MAX_VALUE;
        xLeft[nLeft] = Float.MAX_VALUE;
        aLeft[nLeft] = Integer.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;
        aRight[nRight] = Integer.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = indexLo; k <= indexHi; k++) {
            if (yLeft[i] <= yRight[j]) {
                y[k] = yLeft[i];
                x[k] = xLeft[i];
                a[k] = aLeft[i];
                i += 1;
            } else {
                y[k] = yRight[j];
                x[k] = xRight[j];
                a[k] = aRight[j];
                j += 1;
            }
        }
    }

    /**
     * @param x array
     * @param y array
     * @param indexLo first index of first subarray for sorting
     * @param indexMid first index of second subarray for sorting
     * @param indexHi last index of second subarray for sorting
     * @param incrX sort by increasing x.  if false, will use decreasing x.
     */
    private static void mergeYThenX( float[] x, float[] y, int indexLo, int indexMid, int indexHi, boolean incrX) {

                                                                // cost  times          where n'=(indexHi - indexLo)
        int nLeft = indexMid - indexLo + 1;                     // c01   1
        int nRight = indexHi - indexMid;                        // c02   1

        float[] yLeft = new float[nLeft + 1];                 // c03   4 * n'/2       for assignment to the subsequent 4 arrays
        float[] xLeft = new float[nLeft + 1];

        float[] yRight = new float[nRight + 1];
        float[] xRight = new float[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            yLeft[i] = y[index];
            xLeft[i] = x[index];
        }
        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            yRight[j] = y[index];
            xRight[j] = x[index];
        }

        yLeft[nLeft] = Float.MAX_VALUE;
        xLeft[nLeft] = Float.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;

        i = 0;
        j = 0;

                                                                                     // where M' + P' <= N'; P' is when L==R; M'<N'; P'<N'
        for (int k = indexLo; k <= indexHi; k++) {                                   // c04   n'
            float l = yLeft[i];                                                     // c04   n' * 1
            float r = yRight[j];                                                    // c04   n' * 1

            if (l == r) {                                                            // c04   n' * 1
                // compare the x values.  lowest x value should be moved to the left
                float lx = xLeft[i];                                                // c04   p'
                float rx = xRight[j];                                               // c04   p'

                // if sort by decreasing, we can reverse the lx and rx
                if (!incrX) {                                                        // c04   p'
                    lx = rx;                                                         // c04   p'
                    rx = xLeft[i];                                                   // c04   p'
                }

                if (lx <= rx) {                                                      // 4*c04   p' for total of either block
                    y[k] = yLeft[i];
                    x[k] = xLeft[i];
                    i += 1;
                } else {
                    y[k] = yRight[j];
                    x[k] = xRight[j];
                    j += 1;
                }
            } else if (l < r) {                                                      // 4*c04   m'-n' for total of either block
                y[k] = yLeft[i];
                x[k] = xLeft[i];
                i += 1;
            } else {
                y[k] = yRight[j];
                x[k] = xRight[j];
                j += 1;
            }
        }
    }
}
