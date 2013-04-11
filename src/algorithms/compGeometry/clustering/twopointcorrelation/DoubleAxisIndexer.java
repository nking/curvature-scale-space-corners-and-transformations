package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.sorting.MultiArrayMergeSort;
import java.util.Arrays;

/**
 * Class to hold x and y array of points whose values have been
 * indexed by increasing x and y.  The class is useful for finding
 * the intersection of regions of points.
 *
 * The error arrays are also stored as member variables to keep them
 * ordered parallel to the x and y arrays, though they are not part of the
 * 'behavior' of this class.

   Usage:
      DoubleAxisIndexer indexer = new DoubleAxisIndexer();

      // sort and index if errors do not need to be stored:
      indexer.sortAndIndexXThenY(xPoints, yPoints, nXYPoints);

      // else sort and index using method with errors
      indexer.sortAndIndexXThenY(xPoints, yPoints, xPointErrors, yPointErrors, nXYPoints);

      // find points within a rectangle's bounds:
      int[] regionIndexes = indexer.findIntersectingRegionIndexes(xIndexLo, xIndexHi, yIndexLo, yIndexHi);


 * @author nichole
 */
public class DoubleAxisIndexer {

    /**
     * a copy of the original xPoints array given to the code before sorting.
     * it is used to find values using indexes given by sortedXIndexes
     */
    protected float[] x = null;

    /**
     * a copy of the original yPoints array given to the code before sorting.
     * it is used to find values using indexes given by sortedYIndexes
     */
    protected float[] y = null;

    /**
     * a copy of the original xPointErrors array given to the code.  these are
     * parallel to the array xPoints
     */
    private float[] xErrors = null;

    /**
     * a copy of the original yPointErrors array given to the code.  these are
     * parallel to the array yPoints
     */
    private float[] yErrors = null;

    /**
     * array whose values hold indexes of x.  the order of items in sortedXIndexes
     * is order of increasing value of x, in other words minimum x is in sortedXIndexes[0]
     * and maximum x is in sortedXIndexes[x.length - 1]
     */
    protected int[] sortedXIndexes = null;

    /**
     * array whose values hold indexes of y.  the order of items in sortedYIndexes
     * is order of increasing value of y, in other words minimum y is in sortedYIndexes[0]
     * and maximum y is in sortedYIndexes[y.length - 1]
     */
    protected int[] sortedYIndexes = null;

    /**
     * number of points in arrays x or y
     */
    protected int nXY;

    /**
     * a copy of the reference to xPoints which contains data which have now been altered by
     * sorting and are in the state of 'sorted in increasing order by Y' along with array ySortedByY
     */
    protected float[] xSortedByY = null;

    /**
     * a copy of the reference to yPoints which contains data which have now been altered by
     * sorting and are in the state of 'sorted in increasing order by Y' along with array xSortedByY
     */
    protected float[] ySortedByY = null;

    public DoubleAxisIndexer() {
    }

    public long approximateMemoryUsed() {

        int n = (x == null) ? 0 : x.length;

        long sumBytes = 6*16;

        long sumBits = (6*n*32) + 32;

        if (xErrors != null) {
            sumBits += (2*xErrors.length * 32);
        }

        sumBytes += (sumBits/8);

        // amount of padding needed to make it a round 8 bytes
        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    /**
     * sort the pair (xPoints, yPoints) by x and index them then sort by y and
     * index them for both sorts.
     *
     * The original unsorted arrays are available as copies to use the indexes
     * with as this.x and this.y, but xPoints and yPoints are altered by this method.
     *
     * This method is usually called once after construction of the indexer and
     * thereafter the member variables are used.
     *
     * @param xPoints
     * @param yPoints
     * @param xPointErrors
     * @param yPointErrors
     * @param nPoints
     */
    public void sortAndIndexXThenY(float[] xPoints, float[] yPoints, float[] xPointErrors, float[] yPointErrors, int nPoints) {

        this.x = Arrays.copyOf(xPoints, nPoints);
        this.y = Arrays.copyOf(yPoints, nPoints);
        this.xErrors = Arrays.copyOf(xPointErrors, nPoints);
        this.yErrors = Arrays.copyOf(yPointErrors, nPoints);

        this.nXY = nPoints;

        this.sortedXIndexes = new int[nXY];
        for (int i = 0; i < nXY; i++) {
            sortedXIndexes[i] = i;
        }

        MultiArrayMergeSort.sortByY(yPoints, xPoints, sortedXIndexes, nXY);

        this.sortedYIndexes = Arrays.copyOf(this.sortedXIndexes, nXY);

        MultiArrayMergeSort.sortByY(xPoints, yPoints, sortedYIndexes, nXY);

        this.xSortedByY = xPoints;
        this.ySortedByY = yPoints;
    }

    /**
     * sort the pair (xPoints, yPoints) by x and index them then sort by y and
     * index them for both sorts.
     *
     * The original unsorted arrays are available as copies to use the indexes
     * with as this.x and this.y, but xPoints and yPoints are altered by this method.
     *
     * @param xPoints
     * @param yPoints
     * @param nPoints
     */
    public void sortAndIndexXThenY(float[] xPoints, float[] yPoints, int nPoints) {

        this.x = Arrays.copyOf(xPoints, nPoints);
        this.y = Arrays.copyOf(yPoints, nPoints);
        this.xErrors = null;
        this.yErrors = null;

        this.nXY = nPoints;

        this.sortedXIndexes = new int[nXY];
        for (int i = 0; i < nXY; i++) {
            sortedXIndexes[i] = i;
        }

        MultiArrayMergeSort.sortByY(yPoints, xPoints, sortedXIndexes, nXY);

        this.sortedYIndexes = Arrays.copyOf(this.sortedXIndexes, nXY);

        MultiArrayMergeSort.sortByY(xPoints, yPoints, sortedYIndexes, nXY);

        this.xSortedByY = xPoints;
        this.ySortedByY = yPoints;
    }

    /**
     * find the intersection of regions specified by xIndexLo:xIndexHi and yIndexLo:yIndexHi
     *
     * @param xIndexLo x index w.r.t the sortedXIndexes. used to lookup the
     * index w.r.t. the original unsorted array
     * @param xIndexHi
     * @param yIndexLo y index w.r.t the sortedYIndexes. used to lookup the
     * index w.r.t. the original unsorted array
     * @param yIndexHi
     * @return
     */
    protected int[] findIntersectingRegionIndexes(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi) {

        // since these points are sorted by Y can pick out the x first, and use binary search to find y points
        int n = (xIndexHi - xIndexLo) + (yIndexHi - yIndexLo) + 2;
        if (n < 0) {
            int  z = 1;
        }
        int[] indexes = new int[n];
        int nIndexes = 0;

        for (int i = xIndexLo; i <= xIndexHi; i++) {

            int xIndex = sortedXIndexes[i];

            float xp = x[xIndex];
            float yp = y[xIndex];

            // find the corresponding y index and see if it is in range yIndexLo to yIndexHi
            // where is xIndex in ySortedIndexes?
            int yIndexSorted = Arrays.binarySearch(ySortedByY, yp);

            if ((yIndexSorted >= yIndexLo) && (yIndexSorted <= yIndexHi)) {
                indexes[nIndexes] = xIndex;
                nIndexes++;
            }
        }

        return Arrays.copyOf(indexes, nIndexes);
    }

    /**
     * find the intersection of regions specified by xIndexLo:xIndexHi and yIndexLo:yIndexHi
     * if there are only 2 points, else return an empty array
     *
     * @param xIndexLo x index w.r.t the sortedXIndexes. used to lookup the
     * index w.r.t. the original unsorted array
     * @param xIndexHi
     * @param yIndexLo y index w.r.t the sortedYIndexes. used to lookup the
     * index w.r.t. the original unsorted array
     * @param yIndexHi
     * @param useCompleteSampling if set to true, an exception to returning only
     *     2 point intersecting regions is made and the method might return 3 points
     *     if 2 of the points have the same y.
     * @return
     */
    protected int[] findIntersectingRegionIndexesIfOnlyTwo(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, boolean useCompleteSampling) {

        // since these points are sorted by Y can pick out the x first, and use binary search to find y points
        int n = (xIndexHi - xIndexLo) + (yIndexHi - yIndexLo) + 2;

        int[] indexes = new int[n];
        int nIndexes = 0;

        for (int i = xIndexLo; i <= xIndexHi; i++) {

            int xIndex = sortedXIndexes[i];

            float xp = x[xIndex];
            float yp = y[xIndex];

            // find the corresponding y index and see if it is in range yIndexLo to yIndexHi
            // where is xIndex in ySortedIndexes.  assumes unique y values.
            int yIndexSorted = Arrays.binarySearch(ySortedByY, yp);

            if ((yIndexSorted >= yIndexLo) && (yIndexSorted <= yIndexHi)) {
                indexes[nIndexes] = xIndex;
                nIndexes++;

                if (nIndexes > 2) {

                    // more than one of same y value can be a problem here.
                    // choices:
                    //   (1) ignore that and toss the region as containing more than 2 points.
                    //   (2) search here for the identical y in indexes, and if found, allow
                    //       algorithm to continue.  then discard if more than 3 are found in the region.
                    //       the code that uses this would have to then try to store both pair surface densities
                    //       (this later adds alot to each iteration.  will use it only if complete sampling is selected)
                    //   ==> choosing to implement (2) if useCompleteSampling is true.

                    if (useCompleteSampling) {
                        // if there is more than one of same y value, there can only be that many + 1 total, in other
                        //  words, only one other different y valued point in there,
                        // OR, could be a rectangle, with 2 points having same x's and 2 having same y's
                        boolean isATriangleWithTwoSameYs = isATriangleWithTwoSameYs(indexes, nIndexes, yp);
                        if (isATriangleWithTwoSameYs) {
                            // let algorithm continue
                        } else if (isARectangle(indexes, nIndexes)) {
                            // let algorithm continue
                        } else {
                            return new int[0];
                        }
                    } else {
                        return new int[0];
                    }
                }
            }
        }

        return Arrays.copyOf(indexes, nIndexes);
    }

    /**
     * test whether the points found by the indexes are 3 in number with 2 having the
     * same value as yp (which is the latest entered, so faster to use that).
     *
     * @param indexes
     * @param nIndexes
     * @param yp
     * @return
     */
    protected boolean isATriangleWithTwoSameYs(int[] indexes, int nIndexes, float yp) {
        int countSameYValue = 0;
        for (int ii = 0; ii < nIndexes; ii++) {
            float yt = y[ indexes[ii]];
            if (yt == yp) {
                countSameYValue++;
            }
        }
        if (countSameYValue == (nIndexes - 1)) {
            return true;
        }
        return false;
    }
    protected boolean isARectangle(int[] indexes, int nIndexes) {
        // quickest way to rule out rectangle?  count points with same x's and same y's and should have 2 each
        if (nIndexes != 4) {
            return false;
        }
        int nSameX = 0;
        int nSameY = 0;
        int[] sameXIndexes0 = new int[2]; // a value in sameXIndexes0 is paired with same index in sameXIndexes1
        int[] sameXIndexes1 = new int[2];
        int[] sameYIndexes0 = new int[2];
        int[] sameYIndexes1 = new int[2];
        for (int i = 0; i < nIndexes; i++) {
            float xti = x[ indexes[i] ];
            float yti = y[ indexes[i] ];
            for (int j = (i + 1); j < nIndexes; j++) {
                float xtj = x[ indexes[j] ];
                float ytj = y[ indexes[j] ];
                float diffx = xti - xtj;
                float diffy = yti - ytj;
                if (Math.abs(diffx) <= 0.001*xti) {
                    // [0]  nSameX->0:1
                    // [1]  nSameX->1:2
                    if (nSameX > 1) {
                        return false;
                    } else {
                        sameXIndexes0[nSameX] = i;
                        sameXIndexes1[nSameX] = j;
                        nSameX++;
                    }
                } else if (Math.abs(diffy) <= 0.001*yti) {
                    if (nSameY > 1) {
                        return false;
                    } else {
                        sameYIndexes0[nSameY] = i;
                        sameYIndexes1[nSameY] = j;
                        nSameY++;
                    }
                }
            }
        }
        if (nSameX != 2) {
            return false;
        } else if (nSameY != 2) {
            return false;
        } else {
            // there were 4 points, and we found 2 with same x and 2 with same y and none
            // are identical points.
            // now need to assert that the sets include pattern of points.
            int topYIndex, bottomYIndex, leftXIndex, rightXIndex;
            if (y[sameYIndexes0[0]] > y[sameYIndexes0[1]]) {
                topYIndex = 0;
                bottomYIndex = 1;
            } else {
                topYIndex = 1;
                bottomYIndex = 0;
            }
            if (x[sameXIndexes0[0]] < x[sameXIndexes0[1]]) {
                leftXIndex = 0;
                rightXIndex = 1;
            } else {
                leftXIndex = 1;
                rightXIndex = 0;
            }

            // 1 |1   3      sameYIndexes0[1] = 1   sameYIndexes1[1] = 3
            //   |
            // 0 |0   2      sameYIndexes0[0] = 0   sameYIndexes1[0] = 2
            //   ------
            //    0   1
            //
            //   /\
            //    |
            // sameXIndexes0[0] = 0   sameXIndexes1[0] = 1
            // sameXIndexes0[1] = 2   sameXIndexes1[0] = 3

            // do the top and left contain same index?
            if (sameYIndexes0[topYIndex] == sameXIndexes0[leftXIndex]) {
                // found top left
                // is there a top right with same index?
                if (sameYIndexes1[topYIndex] == sameXIndexes0[rightXIndex]) {
                    // found top left  and top right  since we know 2x's are same and 2y's are same, we have a rectangle now
                    return true;
                } else if (sameYIndexes1[topYIndex] == sameXIndexes1[rightXIndex]) {
                    return true;
                }
            } else if (sameYIndexes1[topYIndex] == sameXIndexes0[leftXIndex]) {
                // found top left
                // is there a top right with same index?
                if (sameYIndexes0[topYIndex] == sameXIndexes0[rightXIndex]) {
                    // found top left  and top right  since we know 2x's are same and 2y's are same, we have a rectangle now
                    return true;
                } else if (sameYIndexes0[topYIndex] == sameXIndexes1[rightXIndex]) {
                    return true;
                }
            } else if (sameYIndexes0[topYIndex] == sameXIndexes1[leftXIndex]) {
                // found top left
                // is there a top right with same index?
                if (sameYIndexes1[topYIndex] == sameXIndexes0[rightXIndex]) {
                    // found top left  and top right  since we know 2x's are same and 2y's are same, we have a rectangle now
                    return true;
                } else if (sameYIndexes1[topYIndex] == sameXIndexes1[rightXIndex]) {
                    return true;
                }
            } else if (sameYIndexes1[topYIndex] == sameXIndexes1[leftXIndex]) {
                // found top left
                // is there a top right with same index?
                if (sameYIndexes0[topYIndex] == sameXIndexes0[rightXIndex]) {
                    // found top left  and top right  since we know 2x's are same and 2y's are same, we have a rectangle now
                    return true;
                } else if (sameYIndexes0[topYIndex] == sameXIndexes1[rightXIndex]) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * find the min and max of the indexes in xSortedByY and ySortedByY specified
     * by xyIndexes
     *
     * @param xyIndexes
     * @return
     */
    public float[] findXYMinMax(int[] xyIndexes) {

        float xmin = Float.MAX_VALUE;
        float ymin = xmin;
        float xmax = Float.MIN_VALUE;
        float ymax = xmax;

        for (int index : xyIndexes) {
            float xp = xSortedByY[index];
            float yp = ySortedByY[index];

            if (xp < xmin) {
                xmin = xp;
            }
            if (yp < ymin) {
                ymin = yp;
            }
            if (xp > xmax) {
                xmax = xp;
            }
            if (yp > ymax) {
                ymax = yp;
            }
        }

        return new float[]{xmin, xmax, ymin, ymax};
    }

    public float[] findXYMinMax() {

        int index = sortedXIndexes[0];
        float xmin = x[index];

        index = sortedXIndexes[x.length - 1];
        float xmax = x[index];

        index = sortedYIndexes[0];
        float ymin = y[index];

        index = sortedYIndexes[y.length - 1];
        float ymax = y[index];

        return new float[]{xmin, xmax, ymin, ymax};
    }

    /**
     * @return the sortedXIndexes
     */
    public int[] getSortedXIndexes() {
        return sortedXIndexes;
    }

    /**
     * @return the sortedYIndexes
     */
    public int[] getSortedYIndexes() {
        return sortedYIndexes;
    }

    /**
     * return the unsorted x array that the indexes refer to
     *
     * @return the x
     */
    public float[] getX() {
        return x;
    }

    /**
     * return the unsorted y array that the indexes refer to
     *
     * @return the y
     */
    public float[] getY() {
        return y;
    }

    /**
     * return reference to the original array which is now sorted by y
     *
     * @return the x
     */
    public float[] getXSortedByY() {
        return xSortedByY;
    }

    /**
     * return reference to the original array which is now sorted by y
     *
     * @return the y
     */
    public float[] getYSortedByY() {
        return ySortedByY;
    }

    protected int findIndexForValue(int[] array, int value) {
        for (int i = 0; i < array.length; i++) {
            if (array[i] == value) {
                return i;
            }
        }
        return -1;
    }

    /**
     * @return the xErrors
     */
    public float[] getXErrors() {
        return xErrors;
    }

    /**
     * @return the yErrors
     */
    public float[] getYErrors() {
        return yErrors;
    }

    public int getNumberOfPoints() {
        return nXY;
    }
}
