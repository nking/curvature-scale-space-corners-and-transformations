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
     * and maximum x is in sortedXIndexes[x.length - 1].  
     * the values in sortedXIndexes can be used as indexes in 
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

    public int getNXY() {
        return nXY;
    }

    public long approximateMemoryUsed() {

        /*
        stack word size is 32 or 64 bits.
        Stack holds:  local variables

        Heap holds:  object references (32 bits) and arrays (32 bits x nItems X itemSize)

        float[] x
        float[] y
        float[] xErrors
        float[] yErrors
        int[] sortedXIndexes
        int[] sortedYIndexes
        float[] xSortedByY
        float[] ySortedByY
        int nXY;
                                                     32-bit platform          64-bit platform
                                                Field   Size on          Field   Size on
           Java types                           size    stack            size    stack
           boolean                              32       32              32      64
           byte                                 32       32              32      64
           char                                 32       32              32      64
           short                                32       32              32      64
           int                                  32       32              32      64
           ﬂoat                                 32       32              32      64
           reference                            32       32              64      64
           array reference                      32       32              32      32
           returnAddress                        32       32              64      64
           long                                 64       64              64      128
           double                               64       64              64      128
           from Table II of http://users.elis.ugent.be/~leeckhou/papers/SPE06.pdf

        */

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int n = (x == null) ? 0 : x.length;

        int nbits = (is32Bit) ? 32 : 64;

        int arrayRefBits = 32;

        int overheadBytes = 16;

        // 6 float and int arrays on the heap
        long sumBits = 6*(arrayRefBits + (n*32));

        if (xErrors != null) {
            // 2 float arrays on the heap
            sumBits += (2*(arrayRefBits + (xErrors.length * 32)));
        }

        long sumBytes = (sumBits/8) + overheadBytes;

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

        MultiArrayMergeSort.sortBy1stArg(xPoints, yPoints, sortedXIndexes, nXY);

        this.sortedYIndexes = Arrays.copyOf(this.sortedXIndexes, nXY);

        MultiArrayMergeSort.sortBy1stArg(yPoints, xPoints, sortedYIndexes, nXY);
        
        //sortedYIndexes hold values that are indexes w.r.t the y array

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

        MultiArrayMergeSort.sortBy1stArg(xPoints, yPoints, sortedXIndexes, nXY);

        this.sortedYIndexes = Arrays.copyOf(this.sortedXIndexes, nXY);

        MultiArrayMergeSort.sortBy1stArg(yPoints, xPoints, sortedYIndexes, nXY);

        this.xSortedByY = xPoints;
        this.ySortedByY = yPoints;
    }
     
    /**
     * given 2 points defined by their indexes with respect to the
     * array sortedXIndexes, return true if there are no other points
     * within their bounds.
     *
     * @param xSortedIndex1 index of first point's location in the array sortedXIndexes
     * @param xSortedIndex2 index of second point's location in the array sortedXIndexes
     *
     * @return true if there are no other points within the bounds of the 2 referenced
     * by their indexes
     */
    protected boolean hasNoOtherPointsWithinBounds(int xSortedIndex1, int xSortedIndex2) {
    
        if (xSortedIndex1 > xSortedIndex2) {
            int tmp = xSortedIndex1;
            xSortedIndex1 = xSortedIndex2;
            xSortedIndex2 = tmp;
        }
        
        int idx1 = sortedXIndexes[xSortedIndex1];
        int idx2 = sortedXIndexes[xSortedIndex2];
        
        // x's are already ordered by increasing x
        float y0 = y[idx1];
        float y1 = y[idx2];
       
        if (y0 > y1) {
            float tmp = y0;
            y0 = y1;
            y1 = tmp;
        }
        
        for (int i = (xSortedIndex1 + 1); i < xSortedIndex2; i++) {
            
            int idx = sortedXIndexes[i];

            // we already know that x is within bounds, just need to test yt and exit quickly if it is within bounds
            float yt = y[idx];
            
            if ( !((yt < y0) || (yt > y1)) ) {
                // it's not outside of bounds
                return false;
            }
        }
        
        return true;
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
     * return the index w.r.t. indexer.sortedXIndexes for the yValue.
     * the index frame of reference is needed for a subsequent use.
     * 
     * @param yValue
     * @param sortedX
     * @return
     */
    public int findSortedXIndexesForY(float yValue) {
        
        int idx = Arrays.binarySearch(ySortedByY, yValue);
        
        if (idx < 0) {
            idx = -1*(idx + 1);
        }
        
        if (idx > (nXY - 1)) {
            idx = nXY - 1;
        }
        
        int idx2 = sortedYIndexes[idx];
        
        // now we have idx2 as the index w.r.t indexer.x and indexer.y, 
        // so just need to find idx2 value in indexer.sortedXIndexes and that's the answer
                    
        for (int ii = 0; ii < sortedXIndexes.length; ii++) {
            if (idx2 == sortedXIndexes[ii]) {
                return ii;
            }
        }
        
        throw new IllegalStateException("could not find index " + idx + " in indexer.sortedXIndexes");
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
        
        if (array.equals(this.ySortedByY)) {
            return Arrays.binarySearch(ySortedByY, value);  
        } 
        
        for (int i = 0; i < array.length; i++) {
            if (array[i] == value) {
                return i;
            }
        }
        return -1;
    }

    /**
     * get the errors in x.  this array is parallel to that returned by getX()
     * @return the xErrors
     */
    public float[] getXErrors() {
        return xErrors;
    }

    /**
     * get the errors in y.  this array is parallel to that returned by getY()
     *
     * @return the yErrors
     */
    public float[] getYErrors() {
        return yErrors;
    }

    public int getNumberOfPoints() {
        return nXY;
    }

}
