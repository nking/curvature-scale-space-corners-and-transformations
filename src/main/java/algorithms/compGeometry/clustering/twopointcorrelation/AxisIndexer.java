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
      AxisIndexer indexer = new AxisIndexer();

      // sort and index if errors do not need to be stored:
      indexer.sortAndIndexXThenY(xPoints, yPoints, nXYPoints);

      // else sort and index using method with errors
      indexer.sortAndIndexXThenY(xPoints, yPoints, xPointErrors, yPointErrors, nXYPoints);

 * @author nichole
 */
public class AxisIndexer {

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
     * number of points in arrays x or y
     */
    protected int nXY;

    protected float xmin = Float.MAX_VALUE;
    
    protected float ymin = xmin;
    
    protected float xmax = Float.MIN_VALUE;
    
    protected float ymax = xmax;
    
    public AxisIndexer() {
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

        // 4 float and 1 int arrays on the heap
        long sumBits = 5*(arrayRefBits + (n*nbits));

        if (xErrors != null) {
            // 2 float arrays on the heap
            sumBits += (2*(arrayRefBits + (xErrors.length * nbits)));
        }

        // 5 primitives
        sumBits += 5*nbits;

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
    public void sortAndIndexX(float[] xPoints, float[] yPoints, float[] xPointErrors, float[] yPointErrors, int nPoints) {

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
        
        xmin = xPoints[0];
        
        xmax = xPoints[xPoints.length - 1];
        
        for (float yp : yPoints) {
            if (yp > ymax) {
                ymax = yp;
            }
            if (yp < ymin) {
                ymin = yp;
            }
        }
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
    public void sortAndIndexX(float[] xPoints, float[] yPoints, int nPoints) {

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

        xmin = xPoints[0];
        
        xmax = xPoints[xPoints.length - 1];
        
        for (float yp : yPoints) {
            if (yp > ymax) {
                ymax = yp;
            }
            if (yp < ymin) {
                ymin = yp;
            }
        }
    }
    
    public float[] findXYMinMax() {

        return new float[]{xmin, xmax, ymin, ymax};
    }
 
    /**
     * @return the sortedXIndexes
     */
    public int[] getSortedXIndexes() {
        return sortedXIndexes;
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
