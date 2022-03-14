package algorithms.compGeometry.convexHull;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.sort.CountingSort;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.List;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 * 
 * The algorithm uses a merge sort which has a worse case runtime is O(N * log_2(N)),
 * but also includes an additional set of operations to remove inner points with the
 * same polar angle from P0.
 *
 */
public class PolarAngleQuickSort {

    /**
     * sort the given points in place by polar angle in counterclockwise order 
     * around the first point x[0], y[0].
     * NOTE that the angles are rounded to degree integers for the "reduce to
     * unique" step.
     * The runtime complexity is O( max(x.length, 360) ).
     * @param x array of x points
     * @param y array of y points of same length as x.
     */
    public static int sortCCWBy1stPoint(long[] x, long[] y) {

        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        
        if (x.length < 3) {
            return x.length;
        }
        
        long x0 = x[0];
        long y0 = y[0];
        
        long[] x1 = Arrays.copyOfRange(x, 1, x.length);
        long[] y1 = Arrays.copyOfRange(y, 1, y.length);
        
        double degRadians;
        int[] polarAngleDegree1 = new int[x.length - 1];
        
        for (int i = 0; i < x1.length; i++) {
            degRadians = AngleUtil.polarAngleCCW((double)(x1[i] - x0), (double)(y1[i] - y0));
            polarAngleDegree1[i] = (int)Math.round(degRadians * (180./Math.PI));
        }
        
        int[] sortedIndexes1 = CountingSort.sortAndReturnIndexes(polarAngleDegree1);
        
        // re-order x1 and y1 by sortedIndexes1
        rewriteBySortedIndexes(sortedIndexes1, x1);
        rewriteBySortedIndexes(sortedIndexes1, y1);
                
        // for same polar angles, keep the one which is furthest from x0, p0.
        // assuming an angular resolution of 1 degree and using rounded integers for the angle degrees
        int nUsable = reduceToUniquePolarAngles(x0, y0, x1, y1, polarAngleDegree1);
        
        System.arraycopy(x1, 0, x, 1, nUsable);
        System.arraycopy(y1, 0, y, 1, nUsable);
                
        // this method is different from the others in separating the 1st point from
        //  the arrays x1, y1 in the reduce method, so needs to add 1 to the return value
        //  for the first point put back into the array
        return nUsable + 1;
    }
    
    public static <T extends PairInt> int sort(T p0, T[] points) {

        if (p0 == null) {
        	throw new IllegalArgumentException("p0 cannot be null");
        }
        if (points == null) {
        	throw new IllegalArgumentException("points cannot be null");
        }
        
        if (points.length == 1) {
            return 1;
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[points.length];
        
        for (int i = 1; i < points.length; i++) {
            
            polarAngle[i] = AngleUtil.polarAngleCCW(
                (double)(points[i].getX() - p0.getX()), 
                (double)(points[i].getY() - p0.getY()));
        }
        
        sortByPolarAngle(points, 1, points.length - 1, polarAngle);
        
        int nUsable = reduceToUniquePolarAngles(p0, points, polarAngle);
        
        return nUsable;
    }
   
    /**
     * sort list points by polar angle w.r.t. point p0.
     * Note that this sort does not remove any points
     * for having same angle.
     * @param <T>
     * @param p0
     * @param points
     * @return 
     */
    public static <T extends PairInt> int sort2(T p0, List<T> points) {

        if (p0 == null) {
        	throw new IllegalArgumentException("p0 cannot be null");
        }
        if (points == null) {
        	throw new IllegalArgumentException("points cannot be null");
        }
        
        if (points.size() == 1) {
            return 1;
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[points.size()];
        
        for (int i = 1; i < points.size(); i++) {
            
            polarAngle[i] = AngleUtil.polarAngleCCW(
                (double)(points.get(i).getX() - p0.getX()), 
                (double)(points.get(i).getY() - p0.getY()));
        }
        
        sortByPolarAngle(points, 1, points.size() - 1, polarAngle);
                
        return points.size();
    }
    
    static int reduceToUniquePolarAngles(float xP0, float yP0, float[] x, 
        float[] y, double[] polarAngle) {

        double maxDist = Double.NEGATIVE_INFINITY;
        int iMaxDist;
        double dist;
        int nextI;
        int i2 = 1;
        
        double eps = 0;

        for (int i = 1; i < x.length; i++) {

            // look ahead
            nextI = i + 1;
            iMaxDist = i;
            
            if ( (nextI < x.length)  
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                maxDist = relativeLengthOfLine(xP0, yP0, x[i], y[i]);
            }

            while ( (nextI < x.length) 
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                dist = relativeLengthOfLine(xP0, yP0, x[nextI], y[nextI]);
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }
            
            x[i2] = x[iMaxDist];
            y[i2] = y[iMaxDist];
            i = nextI - 1;
            ++i2;
        }
        
        return i2;
    }
    
    static <T extends PairInt> int reduceToUniquePolarAngles(T p0, T[] points, 
        double[] polarAngle) {

        double maxDist = Double.NEGATIVE_INFINITY;
        int iMaxDist;
        double dist;
        int nextI;
        int i2 = 1;
        
        double eps = 0;

        for (int i = 1; i < points.length; i++) {

            // look ahead
            nextI = i + 1;
            iMaxDist = i;
            
            if ( (nextI < points.length)  
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                maxDist = relativeLengthOfLine(p0, points[i]);
            }

            while ( (nextI < points.length) 
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                dist = relativeLengthOfLine(p0, points[nextI]);
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }
            
            points[i2] = points[iMaxDist];
            i = nextI - 1;
            ++i2;
        }
        
        return i2;
    }

    static <T extends PairInt> int reduceToUniquePolarAngles(
        T p0, List<T> points, double[] polarAngle) {
         
        double maxDist = Double.NEGATIVE_INFINITY;
        int iMaxDist;
        double dist;
        int nextI;
        int i2 = 1;
        
        double eps = 0;

        for (int i = 1; i < points.size(); i++) {

            // look ahead
            nextI = i + 1;
            iMaxDist = i;
            
            if ( (nextI < points.size()) 
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                maxDist = relativeLengthOfLine(p0, points.get(i));
            }

            while ( (nextI < points.size()) 
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                dist = relativeLengthOfLine(p0, points.get(nextI));
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }
            
            points.set(i2, points.get(iMaxDist));
            i = nextI - 1;
            ++i2;
        }
        
        return i2;
    }

    static <T extends PairInt> void sortByPolarAngle(T[] a, int idxLo, 
        int idxHi, double[] polarAngle) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {
            int idxMid = partitionByPolarAngle(a, idxLo, idxHi, polarAngle);
            sortByPolarAngle(a, idxLo, idxMid - 1, polarAngle);
            sortByPolarAngle(a, idxMid + 1, idxHi, polarAngle);
        }
    }
    
    static <T extends PairInt> void sortByPolarAngle(List<T> a, int idxLo, 
        int idxHi, double[] polarAngle) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {
            int idxMid = partitionByPolarAngle(a, idxLo, idxHi, polarAngle);
            sortByPolarAngle(a, idxLo, idxMid - 1, polarAngle);
            sortByPolarAngle(a, idxMid + 1, idxHi, polarAngle);
        }
    }
    
    private static <T extends PairInt> int partitionByPolarAngle(T[] a, int idxLo, 
        int idxHi, double[] polarAngle) {
     
        double x = polarAngle[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (polarAngle[i] < x) {
                doSwap = true;
            }
            if (doSwap) {
                store++;
                T swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                double swap2 = polarAngle[store];
                polarAngle[store] = polarAngle[i];
                polarAngle[i] = swap2;
            }
        }
        store++;
        T swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        double swap2 = polarAngle[store];
        polarAngle[store] = polarAngle[idxHi];
        polarAngle[idxHi] = swap2;
        
        return store;
    }
    
    private static <T extends PairInt> int partitionByPolarAngle(List<T> a, int idxLo, 
        int idxHi, double[] polarAngle) {
     
        double x = polarAngle[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (polarAngle[i] < x) {
                doSwap = true;
            }
            if (doSwap) {
                store++;
                T swap = a.get(store);
                a.set(store, a.get(i));
                a.set(i, swap);
                double swap2 = polarAngle[store];
                polarAngle[store] = polarAngle[i];
                polarAngle[i] = swap2;
            }
        }
        store++;
        T swap = a.get(store);
        a.set(store, a.get(idxHi));
        a.set(idxHi, swap);
        double swap2 = polarAngle[store];
        polarAngle[store] = polarAngle[idxHi];
        polarAngle[idxHi] = swap2;
        
        return store;
    }
    
    protected static void sort(float[] x, float[] y, 
        int indexLo, int indexHi, double[] polarAngle) {
        
        if (indexLo < indexHi) {

            int indexMid = (indexLo + indexHi) >> 1;

            sort(x, y, indexLo, indexMid, polarAngle);
            sort(x, y, indexMid + 1, indexHi, polarAngle);
            merge(x, y, indexLo, indexMid, indexHi, polarAngle);
        }
    }
    
    private static void merge(float[] x, float[] y, 
        int indexLo, int indexMid, int indexHi, final double[] polarAngle) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        float[] xLeft = Arrays.copyOfRange(x, indexLo, indexMid + 2);       // add 1 for sentinel
        float[] yLeft = Arrays.copyOfRange(y, indexLo, indexMid + 2);
        double[] angleLeft = Arrays.copyOfRange(polarAngle, indexLo, indexMid + 2);

        float[] xRight = Arrays.copyOfRange(x, indexMid + 1, indexHi + 2);  // add 1 for sentinel
        float[] yRight = Arrays.copyOfRange(y, indexMid + 1, indexHi + 2);
        double[] angleRight = Arrays.copyOfRange(polarAngle, indexMid + 1, indexHi + 2);

        xLeft[nLeft] = Float.MAX_VALUE;
        yLeft[nLeft] = Float.MAX_VALUE;
        angleLeft[nLeft] = Double.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        angleRight[nRight] = Double.MAX_VALUE;

        int i = 0;
        int j = 0;

        for (int k = indexLo; k <= indexHi; k++) {

            if (angleLeft[i] <= angleRight[j]) {

                y[k] = yLeft[i];
                x[k] = xLeft[i];
                polarAngle[k] = angleLeft[i];
                i++;
            } else {

                y[k] = yRight[j];
                x[k] = xRight[j];
                polarAngle[k] = angleRight[j];
                j++;
            }
        }
    }

    static double relativeLengthOfLine(double x1, double y1, double x2, double y2) {
        double dx2 = (x2 - x1);
        dx2 *= dx2;
        double dy2 = (y2 - y1);
        dy2 *= dy2;
        //double d = Math.sqrt(dx2 + dy2);
        return dx2 + dy2;
    }
    
    static <T extends PairInt> double relativeLengthOfLine(T p1, T p2) {
        double dx2 = (p2.getX() - p1.getX());
        dx2 *= dx2;
        double dy2 = (p2.getY() - p1.getY());
        dy2 *= dy2;
        //double d = Math.sqrt(dx2 + dy2);
        return dx2 + dy2;
    }

    private static void sortCCWBy1stPoint(long x0, long y0, long[] x, long[] y, 
        double[] outPolarAngle) {
        
        if (x.length != outPolarAngle.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and outPolarAngle must be same lengths");
        }
        
        for (int i = 0; i < outPolarAngle.length; i++) {
            outPolarAngle[i] = AngleUtil.polarAngleCCW((double)(x[i] - x0), (double)(y[i] - y0));
        }
        
        // sort x, y by outPolarAngle
        sortCCW(x, y, outPolarAngle, 0, x.length - 1);
    }

    private static void sortCCW(long[] x, long[] y, double[] pA, int iLo, int iHi) {
        if (iLo < iHi) {
            int iMid = partitionCCW(x, y, pA, iLo, iHi);
            sortCCW(x, y, pA, iLo, iMid-1);
            sortCCW(x, y, pA, iMid+1, iHi);
        }
    }

    private static int partitionCCW(long[] x, long[] y, double[] pA, int iLo, int iHi) {
        double pAR = pA[iHi];
        int i = iLo - 1;
        long swap;
        double swap2;
        
        for (int k = iLo; k < iHi; ++k) {
            if (pA[k] <= pAR) {
                ++i;
                
                swap2 = pA[k];
                pA[k] = pA[i];
                pA[i] = swap2;
                
                swap = x[k];
                x[k] = x[i];
                x[i] = swap;
                
                swap = y[k];
                y[k] = y[i];
                y[i] = swap;
            }
        }
        ++i;
        swap2 = pA[iHi];
        pA[iHi] = pA[i];
        pA[i] = swap2;
        
        swap = x[iHi];
        x[iHi] = x[i];
        x[i] = swap;
        
        swap = y[iHi];
        y[iHi] = y[i];
        y[i] = swap;
        
        return i;
    }

    /**
     * traverse the polar angles in pA, convert the angles to rounded integer degrees, and if
     * points have the same polar angles only keep the one furthest from (x0, y0).
     * note that all points (x, y) have been sorted by increasing angle in
     * pA inCCW order.
     * the arrays x, y, and pA are compacted so that the usable values are at the
     * top r indexes where r is returned by this method
     * @param x0
     * @param y0
     * @param x
     * @param y
     * @param pA 
     * @return returns the number of indexes usable in each of x, y, and pA
     * after compacting the arrays to remove redundant polar angle degrees.
     */
    private static int reduceToUniquePolarAngles(long x0, long y0, long[] x, 
        long[] y, double[] pA) {
        
        if (x.length != pA.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and pA must be same lengths");
        }
        
        // traverse the list of pA, converting to rounded degrees, and store each
        // pA in a list by index
        
        // then traverse the degree list and if there is only one point with
        //   that angle, store it in (x2,y2,pa2) else compare the distances form (x0,y0) among
        //   those with the same angle and keep the largest distance.
        
        int i;
        int[] deg = new int[pA.length];
        for (i = 0; i < pA.length; ++i) {
            deg[i] = (int)Math.round(pA[i]*(180./Math.PI));
        }
        
        long maxDist = Long.MIN_VALUE;
        int iMaxDist;
        long dist;
        int nextI;
        int i2 = 0;
        
        long xd;
        long yd;
        
        for (i = 0; i < pA.length; ++i) {

            // look ahead
            nextI = i + 1;
            iMaxDist = i;
            
            if ( (nextI < pA.length) && (deg[i] == deg[nextI]) ) {
                xd = x0 - x[i];
                yd = y0 - y[i];
                maxDist = (xd * xd + yd * yd);
            }
            
            while ( (nextI < pA.length) && (deg[i] == deg[nextI]) ) {
                xd = x0 - x[nextI];
                yd = y0 - y[nextI];
                dist = (xd * xd + yd * yd);
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }
            
            x[i2] = x[iMaxDist];
            y[i2] = y[iMaxDist];
            pA[i2] = pA[iMaxDist];
            i = nextI - 1;
            ++i2;            
        }
        
        return i2;
    }
    
    /**
     * traverse the polar angles in pA, convert the angles to rounded integer degrees, and if
     * points have the same polar angles only keep the one furthest from (x0, y0).
     * note that all points (x, y) have been sorted by increasing angle in
     * pA inCCW order.
     * the arrays x, y, and pA are compacted so that the usable values are at the
     * top r indexes where r is returned by this method
     * @param x0
     * @param y0
     * @param x
     * @param y
     * @param deg polar angle in degrees
     * @return returns the number of indexes usable in each of x, y, and pA
     * after compacting the arrays to remove redundant polar angle degrees.
     */
    private static int reduceToUniquePolarAngles(long x0, long y0, long[] x, 
        long[] y, int[] deg) {
        
        if (x.length != deg.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and deg must be same lengths");
        }
        
        // traverse the list of pA, converting to rounded degrees, and store each
        // pA in a list by index
        
        // then traverse the degree list and if there is only one point with
        //   that angle, store it in (x2,y2,pa2) else compare the distances form (x0,y0) among
        //   those with the same angle and keep the largest distance.
        
        long maxDist = Long.MIN_VALUE;
        int iMaxDist;
        long dist;
        int nextI;
        int i2 = 0;
        
        long xd;
        long yd;
        
        int i;
        
        for (i = 0; i < deg.length; ++i) {

            // look ahead
            nextI = i + 1;
            iMaxDist = i;
            
            if ( (nextI < deg.length) && (deg[i] == deg[nextI]) ) {
                xd = x0 - x[i];
                yd = y0 - y[i];
                maxDist = (xd * xd + yd * yd);
            }
            
            while ( (nextI < deg.length) && (deg[i] == deg[nextI]) ) {
                xd = x0 - x[nextI];
                yd = y0 - y[nextI];
                dist = (xd * xd + yd * yd);
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }
            
            x[i2] = x[iMaxDist];
            y[i2] = y[iMaxDist];
            deg[i2] = deg[iMaxDist];
            i = nextI - 1;
            ++i2;            
        }
        
        return i2;
    }

    private static void rewriteBySortedIndexes(int[] sortedIndexes1, 
        long[] x) {
        
        int n = x.length;
        
        long[] x2 = new long[n];
        for (int i = 0; i < n; ++i) {
            x2[i] = x[sortedIndexes1[i]];
        }
        
        System.arraycopy(x2, 0, x, 0, n);
    }
    
    private static void rewriteBySortedIndexes(int[] sortedIndexes1, 
        int[] x) {
        
        int n = x.length;
        
        int[] x2 = new int[n];
        for (int i = 0; i < n; ++i) {
            x2[i] = x[sortedIndexes1[i]];
        }
        
        System.arraycopy(x2, 0, x, 0, n);
    }
}
