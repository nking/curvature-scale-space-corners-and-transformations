package algorithms.compGeometry.convexHull;

import algorithms.imageProcessing.util.AngleUtil;
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

        int lastKeptIndex = 0;
        
        double eps = 0;

        for (int i = 1; i < x.length; i++) {

            // store
            x[lastKeptIndex + 1] = x[i];
            y[lastKeptIndex + 1] = y[i];

            // look ahead
            int nSkip = 0;
            int nextI = i + 1;
            double maxDistance = relativeLengthOfLine(xP0, yP0, x[i], y[i]);
            int indexMaxDistance = i;

            while ( (nextI < x.length) && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                double dist = relativeLengthOfLine(xP0, yP0, x[nextI], y[nextI]);
                if (maxDistance < dist) {
                    maxDistance = dist;
                    indexMaxDistance = nextI;
                }
                nSkip++;
                nextI++;
            }
            if (nSkip > 0) {
                x[lastKeptIndex + 1] = x[indexMaxDistance];
                y[lastKeptIndex + 1] = y[indexMaxDistance];
                i = nextI - 1;
            }

            lastKeptIndex++;
        }
        return lastKeptIndex + 1;
    }
    
    static <T extends PairInt> int reduceToUniquePolarAngles(T p0, T[] points, 
        double[] polarAngle) {

        int lastKeptIndex = 0;
        
        double eps = 0;

        for (int i = 1; i < points.length; i++) {

            // store
            points[lastKeptIndex + 1] = points[i];

            // look ahead
            int nSkip = 0;
            int nextI = i + 1;
            double maxDistance = relativeLengthOfLine(p0, points[i]);
            int indexMaxDistance = i;

            while ( (nextI < points.length) && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                double dist = relativeLengthOfLine(p0, points[nextI]);
                if (maxDistance < dist) {
                    maxDistance = dist;
                    indexMaxDistance = nextI;
                }
                nSkip++;
                nextI++;
            }
            
            if (nSkip > 0) {
                points[lastKeptIndex + 1] = points[indexMaxDistance];
                i = nextI - 1;
            }

            lastKeptIndex++;
        }
        
        return lastKeptIndex + 1;
    }

    static <T extends PairInt> int reduceToUniquePolarAngles(
        T p0, List<T> points, double[] polarAngle) {

        int lastKeptIndex = 0;
        
        double eps = 0;

        for (int i = 1; i < points.size(); i++) {

            // store
            points.set(lastKeptIndex + 1, points.get(i));

            // look ahead
            int nSkip = 0;
            int nextI = i + 1;
            double maxDistance = relativeLengthOfLine(p0, points.get(i));
            int indexMaxDistance = i;

            while ( (nextI < points.size()) 
                && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= eps) ) {
                double dist = relativeLengthOfLine(p0, points.get(nextI));
                if (maxDistance < dist) {
                    maxDistance = dist;
                    indexMaxDistance = nextI;
                }
                nSkip++;
                nextI++;
            }
            
            if (nSkip > 0) {
                points.set(lastKeptIndex + 1, points.get(indexMaxDistance));
                i = nextI - 1;
            }

            lastKeptIndex++;
        }
        
        return lastKeptIndex + 1;
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
}
