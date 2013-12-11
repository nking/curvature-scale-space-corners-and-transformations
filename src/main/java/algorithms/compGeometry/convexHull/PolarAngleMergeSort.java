package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.LinesAndAngles;
import java.util.Arrays;

/**
 * The algorithm uses a merge sort which has a worse case runtime is O(N * log_2(N)),
 * but also includes an additional set of operations to remove inner points with the
 * same polar angle from P0.

 *
 * @author Nichole King
 */
public class PolarAngleMergeSort {

    /**
     * value allowed in determining difference between 2 angles
     */
    protected static final double eps = 0.0001;

	/**
     *
     * The sort is for the polar angle between (xP0, yP0) and each other point
     * in x and y, excluding the very first point in x and y which is assumed to be
     * xP0 and yP0.
     *
     * Another pass is made through the data afterwards to remove points with the
     * same angle (keeping the furthest point only).
     *
     * The return value is the number of usable items in the sorted arrays
     *    x and y, specifically, x[0] through x[nUsable-1] and y[0] through y[nUsable-1].
     *
     * The caller of this method must handle not using points beyond the return value index.
     *
     * @param xP0 x-point from which to compute the polar angle to each x,y point
     * @param yP0 y-point from which to compute the polar angle to each x,y point
     * @param x x-array to sort (the first point is ignored!)
     * @param y y-array to sort (the first point is ignored!)
     * @return number of usable items in the array after sorting (the points that should
     * be deleted are at the end of the array, so should be ignored).
     */
	public static int sort(double xP0, double yP0, double[] x, double[] y) {

        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
        	throw new IllegalArgumentException("number of items in x must be the same as in y");
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[x.length];

        sort(xP0, yP0, x, y, 1, x.length - 1, polarAngle);

        int nUsable = reduceToUniquePolarAngles(xP0, yP0, x, y, polarAngle);

        return nUsable;
    }

    public static int sort(float xP0, float yP0, float[] x, float[] y) {

        if (x == null) {
        	throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
        	throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
        	throw new IllegalArgumentException("number of items in x must be the same as in y");
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[x.length];

        sort(xP0, yP0, x, y, 1, x.length - 1, polarAngle);
        
System.out.println("in PA sort x=" + Arrays.toString(x));
System.out.println("in PA sort y=" + Arrays.toString(y));

        int nUsable = reduceToUniquePolarAngles(xP0, yP0, x, y, polarAngle);
        
System.out.println("in PA after reduce x=" + Arrays.toString(x));
System.out.println("in PA after reduce y=" + Arrays.toString(y));
        
        return nUsable;
    }

    static int reduceToUniquePolarAngles(double xP0, double yP0, double[] x, double[] y, double[] polarAngle) {

        int lastKeptIndex = 0;

        for (int i = 1; i < x.length; i++) {

            // store
            x[lastKeptIndex + 1] = x[i];
            y[lastKeptIndex + 1] = y[i];

            // look ahead
            int nSkip = 0;
            int nextI = i + 1;
            double maxDistance = relativeLengthOfLine(xP0, yP0, x[i], y[i]);
            int indexMaxDistance = i;

            while ( (nextI < x.length) && (Math.abs( polarAngle[i] - polarAngle[nextI] ) < eps) ) {
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

    static int reduceToUniquePolarAngles(float xP0, float yP0, float[] x, float[] y, double[] polarAngle) {

        int lastKeptIndex = 0;

        for (int i = 1; i < x.length; i++) {

            // store
            x[lastKeptIndex + 1] = x[i];
            y[lastKeptIndex + 1] = y[i];

            // look ahead
            int nSkip = 0;
            int nextI = i + 1;
            double maxDistance = relativeLengthOfLine(xP0, yP0, x[i], y[i]);
            int indexMaxDistance = i;

            while ( (nextI < x.length) && (Math.abs( polarAngle[i] - polarAngle[nextI] ) < eps) ) {
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

    static void sort(double xP0, double yP0, double[] x, double[] y, int indexLo, int indexHi, double[] polarAngle) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi)/2;

            sort(xP0, yP0, x, y, indexLo, indexMid, polarAngle);
            sort(xP0, yP0, x, y, indexMid + 1, indexHi, polarAngle);
            merge(xP0, yP0, x, y, indexLo, indexMid, indexHi, polarAngle);
        }
    }

    private static void sort(float xP0, float yP0, float[] x, float[] y, int indexLo, int indexHi, double[] polarAngle) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi)/2;

            sort(xP0, yP0, x, y, indexLo, indexMid, polarAngle);
            sort(xP0, yP0, x, y, indexMid + 1, indexHi, polarAngle);
            merge(xP0, yP0, x, y, indexLo, indexMid, indexHi, polarAngle);
        }
    }


    /**
     * @param p0 point in which to calculate polar angle between points
     * @param x array
     * @param y array
     * @param indexLo first index of first subarray
     * @param indexMid first index of second subarray
     * @param indexHi last index of second subarray
     * @param remove a boolean array indicating that a point should be ignored and later removed
     */
    private static void merge( double xP0, double yP0, double[] x, double[] y, int indexLo, int indexMid, int indexHi, double[] polarAngle) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        double[] xLeft = Arrays.copyOfRange(x, indexLo, indexMid + 2);       // add 1 for sentinel
        double[] yLeft = Arrays.copyOfRange(y, indexLo, indexMid + 2);
        double[] angleLeft = new double[nLeft + 1];

        double[] xRight = Arrays.copyOfRange(x, indexMid + 1, indexHi + 2);  // add 1 for sentinel
        double[] yRight = Arrays.copyOfRange(y, indexMid + 1, indexHi + 2);
        double[] angleRight = new double[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            angleLeft[i] = LinesAndAngles.calculatePolarSineTheta(xP0, yP0, x[index], y[index]);
        }

        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            angleRight[j] = LinesAndAngles.calculatePolarSineTheta(xP0, yP0, x[index], y[index]);
        }

        xLeft[nLeft] = Double.MAX_VALUE;
        yLeft[nLeft] = Double.MAX_VALUE;
        angleLeft[nLeft] = Double.MAX_VALUE;
        xRight[nRight] = Double.MAX_VALUE;
        yRight[nRight] = Double.MAX_VALUE;
        angleRight[nRight] = Double.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = indexLo; k <= indexHi; k++) {
            
            double angDiff = angleLeft[i] - angleRight[j];
            
            if (angDiff < 0) {
                angDiff *= -1.;
            }
            if ( (angDiff < eps ) || (angleLeft[i] < angleRight[j] ) ) {

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

    private static void merge( float xP0, float yP0, float[] x, float[] y, int indexLo, int indexMid, int indexHi, double[] polarAngle) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        float[] xLeft = Arrays.copyOfRange(x, indexLo, indexMid + 2);       // add 1 for sentinel
        float[] yLeft = Arrays.copyOfRange(y, indexLo, indexMid + 2);
        double[] angleLeft = new double[nLeft + 1];

        float[] xRight = Arrays.copyOfRange(x, indexMid + 1, indexHi + 2);  // add 1 for sentinel
        float[] yRight = Arrays.copyOfRange(y, indexMid + 1, indexHi + 2);
        double[] angleRight = new double[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            angleLeft[i] = LinesAndAngles.calculatePolarSineTheta(xP0, yP0, x[index], y[index]);
        }

        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            angleRight[j] = LinesAndAngles.calculatePolarSineTheta(xP0, yP0, x[index], y[index]);
        }

        xLeft[nLeft] = Float.MAX_VALUE;
        yLeft[nLeft] = Float.MAX_VALUE;
        angleLeft[nLeft] = Double.MAX_VALUE;
        xRight[nRight] = Float.MAX_VALUE;
        yRight[nRight] = Float.MAX_VALUE;
        angleRight[nRight] = Double.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = indexLo; k <= indexHi; k++) {

            double angDiff = angleLeft[i] - angleRight[j];
            if (angDiff < 0) {
                angDiff *= -1.;
            }
            if ( (angDiff < eps ) || (angleLeft[i] < angleRight[j] ) ) {

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
}
