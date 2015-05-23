package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairIntArray;
import java.util.Arrays;

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
public class PolarAngleMergeSort {

    /**
     * value allowed in determining difference between 2 angles
     */
    protected static final double eps = 0.0000000001;//Double.MAX_VALUE;

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
        
        if (x.length == 1) {
            return 1;
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[x.length];

        sort(xP0, yP0, x, y, 1, x.length - 1, polarAngle);
        
        int nUsable = reduceToUniquePolarAngles(xP0, yP0, x, y, polarAngle);
        
        return nUsable;
    }
    
    public static int sort(int xP0, int yP0, PairIntArray xy) {

        if (xy == null) {
        	throw new IllegalArgumentException("xy cannot be null");
        }

        if (xy.getN() == 1) {
            return 1;
        }

        // for angles which are same, a delete operation is needed after all processing
        //    and ability to ignore the point to be deleted.
        double[] polarAngle = new double[xy.getN()];

        sort(xP0, yP0, xy, 1, xy.getN() - 1, polarAngle);
        
        int nUsable = reduceToUniquePolarAngles(xP0, yP0, xy, polarAngle);
        
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

    static int reduceToUniquePolarAngles(int xP0, int yP0, PairIntArray xy, 
        double[] polarAngle) {

        int lastKeptIndex = 0;

        for (int i = 1; i < xy.getN(); i++) {

            // store
            xy.set(lastKeptIndex + 1, xy.getX(i), xy.getY(i));

            // look ahead
            int nSkip = 0;
            int nextI = i + 1;
            double maxDistance = relativeLengthOfLine(xP0, yP0, 
                xy.getX(i), xy.getY(i));
            int indexMaxDistance = i;

            while ( (nextI < xy.getN()) && (Math.abs( polarAngle[i] - polarAngle[nextI] ) < eps) ) {
                double dist = relativeLengthOfLine(xP0, yP0, 
                    xy.getX(nextI), xy.getY(nextI));
                if (maxDistance < dist) {
                    maxDistance = dist;
                    indexMaxDistance = nextI;
                }
                nSkip++;
                nextI++;
            }
            
            if (nSkip > 0) {
                
                xy.set(lastKeptIndex + 1, xy.getX(indexMaxDistance), 
                    xy.getY(indexMaxDistance));
                
                i = nextI - 1;
            }

            lastKeptIndex++;
        }
        return lastKeptIndex + 1;
    }
    
    protected static void sort(float xP0, float yP0, float[] x, float[] y, int indexLo, int indexHi, double[] polarAngle) {

        int indexMid = -1;

        if (indexLo < indexHi) {

            indexMid = (indexLo + indexHi) >> 1;

            sort(xP0, yP0, x, y, indexLo, indexMid, polarAngle);
            sort(xP0, yP0, x, y, indexMid + 1, indexHi, polarAngle);
            merge(xP0, yP0, x, y, indexLo, indexMid, indexHi, polarAngle);
        }
    }
    
    static void sort(int xP0, int yP0, PairIntArray xy, 
        int indexLo, int indexHi, double[] polarAngle) {

        if (indexLo < indexHi) {

            int indexMid = (indexLo + indexHi) >> 1;

            sort(xP0, yP0, xy, indexLo, indexMid, polarAngle);
            sort(xP0, yP0, xy, indexMid + 1, indexHi, polarAngle);
            merge(xP0, yP0, xy, indexLo, indexMid, indexHi, polarAngle);
        }
    }

    private static void merge(float xP0, float yP0, float[] x, float[] y, int indexLo, int indexMid, int indexHi, double[] polarAngle) {

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

    private static void merge(int xP0, int yP0, PairIntArray xy, 
        int indexLo, int indexMid, int indexHi, double[] polarAngle) {

        int nLeft = indexMid - indexLo + 1;
        int nRight = indexHi - indexMid;

        int[] xLeft = Arrays.copyOfRange(xy.getX(), indexLo, indexMid + 2);       // add 1 for sentinel
        int[] yLeft = Arrays.copyOfRange(xy.getY(), indexLo, indexMid + 2);
        double[] angleLeft = new double[nLeft + 1];

        int[] xRight = Arrays.copyOfRange(xy.getX(), indexMid + 1, indexHi + 2);  // add 1 for sentinel
        int[] yRight = Arrays.copyOfRange(xy.getY(), indexMid + 1, indexHi + 2);
        double[] angleRight = new double[nRight + 1];

        int i, j, index;

        for (i = 0; i < nLeft; i++) {
            index = indexLo + i;
            angleLeft[i] = LinesAndAngles.calculatePolarSineTheta(xP0, yP0, 
                xy.getX(index), xy.getY(index));
        }

        for (j = 0; j < nRight; j++) {
            index = indexMid + j + 1;
            angleRight[j] = LinesAndAngles.calculatePolarSineTheta(xP0, yP0, 
                xy.getX(index), xy.getY(index));
        }

        xLeft[nLeft] = Integer.MAX_VALUE;
        yLeft[nLeft] = Integer.MAX_VALUE;
        angleLeft[nLeft] = Double.MAX_VALUE;
        xRight[nRight] = Integer.MAX_VALUE;
        yRight[nRight] = Integer.MAX_VALUE;
        angleRight[nRight] = Double.MAX_VALUE;

        i = 0;
        j = 0;

        for (int k = indexLo; k <= indexHi; k++) {

            double angDiff = angleLeft[i] - angleRight[j];
            if (angDiff < 0) {
                angDiff *= -1.;
            }
            if ( (angDiff < eps ) || (angleLeft[i] < angleRight[j] ) ) {

                xy.set(k, xLeft[i], yLeft[i]);
                polarAngle[k] = angleLeft[i];
                i++;
            } else {
                xy.set(k, xRight[j], yRight[j]);
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
