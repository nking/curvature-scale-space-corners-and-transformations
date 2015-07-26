package algorithms.compGeometry.convexHull;

import algorithms.imageProcessing.util.AngleUtil;
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
        
        for (int i = 1; i < x.length; i++) {
            polarAngle[i] = AngleUtil.polarAngleCCW((double)(x[i] - xP0), 
                (double)(y[i] - yP0));
        }
        
        sort(xP0, yP0, x, y, 1, x.length - 1, polarAngle);
        
        int nUsable = reduceToUniquePolarAngles(xP0, yP0, x, y, polarAngle);
        
        return nUsable;
    }
    
    static int reduceToUniquePolarAngles(float xP0, float yP0, float[] x, 
        float[] y, double[] polarAngle) {

        int lastKeptIndex = 0;
        
        //TODO: consider making this zero
        double eps = 0;//0.000000001;

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

    protected static void sort(float xP0, float yP0, float[] x, float[] y, 
        int indexLo, int indexHi, double[] polarAngle) {

        if (indexLo < indexHi) {

            int indexMid = (indexLo + indexHi) >> 1;

            sort(xP0, yP0, x, y, indexLo, indexMid, polarAngle);
            sort(xP0, yP0, x, y, indexMid + 1, indexHi, polarAngle);
            merge(xP0, yP0, x, y, indexLo, indexMid, indexHi, polarAngle);
        }
    }
    
    private static void merge(float xP0, float yP0, float[] x, float[] y, 
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
}
