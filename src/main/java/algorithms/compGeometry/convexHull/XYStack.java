package algorithms.compGeometry.convexHull;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class XYStack {

    /**
     *
     */
    public float[] x = null;

    /**
     *
     */
    public float[] y = null;

    /**
     *
     */
    protected int nXY = 0;

    /**
     *
     * @param capacity
     */
    public XYStack(int capacity) {
        x = new float[capacity];
        y = new float[capacity];
        nXY = 0;
    }

    /**
     *
     */
    public XYStack() {
        x = new float[10];
        y = new float[10];
        nXY = 0;
    }

    /**
     *
     * @param xPoint
     * @param yPoint
     */
    public void push(float xPoint, float yPoint) {
        if (nXY + 1 > (x.length - 1)) {
            int n = (int) (2.5*x.length);
            x = Arrays.copyOf(x, n);
            y = Arrays.copyOf(y, n);
        }

        x[nXY] = xPoint;
        y[nXY] = yPoint;
        nXY++;
    }

    /**
     *
     * @return
     */
    public float peekTopX() {
        if (nXY == 0) {
            return Float.NEGATIVE_INFINITY;
        }
        return x[nXY-1];
    }

    /**
     *
     * @return
     */
    public float peekTopY() {
        if (nXY == 0) {
            return Float.NEGATIVE_INFINITY;
        }
        return y[nXY-1];
    }

    /**
     *
     * @return
     */
    public float peekNextToTopX() {
        if (nXY < 2) {
            return Float.NEGATIVE_INFINITY;
        }
        return x[nXY-2];
    }

    /**
     *
     * @return
     */
    public float peekNextToTopY() {
        if (nXY < 2) {
            return Float.NEGATIVE_INFINITY;
        }
        return y[nXY-2];
    }

    /**
     *
     * @return
     */
    public float[] pop() {
        if (nXY == 0) {
            return null;
        }
        float[] a = new float[]{x[nXY - 1], y[nXY - 1]};
        nXY--;
        return a;
    }

    /**
     *
     * @return
     */
    public boolean isEmpty() {
        return (nXY == 0);
    }

    /**
     *
     */
    public void compressArrays() {
        if (nXY < x.length) {
            x = Arrays.copyOf(x, nXY);
            y = Arrays.copyOf(y, nXY);
        }
    }

    /**
     *
     * @return
     */
    public int getNPoints() {
        return nXY;
    }

}
