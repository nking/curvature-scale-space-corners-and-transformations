package algorithms.compGeometry.convexHull;

import algorithms.util.PairFloatArray;
import java.util.Arrays;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 * 
 */
public class XYStack extends PairFloatArray {

    public XYStack(int capacity) {
        super(capacity);
    }
    
    public XYStack() {
        super();
    }

    /**
     * push a value onto the stack.
     * @param xPoint
     * @param yPoint 
     */
    public void push(float xPoint, float yPoint) {
        /*
        value is physically at last used position in array
        */
        
        add(xPoint, yPoint);
    }

    public float peekTopX() {
        if (this.n == 0) {
            return Float.NEGATIVE_INFINITY;
        }
        return x[n - 1];
    }

    public float peekTopY() {
        if (this.n == 0) {
            return Float.NEGATIVE_INFINITY;
        }
        return y[n - 1];
    }
    public float peekNextToTopX() {
        if (this.n < 2) {
            return Float.NEGATIVE_INFINITY;
        }
        return x[n - 2];
    }

    public float peekNextToTopY() {
        if (n < 2) {
            return Float.NEGATIVE_INFINITY;
        }
        return y[n - 2];
    }

    public float[] pop() {
        if (n == 0) {
            return null;
        }
        
        float[] a = new float[]{x[n - 1], y[n - 1]};
        
        removeRange(n - 1, n - 1);
        
        return a;
    }
    
    public void popWithoutReturn() {
        if (n == 0) {
            return;
        }
                
        removeRange(n - 1, n - 1);        
    }

    public boolean isEmpty() {
        return (n == 0);
    }

    public void compressArrays() {
        if (n < x.length) {
            x = Arrays.copyOf(x, n);
            y = Arrays.copyOf(y, n);
        }
    }

}
