package algorithms.util;

import java.util.Arrays;

/**
 * class to hold x and y arrays of points
 * 
 * @author nichole
 */
public class PairFloatArrayUnmodifiable {
    
    private final float[] x;
    
    private final float[] y;
    
    protected int n = 0;
    
    public PairFloatArrayUnmodifiable(float[] xPoints, float[] yPoints) {
        
        x = Arrays.copyOf(xPoints, xPoints.length);
        
        y = Arrays.copyOf(yPoints, yPoints.length);
        
        n = x.length;
    }
    
    public PairFloatArrayUnmodifiable(float[] xPoints, float[] yPoints, int length) {
        
        if (xPoints == null) {
            throw new IllegalArgumentException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalArgumentException("yPoints cannot be null");
        }
        if (xPoints.length != yPoints.length) {
            throw new IllegalArgumentException(
            "xPoints.length must equal to yPoints.length");
        }
        if (length > xPoints.length) {
            throw new IllegalArgumentException(
            "length must be less than or equal to xPoints.length");
        }
       
        x = Arrays.copyOfRange(xPoints, 0, length);
        
        y = Arrays.copyOfRange(yPoints, 0, length);
        
        n = x.length;
    }
    
    public int getN() {
        return n;
    }
  
    public float getX(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return x[index];
    }
    
    public float getY(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return y[index];
    }
    
    public float[] getX() {
        return Arrays.copyOf(x, x.length);
    }
    
    public float[] getY() {
        return Arrays.copyOf(y, y.length);
    }
    
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append("x=").append(x[i]).append("y=").append(y[i]).append("\n");
        }
        return sb.toString();
    }
}
