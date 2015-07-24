package algorithms.util;

import java.util.Arrays;

/**
 * class to hold x and y arrays of points
 * 
 * @author nichole
 */
public class PairFloatArray {
    
    protected float[] x = null;
    
    protected float[] y = null;
    
    protected int n = 0;
    
    public PairFloatArray(int capacity) {
        
        x = new float[capacity];
        
        y = new float[capacity];
    }
    
    public PairFloatArray() {
        
        x = new float[10];
        
        y = new float[10];
    }
    
    public int getN() {
        return n;
    }
    
    public void add(float xPoint, float yPoint) {
        
        expandIfNeeded(n + 1);
        
        x[n] = xPoint;
        y[n] = yPoint;
        
        n++;
    }
    
    public void set(int index, float xPoint, float yPoint) {
        
        if (index < 0) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(index + 1);
        
        x[index] = xPoint;
        y[index] = yPoint;
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
        return x;
    }
    
    public float[] getY() {
        return y;
    }
    
    protected void expandIfNeeded(int nTotal) {
        
        if (nTotal > x.length) {
            
            int n2 = x.length + 10;
            
            if (nTotal > n2) {
                n2 = nTotal;
            }
            
            x = Arrays.copyOf(x, n2);
            
            y = Arrays.copyOf(y, n2);
        }
    }
    
    public PairFloatArray copy() {
        
        PairFloatArray clone = new PairFloatArray(n);
        
        System.arraycopy(x, 0, clone.x, 0, n);
        System.arraycopy(y, 0, clone.y, 0, n);
        
        clone.n = n;
        
        return clone;
    }
    
    /**
     * remove indexes from idxLo to idxHi, inclusive
     * @param idxLo first index to be removed, inclusive
     * @param idxHi last index to be removed, inclusive
     */
    public void removeRange(int idxLo, int idxHi) {
        
        if ((idxLo < 0) || (idxLo > (n - 1))) {
            throw new IllegalArgumentException("idxLo is out of range");
        }
        if ((idxHi < 0) || (idxHi > (n - 1))) {
            throw new IllegalArgumentException("idxHi is out of range");
        }
        if (idxHi < idxLo) {
            throw new IllegalArgumentException("idxHi has to be >= idxLo");
        }
        
        int nRemove = idxHi - idxLo + 1;

        int moveIdx = idxHi + 1;
        if (moveIdx <= (n - 1)) {
            for (int moveToIdx = idxLo; moveToIdx < (n - nRemove); moveToIdx++) {
                x[moveToIdx] = x[moveIdx];
                y[moveToIdx] = y[moveIdx];
                moveIdx++;
            }
        }
        
        // not necessary, but easier debugging to reset the last nRemove to 0
        for (int i = (n - nRemove); i < n; i++) {
            x[i] = 0;
            y[i] = 0;
        }
        
        n -= nRemove;
    }
    
    public PairIntArray toPairIntArray() {
        PairIntArray out = new PairIntArray();
        for (int i = 0; i < n; i++) {
            out.add(Math.round(x[i]), Math.round(y[i]));
        }
        return out;
    }
    
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append("x=").append(x[i]).append("y=").append(y[i]).append("\n");
        }
        return sb.toString();
    }
}
