package algorithms.util;

import algorithms.imageProcessing.SIGMA;
import java.util.Arrays;

/**
 * class to hold x and y arrays of points along with the edge indexes
 * and the SIGMA of the scale space curve that the points were last
 * updated from.  Note that most uses of the class expect that
 * the points are ordered by increasing idx values.
 * 
 * @author nichole
 */
public class CornerArray {
    
    protected float[] x = null;
    
    protected float[] y = null;
    
    protected int[] idx = null;
    
    protected SIGMA[] sigma = null;
    
    protected int n = 0;
    
    public CornerArray(int capacity) {
        
        x = new float[capacity];
        
        y = new float[capacity];
        
        idx = new int[capacity];
        
        sigma = new SIGMA[capacity];
    }
    
    public CornerArray() {
        
        x = new float[10];
        
        y = new float[10];
        
        idx = new int[10];
        
        sigma = new SIGMA[10];
    }
    
    public int getN() {
        return n;
    }
    
    /**
     * add a row of data to the instance.  Note that most uses of this class
     * expect that the instance data are always ordered by increasing idx 
     * (that is, anInt) values.
     * @param xPoint
     * @param yPoint
     * @param anInt
     * @param s 
     */
    public void add(float xPoint, float yPoint, int anInt, SIGMA s) {
        
        expandIfNeeded(n + 1);
        
        x[n] = xPoint;
        y[n] = yPoint;
        idx[n] = anInt;
        sigma[n] = s;
        
        n++;
    }
    
    /**
     * set a row of data in the instance (replacing existing at given index).  
     * Note that most uses of this class
     * expect that the instance data are always ordered by increasing idx 
     * (that is, anInt) values.
     * @param index
     * @param xPoint
     * @param yPoint
     * @param anInt
     * @param s 
     */
    public void set(int index, float xPoint, float yPoint, int anInt, SIGMA s) {
        
        if (index < 0) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(index + 1);
        
        x[index] = xPoint;
        y[index] = yPoint;
        idx[index] = anInt;
        sigma[index] = s;
    }
    
    /**
     * insert a row of data in the instance.  
     * Note that most uses of this class
     * expect that the instance data are always ordered by increasing idx 
     * (that is, anInt) values.
     * @param index
     * @param xPoint
     * @param yPoint
     * @param anInt
     * @param s 
     */
    public void insert(int index, float xPoint, float yPoint, int anInt, SIGMA s) {
        
        if (index < 0 || (index > n)) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(n + 1);
        
        // move everything at index thru n-1 to higher index
        for (int i = n; i > index; i--) {
            x[i] = x[i - 1];
            y[i] = y[i - 1];
            idx[i] = idx[i - 1];
            sigma[i] = sigma[i - 1];
        }
        
        x[index] = xPoint;
        y[index] = yPoint;
        idx[index] = anInt;
        sigma[index] = s;
        
        n++;
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
    
    public int getInt(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return idx[index];
    }
    
    public SIGMA getSIGMA(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return sigma[index];
    }
    
    public float[] getX() {
        return x;
    }
    
    public float[] getY() {
        return y;
    }
    
    public int[] getYInt() {
        return idx;
    }
    
    protected void expandIfNeeded(int nTotal) {
        
        if (nTotal > x.length) {
            
            int n2 = x.length + 10;
            
            if (nTotal > n2) {
                n2 = nTotal;
            }
            
            x = Arrays.copyOf(x, n2);
            
            y = Arrays.copyOf(y, n2);
            
            idx = Arrays.copyOf(idx, n2);
            
            sigma = Arrays.copyOf(sigma, n2);
        }
    }
    
    public CornerArray copy() {
        
        CornerArray clone = new CornerArray(n);
        
        System.arraycopy(x, 0, clone.x, 0, n);
        System.arraycopy(y, 0, clone.y, 0, n);
        System.arraycopy(idx, 0, clone.idx, 0, n);
        System.arraycopy(sigma, 0, clone.sigma, 0, n);
        
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
                idx[moveToIdx] = idx[moveIdx];
                sigma[moveToIdx] = sigma[moveIdx];
                moveIdx++;
            }
        }
        
        // not necessary, but easier debugging to reset the last nRemove to 0
        for (int i = (n - nRemove); i < n; i++) {
            x[i] = 0;
            y[i] = 0;
            idx[i] = 0;
            sigma[i] = null;
        }
        
        n -= nRemove;
    }
    
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append("x=").append(x[i]).append(",y=").append(y[i])
                .append(",idx=").append(idx[i]).append(",sigma=").append(sigma[i])
                .append("\n");
        }
        return sb.toString();
    }
}
