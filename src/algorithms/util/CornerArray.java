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
    
    protected final SIGMA sigma;
    
    protected PairIntArray xy = null;
    
    protected int[] idx = null;
        
    protected float[] x1stDeriv = null;
    
    protected float[] x2ndDeriv = null;
    
    protected float[] y1stDeriv = null;
    
    protected float[] y2ndDeriv = null;
    
    protected boolean isFromClosedCurve = false;
    
    //curvature:
    protected float[] k = null;
    
    protected int n = 0;
    
    public CornerArray(SIGMA theSigma, int capacity) {
        
        this.sigma = theSigma;
        
        initArrays(capacity);
    }
    
    public CornerArray(SIGMA theSigma) {
        
        this.sigma = theSigma;
        
        initArrays(10);
    }
    
    private void initArrays(int capacity) {
        
        xy = new PairIntArray(capacity);
                
        k = new float[capacity];
        
        idx = new int[capacity];
        
        x1stDeriv = new float[capacity];
    
        x2ndDeriv = new float[capacity];
    
        y1stDeriv = new float[capacity];
    
        y2ndDeriv = new float[capacity];
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
     * @param curvature
     * @param x1d
     * @param x2d
     * @param y1d
     * @param y2d
     * @param anInt
     */
    public void add(int xPoint, int yPoint, float curvature,
        float x1d, float x2d, float y1d, float y2d, int anInt) {
        
        expandIfNeeded(n + 1);
        
        assert(xy.getN() == n);
        
        xy.add(xPoint, yPoint);
        k[n] = curvature;
        x1stDeriv[n] = x1d;
        x2ndDeriv[n] = x2d;
        y1stDeriv[n] = y1d;
        y2ndDeriv[n] = y2d;
        idx[n] = anInt;
        
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
     * @param curvature
     * @param x1d
     * @param x2d
     * @param y1d
     * @param y2d
     * @param anInt
     */
    public void set(int index, int xPoint, int yPoint, float curvature,
        float x1d, float x2d, float y1d, float y2d, int anInt) {
        
        if (index < 0) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(index + 1);
        
        xy.set(index, xPoint, yPoint);
        
        k[index] = curvature;
        x1stDeriv[index] = x1d;
        x2ndDeriv[index] = x2d;
        y1stDeriv[index] = y1d;
        y2ndDeriv[index] = y2d;
        idx[index] = anInt;
    }
    
    /**
     * insert a row of data in the instance.  
     * Note that most uses of this class
     * expect that the instance data are always ordered by increasing idx 
     * (that is, anInt) values.
     * @param index
     * @param xPoint
     * @param yPoint
     * @param curvature
     * @param x1d
     * @param x2d
     * @param y1d
     * @param y2d
     * @param anInt
     */
    public void insert(int index, int xPoint, int yPoint, float curvature,
        float x1d, float x2d, float y1d, float y2d, int anInt) {
        
        if (index < 0 || (index > n)) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        assert(xy.getN() == n);
        
        expandIfNeeded(n + 1);
        
        xy.insert(index, xPoint, yPoint);
        
        // move everything at index thru n-1 to higher index
        for (int i = n; i > index; i--) {
            k[i] = k[i - 1];
            x1stDeriv[i] = x1stDeriv[i - 1];
            x2ndDeriv[i] = x2ndDeriv[i - 1];
            y1stDeriv[i] = y1stDeriv[i - 1];
            y2ndDeriv[i] = y2ndDeriv[i - 1];
            idx[i] = idx[i - 1];
        }
        
        k[index] = curvature;
        idx[index] = anInt;
        x1stDeriv[index] = x1d;
        x2ndDeriv[index] = x2d;
        y1stDeriv[index] = y1d;
        y2ndDeriv[index] = y2d;
        
        n++;
        
        assert(xy.getN() == n);
    }
  
    public int getX(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return xy.getX(index);
    }
    
    public float getXFirstDeriv(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return x1stDeriv[index];
    }
    
    public float getXSecondDeriv(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return x2ndDeriv[index];
    }
    
    public int getY(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return xy.getY(index);
    }
    
    public float getYFirstDeriv(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return y1stDeriv[index];
    }
    
    public float getYSecondDeriv(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return y2ndDeriv[index];
    }
    
    public float getCurvature(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return k[index];
    }
    
    public int getInt(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return idx[index];
    }
    
    public SIGMA getSIGMA() {
        return sigma;
    }
    
    public int[] getX() {
        return xy.getX();
    }
    
    public float[] getXFirstDeriv() {
        return x1stDeriv;
    }
    
    public float[] getXSecondDeriv() {
        return x2ndDeriv;
    }
    
    public int[] getY() {
        return xy.getY();
    }
    
    public float[] getYFirstDeriv() {
        return y1stDeriv;
    }
    
    public float[] getYSecondDeriv() {
        return y2ndDeriv;
    }
    
    public float[] getCurvature() {
        return k;
    }
    
    public int[] getYInt() {
        return idx;
    }
    
    protected void expandIfNeeded(int nTotal) {
        
        if (nTotal > k.length) {
            
            int n2 = k.length + 10;
            
            if (nTotal > n2) {
                n2 = nTotal;
            }
            
            k = Arrays.copyOf(k, n2);
            
            x1stDeriv = Arrays.copyOf(x1stDeriv, n2);
            
            x2ndDeriv = Arrays.copyOf(x2ndDeriv, n2);
            
            y1stDeriv = Arrays.copyOf(y1stDeriv, n2);
            
            y2ndDeriv = Arrays.copyOf(y2ndDeriv, n2);
            
            idx = Arrays.copyOf(idx, n2);            
        }
    }
    
    public CornerArray copy() {
        
        CornerArray clone = new CornerArray(sigma, n);
        
        clone.xy = xy.copy();
        
        System.arraycopy(k, 0, clone.k, 0, n);
        System.arraycopy(idx, 0, clone.idx, 0, n);
        System.arraycopy(x1stDeriv, 0, clone.x1stDeriv, 0, n);
        System.arraycopy(x2ndDeriv, 0, clone.x2ndDeriv, 0, n);
        System.arraycopy(y1stDeriv, 0, clone.y1stDeriv, 0, n);
        System.arraycopy(y2ndDeriv, 0, clone.y2ndDeriv, 0, n);
        
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
        
        xy.removeRange(idxLo, idxHi);
        
        int nRemove = idxHi - idxLo + 1;

        int moveIdx = idxHi + 1;
        if (moveIdx <= (n - 1)) {
            for (int moveToIdx = idxLo; moveToIdx < (n - nRemove); moveToIdx++) {
                k[moveToIdx] = k[moveIdx];
                idx[moveToIdx] = idx[moveIdx];
                x1stDeriv[moveToIdx] = x1stDeriv[moveIdx];
                x2ndDeriv[moveToIdx] = x2ndDeriv[moveIdx];
                y1stDeriv[moveToIdx] = y1stDeriv[moveIdx];
                y2ndDeriv[moveToIdx] = y2ndDeriv[moveIdx];
                moveIdx++;
            }
        }
        
        // not necessary, but easier debugging to reset the last nRemove to 0
        for (int i = (n - nRemove); i < n; i++) {
            k[i] = 0;
            idx[i] = 0;
            x1stDeriv[i] = 0;
            x2ndDeriv[i] = 0;
            y1stDeriv[i] = 0;
            y2ndDeriv[i] = 0;
        }
        
        n -= nRemove;
    }

    public void setIsClosedCurve() {
        this.isFromClosedCurve = true;
    }
    
    public boolean isFromAClosedCurve() {
        return isFromClosedCurve;
    }
    
    public PairIntArray getXYCurve() {
        return xy;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("sigma=").append(sigma);
        for (int i = 0; i < n; i++) {
            sb.append(" ").append(xy.toString())
                .append(",curvature=").append(k[i])
                .append(",x1stDeriv=").append(x1stDeriv[i])
                .append(",x2ndDeriv=").append(x2ndDeriv[i])
                .append(",y1stDeriv=").append(y1stDeriv[i])
                .append(",y2ndDeriv=").append(y2ndDeriv[i])
                .append(",idx=").append(idx[i])
                .append("\n");
        }
        return sb.toString();
    }

}
