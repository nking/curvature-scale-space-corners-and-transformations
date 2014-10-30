package algorithms.imageProcessing;

import java.util.Arrays;

/**
 * class to hold x and y arrays of points
 * 
 * @author nichole
 */
public class PairIntArray {
    
    protected int[] x = null;
    
    protected int[] y = null;
    
    protected int n = 0;
    
    public PairIntArray(int capacity) {
        
        x = new int[capacity];
        
        y = new int[capacity];
    }
    
    public PairIntArray() {
        
        x = new int[10];
        
        y = new int[10];
    }
    
    public int getN() {
        return n;
    }
    
    public void add(int xPoint, int yPoint) {
        
        expandIfNeeded(n + 1);
        
        x[n] = xPoint;
        y[n] = yPoint;
        
        n++;
    }
    
    /**
     * remove indexes from idxLo to idxHi, inclusive
     * @param idxLo
     * @param idxHi 
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
        if (moveIdx < (n - 1)) {
            for (int moveToIdx = idxLo; moveToIdx < (n - nRemove); moveToIdx++) {
                x[moveToIdx] = x[moveIdx]; // 1<--3,2<--4,3<--5,4<--6
                y[moveToIdx] = y[moveIdx]; // 
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
    
    public void insertSpaceAtTopOfArrays(int numberOfInserts) {
        
        if (x.length >= (n + numberOfInserts)) {
            
            for (int i = (n - 1); i > -1; i--) {
                x[i + numberOfInserts] = x[i];
                y[i + numberOfInserts] = y[i];
            }
            for (int i = 0; i < numberOfInserts; i++) {
                x[i] = 0;
                y[i] = 0;
            }
            
        } else {
            int[] xx = new int[n + numberOfInserts];
            int[] yy = new int[n + numberOfInserts];
            System.arraycopy(x, 0, xx, numberOfInserts, n);
            System.arraycopy(y, 0, yy, numberOfInserts, n);
            x = xx;
            y = yy;
        }
        
        n += numberOfInserts;
    }
    
    public void set(int index, int xPoint, int yPoint) {
        
        if (index < 0) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(index + 1);
        
        x[index] = xPoint;
        y[index] = yPoint;
    }
    
    public void swapContents(PairIntArray other) {
        
        int[] swap = x;
        x = other.x;
        other.x = swap;
        
        swap = y;
        y = other.y;
        other.y = swap;
        
        int swap2 = n;
        n = other.n;
        other.n = swap2;
    }
    
    public int getX(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return x[index];
    }
    
    public int getY(int index) {
        if (index > (n - 1)) {
            throw new IllegalArgumentException("index is out of range");
        }
        return y[index];
    }
    
    public int[] getX() {
        return x;
    }
    
    public int[] getY() {
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
    
    public void reverse() {
        
        if (n < 2) {
            return;
        }
                
        int end = n >> 1;
        // 0 1 2 3 4
        for (int i = 0; i < end; i++) {
            int idx2 = n - i - 1;
            int swap = x[i];
            x[i] = x[idx2];
            x[idx2] = swap;
            
            swap = y[i];
            y[i] = y[idx2];
            y[idx2] = swap;
        }
    }
    
    public PairIntArray clone() {
        
        PairIntArray clone = new PairIntArray(n);
        
        clone.insertSpaceAtTopOfArrays(n);
        System.arraycopy(x, 0, clone.x, 0, n);
        System.arraycopy(y, 0, clone.y, 0, n);
        clone.n = n;
        
        return clone;
    }
    
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append("(").append(x[i]).append(", ").append(y[i]).append(") ");
        }
        return sb.toString();
    }
}
