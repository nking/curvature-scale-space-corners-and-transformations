package algorithms.misc;

public class Statistic {

    protected float avg = Float.MIN_VALUE;
    
    protected final int n;
    
    protected float stDev = -1;
    
    protected final int[] m;
    
    /**
     * the starting x coordinate for items
     */
    protected final float[] itemsX;
    
    /**
     * the starting y coordinate for items
     */
    protected final float[] itemsY;
    
    /**
     * the size of an x cell.  for example, item[0] starts at (x[0], y[0]) and extends to (x[0] + xSz, y[0] + ySz)
     */
    protected final float xSz;
    
    /**
     * the size of y cell.
     */
    protected final float ySz;
    
    public Statistic(int[] items, float[] xCells, float[] yCells, float xDivSz, float yDivSz) {
        
        if (items == null) {
            throw new IllegalArgumentException("items cannot be null");
        }
        
        this.n = items.length;
        
        this.m = items;
        
        this.xSz = xDivSz;
        
        this.ySz = yDivSz;
        
        this.itemsX = xCells;
        
        this.itemsY = yCells;
        
        calculateStats();
    }
    
    public float getAverage() {
        return avg;
    }
    public float getStandardDeviation() {
        return stDev;
    }
    public int[] getItems() {
        return m;
    }
    public int getNumberOfItems() {
        return n;
    }
    
    public float[] getItemsX() {
        return itemsX;
    }
    public float[] getItemsY() {
        return itemsY;
    }
    public float getXSz() {
        return xSz;
    }
    public float getYSz() {
        return ySz;
    }
    
    protected void calculateStats() {
        if (avg > Float.MIN_VALUE) {
            return;
        }
        
        float s = 0.f;
        for (int item : m) {
            s += item;
        }
        avg = s/n;
        
        s = 0.f;
        for (int i = 0; i < m.length; i++) {
            s += Math.pow((m[i] - avg), 2);
        }
        
        stDev = (float) Math.sqrt(s/(n-1));;
    }
}
