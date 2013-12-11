package algorithms.misc;

public class Statistic {

    protected float avg = Float.MIN_VALUE;
    
    protected final int n;
    
    protected float stDev = -1;
    
    protected final int[] m;
    
    public Statistic(int[] items) {
        
        if (items == null) {
            throw new IllegalArgumentException("items cannot be null");
        }
        
        this.n = items.length;
        
        this.m = items;
        
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
