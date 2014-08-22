package algorithms.util;

public class ArrayPair {

    protected final float[] x;
    
    protected final float[] y;
    
    protected final float[] xe;
    
    protected final float[] ye;
    
    public ArrayPair(float[] a1, float[] a2) {
        
        this.x = a1;
        
        this.y = a2;
        
        this.xe = null;
        
        this.ye = null;
    }
    
    public ArrayPair(float[] xArray, float[] yArray, float[] xErrorArray, 
        float[] yErrorArray) {
        
        this.x = xArray;
        
        this.y = yArray;
        
        this.xe = xErrorArray;
        
        this.ye = yErrorArray;
    }
    
    public float[] getX() {
        return x;
    }
    
    public float[] getY() {
        return y;
    }
    
    public float[] getXErrors() {
        return xe;
    }
    
    public float[] getYErrors() {
        return ye;
    }
}
