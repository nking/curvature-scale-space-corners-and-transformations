package algorithms.util;

public class ArrayPair {

    protected final float[] x;
    
    protected final float[] y;
    
    public ArrayPair(float[] a1, float[] a2) {
        
        this.x = a1;
        
        this.y = a2;
            
    }
    
    public float[] getX() {
        return x;
    }
    
    public float[] getY() {
        return y;
    }
}
