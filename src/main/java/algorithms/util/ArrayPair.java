package algorithms.util;

/**
 *
 * @author nichole
 */
public class ArrayPair {

    /**
     *
     */
    protected final float[] x;
    
    /**
     *
     */
    protected final float[] y;
    
    /**
     *
     */
    protected final float[] xe;
    
    /**
     *
     */
    protected final float[] ye;
    
    /**
     *
     * @param a1
     * @param a2
     */
    public ArrayPair(float[] a1, float[] a2) {
        
        this.x = a1;
        
        this.y = a2;
        
        this.xe = null;
        
        this.ye = null;
    }
    
    /**
     *
     * @param xArray
     * @param yArray
     * @param xErrorArray
     * @param yErrorArray
     */
    public ArrayPair(float[] xArray, float[] yArray, float[] xErrorArray, 
        float[] yErrorArray) {
        
        this.x = xArray;
        
        this.y = yArray;
        
        this.xe = xErrorArray;
        
        this.ye = yErrorArray;
    }
    
    /**
     *
     * @return
     */
    public float[] getX() {
        return x;
    }
    
    /**
     *
     * @return
     */
    public float[] getY() {
        return y;
    }
    
    /**
     *
     * @return
     */
    public float[] getXErrors() {
        return xe;
    }
    
    /**
     *
     * @return
     */
    public float[] getYErrors() {
        return ye;
    }
}
