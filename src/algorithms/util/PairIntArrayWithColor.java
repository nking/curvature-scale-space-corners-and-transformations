package algorithms.util;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class PairIntArrayWithColor extends PairIntArray {
    
    /**
     * a color for the array of points.  can be used to store a value that 
     * is interpreted as a closed curve, for example.
     */
    protected int color = 0;
    
    public PairIntArrayWithColor(PairIntArray obj) {
        this.x = Arrays.copyOf(obj.x, obj.x.length);
        this.y = Arrays.copyOf(obj.y, obj.y.length);
        this.n = obj.n;
    }
    
    public PairIntArrayWithColor(int nPoints) {
        super(nPoints);
    }
    
    public int getColor() {
        return color;
    }
    
    /**
     * set the color for this set of points.  the field can be used to hold
     * a value that's interpreted as the points being a closed curve, for
     * example.
     * 
     * @param clr 
     */
    public void setColor(int clr) {
        color = clr;
    }
     
}
