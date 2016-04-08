package algorithms.util;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class PairIntArrayWithColor extends PairIntArray {
    
    /**
     * a color for the array of points.  can be used to store a value that 
     * holds special meaning.
     * have added logic for "closed curve" and "ordered points"
     * 
     * changing to use set bit and test bit
     * <pre>
     * set color:
     *    color |= (1 << 2)
     * test color:
     *     (color & (1 << x)) != 0
     * </pre>
     */
    private int color = 0;
    
    public PairIntArrayWithColor(PairIntArray obj) {
        this.x = Arrays.copyOf(obj.x, obj.x.length);
        this.y = Arrays.copyOf(obj.y, obj.y.length);
        this.n = obj.n;
    }
    
    public PairIntArrayWithColor(int nPoints) {
        super(nPoints);
    }
    
    public boolean isClosedCurve() {
        return (color & (1 << 1)) != 0;
    }
    
    public void setAsClosedCurve() {
        color |= (1 << 1);
    }
    
    public boolean isOrderedCurve() {
        return (color & (1 << 2)) != 0;
    }
    
    public void setAsOrderedCurve() {
        color |= (1 << 2);
    }
}
