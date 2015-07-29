package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class PixelColors {
    
    private final int red;
    private final int green;
    private final int blue;
    private final float cieX;
    private final float cieY;
    
    public PixelColors(int r, int g, int b, float cieX, float cieY) {
        this.red = r;
        this.green = g;
        this.blue = b;
        this.cieX = cieX;
        this.cieY = cieY;
    }
   
    public int getRed() {
        return red;
    }
    
    public int getGreen() {
        return green;
    }
    
    public int getBlue() {
        return blue;
    }
    
    public float getCIEX() {
        return cieX;
    }
    
    public float getCIEY() {
        return cieY;
    }
}
