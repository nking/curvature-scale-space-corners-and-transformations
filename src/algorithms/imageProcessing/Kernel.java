package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class Kernel {
   
    final int[] a;
    
    private final int width;
    
    private final int height;
    
    private final int nPixels;
    /**
     * @param theWidth
     * @param theHeight
     */
    public Kernel (int theWidth, int theHeight) {
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
     
        a = new int[nPixels];
    }
    
    public void setValue(int col, int row, int value) {
        
        int idx = (row * width) + col;
       
        a[idx] = value;
    }
        
    public int getValue(int col, int row) {
        
        int idx = (row * width) + col;
       
        return a[idx];
    }

    /**
     * @return the width
     */
    public int getWidth() {
        return width;
    }

    /**
     * @return the height
     */
    public int getHeight() {
        return height;
    }

    /**
     * @return the nPixels
     */
    public int getNPixels() {
        return nPixels;
    }
}
