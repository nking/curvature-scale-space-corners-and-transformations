package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class Image {
    
    //TODO:  add alpha when needed
    final int[] r;
    final int[] g;
    final int[] b;
    
    private final int width;
    
    private final int height;
    
    private final int nPixels;
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public Image (int theWidth, int theHeight) {
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        r = new int[nPixels];
        
        g = new int[nPixels];
        
        b = new int[nPixels];
    }
    
    public void setRGB(int col, int row, int rPix, int gPix, int bPix) {
        
        /*if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException(
                "col is out of bounds");
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException(
                "row is out of bounds");
        }*/
        
        int idx = (row * getWidth()) + col;
       
        /*if ((idx < 0) || (idx > (r.length - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }*/
        
        r[idx] = rPix;
        g[idx] = gPix;
        b[idx] = bPix;
    }
    
    public void setRGB(int col, int row, int rgb) {
        
        int idx = (row * getWidth()) + col;
        
        if ((idx < 0) || (idx > (r.length - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
       
        int rPix = (rgb >> 16) & 0xFF;
        int gPix = (rgb >> 8) & 0xFF;
        int bPix = rgb & 0xFF;
        
        r[idx] = rPix;
        g[idx] = gPix;
        b[idx] = bPix;
    }
        
    public int getR(int col, int row) {
        
        int idx = (row * getWidth()) + col;
       
        return r[idx];
    }
    
    public int getB(int col, int row) {
        
        int idx = (row * getWidth()) + col;
       
        return b[idx];
    }
    
    public int getG(int col, int row) {
        
        int idx = (row * getWidth()) + col;
       
        return g[idx];
    }
    
    public int getRGB(int col, int row) {
    
        int idx = (row * getWidth()) + col;
        
        int rgb = (((r[idx] & 0x0ff) << 16) 
            | ((g[idx] & 0x0ff) << 8) | (b[idx] & 0x0ff));
        
        return rgb;
    }
    
    public Image copyImage() {
       
        Image img2 = new Image(getWidth(), getHeight());
        
        System.arraycopy(r, 0, img2.r, 0, nPixels);
        System.arraycopy(g, 0, img2.g, 0, nPixels);
        System.arraycopy(b, 0, img2.b, 0, nPixels);
       
        return img2;
    }
    
    public void resetTo(Image copyThis) {
        
        if (copyThis.getNPixels() != nPixels) {
            throw new IllegalArgumentException("cannot convert this fixed " 
                + "image size to the size of copyThis");
        }
        
        System.arraycopy(copyThis.r, 0, r, 0, nPixels);
        System.arraycopy(copyThis.g, 0, g, 0, nPixels);
        System.arraycopy(copyThis.b, 0, b, 0, nPixels);
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
