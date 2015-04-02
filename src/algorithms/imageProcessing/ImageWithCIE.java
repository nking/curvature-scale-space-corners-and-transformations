package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ImageWithCIE extends Image {
    
    protected float[] cieX;
    protected float[] cieY;
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public ImageWithCIE (int theWidth, int theHeight) {
        
        super(theWidth, theHeight);
    }
    
    protected void init() {
        
        this.cieX = new float[nPixels];
        this.cieY = new float[nPixels];
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int col = 0; col < width; col++) {
            for (int row = 0; row < height; row++) {
                
                int idx = (row * width) + col;
                
                int rPix = r[idx];
                int gPix = g[idx];
                int bPix = b[idx];
                
                float[] xy = cieC.rgbToXYChromaticity(rPix, gPix, bPix);
                cieX[idx] = xy[0];
                cieY[idx] = xy[1];
            }
        }
    }
  
    public Image copyImage() {
       
        ImageWithCIE img2 = new ImageWithCIE(getWidth(), getHeight());
        
        System.arraycopy(r, 0, img2.r, 0, nPixels);
        System.arraycopy(g, 0, img2.g, 0, nPixels);
        System.arraycopy(b, 0, img2.b, 0, nPixels);
        
        System.arraycopy(cieX, 0, img2.cieX, 0, nPixels);
        System.arraycopy(cieY, 0, img2.cieY, 0, nPixels);
       
        return img2;
    }
   
    public void resetTo(Image copyThis) {
        
        if (copyThis.getNPixels() != nPixels) {
            throw new IllegalArgumentException("cannot convert this fixed " 
                + "image size to the size of copyThis");
        }
        
        if (!(copyThis instanceof ImageWithCIE)) {
            throw new IllegalArgumentException(
            "copyThis has to be instance of ImageWithCIE");
        }
        
        System.arraycopy(copyThis.r, 0, r, 0, nPixels);
        System.arraycopy(copyThis.g, 0, g, 0, nPixels);
        System.arraycopy(copyThis.b, 0, b, 0, nPixels);
        
        System.arraycopy(((ImageWithCIE)copyThis).cieX, 0, cieX, 0, nPixels);
        System.arraycopy(((ImageWithCIE)copyThis).cieY, 0, cieY, 0, nPixels);
    }

}
