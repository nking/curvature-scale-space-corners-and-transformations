package algorithms.imageProcessing;

import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;

/**
 *
 * @author nichole
 */
public class Image {
    
    //TODO:  add alpha when needed
    final int[] r;
    final int[] g;
    final int[] b;
    
    protected final int width;
    
    protected final int height;
    
    protected final int nPixels;
    
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
        
        init();
    }
    
    protected void init() {
        // not used for base class
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
    
    public GreyscaleImage copyToGreyscale() {
        
        /*
        to maintain same image conversion as used in ImageIOHelper,
        will use the java methods.  
        TODO: There should be a way to perform
        a pixel by pixel conversion instead of creating a
        BufferedImage.
        */
        
        BufferedImage outputImage = new BufferedImage(width, height, 
            BufferedImage.TYPE_BYTE_GRAY);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                int rgbValue = getRGB(i, j);
                outputImage.setRGB(i, j, rgbValue);
            }
        }
        
        GreyscaleImage out = new GreyscaleImage(width, height);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                // presumably, this is already combined?
                // or does it need separation into rgb and then averaged?
                
                int rgb = outputImage.getRGB(i, j);
                
                // prefer GREEN?
                
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;  
                    
                int v = Math.round((r + g + b)/3.f);
                
                out.setValue(i, j, v);
            }
        }
        
        return out;
    }
    
    public GreyscaleImage copyRedToImage() {
        return copyToImage(0);
    }
    
    public GreyscaleImage copyGreenToImage() {
        return copyToImage(1);
    }
    
    public GreyscaleImage copyBlueToImage() {
        return copyToImage(2);
    }
    
    protected GreyscaleImage copyToImage(int redOrGreenOrBlue) {
        
        /*
        to maintain same image conversion as used in ImageIOHelper,
        will use the java methods.  
        TODO: There should be a way to perform
        a pixel by pixel conversion instead of creating a
        BufferedImage.
        */
        
        BufferedImage outputImage = new BufferedImage(width, height, 
            BufferedImage.TYPE_BYTE_GRAY);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                int rgbValue = getRGB(i, j);
                outputImage.setRGB(i, j, rgbValue);
            }
        }
        
        GreyscaleImage out = new GreyscaleImage(width, height);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                // presumably, this is already combined?
                // or does it need separation into rgb and then averaged?
                
                int rgb = outputImage.getRGB(i, j);
                                
                int v;
                if (redOrGreenOrBlue == 0) {
                    v = (rgb >> 16) & 0xFF;
                } else if (redOrGreenOrBlue == 1) {
                    v = (rgb >> 8) & 0xFF;
                } else {
                    v = rgb & 0xFF;
                }
                                    
                out.setValue(i, j, v);
            }
        }
        
        return out;
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
