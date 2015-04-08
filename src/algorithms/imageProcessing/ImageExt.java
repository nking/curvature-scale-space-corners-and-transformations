package algorithms.imageProcessing;

import java.awt.Color;

/**
 * image specialized to hold CIE 1931 xy chromaticity coordinates (colors) 
 * and HSB colors populated on demand.
 * 
 * @author nichole
 */
public class ImageExt extends Image {
    
    // TODO: consider more compact forms of these numbers.
    // scaled to integers that are small in range?
    // or float converted to bit format IEEE 754 could be a small enough
    // range to use shift and add storage in an integer?  
    // any change would need to be worth the trade-off in space for increasing
    // number of steps to get the original number.
    
    protected float[] cieX;
    protected float[] cieY;
    protected float[] hue;
    protected float[] saturation;
    protected float[] brightness;
    
    //TODO: this could be a bit vector instead
    private boolean[] extPopulated;
    
    private int radiusForPopulateOnDemand = 1;
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public ImageExt (int theWidth, int theHeight) {
        
        super(theWidth, theHeight);
    }
    
    protected void init() {
        
        this.cieX = new float[nPixels];
        this.cieY = new float[nPixels];
        this.hue = new float[nPixels];
        this.saturation = new float[nPixels];
        this.brightness = new float[nPixels];
        this.extPopulated = new boolean[nPixels];
        
    }
    
    /**
     * set the block radius of other pixel color calculations when a pixel's 
     * colors are calculated (populated indirectly from a get method when the 
     * color has not yet been calculated).  In other words, if
     * radius is 0, only the specified pixel's color will be calculated, else if
     * radius is 1 for example, the pixel colors of the 8 surrounding neighbors 
     * will also be calculated.  By default, the instance has a radius of 1.
     * 
     * @param radius 
     */
    public void setRadiusForPopulateOnDemand(int radius) {
        this.radiusForPopulateOnDemand = radius;
    }
    
    public float getCIEX(int col, int row) {
        
        int idx = getInternalIndex(col, row);
        
        return getCIEX(idx);
    }
    
    public float getCIEY(int col, int row) {
        
        int idx = getInternalIndex(col, row);
        
        return getCIEY(idx);
    }
  
    public float getHue(int col, int row) {
        
        int idx = getInternalIndex(col, row);
        
        return getHue(idx);
    }
    
    public float getSaturation(int col, int row) {
        
        int idx = getInternalIndex(col, row);
        
        return getSaturation(idx);
    }
    
    public float getBrightness(int col, int row) {
        
        int idx = getInternalIndex(col, row);
       
        return getBrightness(idx);
    }
    
    public float getCIEX(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (!extPopulated[internalIndex]) {
            calculateColorIncludingNeighbors(internalIndex, 
                radiusForPopulateOnDemand);
        }
       
        return cieX[internalIndex];
    }
    
    public float getCIEY(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (!extPopulated[internalIndex]) {
            calculateColorIncludingNeighbors(internalIndex, 
                radiusForPopulateOnDemand);
        }
       
        return cieY[internalIndex];
    }
    
    public float getHue(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (!extPopulated[internalIndex]) {
            calculateColorIncludingNeighbors(internalIndex, 
                radiusForPopulateOnDemand);
        }
       
        return hue[internalIndex];
    }
    
    public float getSaturation(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (!extPopulated[internalIndex]) {
            calculateColorIncludingNeighbors(internalIndex, 
                radiusForPopulateOnDemand);
        }
       
        return saturation[internalIndex];
    }
    
    public float getBrightness(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (!extPopulated[internalIndex]) {
            calculateColorIncludingNeighbors(internalIndex, 
                radiusForPopulateOnDemand);
        }
       
        return brightness[internalIndex];
    }
  
    protected void calculateColor(int idx) {
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        int rPix = r[idx];
        int gPix = g[idx];
        int bPix = b[idx];

        float[] xy = cieC.rgbToXYChromaticity(rPix, gPix, bPix);
        
        cieX[idx] = xy[0];
        cieY[idx] = xy[1];
        
        float[] hsb = new float[3];
        Color.RGBtoHSB(rPix, gPix, bPix, hsb);
       
        hue[idx] = hsb[0];
        saturation[idx] = hsb[1];
        brightness[idx] = hsb[2];
        
        extPopulated[idx] = true;
    }
    
    protected void calculateColorIncludingNeighbors(int idx, int neighborRadius) {
        
        int col0 = this.getCol(idx);
        int row0 = this.getRow(idx);
        
        for (int col = (col0 - neighborRadius); col <= (col0 + neighborRadius); 
            col++) {
            
            if ((col < 0) || (col > (this.width - 1))) {
                continue;
            }
            
            for (int row = (row0 - neighborRadius); row <= 
                (row0 + neighborRadius); row++) {
                
                if ((row < 0) || (row > (this.height - 1))) {
                    continue;
                }
                
                int index = getInternalIndex(col, row);
                
                calculateColor(index);
            }
        }
    }
  
    public Image copyImage() {
       
        ImageExt img2 = new ImageExt(getWidth(), getHeight());
        
        System.arraycopy(r, 0, img2.r, 0, nPixels);
        System.arraycopy(g, 0, img2.g, 0, nPixels);
        System.arraycopy(b, 0, img2.b, 0, nPixels);
        
        System.arraycopy(cieX, 0, img2.cieX, 0, nPixels);
        System.arraycopy(cieY, 0, img2.cieY, 0, nPixels);
        System.arraycopy(hue, 0, img2.hue, 0, nPixels);
        System.arraycopy(saturation, 0, img2.saturation, 0, nPixels);
        System.arraycopy(brightness, 0, img2.brightness, 0, nPixels);
        System.arraycopy(extPopulated, 0, img2.extPopulated, 0, nPixels);
       
        return img2;
    }
   
    public void resetTo(Image copyThis) {
        
        if (copyThis.getNPixels() != nPixels) {
            throw new IllegalArgumentException("cannot convert this fixed " 
                + "image size to the size of copyThis");
        }
        
        if (!(copyThis instanceof ImageExt)) {
            throw new IllegalArgumentException(
            "copyThis has to be instance of ImageWithCIE");
        }
        
        System.arraycopy(copyThis.r, 0, r, 0, nPixels);
        System.arraycopy(copyThis.g, 0, g, 0, nPixels);
        System.arraycopy(copyThis.b, 0, b, 0, nPixels);
        
        System.arraycopy(((ImageExt)copyThis).cieX, 0, cieX, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).cieY, 0, cieY, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).hue, 0, hue, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).saturation, 0, saturation, 0, 
            nPixels);
        System.arraycopy(((ImageExt)copyThis).brightness, 0, brightness, 0, 
            nPixels);
        System.arraycopy(((ImageExt)copyThis).extPopulated, 0, extPopulated, 0, 
            nPixels);
    }

}
