package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
import java.awt.Color;

/**
 * image specialized to hold CIE 1931 xy chromaticity coordinates (colors) 
 * and HSB colors populated on demand.
 * 
 * @author nichole
 */
public class ImageLtExt extends ImageLt {
    
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
    
    /**
     * luma is the y of yuv.
     */
    protected float[] luma;
    
    //TODO: this could be a bit vector instead
    private boolean[] extPopulated;
    
    private int radiusForPopulateOnDemand = 1;
    
    static public final double[][] rgbToLumaMatrix = new double[3][];
    static {
        rgbToLumaMatrix[0] = new double[]{0.256, 0.504, 0.098};
        rgbToLumaMatrix[1] = new double[]{-0.148, -0.291, 0.439};
        rgbToLumaMatrix[2] = new double[]{0.439, -0.368, -0.072};
    }
    
    protected CIEChromaticity cieC = new CIEChromaticity();
        
    /**
     * @param theWidth
     * @param theHeight
     */
    public ImageLtExt (int theWidth, int theHeight) {
        
        super(theWidth, theHeight);
    }
    
    protected void init() {
        
        this.cieX = new float[nPixels];
        this.cieY = new float[nPixels];
        this.hue = new float[nPixels];
        this.saturation = new float[nPixels];
        this.brightness = new float[nPixels];
        this.luma = new float[nPixels];
        
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
    
    public float getLuma(int col, int row) {
        
        int idx = getInternalIndex(col, row);
       
        return getLuma(idx);
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
    
    public float getLuma(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (!extPopulated[internalIndex]) {
            calculateColorIncludingNeighbors(internalIndex, 
                radiusForPopulateOnDemand);
        }
       
        return luma[internalIndex];
    }
  
    protected void calculateColor(int idx) {
        
        if (extPopulated[idx]) {
            return;
        }
        
        int rPix = getR(idx);
        int gPix = getG(idx);
        int bPix = getB(idx);

        float[] xy = cieC.rgbToXYChromaticity(rPix, gPix, bPix);
        
        cieX[idx] = xy[0];
        cieY[idx] = xy[1];
        
        float[] hsb = new float[3];
        Color.RGBtoHSB(rPix, gPix, bPix, hsb);
       
        hue[idx] = hsb[0];
        saturation[idx] = hsb[1];
        brightness[idx] = hsb[2];
        
        double[] yuv = MatrixUtil.multiply(rgbToLumaMatrix, 
            new double[]{rPix, gPix, bPix});
        
        luma[idx] = (float)yuv[0];
        
        extPopulated[idx] = true;
    }
    
    public void calculateColorIncludingNeighbors(int idx, int neighborRadius) {
        
        if (neighborRadius == 0) {
            calculateColor(idx);
            return;
        }
        
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

    @Override
    public ImageLt copyImage() {
    
        ImageLtExt img2 = copyToImageExt();
        
        System.arraycopy(cieX, 0, img2.cieX, 0, nPixels);
        System.arraycopy(cieY, 0, img2.cieY, 0, nPixels);
        System.arraycopy(hue, 0, img2.hue, 0, nPixels);
        System.arraycopy(saturation, 0, img2.saturation, 0, nPixels);
        System.arraycopy(brightness, 0, img2.brightness, 0, nPixels);
        System.arraycopy(luma, 0, img2.luma, 0, nPixels);
        System.arraycopy(extPopulated, 0, img2.extPopulated, 0, nPixels);
       
        return img2;
    }
   
    public void resetTo(ImageLt copyThis) {
        
        if (copyThis.getNPixels() != nPixels) {
            throw new IllegalArgumentException("cannot convert this fixed " 
                + "image size to the size of copyThis");
        }
        
        if (!(copyThis instanceof ImageLtExt)) {
            throw new IllegalArgumentException(
            "copyThis has to be instance of ImageWithCIE");
        }
        
        super.resetTo(copyThis);
        
        System.arraycopy(((ImageLtExt)copyThis).cieX, 0, cieX, 0, nPixels);
        System.arraycopy(((ImageLtExt)copyThis).cieY, 0, cieY, 0, nPixels);
        System.arraycopy(((ImageLtExt)copyThis).hue, 0, hue, 0, nPixels);
        System.arraycopy(((ImageLtExt)copyThis).saturation, 0, saturation, 0, 
            nPixels);
        System.arraycopy(((ImageLtExt)copyThis).brightness, 0, brightness, 0, 
            nPixels);
        System.arraycopy(((ImageLtExt)copyThis).luma, 0, luma, 0, nPixels);
        System.arraycopy(((ImageLtExt)copyThis).extPopulated, 0, extPopulated, 0, 
            nPixels);
    }

}
