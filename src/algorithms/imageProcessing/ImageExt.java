package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
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
    
    //CIE LAB 1994
    protected float[][] lab;
   
    private enum CIE {
        DEFAULT, LAB1931, LUV1976
    }
    private CIE cieType = CIE.DEFAULT;
    
    // a guard to prevent mixing of cie 1941 and 1931
    private boolean aCIECalcOccurred = false;
    
    public void overrideToUseCIELAB1931() {
        if (cieType.equals(CIE.LAB1931)) {
            return;
        }
        if (aCIECalcOccurred) {
            //TODO: may change this to a warning in the future
            // and wipe out existing CIE values.
            throw new IllegalStateException(
                "values have already been stored"
                + " as cie 1994, so use reset CIELAB first.");
        }
        cieType = CIE.LAB1931;
    }
    
    public void overrideToUseCIELUV1976() {
        if (cieType.equals(CIE.LUV1976)) {
            return;
        }
        if (aCIECalcOccurred) {
            //TODO: may change this to a warning in the future
            // and wipe out existing CIE values.
            throw new IllegalStateException(
                "values have already been stored"
                + " as cie 1994, so use reset CIELAB first.");
        }
        cieType = CIE.LUV1976;
    }
    
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
    public ImageExt (int theWidth, int theHeight) {
        
        super(theWidth, theHeight);
    }
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public ImageExt (int theWidth, int theHeight, boolean use32Bit) {
        
        super(theWidth, theHeight, use32Bit);
    }
    
    protected void init() {
        
        this.cieX = new float[nPixels];
        this.cieY = new float[nPixels];
        this.hue = new float[nPixels];
        this.saturation = new float[nPixels];
        this.brightness = new float[nPixels];
        this.luma = new float[nPixels];
        
        this.extPopulated = new boolean[nPixels];
        
        this.lab = new float[nPixels][];
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
    
    public float[] getCIELAB(int col, int row) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException(
                "col is out of bounds:");
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException(
                "row is out of bounds:");
        }
        
        int idx = getInternalIndex(col, row);
        
        return getCIELAB(idx);
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
    
    public float[] getCIEXY_(int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        int rPix = getR(internalIndex);
        int gPix = getG(internalIndex);
        int bPix = getB(internalIndex);

        float[] xy = cieC._rgbToXYChromaticity(rPix, gPix, bPix);
        
        return xy;
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
  
    public float[] getCIELAB(int internalIndex) {
                
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (lab[internalIndex] == null) {
            int rPix = getR(internalIndex);
            int gPix = getG(internalIndex);
            int bPix = getB(internalIndex);
            
            float[] cieLAB = null;
            if (cieType.equals(CIE.DEFAULT)) {
                cieLAB = cieC.rgbToCIELAB(rPix, gPix, bPix);
            } else if (cieType.equals(CIE.LAB1931)) {
                // this uses cie lab 1931
                cieLAB = cieC.rgbToCIELAB1931(rPix, gPix, bPix);
            } else if (cieType.equals(CIE.LUV1976)) {
                cieLAB = cieC.rgbToCIELUV(rPix, gPix, bPix);
            }            
            lab[internalIndex] = cieLAB;
        
            aCIECalcOccurred = true;
        }
       
        return lab[internalIndex];
    }
    
    protected void calculateColor(int idx) {
        
        if (extPopulated[idx]) {
            return;
        }
        
        int rPix = getR(idx);
        int gPix = getG(idx);
        int bPix = getB(idx);

        // TODO: revisit this... should be using the normalized r,g,b method
        // possibly the same for the HSB and YUV...
        // _rgbToXYChromaticity
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
    
    @Override
    protected Image createNewImage(int width, int height) {
        return new ImageExt(width, height, !is64Bit);
    }
    
    @Override
    public ImageExt createWithDimensions() {
       
        ImageExt img2 = new ImageExt(width, height, is64Bit);
        
        return img2;
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
    public ImageExt createWithDimensions(int theWidth, int theHeight) {
       
        ImageExt img2 = new ImageExt(theWidth, theHeight, is64Bit);
        
        return img2;
    }
    
    @Override
    public Image copyImage() {
    
        ImageExt img2 = copyToImageExt();
        
        System.arraycopy(cieX, 0, img2.cieX, 0, nPixels);
        System.arraycopy(cieY, 0, img2.cieY, 0, nPixels);
        System.arraycopy(hue, 0, img2.hue, 0, nPixels);
        System.arraycopy(saturation, 0, img2.saturation, 0, nPixels);
        System.arraycopy(brightness, 0, img2.brightness, 0, nPixels);
        System.arraycopy(luma, 0, img2.luma, 0, nPixels);
        System.arraycopy(extPopulated, 0, img2.extPopulated, 0, nPixels);
       
        return img2;
    }
   
    /**
     * copy the image to another instance
     * @param x0 inclusive start x coordinate of subimage
     * @param x1 exclusive stop x coordinate of sub image
     * @param y0 inclusive
     * @param y1 exclusive
     * @return 
     */
    @Override
    public Image copySubImage(int x0, int x1, int y0, int y1) {
        
        ImageExt img2 = (ImageExt) super.copySubImage(x0, x1, y0, y1);
               
        for (int i = x0; i < x1; ++i) {
            for (int j = y0; j < y1; ++j) {
                int pixIdx = getInternalIndex(i, j);
                int pixIdx2 = img2.getInternalIndex(i - x0, j - y0);
                img2.cieX[pixIdx2] = cieX[pixIdx];
                img2.cieY[pixIdx2] = cieY[pixIdx];
                img2.hue[pixIdx2] = hue[pixIdx];
                img2.saturation[pixIdx2] = saturation[pixIdx];
                img2.brightness[pixIdx2] = brightness[pixIdx];
                img2.luma[pixIdx2] = luma[pixIdx];
                img2.extPopulated[pixIdx2] = extPopulated[pixIdx];
            }
        }
       
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
        
        super.resetTo(copyThis);
        
        System.arraycopy(((ImageExt)copyThis).cieX, 0, cieX, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).cieY, 0, cieY, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).hue, 0, hue, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).saturation, 0, saturation, 0, 
            nPixels);
        System.arraycopy(((ImageExt)copyThis).brightness, 0, brightness, 0, 
            nPixels);
        System.arraycopy(((ImageExt)copyThis).luma, 0, luma, 0, nPixels);
        System.arraycopy(((ImageExt)copyThis).extPopulated, 0, extPopulated, 0, 
            nPixels);
    }

}
