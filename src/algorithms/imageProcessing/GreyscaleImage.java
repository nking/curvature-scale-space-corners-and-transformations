package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class GreyscaleImage {
    
    private int[] a;
    
    private int width;
    
    private int height;
    
    private int nPixels;
    
    private int xRelativeOffset = 0;
    
    private int yRelativeOffset = 0;
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public GreyscaleImage (int theWidth, int theHeight) {
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        a = new int[nPixels];
    }
    
    /**
     * convenience method to fill image with the value. 
     * @param value 
     */
    public void fill(int value) {
        Arrays.fill(a, value);
    }
    
    public void multiply(int factor) {
        
        for (int i = 0; i < nPixels; i++) {
            a[i] *= factor;
        }
    }
    
    public void divide(int number) {
        
        for (int i = 0; i < nPixels; i++) {
            a[i] /= number;
        }
    }
    
    public void normalizeToMax255() {
        
        int max = MiscMath.findMax(a);
        if (max == 0) {
            return;
        }
        int maxIdx = MiscMath.findYMaxIndex(a);
        int c = getCol(maxIdx);
        int r = getRow(maxIdx);
        
        double factor = 255./max;
        
        for (int i = 0; i < nPixels; i++) {
            double v = a[i] * factor;
            a[i] = (int)v;
        }
        
    }
    
    public void setValue(int col, int row, int value) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException("col is out of bounds: col=" 
                + col + " width=" + width);
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException("row is out of bounds: row=" 
                + row + " height=" + height);
        }
        
        int idx = (row * width) + col;
        
        if ((idx < 0) || (idx > (a.length - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
        
        /*
        int rPix = (value >> 16) & 0xFF;
        int gPix = (value >> 8) & 0xFF;
        int bPix = value & 0xFF;
        */
       
        a[idx] = value;
    }
    
    public void getXY(PairIntArray array, int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        // int idx = (row * width) + col;
        //     idx/w = row  + 0;
        int row = internalIndex/width;
        
        int col = internalIndex - (row * width);
        
        array.add(col, row);
    }
    
    public int getRow(int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        // int idx = (row * width) + col;
        //     idx/w = row  + 0;
        int row = internalIndex/width;
        
        return row;
    }
    
    public int getCol(int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        // int idx = (row * width) + col;
        //     idx/w = row  + 0;
        int row = internalIndex/width;
        
        int col = internalIndex - (row * width);
        
        return col;
    }
    
    public void setValue(int internalIndex, int value) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        a[internalIndex] = value;
    }
    
    public int getValue(int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        return a[internalIndex];
    }
    
    public int getIndex(int col, int row) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException("col is out of bounds: col=" 
                + col + " width=" + width);
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException("row is out of bounds: row=" 
                + row + " height=" + height);
        }
        
        int idx = (row * width) + col;
        
        return idx;
    }
        
    public int getValue(int col, int row) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException("col is out of bounds: col=" 
                + col + " width=" + width);
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException("row is out of bounds: row=" 
                + row + " height=" + height);
        }
        
        int idx = (row * width) + col;
       
        if (idx > a.length) {
            throw new IllegalArgumentException("coords are out of bounds");
        }
        /*
        int rgb = (int)(((r[idx] & 0x0ff) << 16) 
            | ((g[idx] & 0x0ff) << 8) | (b[idx] & 0x0ff));
        */
        return a[idx];
    }
    
    public int[] getValues() {
        return a;
    }
 
    public int getWidth() {
        return width;
    }
    
    public int getHeight() {
        return height;
    }
    
    public float[] getFloatValues() {
        float[] t = new float[nPixels];
        for (int i = 0; i < nPixels; ++i) {
            t[i] = a[i];
        }
        return t;
    }
    
    public GreyscaleImage copyImage() {
       
        GreyscaleImage img2 = createWithDimensions();
                
        System.arraycopy(a, 0, img2.a, 0, nPixels);
        
        img2.nPixels = nPixels;
        
        return img2;
    }
    
    public GreyscaleImage createWithDimensions() {
       
        GreyscaleImage img2 = new GreyscaleImage(width, height);
                
        img2.xRelativeOffset = xRelativeOffset;
        img2.yRelativeOffset = yRelativeOffset;
        
        return img2;
    }
    
    public GreyscaleImage subImage(int xCenter, int yCenter, int subWidth, 
        int subHeight) {
       
        GreyscaleImage img2 = new GreyscaleImage(subWidth, subHeight);
                
        int col2 = 0;
        for (int col = (xCenter - (subWidth/2)); col < (xCenter + (subWidth/2));
            ++col) {
            
            int row2 = 0;
            for (int row = (yCenter - (subHeight/2)); row < (yCenter + (subHeight/2));
            ++row) {
                
                int v = getValue(col, row);
                
                img2.setValue(col2, row2, v);
                
                row2++;
            }
            
            col2++;
        }
                
        return img2;
    }
    
    public Image copyImageToGreen() {
        
        Image img2 = new Image(width, height);
        
        System.arraycopy(a, 0, img2.g, 0, nPixels);
        
        return img2;
    }
    
    public ImageExt createColorGreyscaleExt() {
        
        ImageExt img2 = new ImageExt(width, height);
        
        System.arraycopy(a, 0, img2.r, 0, nPixels);
        System.arraycopy(a, 0, img2.g, 0, nPixels);
        System.arraycopy(a, 0, img2.b, 0, nPixels);
        
        return img2;
    }
    
    public void resetTo(final GreyscaleImage copyThis) {
        
        if (copyThis.nPixels == nPixels) {
            
            System.arraycopy(copyThis.a, 0, a, 0, nPixels);
            
        } else {
            
            a = copyThis.a;
    
            width = copyThis.width;
    
            height = copyThis.height;
    
            nPixels = copyThis.nPixels;
            
            xRelativeOffset = copyThis.xRelativeOffset;
            
            yRelativeOffset = copyThis.yRelativeOffset;
        }
        
    }

    /**
     * @return the nPixels
     */
    public int getNPixels() {
        return nPixels;
    }

    /**
     * @return the xRelativeOffset
     */
    public int getXRelativeOffset() {
        return xRelativeOffset;
    }

    /**
     * @param xRelativeOffset the xRelativeOffset to set
     */
    public void setXRelativeOffset(int xRelativeOffset) {
        this.xRelativeOffset = xRelativeOffset;
    }

    /**
     * @return the yRelativeOffset
     */
    public int getYRelativeOffset() {
        return yRelativeOffset;
    }

    /**
     * @param yRelativeOffset the yRelativeOffset to set
     */
    public void setYRelativeOffset(int yRelativeOffset) {
        this.yRelativeOffset = yRelativeOffset;
    }

    public boolean isThisSize(int[] offsetsAndDimensions) {
        
        int w = offsetsAndDimensions[2];
        int h = offsetsAndDimensions[3];
        
        int offsetX = offsetsAndDimensions[0];
        int offsetY = offsetsAndDimensions[1];
        
        if (this.width != w || this.height != h || 
            this.xRelativeOffset != offsetX || this.yRelativeOffset != offsetY) {
            return false;
        }
        
        return true;
    }
}
