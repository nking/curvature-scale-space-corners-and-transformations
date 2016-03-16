package algorithms.imageProcessing;

import java.awt.image.BufferedImage;

/**
 * space conserving image data holder.
 * 
 * @author nichole
 */
public class Image {
    
    //TODO:  add alpha when needed
    
    public final boolean is64Bit;
    
    /**
     * array holding red image values in range 0-255 and 4 values stored in
     * one int.
     */
    protected final int[] r;
    
    /**
     * 64 bit version of r array
     */
    protected final long[] rL;
    
    /**
     * array holding green image values in range 0-255 and 4 values stored in
     * one int.
     */
    protected final int[] g;
    
    /**
     * 64 bit version of g array
     */
    protected final long[] gL;
    
    /**
     * array holding blue image values in range 0-255 and 4 values stored in
     * one int.
     */
    protected final int[] b;
    
    /**
     * 64 bit version of b array
     */
    protected final long[] bL;
    
    protected final int width;
    
    protected final int height;
    
    protected final int nPixels;
    
    protected final int itemByteLength;
    
    protected final int len;
            
    /**
     * @param theWidth
     * @param theHeight
     */
    public Image (int theWidth, int theHeight) {
        
        String arch = System.getProperty("sun.arch.data.model");
        
        is64Bit = ((arch != null) && arch.equals("64")) ? true : false;
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        if (is64Bit) {
            
            itemByteLength = 8;
            
            len = (nPixels/itemByteLength) + 1;
            
            rL = new long[len];

            gL = new long[len];

            bL = new long[len];
            
            r = null;
            g = null;
            b = null;
            
        } else {
            
            itemByteLength = 4;
            
            len = (nPixels/itemByteLength) + 1;
                
            r = new int[len];

            g = new int[len];

            b = new int[len];
            
            rL = null;
            gL = null;
            bL = null;
        }
        
        init();
    }
    
    /**
     * package protected constructor for use in testing.
     * 
     * @param theWidth
     * @param theHeight
     */
    Image (int theWidth, int theHeight, boolean use32Bit) {
                
        is64Bit = !use32Bit;
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        if (is64Bit) {
            
            itemByteLength = 8;
            
            len = (nPixels/itemByteLength) + 1;
            
            rL = new long[len];

            gL = new long[len];

            bL = new long[len];
            
            r = null;
            g = null;
            b = null;
            
        } else {
            
            itemByteLength = 4;
            
            len = (nPixels/itemByteLength) + 1;
                
            r = new int[len];

            g = new int[len];

            b = new int[len];
            
            rL = null;
            gL = null;
            bL = null;
        }
        
        init();
    }
    
    protected void init() {
        // not used for base class
    }
    
    /**
     * given the pixel index which is a number representing
     * (row * width) + col for an uncondensed array,
     * return the internal array element index (for example the index
     * for the r array holding the value for pixelIdx would be at
     * r[elementIdx] which itself holds 4 values on a 32 bit system
     * or 8 on a 64 bit system.
     * 
     * @param pixelIdx
     * @return 
     */
    protected int getRowNumber(int pixelIdx) {
        
        int nthElement = (pixelIdx/itemByteLength);
        
        return nthElement;
    }
    
    /**
     * given the image column and row, return a pixel number which can be used
     * to retrieve data.
     * @param col
     * @param row
     * @return 
     */
    public int getInternalIndex(int col, int row) {
        return (row * width) + col;
    }
  
    /**
     * set the r, g, and b values for the pixel at location with col and row
     * @param col
     * @param row
     * @param rPix
     * @param gPix
     * @param bPix 
     */
    public void setRGB(int col, int row, int rPix, int gPix, int bPix) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException(
                "col is out of bounds");
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException(
                "row is out of bounds");
        }
        
        int idx = getInternalIndex(col, row);
       
        setRGB(idx, rPix, gPix, bPix);
    }
    
    /**
     * set the r, g, and b values for the pixel at location with col and row
     * @param pixIdx
     * @param rPix
     * @param gPix
     * @param bPix 
     */
    public void setRGB(int pixIdx, int rPix, int gPix, int bPix) {
        
        if ((pixIdx < 0) || (pixIdx > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
        
        if ((rPix > 255) || (rPix < 0)){
            throw new IllegalArgumentException(
                "values must be between 0 and 255, inclusive. rPix=" + rPix);
        }
        if ((gPix > 255) || (gPix < 0)){
            throw new IllegalArgumentException(
                "values must be between 0 and 255, inclusive. gPix=" + gPix);
        }
        if ((bPix > 255) || (bPix < 0)){
            throw new IllegalArgumentException(
                "values must be between 0 and 255, inclusive. bPix=" + bPix);
        }
        
        int elementIdx = pixIdx/itemByteLength;
        
        int byteNumber = pixIdx - (elementIdx*itemByteLength);
                
        if (is64Bit) {
            rL[elementIdx] = set64BitValue(rL[elementIdx], rPix, byteNumber);
            gL[elementIdx] = set64BitValue(gL[elementIdx], gPix, byteNumber);
            bL[elementIdx] = set64BitValue(bL[elementIdx], bPix, byteNumber);            
        } else {
            r[elementIdx] = set32BitValue(r[elementIdx], rPix, byteNumber);
            g[elementIdx] = set32BitValue(g[elementIdx], gPix, byteNumber);
            b[elementIdx] = set32BitValue(b[elementIdx], bPix, byteNumber);  
        }
        
    }
    
    protected long set64BitValue(long rowValue, int setValue, int byteNumber) {
        
        long shift = 8L * (long)byteNumber;
        
        long prevValue = (rowValue >> shift) & 255L;
        
        long shifted = (prevValue & 255L) << shift;
        
        rowValue -= shifted;
        
        shifted = (setValue & 255L) << shift;
                
        rowValue += shifted;
                
        assert(((rowValue >> shift) & 255L) == setValue);
        
        return rowValue;
    }
    
    protected int set32BitValue(int rowValue, int setValue, int byteNumber) {
            
        int shift = 8 * byteNumber;
        
        int prevValue = (rowValue >> shift) & 255;
        
        int shifted = (prevValue & 255) << shift;
        
        rowValue -= shifted;
        
        shifted = (setValue & 255) << shift;
                
        rowValue += shifted;
                
        assert(((rowValue >> shift) & 255) == setValue);
        
        return rowValue;
    }
    
    /**
     * set the r, g, and b values for the pixel at location with col and row
     * @param col
     * @param row
     * @param rgb 
     */
    public void setRGB(int col, int row, int rgb) {
        
        int idx = getInternalIndex(col, row);
        
        if ((idx < 0) || (idx > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
       
        int rPix = (rgb >> 16) & 0xFF;
        int gPix = (rgb >> 8) & 0xFF;
        int bPix = rgb & 0xFF;
        
        setRGB(col, row, rPix, gPix, bPix);
        
    }
        
    protected int get32BitValue(int rowValue, int byteNumber) {
                
        return (rowValue >> (byteNumber * 8)) & 255;
    }
    
    protected int get64BitValue(long rowValue, int byteNumber) {
                
        long shift = 8L * (long)byteNumber;
        
        long byteValue = (rowValue >> shift) & 255L;
        
        return (int)byteValue;
    }
    
    /**
     * get r value at image location col, row
     * @param col
     * @param row
     * @return 
     */
    public int getR(int col, int row) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException(
                "col is out of bounds");
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException(
                "row is out of bounds");
        }
        
        int idx = getInternalIndex(col, row);
       
        if ((idx < 0) || (idx > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
        
        int elementIdx = idx/itemByteLength;
        
        int byteNumber = idx - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(rL[elementIdx], byteNumber);
        } else {
            return get32BitValue(r[elementIdx], byteNumber);
        }
        
    }
    
    /**
     * get b value at image location col, row
     * @param col
     * @param row
     * @return 
     */
    public int getB(int col, int row) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException(
                "col is out of bounds");
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException(
                "row is out of bounds");
        }
        
        int idx = getInternalIndex(col, row);
       
        if ((idx < 0) || (idx > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
        
        int elementIdx = idx/itemByteLength;
        
        int byteNumber = idx - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(bL[elementIdx], byteNumber);
        } else {
            return get32BitValue(b[elementIdx], byteNumber);
        }
    }
    
    /**
     * get g value at image location col, row
     * @param col
     * @param row
     * @return 
     */
    public int getG(int col, int row) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException(
                "col is out of bounds");
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException(
                "row is out of bounds");
        }
        
        int idx = getInternalIndex(col, row);
       
        if ((idx < 0) || (idx > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "col and/or are out of bounds");
        }
        
        int elementIdx = idx/itemByteLength;
        
        int byteNumber = idx - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(gL[elementIdx], byteNumber);
        } else {
            return get32BitValue(g[elementIdx], byteNumber);
        }
    }
    
    /**
     * get r value at image pixel index derived from location col, row
     * @param pixelIndex
     * @return 
     */
    public int getR(int pixelIndex) {
               
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "pixelIndex is out of bounds");
        }
        
        int elementIdx = pixelIndex/itemByteLength;
        
        int byteNumber = pixelIndex - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(rL[elementIdx], byteNumber);
        } else {
            return get32BitValue(r[elementIdx], byteNumber);
        }
    }
    
    /**
     * get b value at image pixel index derived from location col, row
     * @param pixelIndex
     * @return 
     */
    public int getB(int pixelIndex) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "pixelIndex is out of bounds");
        }
        
        int elementIdx = pixelIndex/itemByteLength;
        
        int byteNumber = pixelIndex - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(bL[elementIdx], byteNumber);
        } else {
            return get32BitValue(b[elementIdx], byteNumber);
        }
    }
    
    /**
     * get g value at image pixel index derived from location col, row
     * @param pixelIndex
     * @return 
     */
    public int getG(int pixelIndex) {

        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "pixelIndex is out of bounds");
        }
        
        int elementIdx = pixelIndex/itemByteLength;
        
        int byteNumber = pixelIndex - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(gL[elementIdx], byteNumber);
        } else {
            return get32BitValue(g[elementIdx], byteNumber);
        }
    }
    
    /**
     * get the combined rgb value at image location col, row
     * @param col
     * @param row
     * @return 
     */
    public int getRGB(int col, int row) {
    
        int pixelIndex = getInternalIndex(col, row);
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        int rPix = getR(pixelIndex);
        int gPix = getG(pixelIndex);
        int bPix = getB(pixelIndex);
        
        int rgb = (((rPix & 0x0ff) << 16) 
            | ((gPix & 0x0ff) << 8) | (bPix & 0x0ff));
        
        return rgb;
    }
    
    /**
     * copy the image to another instance
     * @return 
     */
    public Image copyImage() {
       
        Image img2 = new Image(width, height, !is64Bit);
        
        if (is64Bit) {
            System.arraycopy(rL, 0, img2.rL, 0, len);
            System.arraycopy(gL, 0, img2.gL, 0, len);
            System.arraycopy(bL, 0, img2.bL, 0, len);
        } else {
            System.arraycopy(r, 0, img2.r, 0, len);
            System.arraycopy(g, 0, img2.g, 0, len);
            System.arraycopy(b, 0, img2.b, 0, len);
        }
       
        return img2;
    }
    
    public int[] getRValues() {
        
        int[] values = new int[nPixels];
        
        for (int i = 0; i < nPixels; ++i) {
            values[i] = getR(i);
        }
        
        return values;
    }
    
    public int[] getGValues() {
        
        int[] values = new int[nPixels];
        
        for (int i = 0; i < nPixels; ++i) {
            values[i] = getG(i);
        }
        
        return values;
    }
    
    public int[] getBValues() {
        
        int[] values = new int[nPixels];
        
        for (int i = 0; i < nPixels; ++i) {
            values[i] = getB(i);
        }
        
        return values;
    }
    
    /**
     * copy the image to another instance, but of subclass type ImageExt
     * @return 
     */
    public ImageExt copyToImageExt() {
       
        ImageExt img2 = new ImageExt(width, height, !is64Bit);
        
        if (is64Bit) {
            System.arraycopy(rL, 0, img2.rL, 0, len);
            System.arraycopy(gL, 0, img2.gL, 0, len);
            System.arraycopy(bL, 0, img2.bL, 0, len);
        } else {
            System.arraycopy(r, 0, img2.r, 0, len);
            System.arraycopy(g, 0, img2.g, 0, len);
            System.arraycopy(b, 0, img2.b, 0, len);
        }
       
        return img2;
    }
    
    /**
     * using the color tables of awt and BufferedImage, convert this image
     * to greyscale.
     * @return 
     */
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
        
        GreyscaleImage.Type type = is64Bit ? GreyscaleImage.Type.Bits64 :
            GreyscaleImage.Type.Bits32;
        
        GreyscaleImage out = new GreyscaleImage(width, height, type);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                // presumably, this is already combined?
                // or does it need separation into rgb and then averaged?
                
                int rgb = outputImage.getRGB(i, j);
                                
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;  
                    
                int v = Math.round((r + g + b)/3.f);
                
                assert(v <= 255);
                assert(v >= 0);
                
                out.setValue(i, j, v);
            }
        }
        
        return out;
    }
    
    /**
     * using the color tables of awt and BufferedImage, convert this image
     * to greyscale.
     * @return 
     */
    public GreyscaleImage copyRedToGreyscale() {
        
        GreyscaleImage.Type type = is64Bit ? GreyscaleImage.Type.Bits64 :
            GreyscaleImage.Type.Bits32;
        
        GreyscaleImage out = new GreyscaleImage(width, height, type);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                int rgb = getRGB(i, j);
                                
                int r = (rgb >> 16) & 0xFF;
                                    
                out.setValue(i, j, r);
            }
        }
        
        return out;
    }
    
    /**
     * using the color tables of awt and BufferedImage, convert this image
     * to greyscale.
     * @return 
     */
    public GreyscaleImage copyGreenToGreyscale() {
        
        GreyscaleImage.Type type = is64Bit ? GreyscaleImage.Type.Bits64 :
            GreyscaleImage.Type.Bits32;
        
        GreyscaleImage out = new GreyscaleImage(width, height, type);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                int rgb = getRGB(i, j);
                                
                int g = (rgb >> 8) & 0xFF;
                                    
                out.setValue(i, j, g);
            }
        }
        
        return out;
    }
    
    /**
     * using the color tables of awt and BufferedImage, convert this image
     * to greyscale.
     * @return 
     */
    public GreyscaleImage copyBlueToGreyscale() {
        
        GreyscaleImage.Type type = is64Bit ? GreyscaleImage.Type.Bits64 :
            GreyscaleImage.Type.Bits32;
        
        GreyscaleImage out = new GreyscaleImage(width, height, type);
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                int rgb = getRGB(i, j);
                                
                int b = rgb & 0xFF;  
                                    
                out.setValue(i, j, b);
            }
        }
        
        return out;
    }
    
    /**
     * reset the contents of this image to the contents of copyThis, but
     * only if the dimensions are the same.
     * 
     * @param copyThis 
     */
    public void resetTo(Image copyThis) {
        
        if (copyThis.getNPixels() != nPixels) {
            throw new IllegalArgumentException("cannot convert this fixed " 
                + "image size to the size of copyThis");
        }
        
        if (copyThis.width != width) {
            throw new IllegalArgumentException(
                "copyThis.width must be same as this width");
        }
        if (copyThis.height != height) {
            throw new IllegalArgumentException(
                "copyThis.height must be same as this width");
        }
        
        if (is64Bit) {
            System.arraycopy(copyThis.rL, 0, rL, 0, copyThis.len);
            System.arraycopy(copyThis.gL, 0, gL, 0, copyThis.len);
            System.arraycopy(copyThis.bL, 0, bL, 0, copyThis.len);
        } else {
            System.arraycopy(copyThis.r, 0, r, 0, copyThis.len);
            System.arraycopy(copyThis.g, 0, g, 0, copyThis.len);
            System.arraycopy(copyThis.b, 0, b, 0, copyThis.len);
        }
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
    
    /**
     * get the image pixel col location from the pixel index
     * @param internalIndex
     * @return 
     */
    public int getCol(int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        int row = internalIndex/width;
        
        int col = internalIndex - (row * width);
        
        return col;
    }
    
    /**
     * get the image pixel row location from the pixel index
     * @param internalIndex
     * @return 
     */
    public int getRow(int internalIndex) {
        
        if ((internalIndex < 0) || (internalIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        int row = internalIndex/width;
        
        return row;
    }
    
    public void debugPrint() {
        StringBuilder sb = new StringBuilder();
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                int r = getR(col, row);
                int g = getG(col, row);
                int b = getB(col, row);
                String str = String.format("(%3d %3d %3d) ", r, g, b);
                sb.append(str);
            }
            sb.append("\n");
        }
        System.out.println(sb.toString());
    }
}
