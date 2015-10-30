package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class GreyscaleImageLt {
    
    public final boolean is64Bit;
    
    /**
     * array holding image values in range 0-255 and 4 values stored in
     * one int.
     */
    final int[] a;
    
    /**
     * 64 bit version of a array
     */
    final long[] aL;
    
    protected final int width;
    
    protected final int height;
    
    protected final int nPixels;
    
    protected final int itemByteLength;
    
    protected final int len;
    
    private int xRelativeOffset = 0;
    
    private int yRelativeOffset = 0;
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public GreyscaleImageLt (int theWidth, int theHeight) {
        
        String arch = System.getProperty("sun.arch.data.model");
        
        is64Bit = ((arch != null) && arch.equals("64")) ? true : false;
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        if (is64Bit) {
            
            itemByteLength = 8;
            
            len = (nPixels/itemByteLength) + 1;
            
            aL = new long[len];
            
            a = null;
            
        } else {
            
            itemByteLength = 4;
            
            len = (nPixels/itemByteLength) + 1;
                
            a = new int[len];
            
            aL = null;
        }
        
    }
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public GreyscaleImageLt (int theWidth, int theHeight, boolean use32Bit) {
                
        is64Bit = !use32Bit;
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        if (is64Bit) {
            
            itemByteLength = 8;
            
            len = (nPixels/itemByteLength) + 1;
            
            aL = new long[len];
            
            a = null;
            
        } else {
            
            itemByteLength = 4;
            
            len = (nPixels/itemByteLength) + 1;
                
            a = new int[len];
            
            aL = null;
        }
        
    }
    
    /**
     * convenience method to fill image with the value. 
     * @param value 
     */
    public void fill(int value) {
        
        if (value < 0 || (value > 255)) {
            throw new IllegalArgumentException(
                "value must be between 0 and 255 inclusive");
        }
        
        if (is64Bit) {
            
            long total = 0;
            for (int i = 0; i < 8; i++) {
                total += ((value & 255L) << (i * 8L));
            }
            Arrays.fill(aL, total);
            
        } else {
            
            int total = 0;
            for (int i = 0; i < 4; i++) {
                total += ((value & 255) << (i * 8));
            }
            Arrays.fill(a, total);
            
        }
        
    }
    
    /**
     * convenience method to multiply entire image by factor.  Note that
     * if the result for a pixel is > 255, it is set to 255.
     * @param factor a positive number that results in pixel values in range
     * 0 to 255;
     */
    public void multiply(float factor) {
        
        if (factor < 0) {
            throw new IllegalArgumentException(
                "factor has to be a positive number");
        }
        
        if (is64Bit) {
            
            for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
                long total = aL[elementIdx];
                long total2 = 0;
                for (int i = 0; i < 8; i++) {
                    long v = (total >> (i * 8L)) & 255L;
                    v = (long)(v * factor);
                    total2 += ((v & 255L) << (i * 8L));
                }
                aL[elementIdx] = total2;
            }
            
        } else {
            
            for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
                int total = a[elementIdx];
                int total2 = 0;
                for (int i = 0; i < 8; i++) {
                    long v = (total >> (i * 8)) & 255;
                    v = (long)(v * factor);
                    total2 += ((v & 255) << (i * 8));
                }
                a[elementIdx] = total2;
            }
        }
        
    }
    
    public void divide(int number) {
        
        if (number < 0) {
            throw new IllegalArgumentException(
                "number has to be a positive number");
        }
        
        if (is64Bit) {
            
            long f = number;
            for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
                long total = aL[elementIdx];
                long total2 = 0;
                for (int i = 0; i < 8; i++) {
                    long v = (total >> (i * 8L)) & 255L;
                    v = v / f;
                    total2 += ((v & 255L) << (i * 8L));
                }
                aL[elementIdx] = total2;
            }
            
        } else {
            
            for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
                int total = a[elementIdx];
                int total2 = 0;
                for (int i = 0; i < 8; i++) {
                    long v = (total >> (i * 8)) & 255;
                    v = v / number;
                    total2 += ((v & 255) << (i * 8));
                }
                a[elementIdx] = total2;
            }
        }
    }
    
    public void normalizeToMax255() {
        
        int max;
        
        if (is64Bit) {
            max = (int)MiscMath.findMaxForByteCompressed(aL, len, 8);
        } else {
            max = MiscMath.findMaxForByteCompressed(a, len, 4);
        }
        
        if (max == 0) {
            return;
        }
        
        float factor = 255.f/(float)max;
        
        multiply(factor);
        
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
    
    public void setValue(int col, int row, int value) {
        
        if ((col < 0) || (col > (width - 1))) {
            throw new IllegalArgumentException("col is out of bounds: col=" 
                + col + " width=" + width);
        }
        if ((row < 0) || (row > (height - 1))) {
            throw new IllegalArgumentException("row is out of bounds: row=" 
                + row + " height=" + height);
        }
        
        int idx = getInternalIndex(col, row);
        
        setValue(idx, value);
       
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
    
    public void setValue(int pixelIndex, int value) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if ((value > 255) || (value < 0)) {
            throw new IllegalArgumentException(
                "value must be between 0 and 255, inclusive");
        }
        
        int elementIdx = pixelIndex/itemByteLength;
        
        int byteNumber = pixelIndex - (elementIdx*itemByteLength);
                
        if (is64Bit) {
            aL[elementIdx] = set64BitValue(aL[elementIdx], value, byteNumber);
        } else {
            a[elementIdx] = set32BitValue(a[elementIdx], value, byteNumber);
        }
    }
    
    protected int get32BitValue(int rowValue, int byteNumber) {
                
        return (rowValue >> (byteNumber * 8)) & 255;
    }
    
    protected int get64BitValue(long rowValue, int byteNumber) {
                
        return (int)((rowValue >> (byteNumber * 8L)) & 255L);
    }
    
    public int getValue(int pixelIndex) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "pixelIndex is out of bounds:");
        }
        
        int elementIdx = pixelIndex/itemByteLength;
        
        int byteNumber = pixelIndex - (elementIdx*itemByteLength);
        
        if (is64Bit) {
            return get64BitValue(aL[elementIdx], byteNumber);
        } else {
            return get32BitValue(a[elementIdx], byteNumber);
        }
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
        
        int idx = getInternalIndex(col, row);
        
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
        
        int idx = getInternalIndex(col, row);
       
        return getValue(idx);
    }
    
    public int getWidth() {
        return width;
    }
    
    public int getHeight() {
        return height;
    }
    
    public float[] getFloatValues() {
        
        float[] t = new float[nPixels];
        
        int count = 0;
        
        if (is64Bit) {
            for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
                long total = aL[elementIdx];
                for (int i = 0; i < 8; i++) {
                    long v = (total >> (i * 8L)) & 255L;
                    t[count] = (int)v;
                    count++;
                }
            }
        } else {
            for (int elementIdx = 0; elementIdx < len; ++elementIdx) {
                int total = a[elementIdx];
                for (int i = 0; i < 4; i++) {
                    int v = (total >> (i * 8)) & 255;
                    t[count] = v;
                    count++;
                }
            }
        }
        
        return t;
    }
    
    public GreyscaleImageLt copyImage() {
       
        GreyscaleImageLt img2 = createWithDimensions();
                
        if (is64Bit) {
            System.arraycopy(aL, 0, img2.aL, 0, len);
        } else {
            System.arraycopy(a, 0, img2.a, 0, len);
        }
        
        return img2;
    }
    
    public GreyscaleImageLt createWithDimensions() {
       
        GreyscaleImageLt img2 = new GreyscaleImageLt(width, height, !is64Bit);
                
        img2.xRelativeOffset = xRelativeOffset;
        img2.yRelativeOffset = yRelativeOffset;
        
        return img2;
    }
    
    public GreyscaleImageLt subImage(int xCenter, int yCenter, int subWidth, 
        int subHeight) {
       
        GreyscaleImageLt img2 = new GreyscaleImageLt(subWidth, subHeight, !is64Bit);
                
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
    
    public ImageLtExt createColorGreyscaleLtExt() {
        
        ImageLtExt img2 = new ImageLtExt(width, height);
        
        if (is64Bit) {
            System.arraycopy(aL, 0, img2.rL, 0, len);
            System.arraycopy(aL, 0, img2.gL, 0, len);
            System.arraycopy(aL, 0, img2.bL, 0, len);
        } else {
            System.arraycopy(a, 0, img2.r, 0, len);
            System.arraycopy(a, 0, img2.g, 0, len);
            System.arraycopy(a, 0, img2.b, 0, len);
        }
        
        return img2;
    }
    
    public ImageLt createColorGreyscaleLt() {
        
        ImageLt img2 = new ImageLt(width, height);
        
        if (is64Bit) {
            System.arraycopy(aL, 0, img2.rL, 0, len);
            System.arraycopy(aL, 0, img2.gL, 0, len);
            System.arraycopy(aL, 0, img2.bL, 0, len);
        } else {
            System.arraycopy(a, 0, img2.r, 0, len);
            System.arraycopy(a, 0, img2.g, 0, len);
            System.arraycopy(a, 0, img2.b, 0, len);
        }
        
        return img2;
    }
    
    public void resetTo(final GreyscaleImageLt copyThis) {
        
        if (copyThis.getNPixels() != nPixels) {
            throw new IllegalArgumentException("cannot convert this fixed " 
                + "image size to the size of copyThis.");
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
            System.arraycopy(copyThis.aL, 0, aL, 0, copyThis.len);
        } else {
            System.arraycopy(copyThis.a, 0, a, 0, copyThis.len);
        }
        
        xRelativeOffset = copyThis.xRelativeOffset;
            
        yRelativeOffset = copyThis.yRelativeOffset;
        
    }
    
    
    public void debugPrint() {
        StringBuilder sb = new StringBuilder();
        for (int row = 0; row < height; ++row) {
            for (int col = 0; col < width; ++col) {
                int v = getValue(col, row);
                String str = String.format("(%3d) ", v);
                sb.append(str);
            }
            sb.append("\n");
        }
        System.out.println(sb.toString());
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
