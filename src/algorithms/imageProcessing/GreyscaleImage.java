package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;

import java.util.Arrays;

/**
 * class to hold data with values in range 0-255 or -255 to 255 depending
 * upon constructor used (the later is needed for gradient images).
 * 
 * @author nichole
 */
public class GreyscaleImage {

    /**
     * type for how data is stored which is by default Bits32 for a 32 bit
     * platform and Bits64 for a 64 bit platform each holding values in 
     * range 0 to 255, inclusive, else user can 
     * specify Bits32Signed or Bits64Signed which is the Bits32 or Bits64
     * storage with a range of -255 to 255 (needed for gradient images for
     * example) else the Bits32FullRangeInt type can hold signed integers.
     */
    public static enum Type {
        Bits32, Bits64, Bits32Signed, Bits64Signed, Bits32FullRangeInt;
        public Type copy(Type t) {
            return Type.values()[t.ordinal()];
        }
    }
    
    private Type type;
    
    public boolean is64Bit;
    
    /**
     * array holding image values in range 0-255 (default) 
     * or -255 to 255 depending upon constructor used.
     */
    private int[] a;
    
    /**
     * 64 bit version of a array
     */
    private long[] aL;
    
    private int width;
    
    private int height;
    
    private int minAllowed;
    
    private int maxAllowed;
    
    private int nPixels;
    
    // number of datum in one item of the array a or aL
    private int itemNDatum;
    
    private int datumNBits;
    
    private int len;
    
    private int xRelativeOffset = 0;
    
    private int yRelativeOffset = 0;
    
    /**
     * @param theWidth
     * @param theHeight
     */
    public GreyscaleImage (int theWidth, int theHeight) {
        
        if (theWidth > 46340 || theWidth < 0) {
            throw new IllegalArgumentException("theWidth must be between 0 "
            + " and 46340");
        }
        if (theHeight > 46340 || theHeight < 0) {
            throw new IllegalArgumentException("theHeight must be between 0 "
            + " and 46340");
        }
        
        String arch = System.getProperty("sun.arch.data.model");
        
        is64Bit = ((arch != null) && arch.equals("64")) ? true : false;
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        int wordSize;
        if (is64Bit) {
            wordSize = 64;
            type = Type.Bits64;
        } else {
            wordSize = 32;
            type = Type.Bits32;
        }
        
        datumNBits = 8;
        
        minAllowed = 0;
            
        maxAllowed = 255;
        
        itemNDatum = wordSize/datumNBits;
                        
        len = (nPixels/itemNDatum) + 1;
            
        if (is64Bit) {
                        
            aL = new long[len];
            
            a = null;
            
        } else {
            
            a = new int[len];
            
            aL = null;
        }
       
    }
    
    /**
     * @param theWidth
     * @param theHeight
     * @param theType image type
     */
    public GreyscaleImage (int theWidth, int theHeight, Type theType) {
        
        if (theWidth > 46340 || theWidth < 0) {
            throw new IllegalArgumentException("theWidth must be between 0 "
            + " and 46340");
        }
        if (theHeight > 46340 || theHeight < 0) {
            throw new IllegalArgumentException("theHeight must be between 0 "
            + " and 46340");
        }
        
        if (theType == null) {
            String arch = System.getProperty("sun.arch.data.model");        
            if ((arch != null) && arch.equals("64")) {
                type = Type.Bits64;
            } else {
                type = Type.Bits32;
            }
        } else {
            type = theType;
        }
        
        int wordSize;
        
        if (type.equals(Type.Bits64) || type.equals(Type.Bits64Signed)) {
            is64Bit = true;
            wordSize = 64;
        } else {
            is64Bit = false;
            wordSize = 32;
        }
        
        nPixels = theWidth * theHeight;
        
        width = theWidth;
        
        height = theHeight;
        
        if (type.equals(Type.Bits32FullRangeInt)) {
            maxAllowed = Integer.MAX_VALUE;
        } else {
            maxAllowed = 255;
        }
        
        if (type.equals(Type.Bits32) || type.equals(Type.Bits64)) {
            
            // 255 is a datum and that's 8 bits
            datumNBits = 8;
            
            minAllowed = 0;
            
        } else if (type.equals(Type.Bits32FullRangeInt)) {
            
            minAllowed = Integer.MIN_VALUE;
            
            datumNBits = 32;
            
        } else {//if (type.equals(Type.Bits32Signed) || type.equals(Type.Bits64Signed)) {
            
            // 512 is a datum and that's 9 bits
            datumNBits = 9;
            
            minAllowed = -256;
        }
        
        itemNDatum = wordSize/datumNBits;
                        
        len = (nPixels/itemNDatum) + 1;
        
        if (is64Bit) {

            aL = new long[len];

            a = null;

        } else {

            a = new int[len];

            aL = null;
        }
        
        if (minAllowed < 0) {
            fill(0);
        }
    }
    
    /**
     * convenience method to fill image with the value. 
     * @param value 
     */
    public void fill(int value) {
        
        if (value < minAllowed || (value > maxAllowed)) {
            throw new IllegalArgumentException(
                "value must be between " + minAllowed + " and " + maxAllowed + 
                 ", inclusive");
        }
       
        if (!type.equals(Type.Bits32FullRangeInt)) {
            value -= minAllowed;
        }
                
        if (is64Bit) {
            long mask = (1L << datumNBits) - 1L;
            long total = 0;
            for (int i = 0; i < itemNDatum; i++) {
                total += ((value & mask) << (long)(i * datumNBits));
            }
            Arrays.fill(aL, total);
        } else if (type.equals(Type.Bits32FullRangeInt)) {
            Arrays.fill(a, value);
        } else {
            int mask = (1 << datumNBits) - 1;
            int total = 0;
            for (int i = 0; i < itemNDatum; i++) {
                total += ((value & mask) << (i * datumNBits));
            }
            Arrays.fill(a, total);
        }
              
    }
    
    /**
     * convenience method to multiply entire image by factor.  Note that
     * if the result for a pixel is outside of minAllowed and maxAllowed,
     * it is clipped to that end of the allowed range.
     * @param factor a positive number that results in pixel values in range
     * 0 to 255;
     */
    public void multiply(float factor) {
        
        if (factor < 0) {
            throw new IllegalArgumentException(
                "factor has to be a positive number");
        }
        
        for (int i = 0; i < getNPixels(); ++i) {
            int v = Math.round(getValue(i)*factor);
            if (v < minAllowed) {
                v = minAllowed;
            } else if (v > maxAllowed) {
                v = maxAllowed;
            }
            setValue(i, v);
        }
        
    }
    
    /**
     * divide each pixel value by number
     * @param number 
     */
    public void divide(int number) {
        
        if (number < 0) {
            throw new IllegalArgumentException(
                "number has to be a positive number");
        }
        
        for (int i = 0; i < getNPixels(); ++i) {
            int v = getValue(i)/number;
            setValue(i, v);
        }
    }
    
    public int getMin() {
        int min;
    
        if (is64Bit) {
            min = (int)MiscMath.findMinForByteCompressed(aL, len, itemNDatum,
                datumNBits, minAllowed);
        } else if (type.equals(Type.Bits32FullRangeInt)) {
            min = MiscMath.findMin(a, len);
        } else {
            min = MiscMath.findMinForByteCompressed(a, len, itemNDatum,
                datumNBits, minAllowed);
        }
        return min;
    }
    
    public int getMax() {
        int max;
       
        if (is64Bit) {
            max = (int)MiscMath.findMaxForByteCompressed(aL, len, itemNDatum,
                datumNBits, minAllowed);
        } else if (type.equals(Type.Bits32FullRangeInt)) {
            max = MiscMath.findMax(a, len);
        } else {
            max = MiscMath.findMaxForByteCompressed(a, len, itemNDatum,
                datumNBits, minAllowed);
        }
        return max;
    }
    
    /**
     * scale pixel values up to 255, inclusive by applying a normalization
     * factor.
     */
    public void normalizeToMax255() {
        
        int max = getMax();
       
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
    
    /**
     * set value of pixel at location col, row
     * @param col
     * @param row
     * @param value 
     */
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
         
    private long set64BitValue(long rowValue, int setValue, int datumNumber) {
               
        setValue -= minAllowed;
        
        long mask = (1L << datumNBits) - 1L;
              
        long shift = (long)(datumNBits * datumNumber);
        
        long prevValue = (rowValue >> shift) & mask;
        
        long shifted = prevValue << shift;
        
        rowValue -= shifted;
        
        shifted = (setValue & mask) << shift;
                
        rowValue += shifted;
        
        assert(((rowValue >> shift) & mask) == setValue);
        
        return rowValue;
    }
    
    private int set32BitValue(int rowValue, int setValue, int datumNumber) {
          
        if (type.equals(Type.Bits32FullRangeInt)) {
            return setValue;
        }
        
        int mask = (1 << datumNBits) - 1;
        
        setValue -= minAllowed;
        
        int shift = datumNBits * datumNumber;
        
        int prevValue = (rowValue >> shift) & mask;
        
        int shifted = (prevValue & mask) << shift;
        
        rowValue -= shifted;
        
        shifted = (setValue & mask) << shift;
                
        rowValue += shifted;
                
        assert(((rowValue >> shift) & mask) == setValue);
        
        return rowValue;
    }
    
    /**
     * add the col and row for pixel index to the given array 
     * @param array
     * @param internalIndex 
     */
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
    
    /**
     * get the image row number for the pixel index
     * @param pixelIndex
     * @return 
     */
    public int getRow(int pixelIndex) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        // int idx = (row * width) + col;
        //     idx/w = row  + 0;
        int row = pixelIndex/width;
        
        return row;
    }
    
    /**
     * get the image col number for the pixel index
     * @param pixelIndex
     * @return 
     */
    public int getCol(int pixelIndex) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
   
        // int idx = (row * width) + col;
        //     idx/w = row  + 0;
        //int row = pixelIndex/width;
        
        //int col = pixelIndex - (row * width);
        
        return pixelIndex % width;
    }
    
    /**
     * set pixel at pixelIndex to value
     * @param pixelIndex
     * @param value 
     */
    public void setValue(int pixelIndex, int value) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "internalIndex is out of bounds:");
        }
        
        if (value < minAllowed || (value > maxAllowed)) {
            throw new IllegalArgumentException(
                "value must be between " + minAllowed + " and " + maxAllowed + 
                 ", inclusive.  value=" + value);
        }
        
        int elementIdx = pixelIndex/itemNDatum;
        
        int datumNumber = pixelIndex - (elementIdx*itemNDatum);
                
        if (is64Bit) {
            aL[elementIdx] = set64BitValue(aL[elementIdx], value, datumNumber);
        } else {
            a[elementIdx] = set32BitValue(a[elementIdx], value, datumNumber);
        }
    }
    
    private int get32BitValue(int rowValue, int datumNumber) {
        
        if (type.equals(Type.Bits32FullRangeInt)) {
            return rowValue;
        }
        
        int mask = (1 << datumNBits) - 1;
        
        int v = (rowValue >> (datumNumber * datumNBits)) & mask;
        
        if (!type.equals(Type.Bits32FullRangeInt)) {
            v += minAllowed;
        }
        
        return v;
    }
    
    private int get64BitValue(long rowValue, int datumNumber) {
        
        long mask = (1L << datumNBits) - 1L;
        
        int v = (int)((rowValue >> (datumNumber * datumNBits)) & mask);
        
        if (!type.equals(Type.Bits32FullRangeInt)) {
            v += minAllowed;
        }
        
        return v;
    }
    
    /**
     * get value at pixelIndex
     * @param pixelIndex
     * @return 
     */
    public int getValue(int pixelIndex) {
        
        if ((pixelIndex < 0) || (pixelIndex > (nPixels - 1))) {
            throw new IllegalArgumentException(
                "pixelIndex is out of bounds:");
        }
        
        int elementIdx = pixelIndex/itemNDatum;
        
        int datumNumber = pixelIndex - (elementIdx*itemNDatum);
        
        if (is64Bit) {
            return get64BitValue(aL[elementIdx], datumNumber);
        } else {
            return get32BitValue(a[elementIdx], datumNumber);
        }
    }
    
    /**
     * get the pixel index for the given col, row
     * @param col
     * @param row
     * @return 
     */
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
        
    /**
     * get the pixel value at image location col, row
     * @param col
     * @param row
     * @return 
     */
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
    
    /**
     * get the pixel value at image location col, row
     * @param p pairint of x and y coordinates
     * @return 
     */
    public int getValue(PairInt p) {
        
        int col = p.getX();
        int row = p.getY();
        
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
        
        for (int i = 0; i < nPixels; ++i) {
            t[i] = getValue(i);
        }
        
        return t;
    }
    
    public int max() {
        
        Integer maxV = Integer.MIN_VALUE;
        
        for (int i = 0; i < nPixels; ++i) {
            int v = getValue(i);
            if (v > maxV) {
                maxV = v;
            }
        }
        
        return maxV;
    }
    
    public int min() {
        
        Integer minV = Integer.MAX_VALUE;
        
        for (int i = 0; i < nPixels; ++i) {
            int v = getValue(i);
            if (v < minV) {
                minV = v;
            }
        }
        
        return minV;
    }
    
    public GreyscaleImage copyImage() {
       
        GreyscaleImage img2 = createWithDimensions();
                
        if (is64Bit) {
            System.arraycopy(aL, 0, img2.aL, 0, len);
        } else {
            System.arraycopy(a, 0, img2.a, 0, len);
        }
        
        return img2;
    }
    
    public GreyscaleImage copyToSignedImage() {

        GreyscaleImage img2;
        
        if (type.equals(Type.Bits32) || type.equals(Type.Bits64)) {
            if (type.equals(Type.Bits32)) {
                img2 = new GreyscaleImage(width, height, GreyscaleImage.Type.Bits32Signed);
            } else {
                img2 = new GreyscaleImage(width, height, GreyscaleImage.Type.Bits64Signed);
            }
            for (int i = 0; i < nPixels; ++i) {
                int v = getValue(i);
                img2.setValue(i, v);
            }
        } else {
            img2 = new GreyscaleImage(width, height, type);
            if (is64Bit) {
                System.arraycopy(aL, 0, img2.aL, 0, len);
            } else {
                System.arraycopy(a, 0, img2.a, 0, len);
            }
        }
                
        return img2;
    }
    
    public GreyscaleImage copyToFullRangeIntImage() {

        GreyscaleImage img2;
        
        if (type.equals(Type.Bits32FullRangeInt)) {
            img2 = new GreyscaleImage(width, height, type);
            System.arraycopy(a, 0, img2.a, 0, len);
        } else {
            img2 = new GreyscaleImage(width, height, Type.Bits32FullRangeInt);
            for (int i = 0; i < nPixels; ++i) {
                int v = getValue(i);
                img2.setValue(i, v);
            }
        }
           
        return img2;
    }
    
    public GreyscaleImage createWithDimensions() {
       
        GreyscaleImage img2 = new GreyscaleImage(width, height, type);
                
        img2.xRelativeOffset = xRelativeOffset;
        img2.yRelativeOffset = yRelativeOffset;
        
        if (img2.minAllowed < 0) {
            img2.fill(0);
        }
        
        return img2;
    }
    
    /**
     * create new image of same Type and having same xOffset and yOffset.
     * @param width2
     * @param height2
     * @return 
     */
    public GreyscaleImage createWithDimensions(int width2, int height2) {
       
        GreyscaleImage img2 = new GreyscaleImage(width2, height2, type);
                
        img2.xRelativeOffset = xRelativeOffset;
        img2.yRelativeOffset = yRelativeOffset;
        
        if (img2.minAllowed < 0) {
            img2.fill(0);
        }
        
        return img2;
    }
    
    public GreyscaleImage createSignedWithDimensions() {
        
        GreyscaleImage img2;
        
        if (type.equals(Type.Bits32) || type.equals(Type.Bits64)) {
            if (type.equals(Type.Bits32)) {
                img2 = new GreyscaleImage(width, height, GreyscaleImage.Type.Bits32Signed);
            } else {
                img2 = new GreyscaleImage(width, height, GreyscaleImage.Type.Bits64Signed);
            }
        } else {
            img2 = new GreyscaleImage(width, height, type);
        }
        
        if (img2.minAllowed < 0) {
            img2.fill(0);
        }
        
        return img2;
    }
    
    public GreyscaleImage createFullRangeIntWithDimensions() {
        
        GreyscaleImage img2 = new GreyscaleImage(width, height, 
            GreyscaleImage.Type.Bits32FullRangeInt);
        
        return img2;
    }
    
    /**
     * create a sub image centered on (xCenter, yCenter) of width subWidth
     * and height subHeight
     * @param xCenter
     * @param yCenter
     * @param subWidth
     * @param subHeight
     * @return 
     */
    public GreyscaleImage subImage(int xCenter, int yCenter, int subWidth, 
        int subHeight) {
       
        GreyscaleImage img2 = new GreyscaleImage(subWidth, subHeight, type);
                
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
    
    /**
     * create a new sub-image for the given start and stop, inclusive, pixels.
     * 
     * @param startX
     * @param stopX
     * @param startY
     * @param stopY
     * @return 
     */
    public GreyscaleImage subImage2(int startX, int stopX, 
        int startY, int stopY) {
        
        if (startX < 0) {
            throw new IllegalArgumentException("startX must be positive non-zero");
        }
        if (startY < 0) {
            throw new IllegalArgumentException("startY must be positive non-zero");
        }
        if (stopX >= this.width) {
            throw new IllegalArgumentException("stopX must be < imageWidth");
        }
        if (stopY >= this.height) {
            throw new IllegalArgumentException("stopY must be < imageHeight");
        }
        
        int w2 = stopX - startX + 1;
        int h2 = stopY - startY + 1;
        
        GreyscaleImage img2 = new GreyscaleImage(w2, h2, type);
                
        int col2 = 0;
        for (int col = startX; col <= stopX; ++col) {
            
            int row2 = 0;
            for (int row = startY; row <= stopY; ++row) {
                
                int v = getValue(col, row);
                
                img2.setValue(col2, row2, v);
                
                row2++;
            }
            
            col2++;
        }
            
        return img2;
    }
    
    /**
     * note, can only use this if image has values between 0 and 255, inclusive
     * @return 
     */
    public ImageExt copyToColorGreyscaleExt() {
        
        /*
        TODO:
        consider changing this conversion to use
        https://en.wikipedia.org/wiki/Grayscale
            grey = red * .2126 + green * .07152 + blue * .0722;
        */
        ImageExt img2 = new ImageExt(width, height, !is64Bit);
        
        for (int i = 0; i < nPixels; ++i) {
            int v = getValue(i);
            img2.setRGB(i, v, v, v);
        }
        
        return img2;
    }
    
    /**
     * note, can only use this if image has values between 0 and 255, inclusive
     * @return 
     */
    public Image copyToColorGreyscale() {
        
        /*
        TODO:
        consider changing this conversion to use
        https://en.wikipedia.org/wiki/Grayscale
            grey = red * .2126 + green * .07152 + blue * .0722;
        */
        Image img2 = new Image(width, height, !is64Bit);
        
        for (int i = 0; i < nPixels; ++i) {
            int v = getValue(i);
            img2.setRGB(i, v, v, v);
        }
        
        return img2;
    }
    
    /**
     * copy to an instance of Image and clip values below 0 to 0 and clip
     * values above 255 to 255.
     * @return 
     */
    public Image copyToColorGreyscaleWithClipping() {
        
        Image img2 = new Image(width, height, !is64Bit);
        
        for (int i = 0; i < nPixels; ++i) {
            int v = getValue(i);
            if (v < 0) { 
                v = 0;
            } else if (v > 255) {
                v = 255;
            }
            img2.setRGB(i, v, v, v);
        }

        return img2;
    }
    
    public void resetTo(final GreyscaleImage copyThis) {
        
        type = copyThis.getType().copy(copyThis.type);
        is64Bit = copyThis.is64Bit;
        width = copyThis.width;
        height = copyThis.height;
        minAllowed = copyThis.minAllowed;
        maxAllowed = copyThis.maxAllowed;
        nPixels = copyThis.getNPixels();
        itemNDatum = copyThis.itemNDatum;
        datumNBits = copyThis.datumNBits;
        len = copyThis.len;
        xRelativeOffset = copyThis.xRelativeOffset;
        yRelativeOffset = copyThis.yRelativeOffset;
       
        if (a == null && copyThis.a != null) {
            a = new int[copyThis.a.length];
        } else if ((a != null) && copyThis.a == null) {
            a = null;
        } else if ((a != null) && (a.length != copyThis.a.length)) {
            a = new int[copyThis.a.length];
        }
        
        if (aL == null && copyThis.aL != null) {
            aL = new long[copyThis.aL.length];
        } else if ((aL != null) && copyThis.aL == null) {
            aL = null;
        } else if ((aL != null) && (aL.length != copyThis.aL.length)) {
            aL = new long[copyThis.aL.length];
        }
        
        if (a != null) {
            System.arraycopy(copyThis.a, 0, a, 0, copyThis.len);
        } else {
            System.arraycopy(copyThis.aL, 0, aL, 0, copyThis.len);
        }
        
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
    
    /**
     * subtract img from this image and return result in a new image of same
     * type as this image.
     * @param img
     * @return 
     */
    public GreyscaleImage subtract(GreyscaleImage img) {
        
        if (img.getWidth() != width || img.getHeight() != height) {
            throw new IllegalArgumentException("images must be same dimensions");
        }
                
        GreyscaleImage output = createWithDimensions();
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v0 = getValue(i);
            int v1 = img.getValue(i);
            output.setValue(i, (v0 - v1));
        }
        
        return output;
    }
    
    /**
     * add img to this image and return result in a new image of same
     * type as this image.
     * NOTE if the result for a pixel is outside of minAllowed and maxAllowed,
     * it is clipped to that end of the allowed range.
     * @param img
     * @return 
     */
    public GreyscaleImage add(GreyscaleImage img) {
        
        if (img.getWidth() != width || img.getHeight() != height) {
            throw new IllegalArgumentException("images must be same dimensions");
        }
        
        GreyscaleImage output = createWithDimensions();
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v0 = getValue(i);
            int v1 = img.getValue(i);
            v0 += v1;
            if (v0 < minAllowed) {
                v0 = minAllowed;
            } else if (v0 > maxAllowed) {
                v0 = maxAllowed;
            }
            output.setValue(i, v0);
        }
        
        return output;
    }
    
    public int[][] extractAsDoubleArray() {
        
        int[][] out = new int[width][];
        for (int i = 0; i < width; ++i) {
            out[i] = new int[height];
            for (int j = 0; j < height; ++j) {
                out[i][j] = getValue(i, j);
            }
        }
        
        return out;
    }

    public int[] getValues() {
        int[] v = new int[nPixels];
        for (int i = 0; i < nPixels; ++i) {
            v[i] = getValue(i);
        }
        return v;
    }
    
    /**
     * @return the minAllowed
     */
    public int getMinAllowed() {
        return minAllowed;
    }

    /**
     * @return the maxAllowed
     */
    public int getMaxAllowed() {
        return maxAllowed;
    }

    /**
     * @return the itemNDatum
     */
    public int getItemNDatum() {
        return itemNDatum;
    }

    /**
     * @return the datumNBits
     */
    public int getDatumNBits() {
        return datumNBits;
    }
    
    public Type getType() {
        return type;
    }

    /**
     * assuming that max value is 255, return an array normalized by 255
     * so that all values are in range 0 to 1.0.
     * @return 
     */
    public float[][] toNormalizedRowMajorArray() {
        
        float[][] a = new float[height][];
        for (int i = 0; i < height; ++i) {
            a[i] = new float[width];
            for (int j = 0; j < width; ++j) {
                a[i][j] = ((float)getValue(j, i))/255.f;
            }
        }
        
        return a;
    }

    /**
     * export the image in format row major, that is img[i][j] is img[y][x]
     * where x and y are the pixel coordinates.
     * @return
     */
    public double[][] exportRowMajor() {
        double[][] out = new double[height][width];
        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                out[y][x] = getValue(x, y);
            }
        }
        return out;
    }


    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        for (int row = 0; row < height; ++row) {
            sb.append("row=").append(Integer.valueOf(row)).append(":");
            for (int col = 0; col < width; ++col) {
                sb.append(getValue(col, row)).append(" ");
            }
            sb.append("\n");
        }
        
        return sb.toString();
    }
    
}
