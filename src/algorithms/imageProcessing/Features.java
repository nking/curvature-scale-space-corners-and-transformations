package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;

/**
 * NOT READY FOR USE.  NOT TESTED YET.
 * class to help extract and reuse feature descriptors.
 * @author nichole
 */
public class Features {
    
    protected final GreyscaleImage gsImg;
    
    protected final Image clrImg;
    
    /**
     * the half width of a block.  For example, to extract the 25 pixels
     * centered on (xc, yc), bHalfW would be '2'
     */
    protected final int bHalfW;
    
    protected final boolean useNormalizedIntensities;
    
    protected final float[][] xyOffsets;
    
    /**
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
    *     extracted descriptor
    */
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> intensityBlocks =
        new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    /**
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
    *     extracted descriptor
    */
    protected Map<PairInt, Map<Integer, GradientDescriptor>> gradientBlocks = 
        new HashMap<PairInt, Map<Integer, GradientDescriptor>>();
    
    /**
     * 
     * @param image
     * @param blockHalfWidths the half width of a block.  For example, to 
     * extract the 25 pixels centered on (xc, yc), bHalfW would be '2'
     * @param useNormalizedIntensities if true, the intensity descriptors
     * for each block are normalized by mean and standard deviation
     * I_normalized(pixel) = (I(pixel)-I_mean(block))/I_stdev(block)).
     */
    public Features(GreyscaleImage image, int blockHalfWidths, 
        boolean useNormalizedIntensities) {
        this.gsImg = image;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
    }
    
    /**
     * 
     * @param image
     * @param blockHalfWidths the half width of a block.  For example, to 
     * extract the 25 pixels centered on (xc, yc), bHalfW would be '2'
     * @param useNormalizedIntensities if true, the intensity descriptors
     * for each block are normalized by mean and standard deviation
     * I_normalized(pixel) = (I(pixel)-I_mean(block))/I_stdev(block)).
     */
    public Features(Image image, int blockHalfWidths, boolean 
        useNormalizedIntensities) {
        this.gsImg = null;
        this.clrImg = image;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
    }
    
    /**
     * extract the intensity from the image for the given block center and
     * return it in a descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at (xCenter, yCenter)
     * @return 
     */
    public IntensityDescriptor extractIntensity(int xCenter, int yCenter, 
        int rotation) {
        
        checkBounds(xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotKey = Integer.valueOf(rotation);
        
        Map<Integer, IntensityDescriptor> descriptors = intensityBlocks.get(p);
        
        IntensityDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, IntensityDescriptor>();
            intensityBlocks.put(p, descriptors);
        }
        
        descriptor = extractIntensityForBlock(xCenter, yCenter, rotation);
        
        assert(descriptor != null);
        
        descriptors.put(rotKey, descriptor);
        
        return descriptor;
    }
    
    public GradientDescriptor extractGradient(int xCenter, int yCenter, 
        float rotation) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private IntensityDescriptor extractIntensityForBlock(int xCenter, 
        int yCenter, int rotation) {
        
        Transformer transformer = new Transformer();
        
        float[][] frameOffsets = transformer.transformXY(rotation, xyOffsets);
        
        IntensityDescriptor descriptor;
        
        if (gsImg != null) {
            descriptor = extractGsIntensityForBlock(xCenter, yCenter, 
                frameOffsets);
        } else {
            descriptor = extractClrIntensityForBlock(xCenter, yCenter, 
                frameOffsets);
        }
        
        if (useNormalizedIntensities) {
            descriptor.applyNormalization();
        }
        
        return descriptor;
    }

    private void checkBounds(int xCenter, int yCenter) {
        
        if (xCenter < 0 || yCenter < 0) {
            throw new IllegalArgumentException(
            "xCenter and yCenter must be > -1");
        }
        
        if (gsImg != null) {
            if (xCenter > (gsImg.getWidth() - 1)) {
                throw new IllegalArgumentException(
                "xCenter must be less than image width");
            }
            if (yCenter > (gsImg.getHeight() - 1)) {
                throw new IllegalArgumentException(
                "yCenter must be less than image height");
            }
        } else {
            if (xCenter > (clrImg.getWidth() - 1)) {
                throw new IllegalArgumentException(
                "xCenter must be less than image width");
            }
            if (yCenter > (clrImg.getHeight() - 1)) {
                throw new IllegalArgumentException(
                "yCenter must be less than image height");
            }
        }
    }

    /**
     * extract the intensity from the image and place in the descriptor.
     * Note that if the transformed pixel is out of bounds of the image,
     * a sentinel is the value for that location in the descriptor.
     * @param xCenter
     * @param yCenter
     * @param offsets
     * @return 
     */
    private IntensityDescriptor extractGsIntensityForBlock(int xCenter, 
        int yCenter, float[][] offsets) {
        
        int[] output = new int[offsets.length];
        
        int sentinel = GsIntensityDescriptor.sentinel;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int count = 0;
        for (int i = 0; i < offsets.length; ++i) {
            
            float x1P = xCenter + offsets[i][0];
            
            float y1P = yCenter + offsets[i][1];
            
            if ((x1P < 0) || (Math.ceil(x1P) > (gsImg.getWidth() - 1)) || 
                (y1P < 0) || (Math.ceil(y1P) > (gsImg.getHeight() - 1))) {
                
                output[count] = sentinel;
                
            } else {
                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                double v = imageProcessor.biLinearInterpolation(gsImg, x1P, y1P);
                
                output[count] = (int)Math.round(v);
            }
        }
        
        IntensityDescriptor desc = new GsIntensityDescriptor(output);
        
        return desc;
    }
   
    /**
     * extract the intensity from the image and place in the descriptor.
     * Note that if the transformed pixel is out of bounds of the image,
     * a sentinel is the value for that location in the descriptor.
     * @param xCenter
     * @param yCenter
     * @param offsets
     * @return 
     */
    private IntensityDescriptor extractClrIntensityForBlock(int xCenter, 
        int yCenter, float[][] offsets) {
        
        int[] outputR = new int[offsets.length];
        int[] outputG = new int[offsets.length];
        int[] outputB = new int[offsets.length];
        
        int sentinel = ClrIntensityDescriptor.sentinel;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int count = 0;
        for (int i = 0; i < offsets.length; ++i) {
            
            float x1P = xCenter + offsets[i][0];
            
            float y1P = yCenter + offsets[i][1];
            
            if ((x1P < 0) || (Math.ceil(x1P) > (clrImg.getWidth() - 1)) || 
                (y1P < 0) || (Math.ceil(y1P) > (clrImg.getHeight() - 1))) {
                
                outputR[count] = sentinel;
                outputG[count] = sentinel;
                outputB[count] = sentinel;
                
            } else {
                
                double[] v = imageProcessor.biLinearInterpolation(clrImg, x1P, y1P);
                
                outputR[count] = (int)Math.round(v[0]);
                outputG[count] = (int)Math.round(v[1]);
                outputB[count] = (int)Math.round(v[2]);
            }
            
        }
        
        IntensityDescriptor desc = new ClrIntensityDescriptor(outputR, outputG,
            outputB);
        
        return desc;
    }

}
