package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;

/**
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
    }
    
    public IntensityDescriptor extractIntensity(int xCenter, int yCenter, 
        float rotation) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    public GradientDescriptor extractGradient(int xCenter, int yCenter, 
        float rotation) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private IntensityDescriptor extractNormalizedIntensity(int xCenter, 
        int yCenter, float rotation) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private void extractNormalizedIntensity(int xCenter, 
        int yCenter, float rotation, float[] output) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    private void normalize(float[] output) {
        throw new UnsupportedOperationException("not yet implemented");
    }
   
}
