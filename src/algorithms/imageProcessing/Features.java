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
    
    // TODO: consider adding the color gradients when the use of gradients
    //     is tested and working
    protected final GreyscaleImage thetaImg;
    
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
    public Features(GreyscaleImage image, GreyscaleImage gradientThetaImg,
        int blockHalfWidths, 
        boolean useNormalizedIntensities) {
        this.gsImg = image;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.thetaImg = gradientThetaImg;
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
    public Features(Image image, GreyscaleImage gradientThetaImg,
        int blockHalfWidths, boolean useNormalizedIntensities) {
        this.gsImg = null;
        this.clrImg = image;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.thetaImg = gradientThetaImg;
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

    private void checkBounds(int x, int y) {
        
        if (isWithinXBounds(y)) {
            throw new IllegalArgumentException("y is out of bounds of image");
        }
        if (isWithinYBounds(y)) {
            throw new IllegalArgumentException("y is out of bounds of image");
        }
    }
    
    public boolean isWithinXBounds(int x) {
        
        if (x < 0) {
            return false;
        }
        
        if (gsImg != null) {
            if (x > (gsImg.getWidth() - 1)) {
                return false;
            }
        } else {
            if (x > (clrImg.getWidth() - 1)) {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean isWithinYBounds(int y) {
        
        if (y < 0) {
            return false;
        }
        
        if (gsImg != null) {
            if (y > (gsImg.getHeight() - 1)) {
                return false;
            }
        } else {
            if (y > (clrImg.getHeight() - 1)) {
                return false;
            }
        }
        
        return true;
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

    /**
     * calculate the intensity based statistic (SSD) between the two descriptors
     * with the caveat that the 2nd descriptor is the one used to calculate
     * the error (so make sure that pattern is consistently used by invoker).
     * 
     * @param desc1
     * @param x1
     * @param y1
     * @param desc2
     * @param x2
     * @param y2
     * @return 
     */
    public static FeatureComparisonStat calculateIntensityStat(
        IntensityDescriptor desc1, final int x1, final int y1,
        IntensityDescriptor desc2, final int x2, final int y2) {
    
        if (desc1 == null) {
            throw new IllegalArgumentException("desc1 cannot be null");
        }
        if (desc2 == null) {
            throw new IllegalArgumentException("desc2 cannot be null");
        }
        
        float err2Sq = desc2.sumSquaredError();
        
        float ssd = desc1.calculateSSD(desc2);
        
        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumSqDiff(ssd);
        stat.setImg2PointErr(err2Sq);
        
        return stat;
    }

}
