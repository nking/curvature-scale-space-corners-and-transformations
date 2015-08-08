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
    
    protected final GreyscaleImage gradientImg;
    
    //TODO: add this to a setter:
    protected final Image clrGradientImg = null;
    
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
    public Features(GreyscaleImage image, GreyscaleImage theGradientImg,
        int blockHalfWidths, 
        boolean useNormalizedIntensities) {
        this.gsImg = image;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.gradientImg = theGradientImg;
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
    public Features(Image image, GreyscaleImage theGradientImg,
        int blockHalfWidths, boolean useNormalizedIntensities) {
        this.gsImg = null;
        this.clrImg = image;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.gradientImg = theGradientImg;
    }
    
    /**
     * extract the intensity from the image for the given block center and
     * return it in a descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at (xCenter, yCenter)
     * @return 
     */
    public IntensityDescriptor extractIntensity(final int xCenter, 
        final int yCenter, final int rotation) {
        
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
    
    /**
     * To reduce the importance of exact rotation alignment, and to allow more
     * tolerance for projection, this is not use single pixel
     * data in a block from the constructor's half width, instead, it is
     * binning the data by 2 in x and y around xCenter and yCenter to 
     * make 16 such binned regions.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at (xCenter, yCenter)
     * @return 
     */
    public GradientDescriptor extractGradient(int xCenter, int yCenter, 
        int rotation) {
                
        checkBounds(xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotKey = Integer.valueOf(rotation);
        
        Map<Integer, GradientDescriptor> descriptors = gradientBlocks.get(p);
        
        GradientDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, GradientDescriptor>();
            gradientBlocks.put(p, descriptors);
        }
        
        if (gradientImg != null) {
            descriptor = extractGsGradientForCells(xCenter, yCenter);
        } else {
            descriptor = extractClrGradientForCells(xCenter, yCenter);
        }
                
        assert(descriptor != null);
        
        descriptors.put(rotKey, descriptor);
        
        return descriptor;
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
        
        if (!isWithinXBounds(x)) {
            throw new IllegalArgumentException("x is out of bounds of image");
        }
        if (!isWithinYBounds(y)) {
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
            
            count++;
        }
        
        IntensityDescriptor desc = new GsIntensityDescriptor(output);
        
        return desc;
    }
    
    /**
     * extract the gradient from the image in 2X2 cells surrounding 
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param offsets
     * @return 
     */
    private GradientDescriptor extractGsGradientForCells(int xCenter, 
        int yCenter) {
        
        /*
          3 [-][-][.][.][-][-][.][.]
          2 [-][-][.][.][-][-][.][.]
          1 [.][.][-][-][.][.][-][-]
          0 [.][.][-][-] @ [.][-][-]         
         -1 [-][-][.][.][-][-][.][.]
         -2 [-][-][.][.][-][-][.][.]
         -3 [.][.][-][-][.][.][-][-]
         -4 [.][.][-][-][.][.][-][-]
            -4 -3 -2 -1  0  1  2  3  
        
        */
        int[] output = new int[16];
        
        int sentinel = GsGradientDescriptor.sentinel;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int count = 0;
        for (int dx = -4; dx < 4; dx+=2) {
            
            float x1P = xCenter + dx;
            float x2P = x1P + 1;
            
            if ((x1P < 0) || (Math.ceil(x1P) > (gradientImg.getWidth() - 1)) ||
                (x2P < 0) || (Math.ceil(x2P) > (gradientImg.getWidth() - 1))
                ) {
                
                output[count] = sentinel;
                count++;
                continue;
            }
            
            for (int dy = -4; dy < 4; dy+=2) {
            
                float y1P = yCenter + dy;
                float y2P = y1P + 1;
            
                if ((y1P < 0) || (Math.ceil(y1P) > (gradientImg.getWidth() - 1)) ||
                    (y2P < 0) || (Math.ceil(y2P) > (gradientImg.getWidth() - 1))) {
                
                    output[count] = sentinel;
                    count++;
                    continue;
                }
                                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                double v1 = imageProcessor.biLinearInterpolation(gradientImg, x1P, y1P);
                double v2 = imageProcessor.biLinearInterpolation(gradientImg, x1P, y2P);
                double v3 = imageProcessor.biLinearInterpolation(gradientImg, x2P, y1P);
                double v4 = imageProcessor.biLinearInterpolation(gradientImg, x2P, y2P);
                
                output[count] = (int)Math.round((v1 + v2 + v3 + v4)/4.);
            }
            
            count++;
        }
        
        GradientDescriptor desc = new GsGradientDescriptor(output);
        
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
            
            count++;
        }
        
        IntensityDescriptor desc = new ClrIntensityDescriptor(outputR, outputG,
            outputB);
        
        return desc;
    }
    
    /**
     * extract the gradient from the image in 2X2 cells surrounding 
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param offsets
     * @return 
     */
    private GradientDescriptor extractClrGradientForCells(int xCenter, 
        int yCenter) {
        
        throw new UnsupportedOperationException("color gradient not yet implemented");
        
        /*
          3 [-][-][.][.][-][-][.][.]
          2 [-][-][.][.][-][-][.][.]
          1 [.][.][-][-][.][.][-][-]
          0 [.][.][-][-] @ [.][-][-]         
         -1 [-][-][.][.][-][-][.][.]
         -2 [-][-][.][.][-][-][.][.]
         -3 [.][.][-][-][.][.][-][-]
         -4 [.][.][-][-][.][.][-][-]
            -4 -3 -2 -1  0  1  2  3  
        */
        /*
        int[] red = new int[16];
        int[] green = new int[16];
        int[] blue = new int[16];
        
        int sentinel = GsGradientDescriptor.sentinel;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int count = 0;
        for (int dx = -4; dx < 4; dx+=2) {
            
            float x1P = xCenter + dx;
            float x2P = x1P + 1;
            
            if ((x1P < 0) || (Math.ceil(x1P) > (clrGradientImg.getWidth() - 1)) ||
                (x2P < 0) || (Math.ceil(x2P) > (clrGradientImg.getWidth() - 1))
                ) {
                
                red[count] = sentinel;
                green[count] = sentinel;
                blue[count] = sentinel;
                
                count++;
                
                continue;
            }
            
            for (int dy = -4; dy < 4; dy+=2) {
            
                float y1P = yCenter + dy;
                float y2P = y1P + 1;
            
                if ((y1P < 0) || (Math.ceil(y1P) > (clrGradientImg.getWidth() - 1)) ||
                    (y2P < 0) || (Math.ceil(y2P) > (clrGradientImg.getWidth() - 1))) {
                
                    red[count] = sentinel;
                    green[count] = sentinel;
                    blue[count] = sentinel;
                
                    count++;
                    continue;
                }
                                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                double[] v1 = imageProcessor.biLinearInterpolation(clrGradientImg, x1P, y1P);
                double[] v2 = imageProcessor.biLinearInterpolation(clrGradientImg, x1P, y2P);
                double[] v3 = imageProcessor.biLinearInterpolation(clrGradientImg, x2P, y1P);
                double[] v4 = imageProcessor.biLinearInterpolation(clrGradientImg, x2P, y2P);
                
                red[count] = (int)Math.round((v1[0] + v2[0] + v3[0] + v4[0])/4.);
                green[count] = (int)Math.round((v1[1] + v2[1] + v3[1] + v4[1])/4.);
                blue[count] = (int)Math.round((v1[2] + v2[2] + v3[2] + v4[2])/4.);
            }
            
            count++;
        }
        
        GradientDescriptor desc = new ClrGradientDescriptor(red, green, blue);
        
        return desc;
        */
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
        IDescriptor desc1, final int x1, final int y1,
        IDescriptor desc2, final int x2, final int y2) {
    
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
