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
    
    /**
     * either gsImg or clrImg will be set by the constructor and they
     * are the grey scale image or the color image, respectively.
     */
    protected final GreyscaleImage gsImg;
    
    /**
     * either gsImg or clrImg will be set by the constructor and they
     * are the grey scale image or the color image, respectively.
     */
    protected final Image clrImg;
    
    /**
     * the gradient image of gsImg or clrImg.
     */
    protected final GreyscaleImage gradientImg;
    
    /**
     * the gradient orientation image
     */
    protected final GreyscaleImage thetaImg;
    
    //TODO: add this to a setter:
    protected final Image clrGradientImg = null;
    
    /**
     * the half width of a block.  For example, to extract the 25 pixels
     * centered on (xc, yc), bHalfW would be '2'
     */
    protected final int bHalfW;
    
    protected final float[][] xyOffsets;
    
    protected final boolean useNormalizedIntensities;
    
    protected final boolean useBinnedCellGradients = false;
        
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
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
    *     extracted descriptor
    */
    protected Map<PairInt, Map<Integer, ThetaDescriptor>> thetaBlocks = 
        new HashMap<PairInt, Map<Integer, ThetaDescriptor>>();
    
    /**
     * 
     * @param image
     * @param theGradientImg gradient image of the image region (usually
     * from the process of creating corners).
     * @param theThetaImg the gradient orientation angle image
     * @param blockHalfWidths the half width of a block.  For example, to 
     * extract the 25 pixels centered on (xc, yc), bHalfW would be '2'
     * @param useNormalizedIntensities normalize the intensities extracted
     * from image if true
     */
    public Features(GreyscaleImage image, GreyscaleImage theGradientImg,
        GreyscaleImage theThetaImg,
        int blockHalfWidths, boolean useNormalizedIntensities) {
        this.gsImg = image;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.gradientImg = theGradientImg;
        this.thetaImg = theThetaImg;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
    }
    
    /**
     * 
     * @param image
     * @param theGradientImg gradient image of the image region (usually
     * from the process of creating corners).
     * @param theThetaImg the gradient orientation angle image
     * @param blockHalfWidths the half width of a block.  For example, to 
     * extract the 25 pixels centered on (xc, yc), bHalfW would be '2'
     * @param useNormalizedIntensities normalize the intensities extracted
     * from image if true
     */
    public Features(Image image, GreyscaleImage theGradientImg,
        GreyscaleImage theThetaImg,
        int blockHalfWidths, boolean useNormalizedIntensities) {
        this.gsImg = null;
        this.clrImg = image;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.gradientImg = theGradientImg;
        this.thetaImg = theThetaImg;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
    }
    
    /**
     * extract the intensity from the image for the given block center and
     * return it in a descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at 
     * (xCenter, yCenter)
     * @return 
     */
    public IntensityDescriptor extractIntensity(final int xCenter, 
        final int yCenter, final int rotation) {
        
        checkBounds(xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        Map<Integer, IntensityDescriptor> descriptors = intensityBlocks.get(p);
        
        IntensityDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, IntensityDescriptor>();
            intensityBlocks.put(p, descriptors);
        }
        
        descriptor = extractIntensityForBlock(xCenter, yCenter, rotation);
        
        assert(descriptor != null);
        
        descriptors.put(rotationKey, descriptor);
        
        return descriptor;
    }
    
    /**
     * 
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at 
     * (xCenter, yCenter)
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
        
        if (useBinnedCellGradients) {
            descriptor = extractGsGradientForCells(xCenter, yCenter, rotation);
        } else {
            descriptor = extractGsGradientForBlock(xCenter, yCenter, rotation);
        }
        
        assert(descriptor != null);
        
        descriptors.put(rotKey, descriptor);
        
        return descriptor;
    }
    
    /**
     * extract theta orientation from the image for (xCenter, yCenter) and
     * return it in a descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at 
     * (xCenter, yCenter)
     * @return 
     */
    public ThetaDescriptor extractTheta(final int xCenter, final int yCenter, 
        final int rotation) {
       
        checkBounds(xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        Map<Integer, ThetaDescriptor> descriptors = thetaBlocks.get(p);
        
        ThetaDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, ThetaDescriptor>();
            thetaBlocks.put(p, descriptors);
        }
         
        /*
        this will be one of 3 types of theta descriptors and the
        choice will be set in constructor
        */
        
        descriptor = extractThetaForBlock(xCenter, yCenter, rotation);
        
        assert(descriptor != null);
        
        descriptors.put(rotationKey, descriptor);
        
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
    
    private GradientDescriptor extractGsGradientForBlock(int xCenter, 
        int yCenter, int rotation) {
        
        Transformer transformer = new Transformer();

        float[][] frameOffsets = transformer.transformXY(rotation, xyOffsets);
        
        GradientDescriptor descriptor = null;
        
        if (gradientImg != null) {
            descriptor = extractGsGradientForBlock(xCenter, yCenter, 
                frameOffsets);
        } /*else {
            descriptor = extractClrGradientForBlock(xCenter, yCenter, 
                frameOffsets);
        }*/
        
        return descriptor;
    }
    
    private ThetaDescriptor extractThetaForBlock(int xCenter, 
        int yCenter, int rotation) {
        
        Transformer transformer = new Transformer();

        float[][] frameOffsets = transformer.transformXY(rotation, xyOffsets);
        
        ThetaDescriptor descriptor = null;
        
        if (thetaImg != null) {
            descriptor = extractThetaForBlock(xCenter, yCenter, 
                frameOffsets);
        } /*else {
            descriptor = extractClrGradientForBlock(xCenter, yCenter, 
                frameOffsets);
        }*/
        
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
     * @return 
     */
    private IntensityDescriptor extractGsIntensityForBlock(int xCenter, 
        int yCenter, float[][] offsets) {
                
        float[] output = new float[offsets.length];
        
        float sentinel = GsIntensityDescriptor.sentinel;
                
        int count = 0;
        
        for (int i = 0; i < offsets.length; ++i) {
            
            float x1P = xCenter + offsets[i][0];
            
            float y1P = yCenter + offsets[i][1];
            
            if ((x1P < 0) || (Math.ceil(x1P) > (gsImg.getWidth() - 1)) || 
                (y1P < 0) || (Math.ceil(y1P) > (gsImg.getHeight() - 1))) {
                
                output[count] = sentinel;
                
            } else {
                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                // bilinear:
                //double v = imageProcessor.biLinearInterpolation(gsImg, x1P, y1P);
                
                // nearest neighbor
                double v = gsImg.getValue(Math.round(x1P), Math.round(y1P));
                
                output[count] = (int)Math.round(v);
            }
            
            count++;
        }
        
        IntensityDescriptor desc = new GsIntensityDescriptor(output);
        
        return desc;
    }
    
    /**
     * extract theta orientation for the (xCenter, yCenter) region 
     * and place in the descriptor.
     * Note that if the transformed pixel is out of bounds of the image,
     * a sentinel is the value for that location in the descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return 
     */
    private ThetaDescriptor extractThetaForBlock(int xCenter, int yCenter, 
        float[][] offsets) {
        
        int[] output = new int[offsets.length];
        
        int sentinel = PixelThetaDescriptor.sentinel;
                
        int count = 0;
        
        for (int i = 0; i < offsets.length; ++i) {
            
            float x1P = xCenter + offsets[i][0];
            
            float y1P = yCenter + offsets[i][1];
            
            if ((x1P < 0) || (Math.ceil(x1P) > (thetaImg.getWidth() - 1)) || 
                (y1P < 0) || (Math.ceil(y1P) > (thetaImg.getHeight() - 1))) {
                
                output[count] = sentinel;
                
            } else {
                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                // bilinear:
                //double v = imageProcessor.biLinearInterpolation(thetaImg, x1P, y1P);
                
                // nearest neighbor
                double v = thetaImg.getValue(Math.round(x1P), Math.round(y1P));
                
                output[count] = (int)Math.round(v);
            }
            
            count++;
        }
        
        ThetaDescriptor desc = new PixelThetaDescriptor(output);
        
        return desc;
    }
    
    /**
     * extract from the image and place in the descriptor.
     * Note that if the transformed pixel is out of bounds of the image,
     * a sentinel is the value for that location in the descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return 
     */
    private GradientDescriptor extractGsGradientForBlock(int xCenter, 
        int yCenter, float[][] offsets) {
        
        int sentinel = GsGradientDescriptor.sentinel;
        
        int[] output = new int[offsets.length];
                       
        int count = 0;
        
        for (int i = 0; i < offsets.length; ++i) {
            
            float x1P = xCenter + offsets[i][0];
            
            float y1P = yCenter + offsets[i][1];
            
            if ((x1P < 0) || (Math.ceil(x1P) > (gradientImg.getWidth() - 1)) || 
                (y1P < 0) || (Math.ceil(y1P) > (gradientImg.getHeight() - 1))) {
                
                output[count] = sentinel;
                
            } else {
                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                // bilinear:
                //double v = imageProcessor.biLinearInterpolation(gradientImg, x1P, y1P);
                
                // nearest neighbor
                double v = gradientImg.getValue(Math.round(x1P), Math.round(y1P));
                
                output[count] = (int)Math.round(v);
            }
            
            count++;
        }
        
        GradientDescriptor desc = new GsGradientDescriptor(output);
        
        return desc;
    }
    
    /**
     * extract the gradient from the image in 2X2 cells surrounding 
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return 
     */
    private GradientDescriptor extractGsGradientForCells(int xCenter, 
        int yCenter, int rotation) {
        
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
        
          3                   
          2  3     7    11     15
          1                    
          0  2     6   10@     14         
         -1                      
         -2  1     5     9     13
         -3                     
         -4  0     4     8     12
            -4 -3 -2 -1  0  1  2  3 
        */
        
        Transformer transformer = new Transformer();
        
        int[] output = new int[16];
        
        int sentinel = GsGradientDescriptor.sentinel;
        
        int count = 0;                
        for (int dx = -4; dx < 4; dx+=2) {
            for (int dy = -4; dy < 4; dy+=2) {
            
                float[] dXY11T = transformer.rotate(rotation, dx, dy); 
                float[] dXY12T = transformer.rotate(rotation, dx + 1, dy);
                float[] dXY21T = transformer.rotate(rotation, dx, dy + 1); 
                float[] dXY22T = transformer.rotate(rotation, dx + 1, dy + 1);
                
                float x11 = xCenter + dXY11T[0];
                float y11 = yCenter + dXY11T[1];
                
                float x12 = xCenter + dXY12T[0];
                float y12 = yCenter + dXY12T[1];
                
                float x21 = xCenter + dXY21T[0];
                float y21 = yCenter + dXY21T[1];
                
                float x22 = xCenter + dXY22T[0];
                float y22 = yCenter + dXY22T[1];

                if ((x11 < 0) || (Math.ceil(x11) > (gradientImg.getWidth() - 1)) ||
                    (y11 < 0) || (Math.ceil(y11) > (gradientImg.getHeight() - 1)) ||

                    (x12 < 0) || (Math.ceil(x12) > (gradientImg.getWidth() - 1)) ||
                    (y12 < 0) || (Math.ceil(y12) > (gradientImg.getHeight() - 1)) ||

                    (x21 < 0) || (Math.ceil(x21) > (gradientImg.getWidth() - 1)) ||
                    (y21 < 0) || (Math.ceil(y21) > (gradientImg.getHeight() - 1)) ||

                    (x22 < 0) || (Math.ceil(x22) > (gradientImg.getWidth() - 1)) ||
                    (y22 < 0) || (Math.ceil(y22) > (gradientImg.getHeight() - 1))
                    ) {
                    output[count] = sentinel;
                    count++;
                    continue;
                }

                int v11 = gradientImg.getValue(Math.round(x11), Math.round(y11));
                int v12 = gradientImg.getValue(Math.round(x12), Math.round(y12));
                int v21 = gradientImg.getValue(Math.round(x21), Math.round(y21));
                int v22 = gradientImg.getValue(Math.round(x22), Math.round(y22));

                output[count] = (int)Math.round((v11 + v12 + v21 + v22)/4.);

                count++;
            }
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
                
        int sentinel = ClrIntensityDescriptor.sentinel;
                
        int[] outputR = new int[offsets.length];
        int[] outputG = new int[outputR.length];
        int[] outputB = new int[outputR.length];
        
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
                
                //non-adaptive algorithms: nearest neighbor or bilinear
                
                // bilinear:
                //double[] v = imageProcessor.biLinearInterpolation(clrImg, x1P, y1P);
                
                // nearest neighbor
                int r = clrImg.getR(Math.round(x1P), Math.round(y1P));
                    
                int g = clrImg.getG(Math.round(x1P), Math.round(y1P));
                    
                int b = clrImg.getB(Math.round(x1P), Math.round(y1P));
                    
                outputR[count] = r;
                outputG[count] = g;
                outputB[count] = b;
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
        stat.setSumIntensitySqDiff(ssd);
        stat.setImg2PointIntensityErr(err2Sq);
        
        return stat;
    }
    
    /**
     * calculate the intensity and gradient based statistics between the two 
     * descriptors with the caveat that the 2nd descriptor is the one used to 
     * calculate the error (so make sure that pattern is consistently used by 
     * invoker).
     * 
     * @param descIntensity1
     * @param descGradient1
     * @param x1
     * @param y1
     * @param descIntensity2
     * @param descGradient2
     * @param x2
     * @param y2
     * @return 
     */
    public static FeatureComparisonStat calculateStats(
        IntensityDescriptor descIntensity1, GradientDescriptor descGradient1, 
        final int x1, final int y1,
        IntensityDescriptor descIntensity2, GradientDescriptor descGradient2, 
        final int x2, final int y2) {
    
        if (descIntensity1 == null) {
            throw new IllegalArgumentException("descIntensity1 cannot be null");
        }
        if (descGradient1 == null) {
            throw new IllegalArgumentException("descGradient1 cannot be null");
        }
        if (descIntensity2 == null) {
            throw new IllegalArgumentException("descIntensity2 cannot be null");
        }
        if (descGradient2 == null) {
            throw new IllegalArgumentException("descGradient2 cannot be null");
        }
        
        float err2SqIntensity = descIntensity2.sumSquaredError();
        
        float ssdIntensity = descIntensity1.calculateSSD(descIntensity2);
        
        float err2SqGradient = descGradient2.sumSquaredError();
        
        float ssdGradient = descGradient1.calculateSSD(descGradient2);
        
        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumIntensitySqDiff(ssdIntensity);
        stat.setImg2PointIntensityErr(err2SqIntensity);
        stat.setSumGradientSqDiff(ssdGradient);
        stat.setImg2PointGradientErr(err2SqGradient);
        
        return stat;
    }

}
