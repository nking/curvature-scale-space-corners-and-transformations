package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author nichole
 */
public class IntensityFeatures {
    
    /**
     * either gsImg or clrImg will be set by the constructor and they
     * are the grey scale image or the color image, respectively.
     */
    protected final Image clrImg;

    /**
     * either gsImg or clrImg will be set by the constructor and they
     * are the grey scale image or the color image, respectively.
     */
    protected final GreyscaleImage gsImg;
    
    /**
     * the half width of a block.  For example, to extract the 25 pixels
     * centered on (xc, yc), bHalfW would be '2'.  a minimum of 5 seems to be
     * necessary.
     */
    protected final int bHalfW;
    
    protected final float[][] xyOffsets;
    protected final boolean useNormalizedIntensities;
    protected final boolean useBinnedCellIntensities = true;
    /**
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
     *     extracted descriptor
     */
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();

    public IntensityFeatures(final GreyscaleImage image, 
        final int blockHalfWidths, final boolean useNormalizedIntensities) {

        this.gsImg = image;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
    }
    
    public IntensityFeatures(final Image image, 
        final int blockHalfWidths, final boolean useNormalizedIntensities) {

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
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensity(final int xCenter, final int yCenter, final int rotation) {
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
        
        if (useBinnedCellIntensities) {
            descriptor = extractIntensityForCells(xCenter, yCenter, rotation);
        } else {
            descriptor = extractIntensityForBlock(xCenter, yCenter, rotation);
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }
        
        return descriptor;
    }

    protected IntensityDescriptor extractIntensityForBlock(int xCenter, int yCenter, int rotation) {
        Transformer transformer = new Transformer();
        
        float[][] frameOffsets = transformer.transformXY(rotation, xyOffsets);
        
        IntensityDescriptor descriptor;
        
        if (gsImg != null) {
            descriptor = extractGsIntensityForBlock(xCenter, yCenter, frameOffsets);
        } else {
            descriptor = extractClrIntensityForBlock(xCenter, yCenter, frameOffsets);
        }
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        return descriptor;
    }

    protected IntensityDescriptor extractIntensityForCells(int xCenter, int yCenter, int rotation) {
        IntensityDescriptor descriptor;
        
        if (gsImg != null) {
            descriptor = extractGsIntensityForCells(xCenter, yCenter, rotation);
        } else {
            descriptor = extractClrIntensityForCells(xCenter, yCenter, rotation);
        }
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        return descriptor;
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
    public static FeatureComparisonStat calculateIntensityStat(IDescriptor desc1, 
        final int x1, final int y1, IDescriptor desc2, final int x2, final int y2) {
        
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
    
    protected void checkBounds(int x, int y) {
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
    protected IntensityDescriptor extractGsIntensityForBlock(int xCenter, int yCenter, float[][] offsets) {
        
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
                output[count] = (int) Math.round(v);
            }
            
            count++;
        }
        
        IntensityDescriptor desc = new GsIntensityDescriptor(output, offsets.length >> 1);
        
        return desc;
    }

    /**
     * extract the intensity from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return
     */
    protected IntensityDescriptor extractGsIntensityForCells(int xCenter, 
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
        float sentinel = GsIntensityDescriptor.sentinel;
        int cellDim = 2;
        int nCellsAcross = 6;
        int range0 = (int) (cellDim * ((float) nCellsAcross / 2.f));
        int nColsHalf = nCellsAcross / 2;
        int centralPixelIndex = (nColsHalf * nCellsAcross) + nColsHalf;
        float[] output = new float[nCellsAcross * nCellsAcross];
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];
        
        int count = 0;
        
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                
                // --- calculate values for the cell ---
                boolean withinBounds = transformCellCoordinates(rotation, 
                    xCenter, yCenter, dx, dy, cellDim, gsImg.getWidth(), 
                    gsImg.getHeight(), xT, yT);
                
                if (!withinBounds) {
                    if (count == centralPixelIndex) {
                        return null;
                    }
                    output[count] = sentinel;
                    count++;
                    continue;
                }
                int cCount = 0;
                int v = 0;
                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int v0 = gsImg.getValue(x, y);
                    v += v0;
                    cCount++; // img changed?  w=480,h=640
                }
                if (cCount == 0) {
                    output[count] = sentinel;
                    count++;
                    continue;
                }
                v /= (float) cCount;
                output[count] = v;
                count++;
            }
        }
        IntensityDescriptor desc = new GsIntensityDescriptor(output, centralPixelIndex);
        return desc;
    }

    /**
     * extract the intensity from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return
     */
    protected IntensityDescriptor extractClrIntensityForCells(int xCenter, 
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
        int sentinel = ClrIntensityDescriptor.sentinel;
        int cellDim = 2;
        int nCellsAcross = 6;
        int range0 = (int) (cellDim * ((float) nCellsAcross / 2.f));
        int nColsHalf = nCellsAcross / 2;
        int centralPixelIndex = (nColsHalf * nCellsAcross) + nColsHalf;
        
        int[] red = new int[nCellsAcross * nCellsAcross];
        int[] green = new int[red.length];
        int[] blue = new int[red.length];
        
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];
        
        int count = 0;
        
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                // --- calculate values for the cell ---
                boolean withinBounds = transformCellCoordinates(rotation, 
                    xCenter, yCenter, dx, dy, cellDim, clrImg.getWidth(), 
                    clrImg.getHeight(), xT, yT);
                
                if (!withinBounds) {
                    if (count == centralPixelIndex) {
                        return null;
                    }
                    red[count] = sentinel;
                    green[count] = sentinel;
                    blue[count] = sentinel;
                    count++;
                    continue;
                }
                int cCount = 0;
                int r = 0;
                int g = 0;
                int b = 0;
                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    r += clrImg.getR(x, y);
                    g += clrImg.getG(x, y);
                    b += clrImg.getB(x, y);
                    cCount++;
                }
                if (cCount == 0) {
                    red[count] = sentinel;
                    green[count] = sentinel;
                    blue[count] = sentinel;
                    count++;
                    continue;
                }
                r /= (float) cCount;
                g /= (float) cCount;
                b /= (float) cCount;
                red[count] = r;
                green[count] = g;
                blue[count] = b;
                count++;
            }
        }
        
        IntensityDescriptor desc = new ClrIntensityDescriptor(red, green, blue, 
            centralPixelIndex);
        
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
    protected IntensityDescriptor extractClrIntensityForBlock(int xCenter, int yCenter, float[][] offsets) {
        int sentinel = ClrIntensityDescriptor.sentinel;
        int[] outputR = new int[offsets.length];
        int[] outputG = new int[outputR.length];
        int[] outputB = new int[outputR.length];
        int count = 0;
        for (int i = 0; i < offsets.length; ++i) {
            float x1P = xCenter + offsets[i][0];
            float y1P = yCenter + offsets[i][1];
            if ((x1P < 0) || (Math.ceil(x1P) > (clrImg.getWidth() - 1)) || (y1P < 0) || (Math.ceil(y1P) > (clrImg.getHeight() - 1))) {
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
        IntensityDescriptor desc = new ClrIntensityDescriptor(outputR, outputG, outputB, offsets.length);
        return desc;
    }

    protected boolean transformCellCoordinates(int rotation, int x0, int y0, 
        int dx, int dy, int cellDim, int imageWidth, int imageHeight, 
        float[] outputX, float[] outputY) {
        
        Transformer transformer = new Transformer();
        
        int cCount = 0;
        
        for (int dxc = 0; dxc < cellDim; ++dxc) {
            for (int dyc = 0; dyc < cellDim; ++dyc) {
                
                transformer.rotate(rotation, dx + dxc, dy + dyc, outputX, 
                    outputY, cCount);
                
                outputX[cCount] += x0;
                outputY[cCount] += y0;
                
                if ((outputX[cCount] < 0) || 
                    (Math.ceil(outputX[cCount]) > (imageWidth - 1)) || 
                    (outputY[cCount] < 0) || 
                    (Math.ceil(outputY[cCount]) > (imageHeight - 1))) {
                    
                    return false;
                }
                
                cCount++;
            }
        }
        return true;
    }
    
    /**
     * calculate the intensity and gradient based statistics between the two
     * descriptors with the caveat that the 2nd descriptor is the one used to
     * calculate the error (so make sure that pattern is consistently used by
     * invoker).
     *
     * @param descIntensity1
     * @param x1
     * @param y1
     * @param descIntensity2
     * @param x2
     * @param y2
     * @return
     */
    public static FeatureComparisonStat calculateStats(
        IntensityDescriptor descIntensity1, final int x1, final int y1,
        IntensityDescriptor descIntensity2, final int x2, final int y2) {

        if (descIntensity1 == null) {
            throw new IllegalArgumentException("descIntensity1 cannot be null");
        }
        if (descIntensity2 == null) {
            throw new IllegalArgumentException("descIntensity2 cannot be null");
        }

        float err2SqIntensity = descIntensity2.sumSquaredError();

        float err1SqIntensity = descIntensity1.sumSquaredError();
        
        //TODO: revisit this:
        //float errSqIntensity = (err2SqIntensity > err1SqIntensity) ?
        //    err2SqIntensity : err1SqIntensity;

        float ssdIntensity = descIntensity1.calculateSSD(descIntensity2);

        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumIntensitySqDiff(ssdIntensity);
        stat.setImg2PointIntensityErr(err2SqIntensity);

        return stat;
    }
}
