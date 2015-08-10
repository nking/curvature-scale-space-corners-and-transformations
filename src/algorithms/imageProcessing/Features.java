package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
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

    /*
    for theta descriptors, 1 of 3 choices:
       (1) PixelThetaDescriptor: each pixel in 5x5 block stored.
           comparion is pixel by pixel angular differences
       (2) PixelThetaDescriptor: pixels binned into 2x2 cells and 16 total
           surrounding the center.
           comparison is item by item angular differences
       (3) HistogramThetaDescriptor: pixels within mxm cells are made into
           histograms of w bin width, and 8 or 16 such cells surrounding the
           center.
           each histogram before compared is normalized by peak?
           corrections for wrap around during comparisons?
           comparisons are then similar to SSD?

       Note that use of (1) does not have S/N correction and it requires
       very good rotation alignment to achieve reasonable results, so the
       "dominant orientation" as a stable relative angle for comparing frames
       would need a refinement such as trying small degrees around it within
       a range to find best match for better use of (1).
    */
    protected final int thetaType = 2;

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
        GreyscaleImage theThetaImg, int blockHalfWidths,
        boolean useNormalizedIntensities) {

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
        GreyscaleImage theThetaImg, int blockHalfWidths,
        boolean useNormalizedIntensities) {

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
     * @return descriptor holding gradient rotated.  NOTE that the method may
     * return null if ant part of the center coordinate cell was was rotated
     * out of the frame.
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

        if (descriptor != null) {
            descriptors.put(rotKey, descriptor);
        }

        return descriptor;
    }

    /**
     * extract theta orientation from the image for (xCenter, yCenter) and
     * return it in a descriptor.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return descriptor for theta extracted from the image.  NOTE that this
     * method may return null if the (xCenter, yCenter) frame is rotated in
     * part out of the frame.
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

        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }

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

        if (thetaType == 1) {
            descriptor = extractThetaForBlock(xCenter, yCenter, rotation,
                frameOffsets);
        } else if (thetaType == 2) {
            descriptor = extractThetaForCells(xCenter, yCenter, rotation);
        } else if (thetaType == 3) {
            throw new UnsupportedOperationException("not yet implemented");
            //descriptor = extractThetaHistograms(xCenter, yCenter, rotation);
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

        IntensityDescriptor desc = new GsIntensityDescriptor(output,
            offsets.length >> 1);

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
        int rotation, float[][] offsets) {

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

                // bilinear: (needs to be adapted for angular addition)
                //double v = imageProcessor.biLinearInterpolation(thetaImg, x1P, y1P);

                // nearest neighbor
                double v = thetaImg.getValue(Math.round(x1P), Math.round(y1P));

                v += rotation;

                if (v >= 360) {
                    v = v - 360;
                }

                output[count] = (int)Math.round(v);
            }

            count++;
        }

        ThetaDescriptor desc = new PixelThetaDescriptor(output, 
            offsets.length >> 1);

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

        GradientDescriptor desc = new GsGradientDescriptor(output,
            offsets.length >> 1);

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

        int sentinel = GsGradientDescriptor.sentinel;
        
        int centralPixelIndex = 10;

        int cellDim = 2;
        int nCellsAcross = 4;
        int range0 = (int)(cellDim * ((float)nCellsAcross/2.f));

        int[] output = new int[nCellsAcross * nCellsAcross];
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];

        int count = 0;
        for (int dx = -range0; dx < range0; dx+=cellDim) {
            for (int dy = -range0; dy < range0; dy+=cellDim) {

                // --- calculate values for the cell ---
                boolean withinBounds = transformCellCoordinates(gradientImg,
                    rotation, xCenter, yCenter, dx, dy, cellDim, xT, yT);

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
                    v += gradientImg.getValue(x, y);
                    cCount++;
                }

                if (cCount == 0) {
                    output[count] = sentinel;
                    count++;
                    continue;
                }

                v /= (float)cCount;

                // subtract the "dominant orientation"
                v -= rotation;

                if (v > 359) {
                    v = v % 360;
                } else if (v < 0) {
                    v = v + 360;
                }

                output[count] = v;

                count++;
            }
        }

        GradientDescriptor desc = new GsGradientDescriptor(output, 
            centralPixelIndex);

        return desc;
    }

    /**
     * extract theta in angular degrees from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells, making the value correction for
     * "dominant orientation" too.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return descriptor for the region given and orientation given.  NOTE
     * that it can return null if the central region is shifted out of the
     * image frame.  If any other cells excepting the center are shifted
     * out, an object is returned with those cells holding the sentinel as a
     * value (further code interprets those as a skip instruction).
     */
    private ThetaDescriptor extractThetaForCells(int xCenter, int yCenter,
        int rotation) {

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

        int sentinel = ThetaDescriptor.sentinel;

        int centralPixelIndex = 10;
        
        int cellDim = 2;
        int nCellsAcross = 4;
        int range0 = (int)(cellDim * ((float)nCellsAcross/2.f));

        int[] output = new int[nCellsAcross * nCellsAcross];
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];

        int count = 0;
        for (int dx = -range0; dx < range0; dx+=cellDim) {
            for (int dy = -range0; dy < range0; dy+=cellDim) {

                // --- calculate values for the cell ---
                boolean withinBounds = transformCellCoordinates(thetaImg,
                    rotation, xCenter, yCenter, dx, dy, cellDim, xT, yT);
                
                if (!withinBounds) {
                    if (count == centralPixelIndex) {
                        return null;
                    }
                    output[count] = sentinel;
                    count++;
                    continue;
                }

                int maxGrd = Integer.MIN_VALUE;
                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int v = gradientImg.getValue(x, y);
                    if (v > maxGrd) {
                        maxGrd = v;
                    }
                }

                float lowLimit = 0.1f*maxGrd;

                /*
                 because these are angles from 0 to 360, need to add them
                 while considering their quadrants.

                 e.g. (0 + 350 + 340)/3. should equal 350, but if
                 the quadrants are not used, the result is 130.

                Can add in pairs, but need to correct the result.

                averaging 4 numbers:
                (a + b + c + d)/4. = (a/4) + (b/4) + (c/4) + (d/4)

                averaging 4 numbers through 3 pair averages:
                1) (a/2) + (b/2)
                2) ((a/2) + (b/2))/2 + (c/2) = (a/4) + (b/4) + (c/2)
                3) ((a/4) + (b/4) + (c/2))/2 + (d/2)
                    = (a/8) + (b/8) + (c/4) + (d/2)

                and that needs to be corrected to (a/4) + (b/4) + (c/4) + (d/4)

                for n=5, the successive pair averages is
                1) (a/2) + (b/2)
                2) ((a/2) + (b/2))/2 + (c/2) = (a/4) + (b/4) + (c/2)
                3) ((a/4) + (b/4) + (c/2))/2 + (d/2)
                    = (a/8) + (b/8) + (c/4) + (d/2)
                4) ((a/8) + (b/8) + (c/4) + (d/2))/2 + (e/2)
                    = (a/16) + (b/16) + (c/8) + (d/4) + (e/2)
                       0        1         2       3       4
                and that needs to be corrected to (a/5) + (b/5) + (c/5) + (d/5) + (e/5)

                for all but the first number:
                    v[i]                 v[i]
                ----------  +  v[i]*x  = ---
                (1<<(n-i))                n
                
                for first variable,
                   total += v[0] * ((1<<(n-1)) - n)/(n*(1<<(n-1)))

                from i=1 to i < n
                   total += v[i] * ((1<<(n-i)) - n)/(n*(1<<(n-i)))

                **(This is probably one reason to prefer histograms).**
                */
                int cCount = 0;
                float avg = 0;

                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int vGradient = gradientImg.getValue(x, y);
                    if (vGradient < lowLimit) {
                        continue;
                    }

                    int v = thetaImg.getValue(x, y);

                    if (cCount == 0) {
                        avg = v;
                    } else {
                        avg += AngleUtil.getAngleAverageInDegrees(avg, v);
                    }

                    cCount++;
                }

                if (cCount == 0) {
                    output[count] = sentinel;
                    count++;
                    continue;
                }

                /*
                correction for adding pairs:
                for first variable,
                   total += v[0] * ((1<<(n-1)) - n)/(n*(1<<(n-1)))

                from i=1 to i < n
                   total += v[i] * ((1<<(n-i)) - n)/(n*(1<<(n-i)))
                */
                int n = cCount;
                cCount = 0;
                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int vGradient = gradientImg.getValue(x, y);
                    if (vGradient < lowLimit) {
                        continue;
                    }

                    //TODO: there's still an error here when the above angleutil
                    //made a quadrant correction.
                    //need to extract the angles used in the angleutil
                    //store them in an array above to iterate over instead of 
                    //fetching the values from the image again here.
                    
                    int v = thetaImg.getValue(x, y);

                    float c0;

                    if (cCount == 0) {
                        c0 = 1 << (n - 1);
                    } else {
                        c0 = 1 << (n - i);
                    }

                    float add = v * (c0 - n)/(n * c0);

                    avg += add;

                    cCount++;
                }

                // subtract the "dominant orientation"
                avg -= rotation;

                if (avg > 359) {
                    avg = avg % 360;
                } else if (avg < 0) {
                    avg = avg + 360;
                }

                output[count] = Math.round(avg);

                count++;
            }
        }

        ThetaDescriptor desc = new PixelThetaDescriptor(output, 
            centralPixelIndex);

        return desc;
    }

    /**
     * extract the gradient from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return an instance of ThetaHistogramDescriptor.  Note that it's
     * internal array may contain null values if a cell was out of bounds
     * so any comparison methods should check for null.
     */
    private ThetaDescriptor extractThetaHistograms(int xCenter, int yCenter,
        int rotation) {

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

        /*
         theta' histogram is made for each group (subdivisions of the block)
         for n=8 intervals of 45 degrees each.
         group is 4x4... or so.  above, still using 2x2 which would be
         too small for a histogram

         S/N corrections:

         consider discarding values that have gradientXY < critical value
         (suggested as 0.1 * maximum of gradientXY in this cell.  reference?)

         create a weight for each pixel within the 4x4 group based upon the gradient value:
         w[pix] = gXY[pix] / sqrt(sum of squares of gXY[pix] in 4x4 group) + eps)
         where eps is a small number to avoid divide by zero errors

         then the contribution of the pixel to the histogram bin is
         weighted by that value, that is histogram increments are fractional
         not integer.
         (the bin the pixel is placed in is the same, just its contribution
         that is weighted by fraction of total intensity).
        */

        int sentinel = ThetaDescriptor.sentinel;

        int centralPixelIndex = 10;
        
        int cellDim = 2;
        int nCellsAcross = 4;
        int range0 = (int)(cellDim * ((float)nCellsAcross/2.f));

        int nHistBins = 36;
        float binSize = 360.f/(float)nHistBins;

        // TODO: change the output here to be an array of histograms
        float[][] output = new float[nCellsAcross * nCellsAcross][];

        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];

        int count = 0;
        for (int dx = -range0; dx < range0; dx+=cellDim) {
            for (int dy = -range0; dy < range0; dy+=cellDim) {

                // --- calculate values for the cell ---
                boolean withinBounds = transformCellCoordinates(thetaImg,
                    rotation, xCenter, yCenter, dx, dy, cellDim, xT, yT);

                if (!withinBounds) {
                    if (count == centralPixelIndex) {
                        return null;
                    }
                    output[count] = null;
                    count++;
                    continue;
                }

                int maxGrd = Integer.MIN_VALUE;
                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int v = gradientImg.getValue(x, y);
                    if (v > maxGrd) {
                        maxGrd = v;
                    }
                }

                float lowLimit = 0.1f*maxGrd;
                int cellGrdSum = 0;
                int cCount = 0;
                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int v = gradientImg.getValue(x, y);
                    if (v < lowLimit) {
                        continue;
                    }
                    cellGrdSum += v;
                    cCount++;
                }

                if (cCount == 0) {
                    output[count] = null;
                    count++;
                    continue;
                }

                cCount = 0;

                output[count] = new float[nHistBins];

                for (int i = 0; i < xT.length; ++i) {
                    int x = Math.round(xT[i]);
                    int y = Math.round(yT[i]);
                    int vGradient = gradientImg.getValue(x, y);
                    if (vGradient < lowLimit) {
                        continue;
                    }

                    // subtract "dominant orientation"
                    int v = thetaImg.getValue(x, y) - rotation;
                    if (v < 0) {
                        v += 360;
                    }

                    float weightedContrib = (float)vGradient/(float)cellGrdSum;

                    //TODO: add weightedContrib to the histogram
                    //  bin the v belongs in

                    cCount++;
                }
            }
        }

        ThetaDescriptor desc = new HistogramThetaDescriptor(output,
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
            outputB, offsets.length);

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

    /**
     * calculate the intensity and gradient based statistics between the two
     * descriptors with the caveat that the 2nd descriptor is the one used to
     * calculate the error (so make sure that pattern is consistently used by
     * invoker).
     *
     * @param descIntensity1
     * @param descGradient1
     * @param descTheta1
     * @param x1
     * @param y1
     * @param descIntensity2
     * @param descGradient2
     * @param descTheta2
     * @param x2
     * @param y2
     * @return
     */
    public static FeatureComparisonStat calculateStats(
        IntensityDescriptor descIntensity1, GradientDescriptor descGradient1,
        ThetaDescriptor descTheta1,
        final int x1, final int y1,
        IntensityDescriptor descIntensity2, GradientDescriptor descGradient2,
        ThetaDescriptor descTheta2,
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
        if (descTheta1 == null) {
            throw new IllegalArgumentException("descTheta1 cannot be null");
        }
        if (descTheta2 == null) {
            throw new IllegalArgumentException("descTheta2 cannot be null");
        }

        float err2SqIntensity = descIntensity2.sumSquaredError();

        float ssdIntensity = descIntensity1.calculateSSD(descIntensity2);

        float err2SqGradient = descGradient2.sumSquaredError();

        float ssdGradient = descGradient1.calculateSSD(descGradient2);

        float compTheta = descTheta1.calculateDifference(descTheta2);
        float compThetaErr = descTheta2.calculateError();

        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumIntensitySqDiff(ssdIntensity);
        stat.setImg2PointIntensityErr(err2SqIntensity);
        stat.setSumGradientSqDiff(ssdGradient);
        stat.setImg2PointGradientErr(err2SqGradient);
        stat.setSumThetaDiff(compTheta);
        stat.setImg2PointThetaErr(compThetaErr);

        return stat;
    }

    protected boolean transformCellCoordinates(GreyscaleImage img,
        int rotation,
        int x0, int y0, int dx, int dy,
        int cellDim, float[] outputX, float[] outputY) {

        Transformer transformer = new Transformer();

        int cCount = 0;
        for (int dxc = 0; dxc < cellDim; ++dxc) {
            for (int dyc = 0; dyc < cellDim; ++dyc) {

                transformer.rotate(rotation, dx + dxc, dy + dyc,
                    outputX, outputY, cCount);

                outputX[cCount] += x0;
                outputY[cCount] += y0;

                if ((outputX[cCount] < 0)
                    || (Math.ceil(outputX[cCount]) > (img.getWidth() - 1))
                    || (outputY[cCount] < 0)
                    || (Math.ceil(outputY[cCount]) > (img.getHeight() - 1))) {

                    return false;
                }

                cCount++;
            }
        }

        return true;
    }

}
