package algorithms.imageProcessing.features;

import algorithms.imageProcessing.transform.Transformer;
import algorithms.QuickSort;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;

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
     * created by method calculateGradientWithGreyscale and kept for use
     * with calculate45DegreeOrientation internally.
     */
    //TODO: since these are used sparsely, consider calculating only for a point and on demand.
    protected GreyscaleImage gXY = null;
    protected GreyscaleImage theta = null;
    
    /**
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
     *     extracted descriptor
     */
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    /**
     * key = pixel coordinates;
     * value = 45 degree resolution orientation angle determined from greyscale
     * intensities.  (NOTE, this is a cache to hold orientation
     * angles for pixels that are corner regions or nearby dithered locations
     * of corner regions.  The CornerRegion.orientation itself is not always
     * w.r.t. the same curves, so cannot always be used.)
     */
    protected Map<PairInt, Integer> mapOf45DegOr = new
        HashMap<PairInt, Integer>();
    
    /**
     * key = pixel coordinates;
     * value = dominant orientation
     */
    protected Map<PairInt, Integer> mapOfOrientation = new HashMap<PairInt, Integer>();
    
    protected final RotatedOffsets rotatedOffsets;
    
    /**
     * Constructor that keeps a reference to the image.  Note, the constructor
     * without the image argument is preferred followed by using methods that
     * accept the image as an argument because that as a model allows the JVM
     * to manage references and garbage collection better.
     * @param image
     * @param blockHalfWidths
     * @param useNormalizedIntensities
     * @param rotatedOffsets 
     */
    public IntensityFeatures(final GreyscaleImage image, 
        final int blockHalfWidths, final boolean useNormalizedIntensities,
        RotatedOffsets rotatedOffsets) {

        this.gsImg = image;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.rotatedOffsets = rotatedOffsets;
    }
    
    /**
     * Constructor that keeps a reference to the image.  Note, the constructor
     * without the image argument is preferred followed by using methods that
     * accept the image as an argument because that as a model allows the JVM
     * to manage references and garbage collection better.
     * @param image
     * @param blockHalfWidths
     * @param useNormalizedIntensities
     * @param rotatedOffsets 
     */
    public IntensityFeatures(final Image image, 
        final int blockHalfWidths, final boolean useNormalizedIntensities,
        RotatedOffsets rotatedOffsets) {

        this.gsImg = null;
        this.clrImg = image;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.rotatedOffsets = rotatedOffsets;
    }
    
    /**
     * constructor for an instance which retains no image references and
     * relies on the user to consistently provide the same image to 
     * instance methods.
     * @param blockHalfWidths
     * @param useNormalizedIntensities 
     * @param rotatedOffsets 
     */
    public IntensityFeatures(final int blockHalfWidths, final boolean 
        useNormalizedIntensities, RotatedOffsets rotatedOffsets) {

        this.gsImg = null;
        this.clrImg = null;
        this.bHalfW = blockHalfWidths;
        this.useNormalizedIntensities = useNormalizedIntensities;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.rotatedOffsets = rotatedOffsets;
    }
    
    public boolean gradientWasCreated() {
        return (gXY != null);
    }
    
    public GreyscaleImage getGradientImage() {
        return gXY;
    }

    /**
     * create the gradient of the given greyscale image in order to use it
     * internally for the method calculate45DegreeOrientation.
     * @param gsImg 
     */
    public void calculateGradientWithGreyscale(GreyscaleImage gsImg) {
        
        if (gXY != null) {
            throw new IllegalArgumentException("the gradient image has already"
            + " been created");
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage gX = gsImg.copyToFullRangeIntImage();
        float[] kX = new float[]{-1, 0, 1};
        imageProcessor.applyKernel1D(gX, kX, true);
        
        GreyscaleImage gY = gsImg.copyToFullRangeIntImage();
        float[] kY = new float[]{1, 0, -1};
        imageProcessor.applyKernel1D(gY, kY, false);
        
        gXY = imageProcessor.combineConvolvedImages(gX, gY);
        
        theta = imageProcessor.computeTheta360_0(gX, gY);
        
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
        
        if (clrImg == null && gsImg == null) {
            throw new IllegalStateException("no images were provided at "
            + "construction so " 
            + "extractIntensity(image, xCenter, yCenter, rotatation) "
            + "must be used instead");
        }
        
        if (clrImg != null) {
            return extractIntensity(clrImg, xCenter, yCenter, rotation);
        } else {
            return extractIntensity(gsImg, xCenter, yCenter, rotation);
        }
    }
        
    /**
     * extract the intensity from the image for the given block center and
     * return it in a descriptor.  Note, this method is a short cut for
     * use when the known internal type extract is by cells.
     * @param img greyscale image which needs a descriptor for given location
     * and orientation.  note that because the instance caches previous
     * calculations, the invoker must supply the same image each time.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    IntensityDescriptor extractIntensityCellDesc(GreyscaleImage img, 
        final int xCenter, final int yCenter, final int rotation) {
        
        if (img == null) {
            throw new IllegalStateException("img cannot be null");
        }
                
        checkBounds(img, xCenter, yCenter);
        
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
        
        descriptor = extractGsIntensityForCells2(img, xCenter, yCenter, 
            rotation);
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }
        
        return descriptor;
    }
    
    /**
     * extract the intensity from the image for the given block center and
     * return it in a descriptor.
     * @param img greyscale image which needs a descriptor for given location
     * and orientation.  note that because the instance caches previous
     * calculations, the invoker must supply the same image each time.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensity(GreyscaleImage img, 
        final int xCenter, final int yCenter, final int rotation) {
        
        if (img == null) {
            throw new IllegalStateException("img cannot be null");
        }
                
        checkBounds(img, xCenter, yCenter);
        
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
            descriptor = extractIntensityForCells(img, xCenter, yCenter, 
                rotation);
        } else {
            descriptor = extractIntensityForBlock(img, xCenter, yCenter, 
                rotation);
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }
        
        return descriptor;
    }
    
    /**
     * extract the intensity from the image for the given block center and
     * return it in a descriptor.
     * @param img color image which needs a descriptor for given location
     * and orientation.  note that because the instance caches previous
     * calculations, the invoker must supply the same image each time.
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensity(Image img, final int xCenter, 
        final int yCenter, final int rotation) {
        
        if (img == null) {
            throw new IllegalStateException("img cannot be null");
        }
                    
        checkBounds(img, xCenter, yCenter);
        
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
            descriptor = extractIntensityForCells(img, xCenter, yCenter, 
                rotation);
        } else {
            descriptor = extractIntensityForBlock(img, xCenter, yCenter, 
                rotation);
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }
        
        return descriptor;
    }
    
    protected IntensityDescriptor extractIntensityForBlock(GreyscaleImage img,
        int xCenter, int yCenter, int rotation) {
        
        Transformer transformer = new Transformer();
        
        float[][] frameOffsets = transformer.transformXY(rotation, xyOffsets);
        
        IntensityDescriptor descriptor;
        
        descriptor = extractGsIntensityForBlock(img, xCenter, yCenter, 
            frameOffsets);
       
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        return descriptor;
    }
    
    protected IntensityDescriptor extractIntensityForBlock(Image img, 
        int xCenter, int yCenter, int rotation) {
        
        Transformer transformer = new Transformer();
        
        float[][] frameOffsets = transformer.transformXY(rotation, xyOffsets);
        
        IntensityDescriptor descriptor = extractClrIntensityForBlock(img, 
            xCenter, yCenter, frameOffsets);
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        return descriptor;
    }

    protected IntensityDescriptor extractIntensityForCells(GreyscaleImage img,
        int xCenter, int yCenter, int rotation) {
        
        IntensityDescriptor descriptor = 
            extractGsIntensityForCells2(img, xCenter, yCenter, rotation);
                
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        return descriptor;
    }
    
    protected IntensityDescriptor extractIntensityForCells(Image img,
        int xCenter, int yCenter, int rotation) {
        
        IntensityDescriptor descriptor = 
            extractClrIntensityForCells(img, xCenter, yCenter, rotation);
        
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
        if (clrImg == null && gsImg == null) {
            throw new IllegalStateException("no images were provided at "
            + "construction so " 
            + "extractIntensity(image, xCenter, yCenter, rotatation) "
            + "must be used instead");
        }
        if (clrImg != null) {
            checkBounds(clrImg, x, y);
        } else {
            checkBounds(gsImg, x, y);
        }
    }
    
    protected void checkBounds(GreyscaleImage img, int x, int y) {
        if (!isWithinXBounds(img, x)) {
            throw new IllegalArgumentException("x is out of bounds of image");
        }
        if (!isWithinYBounds(img, y)) {
            throw new IllegalArgumentException("y is out of bounds of image");
        }
    }
    
    protected void checkBounds(Image img, int x, int y) {
        if (!isWithinXBounds(img, x)) {
            throw new IllegalArgumentException("x is out of bounds of image");
        }
        if (!isWithinYBounds(img, y)) {
            throw new IllegalArgumentException("y is out of bounds of image");
        }
    }

    protected boolean isWithinXBounds(int x) {
        if (clrImg == null && gsImg == null) {
            throw new IllegalStateException("no images were provided at "
            + "construction so " 
            + "extractIntensity(image, xCenter, yCenter, rotatation) "
            + "must be used instead");
        }
        if (clrImg != null) {
            return isWithinXBounds(clrImg, x);
        } else {
            return isWithinXBounds(gsImg, x);
        }
    }
    
    protected boolean isWithinYBounds(int y) {
        if (clrImg == null && gsImg == null) {
            throw new IllegalStateException("no images were provided at "
            + "construction so " 
            + "extractIntensity(image, xCenter, yCenter, rotatation) "
            + "must be used instead");
        }
        if (clrImg != null) {
            return isWithinYBounds(clrImg, y);
        } else {
            return isWithinYBounds(gsImg, y);
        }
    }
    
    public boolean isWithinXBounds(GreyscaleImage img, int x) {
        if (x < 0) {
            return false;
        }
        if (x > (img.getWidth() - 1)) {
            return false;
        }
        return true;
    }
    
    public boolean isWithinYBounds(GreyscaleImage img, int y) {
        if (y < 0) {
            return false;
        }
        if (y > (img.getHeight() - 1)) {
            return false;
        }
        return true;
    }

    public boolean isWithinXBounds(Image img, int x) {
        if (x < 0) {
            return false;
        }
        if (x > (img.getWidth() - 1)) {
            return false;
        }
        return true;
    }
    
    public boolean isWithinYBounds(Image img, int y) {
        if (y < 0) {
            return false;
        }
        if (y > (img.getHeight() - 1)) {
            return false;
        }
        return true;
    }

    /**
     * extract the intensity from the image and place in the descriptor.
     * Note that if the transformed pixel is out of bounds of the image,
     * a sentinel is the value for that location in the descriptor.
     * @param img
     * @param xCenter
     * @param yCenter
     * @return
     */
    protected IntensityDescriptor extractGsIntensityForBlock(GreyscaleImage img,
        int xCenter, int yCenter, float[][] offsets) {
        
        float[] output = new float[offsets.length];
        
        float sentinel = GsIntensityDescriptor.sentinel;
        
        int count = 0;
        
        for (int i = 0; i < offsets.length; ++i) {
            
            float x1P = xCenter + offsets[i][0];
            
            float y1P = yCenter + offsets[i][1];
            
            if ((x1P < 0) || (Math.ceil(x1P) > (img.getWidth() - 1)) || 
                (y1P < 0) || (Math.ceil(y1P) > (img.getHeight() - 1))) {
                output[count] = sentinel;
            } else {
                //non-adaptive algorithms: nearest neighbor or bilinear
                // bilinear:
                //double v = imageProcessor.biLinearInterpolation(gsImg, x1P, y1P);
                // nearest neighbor
                double v = img.getValue(Math.round(x1P), Math.round(y1P));
                output[count] = (int) Math.round(v);
            }
            
            count++;
        }
        
        IntensityDescriptor desc = new GsIntensityDescriptor(output, 
            offsets.length >> 1, (int)Math.sqrt(output.length));
        
        return desc;
    }

    /**
     * extract the intensity from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param img
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return
     */
    protected IntensityDescriptor extractClrIntensityForCells(Image img,
        int xCenter, int yCenter, int rotation) {

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
                    xCenter, yCenter, dx, dy, cellDim, img.getWidth(), 
                    img.getHeight(), xT, yT);
                
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
                    r += img.getR(x, y);
                    g += img.getG(x, y);
                    b += img.getB(x, y);
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
     * @param img
     * @param xCenter
     * @param yCenter
     * @param offsets
     * @return
     */
    protected IntensityDescriptor extractClrIntensityForBlock(Image img,
        int xCenter, int yCenter, float[][] offsets) {
        
        int sentinel = ClrIntensityDescriptor.sentinel;
        int[] outputR = new int[offsets.length];
        int[] outputG = new int[outputR.length];
        int[] outputB = new int[outputR.length];
        int count = 0;
        for (int i = 0; i < offsets.length; ++i) {
            float x1P = xCenter + offsets[i][0];
            float y1P = yCenter + offsets[i][1];
            if ((x1P < 0) || (Math.ceil(x1P) > (img.getWidth() - 1)) 
                || (y1P < 0) || (Math.ceil(y1P) > (img.getHeight() - 1))) {
                outputR[count] = sentinel;
                outputG[count] = sentinel;
                outputB[count] = sentinel;
            } else {
                //non-adaptive algorithms: nearest neighbor or bilinear
                // bilinear:
                //double[] v = imageProcessor.biLinearInterpolation(img, x1P, y1P);
                // nearest neighbor
                int r = img.getR(Math.round(x1P), Math.round(y1P));
                int g = img.getG(Math.round(x1P), Math.round(y1P));
                int b = img.getB(Math.round(x1P), Math.round(y1P));
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

    static boolean transformCellCoordinates(int rotation, int x0, int y0, 
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
    
    public static int getDefaultCellDimForExtract() {
        int cellDim = 2;
        return cellDim;
    }
    public static int getDefaultNCellsAcrossForExtract() {
        int nCellsAcross = 6;
        return nCellsAcross;
    }
    public static int getDefaultLengthForCellExtractOffsets() {
        int nCellsAcross = getDefaultNCellsAcrossForExtract();
        int len = 2*nCellsAcross;
        return len * len;
    }
    
    /**
     * extract the intensity from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param img
     * @param xCenter
     * @param yCenter
     * @param rotation in degrees
     * @return
     */
    public IntensityDescriptor extractGsIntensityForCells2(GreyscaleImage img,
        int xCenter, int yCenter, int rotation) {
        
        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }
     
        float sentinel = GsIntensityDescriptor.sentinel;
        
        int cellDim = 2;        
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        
        int[] xOffsets = rotatedOffsets.getXOffsets(rotation);
        int[] yOffsets = rotatedOffsets.getYOffsets(rotation);
        
        int w = img.getWidth();
        int h = img.getHeight();

        float[] output = new float[nCellsAcross * nCellsAcross];

        int n2 = cellDim * cellDim;
        
        //index of center pixel in array of descriptor:
        int centralPixelIndex = (nCellsAcross*nColsHalf) + nColsHalf;
        
        int count = 0;
        int idx = 0;
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                boolean withinBounds = true;
                int cCount = 0;
                // ---- sum within the cell ----
                int v = 0;
                for (int i = 0; i < n2; ++i) {
                    int xOff = xOffsets[idx];
                    int yOff = yOffsets[idx];
                    idx++;
                    int x = xOff + xCenter;
                    int y = yOff + yCenter;
                    if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                        withinBounds = false;
                        break;
                    }
                    v += img.getValue(x, y);
                    cCount++;
                }
                if (!withinBounds || (cCount == 0)) {
                    if (count == centralPixelIndex) {
                        return null;
                    }
                    output[count] = sentinel;
                    count++;
                    continue;
                }
                v /= (float) cCount;
                output[count] = v;
                count++;
            }
        }

        IntensityDescriptor desc = new GsIntensityDescriptor(output, 
            centralPixelIndex, nCellsAcross);
        
        return desc;
    }
    
    /**
     * Determine whether to remove a feature that is difficult to localize.
     * The method follows Szeliski "Computer Vision: Algorithms and Applications" 
     * equation 4.11, (det(A)/trace(A)) > 10 which is the harmonic mean of
     * the auto-correlation matrix.  references Brown, Szeliski, and Winder, 2005.
     * 
     * The method extracts the same descriptor as used for auto-correlation 
     * to perform the trace and determinant on.
     * @param img
     * @param xCenter
     * @param yCenter
     * @param rotation
     */
    public boolean removeDueToLocalization(GreyscaleImage img,
        int xCenter, int yCenter, int rotation) {
        
        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }
     
        float sentinel = GsIntensityDescriptor.sentinel;
        
        GsIntensityDescriptor descriptor = (GsIntensityDescriptor) 
            extractIntensity(img, xCenter, yCenter, rotation);
        
        if (descriptor == null) {
            return true;
        }
                        
        SimpleMatrix m = createAutoCorrelationMatrix(descriptor);
        
        double det = m.determinant();
        double trace = m.trace();
        
        if (Math.abs(det/trace) < 10) {
            return true;
        }
        return false;
    }
    
    /**
     * calculate the orientation of the pixel at (x, y) using the gradient image
     * created with calculateGradientWithGreyscale.  The method finds the 
     * orientation of the local intensities with respect to the
     * given point in resolution of 45 degrees. The method finds the direction
     * of the largest difference from (x, y) and if negative that direction
     * is returned else the direction opposite of it.
       <pre>
                  90
           135    |    45
                  |
        180 ---------------  0
                  |
           225    |    315
                 270
       </pre>
       Note that this method uses
     * a cache to reuse calculations so it is up to the invoker to make
     * sure that the same image is given to this instance.  
     * (the image instance isn't retained to reduce references to potentially
     * large data structures which would otherwise be garbage collected).
     * @param x
     * @param y
     * @return
     */
    public Integer calculateOrientation(int x, int y) {
        
        if (gXY == null || theta == null) {
            throw new IllegalArgumentException("gradient image is null, so use"
                + " calculateGradientWithGreyscale to create it before using"
                + " this method");
        }
        
        checkBounds(gXY, x, y);
        
        PairInt p = new PairInt(x, y);
        
        Integer orientation = mapOfOrientation.get(p);
        
        if (orientation == null) {
            
            orientation = Integer.valueOf(findDominantOrientation(x, y));
            
            mapOfOrientation.put(p, orientation);
        }
                
        return orientation;
    }

    public RotatedOffsets getRotatedOffsets() {
        return rotatedOffsets;
    }

    private int findDominantOrientation(int xCenter, int yCenter) {
        
        // cell Dim affects angle resolution so cannot make this too small
        int cellDim = 4;//6;        
        int nCellsAcross = 2;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        
        List<Integer> angles = new ArrayList<Integer>();
        List<Integer> magnitudes = new ArrayList<Integer>();
        
        //Logger log = Logger.getLogger(this.getClass().getName());
        
        int idx = 0;
        for (int dy = -range0; dy < range0; dy += cellDim) {
            for (int dx = -range0; dx < range0; dx += cellDim) {
                
                //0 to 360 in bins of 20 degrees
                int[] hist = new int[18];
                Arrays.fill(hist, -1);
                
                int col = xCenter + dx;
                int row = yCenter + dy;
                
                populateOrientationHistogram(col, row, cellDim, hist);
                
                //log.info(String.format("(col, row)=(%5d,%5d)  %s", col, row,
                //    Arrays.toString(hist)));
                
                int[] indexes = new int[hist.length];
                for (int i = 0; i < hist.length; ++i) {
                    indexes[i] = i;
                }
                QuickSort.descendingSort(hist, indexes);
                float limit = (float)hist[0] * 0.8f;
                for (int i = 0; i < hist.length; ++i) {
                    int magnitude = hist[i];
                    if (magnitude >= limit) {
                        angles.add(Integer.valueOf(indexes[i] * 20));
                        magnitudes.add(Integer.valueOf(magnitude));
                    } else {
                        break;
                    }
                }
                
                idx++;
            }
        }
        assert(idx == (nCellsAcross * nCellsAcross));
        
        int[] a = new int[angles.size()];
        int[] w = new int[a.length];
        for (int i = 0; i < angles.size(); ++i) {
            a[i] = angles.get(i).intValue();
            w[i] = magnitudes.get(i).intValue();
            //log.info("angle=" + a[i]);
        }
        
        float avg = AngleUtil.calculateWeightedAverageWithQuadrantCorrections(
            a, w, a.length - 1);
        
        int avgInt = Math.round(avg);
        if (avgInt < 0) {
            avgInt += 360;
        } else if (avgInt > 359) {
            avgInt -= 360;
        }

        return avgInt;
    }
    
    /**
     * @param xStart
     * @param yStart
     * @param cellDimension
     * @param outputHistogram 
     */
    private void populateOrientationHistogram(int xStart, int yStart,
        int cellDimension, int[] outputHistogram) {
        
        int w = gXY.getWidth();
        int h = gXY.getHeight();
        
        int nBins = outputHistogram.length;
        int binSz = 360/nBins;
        
        for (int dxc = 0; dxc < cellDimension; ++dxc) {
            for (int dyc = 0; dyc < cellDimension; ++dyc) {
                int x = xStart + dxc;
                int y = yStart + dyc;
                if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                    continue;
                }
                int t = theta.getValue(x, y);
                int g = gXY.getValue(x, y);
                int hIdx = t / binSz;
                if (hIdx >= nBins) {
                    hIdx = nBins - 1;
                }
                outputHistogram[hIdx] += g;
            }
        }
    }
    
    SimpleMatrix createAutoCorrelationMatrix(IntensityDescriptor desc) {
        
        if (desc == null) {
            throw new IllegalArgumentException("desc cannot be null");
        }
        
        GsIntensityDescriptor descriptor = (GsIntensityDescriptor)desc;
     
        float sentinel = GsIntensityDescriptor.sentinel;
        
        int nCellsAcross = (int)(Math.sqrt(descriptor.grey.length));
        
        float vc = descriptor.grey[descriptor.getCentralIndex()];
        
        SimpleMatrix a = new SimpleMatrix(nCellsAcross, nCellsAcross);
        
        int idx = 0;
        for (int col = 0; col < nCellsAcross; ++col) {
            for (int row = 0; row < nCellsAcross; ++row) {
                float v = descriptor.grey[idx];
                if (v == sentinel) {
                    a.set(row, col, 0);
                } else {
                    a.set(row, col, v - vc);
                }
                idx++;
            }
        }
        
        return a;
    }

}
