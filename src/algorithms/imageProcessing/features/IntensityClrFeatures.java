package algorithms.imageProcessing.features;

import algorithms.imageProcessing.transform.Transformer;
import algorithms.QuickSort;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import static algorithms.imageProcessing.features.GsIntensityDescriptor.sentinel;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class IntensityClrFeatures {

    /**
     * the half width of a block.  For example, to extract the 25 pixels
     * centered on (xc, yc), bHalfW would be '2'.  a minimum of 5 seems to be
     * necessary.
     */
    protected final int bHalfW;
    
    protected final float[][] xyOffsets;
    protected final boolean useNormalizedIntensities = false;
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
        intensity1Blocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensity2Blocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensity3Blocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    /**
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
     *     extracted descriptor
     */
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityLBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityABlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityBBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
        
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
     * @param image greyscale image used to create the gradient and theta images
     * @param blockHalfWidths
     * @param rotatedOffsets 
     */
    public IntensityClrFeatures(final GreyscaleImage image, 
        final int blockHalfWidths, RotatedOffsets rotatedOffsets) {
        this.bHalfW = blockHalfWidths;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.rotatedOffsets = rotatedOffsets;
        calculateGradientWithGreyscale(image);
    }
    
    public GreyscaleImage getGradientImage() {
        return gXY;
    }

    /**
     * create the gradient of the given greyscale image in order to use it
     * internally for the method calculate45DegreeOrientation.
     * @param gsImg 
     */
    private void calculateGradientWithGreyscale(GreyscaleImage gsImg) {
        
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
     * extract the O1 intensity from the image for the given block center and
     * return it in a descriptor.
     * @param redImg
     * @param greenImg
     * @param blueImg
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensityO1(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation) {
     
        /*
         * extract the intensity from the image for the given block center and
         * return it in a descriptor.  Note, this method is a short cut for
         * use when the known internal type extract is by cells.
         * @param img greyscale image which needs a descriptor for given location
         * and orientation.  note that because the instance caches previous
         * calculations, the invoker must supply the same image each time.
         */
        
        if (redImg == null) {
            throw new IllegalStateException("redImg cannot be null");
        }
        if (greenImg == null) {
            throw new IllegalStateException("greenImg cannot be null");
        }
        if (redImg.getWidth() != greenImg.getWidth() || redImg.getHeight() != greenImg.getHeight()) {
            throw new IllegalStateException("redImg and greenImg must be the same size");
        }
                
        checkBounds(redImg, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        Map<Integer, IntensityDescriptor> descriptors = intensity1Blocks.get(p);
        
        IntensityDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, IntensityDescriptor>();
            intensity1Blocks.put(p, descriptors);
        }
        
        descriptor = extractO1IntensityForCells(redImg, greenImg, blueImg, 
            xCenter, yCenter, rotation);
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }
        
        return descriptor;
    }
    
    /**
     * extract the L intensity of CIE LAB color space from the images for the 
     * given block center and return it in a descriptor (note, if not already
     * cached, it computes L, A, and B descriptors and caches them all).
     * @param redImg
     * @param greenImg
     * @param blueImg
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensityLOfCIELAB(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation) {
     
        if (redImg == null) {
            throw new IllegalStateException("redImg cannot be null");
        }
        if (greenImg == null) {
            throw new IllegalStateException("greenImg cannot be null");
        }
        if (redImg.getWidth() != greenImg.getWidth() || redImg.getHeight() != greenImg.getHeight()) {
            throw new IllegalStateException("redImg and greenImg must be the same size");
        }
                
        checkBounds(redImg, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        // if this is present, then L, A, and B are all present, else all missing
        Map<Integer, IntensityDescriptor> descriptors = intensityLBlocks.get(p);
                
        if (descriptors != null) {
            IntensityDescriptor descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        }
        
        calculateAndCacheCIELABForCells(redImg, greenImg, blueImg, xCenter, 
            yCenter, rotation);
        
        descriptors = intensityLBlocks.get(p);
        if (descriptors == null) {
            return null;
        }
        
        // this is possibly null if near the edge of image
        IntensityDescriptor descriptor = descriptors.get(rotationKey);
        
        return descriptor;
    }
    
    /**
     * extract the A intensity of CIE LAB color space from the images for the 
     * given block center and return it in a descriptor (note, if not already
     * cached, it computes L, A, and B descriptors and caches them all).
     * @param redImg
     * @param greenImg
     * @param blueImg
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensityAOfCIELAB(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation) {
     
        if (redImg == null) {
            throw new IllegalStateException("redImg cannot be null");
        }
        if (greenImg == null) {
            throw new IllegalStateException("greenImg cannot be null");
        }
        if (redImg.getWidth() != greenImg.getWidth() || redImg.getHeight() != greenImg.getHeight()) {
            throw new IllegalStateException("redImg and greenImg must be the same size");
        }
                
        checkBounds(redImg, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        // if this is present, then L, A, and B are all present, else all missing
        Map<Integer, IntensityDescriptor> descriptors = intensityABlocks.get(p);
                
        if (descriptors != null) {
            IntensityDescriptor descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        }
        
        calculateAndCacheCIELABForCells(redImg, greenImg, blueImg, xCenter, 
            yCenter, rotation);
        
        descriptors = intensityABlocks.get(p);
        if (descriptors == null) {
            return null;
        }
        
        // this is possibly null if near the edge of image
        IntensityDescriptor descriptor = descriptors.get(rotationKey);
        
        return descriptor;
    }
    
    /**
     * extract the B intensity of CIE LAB color space from the images for the 
     * given block center and return it in a descriptor (note, if not already
     * cached, it computes L, A, and B descriptors and caches them all).
     * @param redImg
     * @param greenImg
     * @param blueImg
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensityBOfCIELAB(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation) {
     
        if (redImg == null) {
            throw new IllegalStateException("redImg cannot be null");
        }
        if (greenImg == null) {
            throw new IllegalStateException("greenImg cannot be null");
        }
        if (redImg.getWidth() != greenImg.getWidth() || redImg.getHeight() != greenImg.getHeight()) {
            throw new IllegalStateException("redImg and greenImg must be the same size");
        }
                
        checkBounds(redImg, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        // if this is present, then L, A, and B are all present, else all missing
        Map<Integer, IntensityDescriptor> descriptors = intensityBBlocks.get(p);
                
        if (descriptors != null) {
            IntensityDescriptor descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        }
        
        calculateAndCacheCIELABForCells(redImg, greenImg, blueImg, xCenter, 
            yCenter, rotation);
        
        descriptors = intensityBBlocks.get(p);
        if (descriptors == null) {
            return null;
        }
        
        // this is possibly null if near the edge of image
        IntensityDescriptor descriptor = descriptors.get(rotationKey);
        
        return descriptor;
    }
    
    /**
     * extract the O2 intensity from the image for the given block center and
     * return it in a descriptor.
     * @param redImg
     * @param greenImg
     * @param blueImg
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensityO2(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation) {
     
        /*
         * extract the intensity from the image for the given block center and
         * return it in a descriptor.  Note, this method is a short cut for
         * use when the known internal type extract is by cells.
         * @param img greyscale image which needs a descriptor for given location
         * and orientation.  note that because the instance caches previous
         * calculations, the invoker must supply the same image each time.
         */
        
        if (redImg == null) {
            throw new IllegalStateException("redImg cannot be null");
        }
        if (greenImg == null) {
            throw new IllegalStateException("greenImg cannot be null");
        }
        if (redImg.getWidth() != greenImg.getWidth() || redImg.getHeight() != greenImg.getHeight()) {
            throw new IllegalStateException("redImg and greenImg must be the same size");
        }
                
        checkBounds(redImg, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        Map<Integer, IntensityDescriptor> descriptors = intensity2Blocks.get(p);
        
        IntensityDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, IntensityDescriptor>();
            intensity2Blocks.put(p, descriptors);
        }
        
        descriptor = extractO2IntensityForCells(redImg, greenImg, blueImg, 
            xCenter, yCenter, rotation);
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
        }
        
        return descriptor;
    }
    
    /**
     * extract the O3 intensity from the image for the given block center and
     * return it in a descriptor.
     * @param redImg
     * @param greenImg
     * @param blueImg
     * @param xCenter
     * @param yCenter
     * @param rotation dominant orientation in degrees for feature at
     * (xCenter, yCenter)
     * @return
     */
    public IntensityDescriptor extractIntensityO3(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg,
        final int xCenter, final int yCenter, final int rotation) {
     
        /*
         * extract the intensity from the image for the given block center and
         * return it in a descriptor.  Note, this method is a short cut for
         * use when the known internal type extract is by cells.
         * @param img greyscale image which needs a descriptor for given location
         * and orientation.  note that because the instance caches previous
         * calculations, the invoker must supply the same image each time.
         */
        
        if (redImg == null) {
            throw new IllegalStateException("redImg cannot be null");
        }
        if (greenImg == null) {
            throw new IllegalStateException("greenImg cannot be null");
        }
        if (redImg.getWidth() != greenImg.getWidth() || redImg.getHeight() != greenImg.getHeight()) {
            throw new IllegalStateException("redImg and greenImg must be the same size");
        }
                
        checkBounds(redImg, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        Map<Integer, IntensityDescriptor> descriptors = intensity3Blocks.get(p);
        
        IntensityDescriptor descriptor = null;
        
        if (descriptors != null) {
            descriptor = descriptors.get(rotationKey);
            if (descriptor != null) {
                return descriptor;
            }
        } else {
            descriptors = new HashMap<Integer, IntensityDescriptor>();
            intensity3Blocks.put(p, descriptors);
        }
        
        descriptor = extractO3IntensityForCells(redImg, greenImg, blueImg, 
            xCenter, yCenter, rotation);
        
        if (useNormalizedIntensities && (descriptor != null)) {
            descriptor.applyNormalization();
        }
        
        if (descriptor != null) {
            descriptors.put(rotationKey, descriptor);
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
     * extract the O3 intensity from the images in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param img
     * @param xCenter
     * @param yCenter
     * @param rotation in degrees
     * @return
     */
    private IntensityDescriptor extractO3IntensityForCells(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg, int xCenter, 
        int yCenter, int rotation) {
        
        return extractOIntensityForCells(redImg, greenImg, blueImg, 
            xCenter, yCenter, rotation, 3);
    }
    
    /**
     * extract the O2 intensity from the images in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param img
     * @param xCenter
     * @param yCenter
     * @param rotation in degrees
     * @return
     */
    private IntensityDescriptor extractO2IntensityForCells(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg, int xCenter, 
        int yCenter, int rotation) {
        
        return extractOIntensityForCells(redImg, greenImg, blueImg, 
            xCenter, yCenter, rotation, 2);
    }
    
    /**
     * extract the O1 intensity from the images in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param img
     * @param xCenter
     * @param yCenter
     * @param rotation in degrees
     * @return
     */
    private IntensityDescriptor extractO1IntensityForCells(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg, int xCenter, 
        int yCenter, int rotation) {
        
        return extractOIntensityForCells(redImg, greenImg, blueImg, 
            xCenter, yCenter, rotation, 1);
    }
    
    /**
     * extract the intensity from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param redImg
     * @param otherImg either green or blue image depending upon oInt
     * @param xCenter
     * @param yCenter
     * @param rotation in degrees
     * @param oInt represents whether to calculate O1, O2, or O3
     * @return
     */
    private IntensityDescriptor extractOIntensityForCells(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg, int xCenter, 
        int yCenter, int rotation, int oInt) {
        
        if ((oInt < 1) || (oInt > 3)) {
            throw new IllegalArgumentException("expecting oInt to be 1, 2, or 3");
        }
        
        float sentinel = GsIntensityDescriptor.sentinel;
        
        int cellDim = 2;        
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        
        int[] xOffsets = rotatedOffsets.getXOffsets(rotation);
        int[] yOffsets = rotatedOffsets.getYOffsets(rotation);
        
        int w = redImg.getWidth();
        int h = redImg.getHeight();

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
                    if (oInt == 1) {
                        // O1 is (R-G)/sqrt(2)
                        v += (redImg.getValue(x, y) - greenImg.getValue(x, y));
                        v /= Math.sqrt(2.);
                    } else if (oInt == 2) {
                        // (R+G-2B)/sqrt(6)
                        v += (redImg.getValue(x, y) + greenImg.getValue(x, y) -
                            2*blueImg.getValue(x, y));
                        v /= Math.sqrt(6.);
                    } else {
                        // (R+G+B)/sqrt(2)
                        v += (redImg.getValue(x, y) + greenImg.getValue(x, y) +
                            blueImg.getValue(x, y));
                        v /= Math.sqrt(3.);
                    }
                    
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
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    public Integer calculateOrientation(int x, int y) throws 
        CornerRegion.CornerRegionDegneracyException {
        
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

    /**
     * Determine whether to remove a feature that is difficult to localize using
     * the O3 color feature.
     * The method follows Szeliski "Computer Vision: Algorithms and Applications" 
     * equation 4.11, (det(A)/trace(A)) > 10 which is the harmonic mean of
     * the auto-correlation matrix.  references Brown, Szeliski, and Winder, 2005.
     * 
     * The method extracts the same descriptor as used for auto-correlation 
     * to perform the trace and determinant on.
     * @param rImg
     * @param gImg
     * @param bImg
     * @param xCenter
     * @param yCenter
     * @param rotation
     * @return 
     */
    public boolean removeDueToLocalization(GreyscaleImage rImg, 
        GreyscaleImage gImg, GreyscaleImage bImg, int xCenter, int yCenter, 
        int rotation) {
        
        IntensityDescriptor descriptor0 = extractIntensityLOfCIELAB(
            rImg, gImg, bImg, xCenter, yCenter, rotation);
        
        if (descriptor0 == null) {
            return true;
        }
             
        SimpleMatrix m = createAutoCorrelationMatrix(descriptor0);
        
        double det = m.determinant();
        double trace = m.trace();
        //Logger.getLogger(this.getClass().getName()).info(String.format("(%d,%d) %.1f",
        //xCenter, yCenter, Math.abs(det/trace)));
        if (Math.abs(det/trace) < 10) {
            return true;
        }
        return false;
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
    
    public int getDefaultFinalDescriptorSize() {
        int nCellsAcross = getDefaultNCellsAcrossForExtract();
        return nCellsAcross * nCellsAcross;
    }

    private void calculateAndCacheCIELABForCells(GreyscaleImage redImg, 
        GreyscaleImage greenImg, GreyscaleImage blueImg, 
        int xCenter, int yCenter, int rotation) {
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float sentinel = GsIntensityDescriptor.sentinel;
        
        int cellDim = 2;        
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        
        int[] xOffsets = rotatedOffsets.getXOffsets(rotation);
        int[] yOffsets = rotatedOffsets.getYOffsets(rotation);
        
        int w = redImg.getWidth();
        int h = redImg.getHeight();

        float[] outputL = new float[nCellsAcross * nCellsAcross];
        float[] outputA = new float[nCellsAcross * nCellsAcross];
        float[] outputB = new float[nCellsAcross * nCellsAcross];

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
                float lSum = 0;
                float aSum = 0;
                float bSum = 0;
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
                    int r = redImg.getValue(x, y);
                    int g = greenImg.getValue(x, y);
                    int b = blueImg.getValue(x, y);
                    
                    float[] lab = cieC.rgbToCIELAB(r, g, b);
                    
                    lSum += lab[0];
                    aSum += lab[1];
                    bSum += lab[2];
                    
                    cCount++;
                }
                if (!withinBounds || (cCount == 0)) {
                    if (count == centralPixelIndex) {
                        return;
                    }
                    outputL[count] = sentinel;
                    outputA[count] = sentinel;
                    outputB[count] = sentinel;
                    count++;
                    continue;
                }
                lSum /= (float) cCount;
                aSum /= (float) cCount;
                bSum /= (float) cCount;
                outputL[count] = lSum;
                outputA[count] = aSum;
                outputB[count] = bSum;
                count++;
            }
        }

        PairInt p = new PairInt(xCenter, yCenter);
        Integer rotationKey = Integer.valueOf(rotation);
        
        IntensityDescriptor descL = new GsIntensityDescriptor(outputL, 
            centralPixelIndex, nCellsAcross);
        IntensityDescriptor descA = new GsIntensityDescriptor(outputA, 
            centralPixelIndex, nCellsAcross);
        IntensityDescriptor descB = new GsIntensityDescriptor(outputB, 
            centralPixelIndex, nCellsAcross);
        
        Map<Integer, IntensityDescriptor> dL = new HashMap<Integer, IntensityDescriptor>();
        dL.put(rotationKey, descL);
        this.intensityLBlocks.put(p, dL);
        
        Map<Integer, IntensityDescriptor> dA = new HashMap<Integer, IntensityDescriptor>();
        dA.put(rotationKey, descA);
        this.intensityABlocks.put(p, dA);
        
        Map<Integer, IntensityDescriptor> dB = new HashMap<Integer, IntensityDescriptor>();
        dB.put(rotationKey, descB);
        this.intensityBBlocks.put(p, dB);
    }
    
    public static FeatureComparisonStat calculateStats(IntensityDescriptor 
        desc1_l, IntensityDescriptor desc1_a, IntensityDescriptor desc1_b, 
        int x1, int y1, 
        IntensityDescriptor desc2_l, IntensityDescriptor desc2_a, 
        IntensityDescriptor desc2_b, int x2, int y2) {
        
        if (desc1_l == null) {
            throw new IllegalArgumentException("desc1_l cannot be null");
        }
        if (desc1_a == null) {
            throw new IllegalArgumentException("desc1_a cannot be null");
        }
        if (desc1_b == null) {
            throw new IllegalArgumentException("desc1_b cannot be null");
        }
        if (desc2_l == null) {
            throw new IllegalArgumentException("desc2_l cannot be null");
        }
        if (desc2_a == null) {
            throw new IllegalArgumentException("desc2_a cannot be null");
        }
        if (desc2_b == null) {
            throw new IllegalArgumentException("desc2_b cannot be null");
        }
        
        float[] l1 = ((GsIntensityDescriptor)desc1_l).grey;
        float[] a1 = ((GsIntensityDescriptor)desc1_a).grey;
        float[] b1 = ((GsIntensityDescriptor)desc1_b).grey;
        
        float[] l2 = ((GsIntensityDescriptor)desc2_l).grey;
        float[] a2 = ((GsIntensityDescriptor)desc2_a).grey;
        float[] b2 = ((GsIntensityDescriptor)desc2_b).grey;
        
        int n = l1.length;
        int centralPixIdx1 = desc1_l.getCentralIndex();
        int centralPixIdx2 = desc2_l.getCentralIndex();
        
        float vcL1 = l1[centralPixIdx1];
        if (vcL1 == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        float vcA1 = a1[centralPixIdx1];
        float vcB1 = b1[centralPixIdx1];
        
        float vcL2 = l2[centralPixIdx2];
        if (vcL2 == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        float vcA2 = a2[centralPixIdx2];
        float vcB2 = b2[centralPixIdx2];
        
        CIEChromaticity cieC = new CIEChromaticity();
                
        int count = 0;
        double autoCorrel = 0;
        double deltaESum = 0;
        
        for (int i = 0; i < n; ++i) {
            
            float vL1 = l1[i];
            if (vL1 == sentinel) {
                continue;
            }
            float vL2 = l2[i];
            if (vL2 == sentinel) {
                continue;
            }
            float vA1 = a1[i];
            float vB1 = b1[i];
            float vA2 = a2[i];
            float vB2 = b2[i];
            
            double deltaE = cieC.calcDeltaECIE94(vL1, vA1, vB1, vL2, vA2, vB2);
            
            deltaESum += (deltaE * deltaE);
            
            double deltaEC = cieC.calcDeltaECIE94(vcL2, vcA2, vcB2, vL2, vA2, vB2);
                    
            autoCorrel += (deltaEC * deltaEC);            
            
            count++;
        }
        
        autoCorrel /= (double)count;
        
        deltaESum /= (double)count;
         
        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumIntensitySqDiff((float)deltaESum);
        stat.setImg2PointIntensityErr((float)autoCorrel);
        
        return stat;
    }
    
    public static FeatureComparisonStat calculateHalfStats(
        IntensityDescriptor desc1_l, IntensityDescriptor desc1_a, 
        IntensityDescriptor desc1_b, int x1, int y1, boolean useTop1, 
        IntensityDescriptor desc2_l, IntensityDescriptor desc2_a, 
        IntensityDescriptor desc2_b, int x2, int y2, boolean useTop2) {
        
        if (desc1_l == null) {
            throw new IllegalArgumentException("desc1_l cannot be null");
        }
        if (desc1_a == null) {
            throw new IllegalArgumentException("desc1_a cannot be null");
        }
        if (desc1_b == null) {
            throw new IllegalArgumentException("desc1_b cannot be null");
        }
        if (desc2_l == null) {
            throw new IllegalArgumentException("desc2_l cannot be null");
        }
        if (desc2_a == null) {
            throw new IllegalArgumentException("desc2_a cannot be null");
        }
        if (desc2_b == null) {
            throw new IllegalArgumentException("desc2_b cannot be null");
        }
        
        int[] indexes1 = useTop1 ? 
            ((GsIntensityDescriptor)desc1_l).getUpperHalfIndexes() :
            ((GsIntensityDescriptor)desc1_l).getLowerHalfIndexes();
        
        int[] indexes2 = useTop2 ? 
            ((GsIntensityDescriptor)desc2_l).getUpperHalfIndexes() :
            ((GsIntensityDescriptor)desc2_l).getLowerHalfIndexes();
        
        if (indexes1.length != indexes2.length) {
            throw new IllegalArgumentException(
            "indexes1 and indexes2 must be same length");
        }
        
        float[] l1 = ((GsIntensityDescriptor)desc1_l).grey;
        float[] a1 = ((GsIntensityDescriptor)desc1_a).grey;
        float[] b1 = ((GsIntensityDescriptor)desc1_b).grey;
        
        float[] l2 = ((GsIntensityDescriptor)desc2_l).grey;
        float[] a2 = ((GsIntensityDescriptor)desc2_a).grey;
        float[] b2 = ((GsIntensityDescriptor)desc2_b).grey;
        
        int n = l1.length;
        int centralPixIdx1 = desc1_l.getCentralIndex();
        int centralPixIdx2 = desc2_l.getCentralIndex();
                
        float vcL1 = l1[centralPixIdx1];
        if (vcL1 == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        //float vcA1 = a1[centralPixIdx1];
        //float vcB1 = b1[centralPixIdx1];
        
        float vcL2 = l2[centralPixIdx2];
        if (vcL2 == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        float vcA2 = a2[centralPixIdx2];
        float vcB2 = b2[centralPixIdx2];
        
        assert(indexes1.length == indexes2.length);
        
        CIEChromaticity cieC = new CIEChromaticity();
                
        int count = 0;
        //TODO: review the math for auto-correlation here since it is estimated from 2 points instead of 1
        double autoCorrel = 0;
        double deltaESum = 0;
        
        n = indexes1.length;
        
        for (int i = 0; i < n; ++i) {
            
            int idx1 = indexes1[i];
            int idx2 = indexes2[i];
            
            float vL1 = l1[idx1];
            if (vL1 == sentinel) {
                continue;
            }
            float vL2 = l2[idx2];
            if (vL2 == sentinel) {
                continue;
            }
            float vA1 = a1[idx1];
            float vB1 = b1[idx1];
            float vA2 = a2[idx2];
            float vB2 = b2[idx2];
            
            double deltaE = cieC.calcDeltaECIE94(vL1, vA1, vB1, vL2, vA2, vB2);
            
            deltaESum += (deltaE * deltaE);
            
            double deltaEC = cieC.calcDeltaECIE94(vcL2, vcA2, vcB2, vL2, vA2, vB2);
                    
            autoCorrel += (deltaEC * deltaEC);            
            
            count++;
        }
        
        autoCorrel /= (double)count;
        
        deltaESum /= (double)count;
         
        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumIntensitySqDiff((float)deltaESum);
        stat.setImg2PointIntensityErr((float)autoCorrel);
        
        return stat;
    }
}
