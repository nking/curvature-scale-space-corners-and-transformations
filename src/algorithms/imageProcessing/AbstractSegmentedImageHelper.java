package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractSegmentedImageHelper implements ISegmentedImageHelper {

    protected int smallestGroupLimit = 100;
    
    protected int largestGroupLimit = 5000;
    
    protected int smallestGroupLimitBinned = smallestGroupLimit;
    
    protected int largestGroupLimitBinned = largestGroupLimit;
    
    protected final ImageExt img;
    
    protected final GreyscaleImage imgGrey;
    
    /**
     * if true, was applied to imgGrey
     */
    protected boolean didApplyHistEq = false;
    
    /**
     * applied to imgGreyBinned and imgSegmentedMap
     */
    protected int binFactor = 1;
    
    protected GreyscaleImage imgGreyBinned = null;
    
    protected Map<SegmentationType, GreyscaleImage> imgSegmentedMap 
        = new HashMap<SegmentationType, GreyscaleImage>();
    
    protected Map<SegmentationType, GreyscaleImage> imgBinnedSegmentedMap 
        = new HashMap<SegmentationType, GreyscaleImage>();
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;
    
    protected String debugTag = "";

    public AbstractSegmentedImageHelper(final ImageExt img) {

        this.img = img;

        this.imgGrey = img.copyToGreyscale();

        this.binFactor = 1;
    }

    /**
     * constructor with a tag for debugging
     *
     * @param img
     * @param debugTag
     */
    public AbstractSegmentedImageHelper(final ImageExt img, final String debugTag) {

        this.img = img;

        this.imgGrey = img.copyToGreyscale();

        this.binFactor = 1;

        debug = true;

        this.debugTag = debugTag;
    }

    @Override
    public void createBinnedGreyscaleImage(int maxDimension) {
        
        if (imgGreyBinned != null) {
            return;
        }
        
        binFactor = (int) Math.ceil(Math.max((float) imgGrey.getWidth() / maxDimension, 
            (float) imgGrey.getHeight() / maxDimension));
        
        smallestGroupLimitBinned = smallestGroupLimit / (binFactor * binFactor);
        
        largestGroupLimitBinned = largestGroupLimit / (binFactor * binFactor);
        
        log.info("binFactor=" + binFactor);
        // prevent from being smaller than needed for a convex hull
        if (smallestGroupLimitBinned < 4) {
            smallestGroupLimitBinned = 4;
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imgGreyBinned = imageProcessor.binImage(imgGrey, binFactor);
    }

    @Override
    public boolean applyEqualizationIfNeededByComparison(ISegmentedImageHelper 
        otherImgHelper) {
        
        // doing this automatically for now, but can use the stats instead to
        // decide
        boolean performHistEq = false;
        boolean useStatsToDecide = true;
        
        if (useStatsToDecide) {
            ImageStatistics stats1 = ImageStatisticsHelper.examineImage(imgGrey, true);
            ImageStatistics stats2 = 
                ImageStatisticsHelper.examineImage(otherImgHelper.getGreyscaleImage(), true);
            double median1DivMedian2 = stats1.getMedian() / stats2.getMedian();
            double meanDivMedian1 = stats1.getMean() / stats1.getMedian();
            double meanDivMedian2 = stats2.getMean() / stats2.getMedian();
            if (((median1DivMedian2 > 1) && ((median1DivMedian2 - 1) > 0.2)) || 
                ((median1DivMedian2 < 1) && (median1DivMedian2 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian1 > 1) && ((meanDivMedian1 - 1) > 0.2)) || 
                ((meanDivMedian1 < 1) && (meanDivMedian1 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian2 > 1) && ((meanDivMedian2 - 1) > 0.2)) || 
                ((meanDivMedian2 < 1) && (meanDivMedian2 < 0.8))) {
                performHistEq = true;
            }
        }
        if (performHistEq) {
            applyEqualization();
            otherImgHelper.applyEqualization();
            didApplyHistEq = true;
        }
        return performHistEq;
    }

    public void applyEqualization() {
        if (!didApplyHistEq) {
            HistogramEqualization hEq = new HistogramEqualization(imgGrey);
            hEq.applyFilter();
            didApplyHistEq = true;
        }
    }

    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.
     *
     * @param type
     * @param applyToBinnedImage
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    @Override
    public void applySegmentation(SegmentationType type, 
        boolean applyToBinnedImage) throws IOException, NoSuchAlgorithmException {
        
        if (applyToBinnedImage) {
            applySegmentationToBinned(type);
        } else {
            applySegmentation(type);
        }
    }

    protected abstract void clearBinnedPointsOfInterestMaps();

    protected abstract void clearUnbinnedPointsOfInterestMaps();

    public abstract int sumPointsOfInterest(SegmentationType segmentationType1, 
        boolean useBinned1);
    
    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.
     *
     * @param type
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    protected void applySegmentationToBinned(SegmentationType type) throws
        IOException, NoSuchAlgorithmException {

        imgBinnedSegmentedMap.clear();

        clearBinnedPointsOfInterestMaps();

        if (type.equals(SegmentationType.BINARY)) {
            applySegmentationToBinned(type, 2);
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            applySegmentationToBinned(type, 2);
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            applySegmentationToBinned(type, -1);
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            applySegmentationToBinned(type, -1);
        }
    }

    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.
     *
     * @param type
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    protected void applySegmentation(SegmentationType type) throws IOException,
        NoSuchAlgorithmException {

        imgSegmentedMap.clear();

        clearUnbinnedPointsOfInterestMaps();

        if (type.equals(SegmentationType.BINARY)) {
            applySegmentation(type, 2);
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            applySegmentation(type, 2);
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            applySegmentation(type, -1);
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            applySegmentation(type, -1);
        }
    }

    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type. k is the number of segmentation
     * bins requested. Note that when type is BINARY, k is ignored. Also note
     * that currently, this class internally will only hold one binned image for
     * the segmentation type, so should be changed if the use changes to need
     * them saved by key k also.
     *
     * @param type
     * @param k
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    protected void applySegmentationToBinned(SegmentationType type, int k)
        throws IOException, NoSuchAlgorithmException {

        GreyscaleImage segImg = getBinnedSegmentationImage(type);
        if (segImg != null) {
            log.warning("segmentation was already applied.  error?");
            //return;
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        GreyscaleImage gsImg = getGreyscaleImageBinned();
        
        if (type.equals(SegmentationType.BINARY)) {
            
            int scl = 1; //20/binFactor;
            if (scl == 0) {
                scl = 1;
            }
            
            segImg = gsImg.copyImage();
            imageSegmentation.applyUsingKMPP(segImg, 2);
            imageProcessor.applyAdaptiveMeanThresholding(segImg, scl);
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            
            segImg = gsImg.copyImage();
            imageSegmentation.applyUsingKMPP(segImg, k);
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            
            // not expecting to use binned color images:
            ImageExt imgBinned = imageProcessor.binImage(img, binFactor);
            segImg = imageSegmentation.applyUsingPolarCIEXYAndFrequency(imgBinned, 0.2f, true);
            imageProcessor.applyAdaptiveMeanThresholding(segImg, 2);
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            
            // not expecting to use binned color images:
            ImageExt imgBinned = imageProcessor.binImage(img, binFactor);
            segImg = imageSegmentation.applyUsingPolarCIEXYAndFrequency(imgBinned, 0.2f, false);
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else {
            throw new UnsupportedOperationException("did not add impl for " 
                + type.name());
        }
        
        if (debug) {
            MiscDebug.writeImageCopy(gsImg, "segmented_" + debugTag + "_" 
                + MiscDebug.getCurrentTimeFormatted() + ".png");
        }
    }

    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type. k is the number of segmentation
     * bins requested. Note that when type is BINARY, k is ignored. Also note
     * that currently, this class internally will only hold one image for the
     * segmentation type, so should be changed if the use changes to need them
     * saved by key k also.
     *
     * @param type
     * @param k
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    protected void applySegmentation(SegmentationType type, final int k) 
        throws IOException, NoSuchAlgorithmException {
        
        GreyscaleImage segImg = getSegmentationImage(type);
        
        if (segImg != null) {
            log.warning("segmentation was already applied.  error?");
            //return;
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        if (type.equals(SegmentationType.BINARY)) {
            
            //TODO: replace w/ resolution based scl
            int scl = 1;
            segImg = imgGrey.copyImage();
            imageSegmentation.applyUsingKMPP(segImg, 2);
            imageProcessor.blur(segImg, SIGMA.ZEROPOINTSEVENONE);
            imageProcessor.applyAdaptiveMeanThresholding(segImg, scl);
            imgSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            
            segImg = imgGrey.copyImage();
            // expecting k=2
            imageSegmentation.applyUsingKMPP(segImg, k);
            imgSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            
            segImg = imageSegmentation.applyUsingPolarCIEXYAndFrequency(img, 0.2f, 
                true);
            imageProcessor.applyAdaptiveMeanThresholding(segImg, 2);
            imgSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            
            segImg = imageSegmentation.applyUsingPolarCIEXYAndFrequency(img, 0.2f, 
                true);
            imgSegmentedMap.put(type, segImg);
            
        } else {
            throw new UnsupportedOperationException("did not add impl for " 
                + type.name());
        }
        if (debug) {
            MiscDebug.writeImageCopy(segImg, "segmented_" + debugTag + "_" 
                + MiscDebug.getCurrentTimeFormatted() + ".png");
        }
    }

    protected abstract void generatePerimeterPointsOfInterestForUnbinned(
        SegmentationType type);

    protected abstract void generatePerimeterPointsOfInterestForBinned(
        SegmentationType type);

    @Override
    public GreyscaleImage getSegmentationImage(SegmentationType type) {
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        GreyscaleImage segImg = imgSegmentedMap.get(type);
        return segImg;
    }

    @Override
    public GreyscaleImage getBinnedSegmentationImage(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        GreyscaleImage segImg = imgBinnedSegmentedMap.get(type);
        
        return segImg;
    }

    @Override
    public GreyscaleImage getGreyscaleImage(boolean getTheBinned) {
        
        if (getTheBinned) {
            return getGreyscaleImageBinned();
        } else {
            return getGreyscaleImage();
        }
    }

    @Override
    public GreyscaleImage getGreyscaleImage() {
        return imgGrey;
    }
    
     @Override
    public ImageExt getImage() {
        return img;
    }

    @Override
    public GreyscaleImage getGreyscaleImageBinned() {
        
        //TODO: may want to return the unbinned image
        //for now throwing exception while testing
        if (imgGreyBinned == null) {
            throw new IllegalStateException(
            "binned image is null.  use createBinnedGreyscaleImage(...)");
        }
        
        return imgGreyBinned;
    }

    @Override
    public int getSmallestGroupLimit() {
        return smallestGroupLimit;
    }

    @Override
    public int getLargestGroupLimit() {
        return largestGroupLimit;
    }

    @Override
    public int getSmallestGroupLimitBinned() {
        return smallestGroupLimitBinned;
    }

    @Override
    public int getLargestGroupLimitBinned() {
        return largestGroupLimitBinned;
    }

    public int getBinFactor() {
        return binFactor;
    }

    /**
     * convenience method that returns '1' if getForBinned is false else returns
     * binFactor
     *
     * @param getForBinned
     * @return
     */
    @Override
    public int getBinFactor(boolean getForBinned) {
        return getForBinned ? binFactor : 1;
    }

}
