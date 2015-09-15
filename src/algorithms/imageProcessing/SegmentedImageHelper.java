package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * a class to maintain the segmented images used by the blob scale finder.
 * it holds methods to perform the segmentation and also methods to 
 * extract blobs and their borders.
 * 
 * @author nichole
 */
public class SegmentedImageHelper {
    
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
    
    protected Map<SegmentationType, BlobsAndContours> imgBlobsAndContoursMap 
        = new HashMap<SegmentationType, BlobsAndContours>();
    
    protected Map<SegmentationType, BlobsAndContours> imgBinnedBlobsAndContoursMap 
        = new HashMap<SegmentationType, BlobsAndContours>();
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;
    
    protected String debugTag = "";
    
    /**
     * <pre>
     * BINARY = KMPP for k=2 followed by adaptive mean thresholding 
     *     (scl=20/binFactor);
     * GREYSCALE_KMPP = KMPP for k=2;
     * COLOR_POLARCIEXY_ADAPT = color to greyscale using polar theta of CIEXY 
     *     scale space followed by adaptive mean thresholding (scl=4?);
     * COLOR_POLARCIEXY = color to greyscale using polar theta of CIEXY.
     * </pre>
     */
    public enum SegmentationType {
        BINARY, 
        GREYSCALE_KMPP, 
        COLOR_POLARCIEXY_ADAPT, 
        COLOR_POLARCIEXY
        ;
    }
    
    public SegmentedImageHelper(final ImageExt img) {
        
        this.img = img;
        
        this.imgGrey = img.copyToGreyscale();
        
        this.binFactor = 1;
    }
    
    /**
     * constructor with a tag for debugging
     * @param img
     * @param debugTag 
     */
    public SegmentedImageHelper(final ImageExt img, final String debugTag) {
        
        this.img = img;
        
        this.imgGrey = img.copyToGreyscale();
        
        this.binFactor = 1;
        
        debug = true;
        
        this.debugTag = debugTag;
    }
    
    public void createBinnedGreyscaleImage(int maxDimension) {
        
        if (imgGreyBinned != null) {
            return;
        }
        
        binFactor = (int) Math.ceil(
            Math.max((float)imgGrey.getWidth()/maxDimension,
            (float)imgGrey.getHeight()/maxDimension));
        
        smallestGroupLimitBinned = smallestGroupLimit/(binFactor*binFactor);
        largestGroupLimitBinned = largestGroupLimit/(binFactor*binFactor);

        log.info("binFactor=" + binFactor);

        // prevent from being smaller than needed for a convex hull
        if (smallestGroupLimitBinned < 4) {
            smallestGroupLimitBinned = 4;
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        imgGreyBinned = imageProcessor.binImage(imgGrey, binFactor);
    }
    
    public boolean applyEqualizationIfNeededByComparison(
    SegmentedImageHelper otherImgHelper) {
        
        // doing this automatically for now, but can use the stats instead to
        // decide
        boolean performHistEq = false;

        boolean useStatsToDecide = true;

        if (useStatsToDecide) {
            ImageStatistics stats1 = ImageStatisticsHelper.examineImage(imgGrey,
                true);
            ImageStatistics stats2 = ImageStatisticsHelper.examineImage(
                otherImgHelper.getGreyscaleImage(), true);
            double median1DivMedian2 = stats1.getMedian() / stats2.getMedian();
            double meanDivMedian1 = stats1.getMean() / stats1.getMedian();
            double meanDivMedian2 = stats2.getMean() / stats2.getMedian();
            if (((median1DivMedian2 > 1) && ((median1DivMedian2 - 1) > 0.2))
                || ((median1DivMedian2 < 1) && (median1DivMedian2 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian1 > 1) && ((meanDivMedian1 - 1) > 0.2))
                || ((meanDivMedian1 < 1) && (meanDivMedian1 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian2 > 1) && ((meanDivMedian2 - 1) > 0.2))
                || ((meanDivMedian2 < 1) && (meanDivMedian2 < 0.8))) {
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
     * @param type
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public void applySegmentationToBinned(SegmentationType type, boolean applyToBinnedImage) 
        throws IOException, NoSuchAlgorithmException {
        
        if (applyToBinnedImage) {
            applySegmentationToBinned(type);
        } else {
            applySegmentation(type);
        }
    }
    
    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.
     * @param type
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    private void applySegmentationToBinned(SegmentationType type) 
        throws IOException, NoSuchAlgorithmException {
        
        GreyscaleImage segImg = getBinnedSegmentationImage(type);
        
        if (segImg != null) {
            log.warning("segmentation was already applied.  error?");
            return;
        }
        
        if (type.equals(SegmentationType.BINARY)) {
            
            applySegmentationToBinned(type, 2);
            
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            
            applySegmentationToBinned(type, 2);
        
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            
            applySegmentationToBinned(type, 4);
        
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            
            applySegmentationToBinned(type, 4);
        
        }
    }
    
    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.
     * @param type
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    private void applySegmentation(SegmentationType type) 
        throws IOException, NoSuchAlgorithmException {
        
        GreyscaleImage segImg = getSegmentationImage(type);
        
        if (segImg != null) {
            log.warning("segmentation was already applied.  error?");
            return;
        }
        
        if (type.equals(SegmentationType.BINARY)) {
            
            applySegmentation(type, 2);
            
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            
            applySegmentation(type, 2);
        
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            
            applySegmentation(type, 4);
        
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            
            applySegmentation(type, 4);
        
        }
    }
    
    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.  k is the number of segmentation
     * bins requested.  Note that when type is BINARY, k is ignored.
     * Also note that currently, this class internally will only hold one 
     * binned image for the segmentation type, so should be changed if the use 
     * changes to need them saved by key k also.
     * @param type
     * @param k
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    private void applySegmentationToBinned(SegmentationType type, int k) 
        throws IOException, NoSuchAlgorithmException {
        
        GreyscaleImage segImg = getBinnedSegmentationImage(type);
        
        if (segImg != null) {
            log.warning("segmentation was already applied.  error?");
            return;
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage gsImg = getGreyscaleImageBinned();

        if (type.equals(SegmentationType.BINARY)) {
                        
            int scl = 20/binFactor;
            
            segImg = gsImg.copyImage();
            
            imageProcessor.applyImageSegmentation(segImg, 2);
            // consider
            //imageProcessor.convertToCIEXYPolarTheta( color binned, 2);
            
            imageProcessor.applyAdaptiveMeanThresholding(segImg, scl);
            
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            
            segImg = gsImg.copyImage();
            
            imageProcessor.applyImageSegmentation(segImg, k);
            
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            
            // not expecting to use binned color images:
            ImageExt imgBinned = imageProcessor.binImage(img, binFactor);
            
            segImg = imageProcessor.createGreyscaleFromColorSegmentation(imgBinned, k);
            
            imageProcessor.applyAdaptiveMeanThresholding(segImg, 2);
            
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            
            // not expecting to use binned color images:
            ImageExt imgBinned = imageProcessor.binImage(img, binFactor);
            
            segImg = imageProcessor.createGreyscaleFromColorSegmentation(imgBinned, k);
            
            // also a possiblity:
            //= imageProcessor.createGreyscaleFromColorSegmentationKMPP(img, k, false);
            
            imgBinnedSegmentedMap.put(type, segImg);
            
        } else {
            throw new UnsupportedOperationException("did not add impl for " + type.name());
        }
        
        if (debug) {
            MiscDebug.writeImageCopy(gsImg, "segmented_" + debugTag + "_" 
            + MiscDebug.getCurrentTimeFormatted() + ".png");
        }
    }
    
    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.  k is the number of segmentation
     * bins requested.  Note that when type is BINARY, k is ignored.
     * Also note that currently, this class internally will only hold one 
     * image for the segmentation type, so should be changed if the use 
     * changes to need them saved by key k also.
     * @param type
     * @param k
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    private void applySegmentation(SegmentationType type, final int k) throws 
        IOException, NoSuchAlgorithmException {
        
        GreyscaleImage segImg = getSegmentationImage(type);
        
        if (segImg != null) {
            log.warning("segmentation was already applied.  error?");
            return;
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        
        if (type.equals(SegmentationType.BINARY)) {
            
            //TODO: replace w/ resolution based scl
            
            int scl = 20;
            
            segImg = imgGrey.copyImage();
            
imageProcessor.blur(segImg, SIGMA.ONE);

            imageProcessor.applyImageSegmentation(segImg, 2);
                        
            imageProcessor.applyAdaptiveMeanThresholding(segImg, scl);
            
            imgSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.GREYSCALE_KMPP)) {
            
            segImg = imgGrey.copyImage();
            
            // expecting k=2
            imageProcessor.applyImageSegmentation(segImg, k);
            
            imgSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT)) {
            
            segImg = imageProcessor.createGreyscaleFromColorSegmentation(img, k);
            
            imageProcessor.applyAdaptiveMeanThresholding(segImg, 2);
            
            imgSegmentedMap.put(type, segImg);
            
        } else if (type.equals(SegmentationType.COLOR_POLARCIEXY)) {
            
            segImg = imageProcessor.createGreyscaleFromColorSegmentation(img, k);
            
            // also a possiblity:
            //= imageProcessor.createGreyscaleFromColorSegmentationKMPP(img, k, false);
            
            imgSegmentedMap.put(type, segImg);
            
        } else {
            throw new UnsupportedOperationException("did not add impl for " + type.name());
        }
        
        if (debug) {
            MiscDebug.writeImageCopy(segImg, "segmented_" + debugTag + "_" 
            + MiscDebug.getCurrentTimeFormatted() + ".png");
        }
    }
    
    public void extractBlobsAndContours(SegmentationType type, boolean applyToBinnedImage) {
        
        if (applyToBinnedImage) {
            extractBlobsAndContoursForBinned(type);
        } else {
            extractBlobsAndContoursForUnbinned(type);
        }
    }
    
    private void extractBlobsAndContoursForUnbinned(SegmentationType type) {
        
        GreyscaleImage segImg = getSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        BlobsAndContours blobsAndContours = imgBlobsAndContoursMap.get(type);

        if (blobsAndContours != null) {
            return;
        }
        
        boolean segmentedToLineDrawing = false;
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT) ||
            type.equals(SegmentationType.BINARY)) {
            segmentedToLineDrawing = true;
        }
        
        if (debug) {
            blobsAndContours = new BlobsAndContours(segImg, smallestGroupLimit, 
                largestGroupLimit, type, segmentedToLineDrawing, debugTag);
        } else {
            blobsAndContours = new BlobsAndContours(segImg, smallestGroupLimit, 
                largestGroupLimit, type, segmentedToLineDrawing);
        }
        
        imgBlobsAndContoursMap.put(type, blobsAndContours);
    }
    
    private void extractBlobsAndContoursForBinned(SegmentationType type) {
        
        GreyscaleImage segImg = getBinnedSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        BlobsAndContours blobsAndContours = imgBinnedBlobsAndContoursMap.get(type);

        if (blobsAndContours != null) {
            return;
        }
        
        boolean segmentedToLineDrawing = false;
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT) ||
            type.equals(SegmentationType.BINARY)) {
            segmentedToLineDrawing = true;
        }
        
        if (debug) {
            blobsAndContours = new BlobsAndContours(segImg, 
                smallestGroupLimitBinned, largestGroupLimitBinned, type,
                segmentedToLineDrawing, debugTag);
        } else {
            blobsAndContours = new BlobsAndContours(segImg, 
                smallestGroupLimitBinned, largestGroupLimitBinned, type,
                segmentedToLineDrawing);
        }
        
        imgBinnedBlobsAndContoursMap.put(type, blobsAndContours);
    }
    
    public GreyscaleImage getSegmentationImage(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        GreyscaleImage segImg = imgSegmentedMap.get(type);
       
        return segImg;
    }
    
    public GreyscaleImage getBinnedSegmentationImage(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        GreyscaleImage segImg = imgBinnedSegmentedMap.get(type);
       
        return segImg;
    }
    
    public BlobsAndContours getBlobsAndContours(SegmentationType type, 
        boolean getTheBinned) {
        
        if (getTheBinned) {
            return getBlobsAndContoursForBinned(type);
        }
               
        return getBlobsAndContoursForUnbinned(type);
    }
    
    private BlobsAndContours getBlobsAndContoursForUnbinned(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        BlobsAndContours bc = imgBlobsAndContoursMap.get(type);
       
        return bc;
    }
    
    private BlobsAndContours getBlobsAndContoursForBinned(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        BlobsAndContours bc = imgBinnedBlobsAndContoursMap.get(type);
       
        return bc;
    }
    
    public GreyscaleImage getGreyscaleImage(boolean getTheBinned) {
        
        if (getTheBinned) {
            return getGreyscaleImageBinned();
        } else {
            return getGreyscaleImage();
        }
        
    }
    
    public GreyscaleImage getGreyscaleImage() {
        return imgGrey;
    }
    
    public GreyscaleImage getGreyscaleImageBinned() {
        
        //TODO: may want to return the unbinned image
        //for now throwing exception while testing
        if (imgGreyBinned == null) {
            throw new IllegalStateException(
            "binned image is null.  use createBinnedGreyscaleImage(...)");
        }
        
        return imgGreyBinned;
    }
    
    public int getSmallestGroupLimit() {
        return smallestGroupLimit;
    }
    
    public int getLargestGroupLimit() {
        return largestGroupLimit;
    }
    
    public int getSmallestGroupLimitBinned() {
        return smallestGroupLimitBinned;
    }
    
    public int getLargestGroupLimitBinned() {
        return largestGroupLimitBinned;
    }
    
    public int getBinFactor() {
        return binFactor;
    }
    
    /**
     * convenience method that returns '1' if getForBinned is false else returns
     * binFactor
     * @param getForBinned
     * @return 
     */
    public int getBinFactor(boolean getForBinned) {
        return getForBinned ? binFactor : 1;
    }
}
