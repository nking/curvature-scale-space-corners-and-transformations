package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class BlobPerimeterHelper {
    
    private final SegmentedImageHelper imgHelper;
    
    private final Map<SegmentationType, List<Set<PairInt>>> segBlobsMap =      
        new HashMap<SegmentationType, List<Set<PairInt>>>();
   
    private final Map<SegmentationType, List<PairIntArray>> 
        segBlobPerimetersMap = new HashMap<SegmentationType, List<PairIntArray>>();
    
    private final Map<SegmentationType, List<Set<PairInt>>> segBinnedBlobsMap =      
        new HashMap<SegmentationType, List<Set<PairInt>>>();
   
    private final Map<SegmentationType, List<PairIntArray>> 
        segBinnedBlobPerimetersMap = new HashMap<SegmentationType, List<PairIntArray>>();
    
    private final Map<SegmentationType, List<Set<PairInt>>> segLargeBlobsMap =      
        new HashMap<SegmentationType, List<Set<PairInt>>>();
   
    private final Map<SegmentationType, List<PairIntArray>> 
        segLargeBlobPerimetersMap = new HashMap<SegmentationType, List<PairIntArray>>();
    
    public BlobPerimeterHelper(final ImageExt img) {
        imgHelper = new SegmentedImageHelper(img);
    }
    
    public BlobPerimeterHelper(final ImageExt img, final String debugTag) {
        imgHelper = new SegmentedImageHelper(img, debugTag);
    }
    
    public void createBinnedGreyscaleImage(int maxDimension) {
        imgHelper.createBinnedGreyscaleImage(maxDimension);
    }
    
    /**
     * change the larger limit for blob sizes.  caveat is that this should be
     * invoked before createBinnedGreyscaleImage for it to be adapted correctly
     * for the binned image sizes.
     * 
     * @param limit 
     */
    public void increaseLargestGroupLimit(int limit) {
        imgHelper.increaseLargestGroupLimit(limit);
    }
    
    /**
     * change the smallest limit for blob sizes.  caveat is that this should be
     * invoked before createBinnedGreyscaleImage for it to be adapted correctly
     * for the binned image sizes.
     * 
     * @param limit 
     */
    public void increaseSmallestGroupLimit(int limit) {
        imgHelper.increaseSmallestGroupLimit(limit);
    }

    public boolean applyEqualizationIfNeededByComparison(SegmentedImageHelper 
        otherImgHelper) {
        
        return imgHelper.applyEqualizationIfNeededByComparison(otherImgHelper);
    }

    public void applyEqualization() {
        imgHelper.applyEqualization();
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
    public void applySegmentation(SegmentationType type, 
        boolean applyToBinnedImage) throws IOException, NoSuchAlgorithmException {
        
        imgHelper.applySegmentation(type, applyToBinnedImage);
    }

    public GreyscaleImage getSegmentationImage(SegmentationType type) {
        return imgHelper.getSegmentationImage(type);
    }

    public GreyscaleImage getBinnedSegmentationImage(SegmentationType type) {
        return imgHelper.getBinnedSegmentationImage(type);
    }

    public GreyscaleImage getGreyscaleImage(boolean getTheBinned) {
        return imgHelper.getGreyscaleImage(getTheBinned);
    }

    public GreyscaleImage getGreyscaleImage() {
        return imgHelper.getGreyscaleImage();
    }
    
    public ImageExt getImage() {
        return imgHelper.getImage();
    }

    public GreyscaleImage getGreyscaleImageBinned() {
        return imgHelper.getGreyscaleImageBinned();
    }

    public int getSmallestGroupLimit() {
        return imgHelper.getSmallestGroupLimit();
    }

    public int getLargestGroupLimit() {
        return imgHelper.getLargestGroupLimit();
    }

    public int getSmallestGroupLimitBinned() {
        return imgHelper.getSmallestGroupLimitBinned();
    }

    public int getLargestGroupLimitBinned() {
        return imgHelper.getLargestGroupLimitBinned();
    }

    public int getBinFactor() {
        return imgHelper.getBinFactor();
    }

    /**
     * convenience method that returns '1' if getForBinned is false else returns
     * binFactor
     *
     * @param getForBinned
     * @return
     */
    public int getBinFactor(boolean getForBinned) {
        return imgHelper.getBinFactor(getForBinned);
    }

    public List<PairIntArray> getBlobPerimeters(SegmentationType type, boolean 
        useBinnedImage) {
                
        if (useBinnedImage) {
            return getBinnedBlobPerimeters(type);
        } else {
            return getUnbinnedBlobPerimeters(type);
        }
       
    }

    private List<PairIntArray> getUnbinnedBlobPerimeters(SegmentationType type) {
                
        List<PairIntArray> blobPerimeters = null;
        
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            blobPerimeters = segLargeBlobPerimetersMap.get(type);
        } else {
            blobPerimeters = segBlobPerimetersMap.get(type);
        }
        
        if (blobPerimeters != null) {
            return blobPerimeters;
        }
        
        List<Set<PairInt>> blobs = null;
        
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            blobs = segLargeBlobsMap.get(type);
        } else {
            blobs = segBlobsMap.get(type);
        }
        
        boolean useBinned = false;
        
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        if (blobs == null) {
            
            long t0 = System.currentTimeMillis();
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned);
            
            long t1 = System.currentTimeMillis();
            long t1Sec = (t1 - t0)/1000;
            Logger.getLogger(this.getClass().getName()).info("blobs(sec)=" 
                + t1Sec);
            
            if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
                segLargeBlobsMap.put(type, blobs);
            } else {
                segBlobsMap.put(type, blobs);
            }
        }
        
        long t0 = System.currentTimeMillis();
        
        blobPerimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
            imgHelper, type, blobs, useBinned,
            discardWhenCavityIsSmallerThanBorder);
        
        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        Logger.getLogger(this.getClass().getName()).info("perimeters(sec)=" 
            + t1Sec);
        
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            segLargeBlobPerimetersMap.put(type, blobPerimeters);
        
            // blobs may have been removed, so re-put them
            segLargeBlobsMap.put(type, blobs);
        } else {
            segBlobPerimetersMap.put(type, blobPerimeters);
        
            // blobs may have been removed, so re-put them
            segBlobsMap.put(type, blobs);
        }
        
        return blobPerimeters;
    }

    private List<PairIntArray> getBinnedBlobPerimeters(SegmentationType type) {
        
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            throw new IllegalStateException(
            "Cannot use segmentation type COLOR_POLARCIEXY_LARGE with binned preference");
        }
        
        List<PairIntArray> blobPerimeters = segBinnedBlobPerimetersMap.get(type);
        
        if (blobPerimeters != null) {
            return blobPerimeters;
        }
        
        List<Set<PairInt>> blobs = segBinnedBlobsMap.get(type);
        
        boolean useBinned = true;
        
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned);
            
            segBinnedBlobsMap.put(type, blobs);
        }
        
        blobPerimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
            imgHelper, type, blobs, useBinned,
            discardWhenCavityIsSmallerThanBorder);
        
        segBinnedBlobPerimetersMap.put(type, blobPerimeters);
        
        // blobs may have been removed, so re-put them
        segBinnedBlobsMap.put(type, blobs);
            
        return blobPerimeters;        
    }
    
    public List<Set<PairInt>> getBlobs(SegmentationType type, boolean 
        useBinnedImage) {
                
        if (useBinnedImage) {
            return getBinnedBlobs(type);
        } else {
            return getUnbinnedBlobs(type);
        }
       
    }

    private List<Set<PairInt>> getBinnedBlobs(SegmentationType type) {
        
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            throw new IllegalStateException(
            "Cannot use segmentation type COLOR_POLARCIEXY_LARGE with binned preference");
        }
        
        List<Set<PairInt>> blobs = segBinnedBlobsMap.get(type);
        
        boolean useBinned = true;
                
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned);
            
            segBinnedBlobsMap.put(type, blobs);
        }
        
        return blobs;            
    }

    private List<Set<PairInt>> getUnbinnedBlobs(SegmentationType type) {
        
        List<Set<PairInt>> blobs = null;
        
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
            blobs = segLargeBlobsMap.get(type);
        } else {
            blobs = segBlobsMap.get(type);
        }
        
        boolean useBinned = false;
                
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned);
            
            if (type.equals(SegmentationType.COLOR_POLARCIEXY_LARGE)) {
                segLargeBlobsMap.put(type, blobs);
            } else {
                segBlobsMap.put(type, blobs);
            }
        }
        
        return blobs;
    }
    
    public boolean isInDebugMode() {
        return imgHelper.isInDebugMode();
    }

    public String getDebugTag() {
        return imgHelper.getDebugTag();
    }
}
