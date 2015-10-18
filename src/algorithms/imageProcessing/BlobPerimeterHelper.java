package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    
    public BlobPerimeterHelper(final ImageExt img) {
        imgHelper = new SegmentedImageHelper(img);
    }
    
    public BlobPerimeterHelper(final ImageExt img, final String debugTag) {
        imgHelper = new SegmentedImageHelper(img, debugTag);
    }
    
    public void createBinnedGreyscaleImage(int maxDimension) {
        imgHelper.createBinnedGreyscaleImage(maxDimension);
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
                
        List<PairIntArray> blobPerimeters = segBlobPerimetersMap.get(type);
        
        if (blobPerimeters != null) {
            return blobPerimeters;
        }
        
        List<Set<PairInt>> blobs = segBlobsMap.get(type);
        
        boolean useBinned = false;
        
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned);
            
            segBlobsMap.put(type, blobs);
        }
        
        blobPerimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
            imgHelper, type, blobs, useBinned,
            discardWhenCavityIsSmallerThanBorder);
        
        segBlobPerimetersMap.put(type, blobPerimeters);
       
        return blobPerimeters;
    }

    private List<PairIntArray> getBinnedBlobPerimeters(SegmentationType type) {
        
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
        
        List<Set<PairInt>> blobs = segBlobsMap.get(type);
        
        boolean useBinned = false;
                
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned);
            
            segBlobsMap.put(type, blobs);
        }
        
        return blobs;
    }
}
