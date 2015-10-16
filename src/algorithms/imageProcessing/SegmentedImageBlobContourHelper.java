package algorithms.imageProcessing;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * a class to maintain the segmented images used by the blob scale finder.
 * it holds methods to perform the segmentation and also methods to 
 * extract blobs and their borders.
 * 
 * @author nichole
 */
public class SegmentedImageBlobContourHelper extends AbstractSegmentedImageHelper {
    
    protected Map<SegmentationType, BlobsAndContours> imgBlobsAndContoursMap 
        = new HashMap<SegmentationType, BlobsAndContours>();
    
    protected Map<SegmentationType, BlobsAndContours> imgBinnedBlobsAndContoursMap 
        = new HashMap<SegmentationType, BlobsAndContours>();
    
    public SegmentedImageBlobContourHelper(final ImageExt img) {
        
        super(img);
    }
    
    /**
     * constructor with a tag for debugging
     * @param img
     * @param debugTag 
     */
    public SegmentedImageBlobContourHelper(final ImageExt img, final String debugTag) {
        
        super(img, debugTag);
    }
    
    protected void clearBinnedPointsOfInterestMaps() {
        imgBinnedBlobsAndContoursMap.clear();
    }

    protected void clearUnbinnedPointsOfInterestMaps() {
        imgBlobsAndContoursMap.clear();
    }
    
    protected void generatePerimeterPointsOfInterestForUnbinned(SegmentationType type) {
        
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
            blobsAndContours = new BlobsAndContours(getGreyscaleImage(),
                segImg, smallestGroupLimit, largestGroupLimit, type, 
                segmentedToLineDrawing, debugTag);
        } else {
            blobsAndContours = new BlobsAndContours(getGreyscaleImage(),
                segImg, smallestGroupLimit, largestGroupLimit, type, 
                segmentedToLineDrawing);
        }
        
        imgBlobsAndContoursMap.put(type, blobsAndContours);
    }
    
    protected void generatePerimeterPointsOfInterestForBinned(SegmentationType type) {
        
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
            blobsAndContours = new BlobsAndContours(getGreyscaleImageBinned(),
                segImg, 
                smallestGroupLimitBinned, largestGroupLimitBinned, type,
                segmentedToLineDrawing, debugTag);
        } else {
            blobsAndContours = new BlobsAndContours(getGreyscaleImageBinned(),
                segImg, 
                smallestGroupLimitBinned, largestGroupLimitBinned, type,
                segmentedToLineDrawing);
        }
        
        imgBinnedBlobsAndContoursMap.put(type, blobsAndContours);
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
    
    @Override
    public int sumPointsOfInterest( 
        SegmentationType segmentationType, boolean useBinned) {
        
        BlobsAndContours bc = getBlobsAndContours(segmentationType, useBinned);
        
        int n = 0;
        
        for (List<CurvatureScaleSpaceContour> list : bc.getContours()) {
            n += list.size();
        }
        
        return n;
    }

}
