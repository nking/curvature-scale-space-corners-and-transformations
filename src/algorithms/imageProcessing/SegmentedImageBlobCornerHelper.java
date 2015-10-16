package algorithms.imageProcessing;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author nichole
 */
public class SegmentedImageBlobCornerHelper extends AbstractSegmentedImageHelper {
    
    protected Map<SegmentationType, BlobsAndCorners> imgBlobsAndCornersMap 
        = new HashMap<SegmentationType, BlobsAndCorners>();
    
    protected Map<SegmentationType, BlobsAndCorners> imgBinnedBlobsAndCornersMap 
        = new HashMap<SegmentationType, BlobsAndCorners>();
    
    public SegmentedImageBlobCornerHelper(ImageExt img) {
        super(img);
    }
    
    public SegmentedImageBlobCornerHelper(ImageExt img, String debugTag) {
        super(img, debugTag);
    }

    protected void clearBinnedPointsOfInterestMaps() {
        imgBinnedBlobsAndCornersMap.clear();
    }

    protected void clearUnbinnedPointsOfInterestMaps() {
        imgBlobsAndCornersMap.clear();
    }

    @Override
    protected void generatePerimeterPointsOfInterestForUnbinned(SegmentationType type) {
        GreyscaleImage segImg = getSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        BlobsAndCorners blobsAndCorners = imgBlobsAndCornersMap.get(type);

        if (blobsAndCorners != null) {
            return;
        }
        
        boolean segmentedToLineDrawing = false;
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT) ||
            type.equals(SegmentationType.BINARY)) {
            segmentedToLineDrawing = true;
        }
        
        if (debug) {
            blobsAndCorners = new BlobsAndCorners(getGreyscaleImage(),
                segImg, smallestGroupLimit, largestGroupLimit, type, 
                segmentedToLineDrawing, debugTag);
        } else {
            blobsAndCorners = new BlobsAndCorners(getGreyscaleImage(),
                segImg, smallestGroupLimit, largestGroupLimit, type, 
                segmentedToLineDrawing);
        }
        
        imgBlobsAndCornersMap.put(type, blobsAndCorners); 
    }

    @Override
    protected void generatePerimeterPointsOfInterestForBinned(SegmentationType type) {
        
        GreyscaleImage segImg = getBinnedSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        BlobsAndCorners blobsAndCorners = imgBinnedBlobsAndCornersMap.get(type);

        if (blobsAndCorners != null) {
            return;
        }
        
        boolean segmentedToLineDrawing = false;
        if (type.equals(SegmentationType.COLOR_POLARCIEXY_ADAPT) ||
            type.equals(SegmentationType.BINARY)) {
            segmentedToLineDrawing = true;
        }
        
        if (debug) {
            blobsAndCorners = new BlobsAndCorners(getGreyscaleImageBinned(),
                segImg, 
                smallestGroupLimitBinned, largestGroupLimitBinned, type,
                segmentedToLineDrawing, debugTag);
        } else {
            blobsAndCorners = new BlobsAndCorners(getGreyscaleImageBinned(),
                segImg, 
                smallestGroupLimitBinned, largestGroupLimitBinned, type,
                segmentedToLineDrawing);
        }
        
        imgBinnedBlobsAndCornersMap.put(type, blobsAndCorners);
    }

    public BlobsAndCorners getBlobsAndCorners(SegmentationType type, 
        boolean getTheBinned) {
        
        if (getTheBinned) {
            return getBlobsAndCornersForBinned(type);
        }
               
        return getBlobsAndCornersForUnbinned(type);
    }
    
    private BlobsAndCorners getBlobsAndCornersForUnbinned(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        BlobsAndCorners bc = imgBlobsAndCornersMap.get(type);
       
        return bc;
    }
    
    private BlobsAndCorners getBlobsAndCornersForBinned(SegmentationType type) {
        
        if (type == null) {
            throw new IllegalArgumentException("type cannot be null");
        }
        
        BlobsAndCorners bc = imgBinnedBlobsAndCornersMap.get(type);
       
        return bc;
    }
    
    @Override
    public int sumPointsOfInterest( 
        SegmentationType segmentationType, boolean useBinned) {
        
        BlobsAndCorners bc = getBlobsAndCorners(segmentationType, useBinned);
        
        int n = 0;
        
        for (List<CornerRegion> list : bc.getCorners()) {
            n += list.size();
        }
        
        return n;
    }

}
