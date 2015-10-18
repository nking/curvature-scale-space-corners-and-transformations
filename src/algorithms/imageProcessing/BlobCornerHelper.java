package algorithms.imageProcessing;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author nichole
 */
public class BlobCornerHelper {
    
    protected Map<SegmentationType, List<List<CornerRegion>>>
        segCornersMap = new HashMap<SegmentationType, 
        List<List<CornerRegion>>>();
    
    protected Map<SegmentationType, List<List<CornerRegion>>>
        segBinnedCornersMap = new HashMap<SegmentationType, 
        List<List<CornerRegion>>>();
    
    protected final BlobPerimeterHelper imgHelper;
    
    protected boolean debug = false;
    
    protected String debugTag = "";
    
    public BlobCornerHelper(final BlobPerimeterHelper imgHelper) {
        
        this.imgHelper = imgHelper;
    }
    
    public BlobCornerHelper(final BlobPerimeterHelper imgHelper, 
        final String debugTag) {
        
        this.imgHelper = imgHelper;
        
        debug = false;
        
        this.debugTag = debugTag;
    }
    
    /**
     * generate contours (or retrieve cached contours) for blob perimeters from 
     * segmented images of type type and binned or not binned.
     * @param type
     * @param applyToBinnedImage 
     * @return  
     */
    public List<List<CornerRegion>> generatePerimeterCorners(
        SegmentationType type, boolean applyToBinnedImage) {
         
        if (applyToBinnedImage) {
            return generatePerimeterCornersForBinned(type);
        } else {
            return generatePerimeterCornersUnbinned(type);
        }
    }
    
    public List<List<CornerRegion>> getPerimeterCorners(
        SegmentationType type, boolean applyToBinnedImage) {
         
        return generatePerimeterCorners(type, applyToBinnedImage);
    }
    
    protected List<List<CornerRegion>> 
        generatePerimeterCornersUnbinned(SegmentationType type) {
                
        GreyscaleImage segImg = imgHelper.getSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException(
            "segmented image hasn't been created yet.  error?");
        }
        
        List<List<CornerRegion>> contours = segCornersMap.get(type);

        if (contours != null) {
            return contours;
        }
        
        boolean useBinned = false;
        
        final boolean outdoorMode = false;
        
        final boolean enableJaggedLineCorrections = false;
        
        final float factorIncreaseForCurvatureMinimum = 1.f;
        
        contours = BlobsAndCorners.populateCorners(imgHelper, type, useBinned,
            outdoorMode, enableJaggedLineCorrections, 
            factorIncreaseForCurvatureMinimum);
       
        segCornersMap.put(type, contours);
        
        return contours;
    }
    
    protected List<List<CornerRegion>> 
        generatePerimeterCornersForBinned(SegmentationType type) {
        
        segBinnedCornersMap.clear();
        
        GreyscaleImage segImg = imgHelper.getBinnedSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException(
            "segmented image hasn't been created yet.  error?");
        }
        
        List<List<CornerRegion>> contours = segBinnedCornersMap.get(type);

        if (contours != null) {
            return contours;
        }
        
        boolean useBinned = true;
        
        final boolean outdoorMode = false;
        
        final boolean enableJaggedLineCorrections = false;
        
        final float factorIncreaseForCurvatureMinimum = 1.f;
        
        contours = BlobsAndCorners.populateCorners(imgHelper, type, useBinned,
            outdoorMode, enableJaggedLineCorrections, 
            factorIncreaseForCurvatureMinimum);
        
        segBinnedCornersMap.put(type, contours);
        
        return contours;        
    }
    
    public int sumPointsOfInterest(SegmentationType segmentationType, 
        boolean useBinned) {
        
        List<List<CornerRegion>> corners = getPerimeterCorners(segmentationType, 
            useBinned);
        
        int n = 0;
        
        for (List<CornerRegion> list : corners) {
            n += list.size();
        }
        
        return n;
    }

}
