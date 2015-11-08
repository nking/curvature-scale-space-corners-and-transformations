package algorithms.imageProcessing;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * 
 * @author nichole
 */
public class BlobContourHelper {
    
    protected Map<SegmentationType, List<List<CurvatureScaleSpaceContour>>> 
        segContoursMap = new HashMap<SegmentationType, 
        List<List<CurvatureScaleSpaceContour>>>();
    
    protected Map<SegmentationType, List<List<CurvatureScaleSpaceContour>>>
        segBinnedContoursMap = new HashMap<SegmentationType, 
        List<List<CurvatureScaleSpaceContour>>>();
    
    protected Map<SegmentationType, List<List<BlobPerimeterRegion>>> 
        segBPRMap = new HashMap<SegmentationType, 
        List<List<BlobPerimeterRegion>>>();
    
    protected Map<SegmentationType, List<List<BlobPerimeterRegion>>>
        segBinnedBPRMap = new HashMap<SegmentationType, 
        List<List<BlobPerimeterRegion>>>();
    
    protected final BlobPerimeterHelper imgHelper;
    
    protected boolean debug = false;
    
    protected String debugTag = "";
    
    public BlobContourHelper(final BlobPerimeterHelper imgHelper) {
        
        this.imgHelper = imgHelper;
    }
    
    public BlobContourHelper(final BlobPerimeterHelper imgHelper, 
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
    public List<List<CurvatureScaleSpaceContour>> generatePerimeterContours(
        SegmentationType type, boolean applyToBinnedImage) {
         
        if (applyToBinnedImage) {
            return generatePerimeterContoursForBinned(type);
        } else {
            return generatePerimeterContoursUnbinned(type);
        }
    }
    
    /**
     * generate contours (or retrieve cached contours) for blob perimeters from 
     * segmented images of type type and binned or not binned.
     * @param type
     * @param applyToBinnedImage 
     * @return  
     */
    public List<List<CurvatureScaleSpaceContour>> getPerimeterContours(
        SegmentationType type, boolean applyToBinnedImage) {
         
        return generatePerimeterContours(type, applyToBinnedImage);
    }
    
    protected List<List<CurvatureScaleSpaceContour>> 
        generatePerimeterContoursUnbinned(SegmentationType type) {
            
        List<List<CurvatureScaleSpaceContour>> contours = segContoursMap.get(type);

        if (contours != null) {
            return contours;
        }
        
        GreyscaleImage segImg = imgHelper.getSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        boolean useBinned = false;
        
        contours = BlobsAndContours.populateContours(imgHelper, type, useBinned);
        
        segContoursMap.put(type, contours);
        
        return contours;
    }
    
    protected List<List<CurvatureScaleSpaceContour>> 
        generatePerimeterContoursForBinned(SegmentationType type) {
            
        List<List<CurvatureScaleSpaceContour>> contours = segBinnedContoursMap.get(type);

        if (contours != null) {
            return contours;
        }
            
        GreyscaleImage segImg = imgHelper.getBinnedSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        boolean useBinned = true;

        contours = BlobsAndContours.populateContours(imgHelper, type, useBinned);
        
        segBinnedContoursMap.put(type, contours);
        
        return contours;        
    }
    
    public int sumPointsOfInterest(SegmentationType segmentationType, 
        boolean useBinned) {
        
        List<List<CurvatureScaleSpaceContour>> contours 
            = getPerimeterContours(segmentationType, useBinned);
        
        int n = 0;
        
        for (List<CurvatureScaleSpaceContour> list : contours) {
            n += list.size();
        }
        
        return n;
    }

    // ---------------------------------
    /**
     * generate (or retrieve cached) for blob perimeter regions from 
     * segmented images of type type and binned or not binned.
     * @param type
     * @param applyToBinnedImage 
     * @return  
     */
    public List<List<BlobPerimeterRegion>> generatePerimeterRegions(
        SegmentationType type, boolean applyToBinnedImage) {
         
        if (applyToBinnedImage) {
            return generatePerimeterRegionsForBinned(type);
        } else {
            return generatePerimeterRegionsUnbinned(type);
        }
    }
    
    /**
     * generate (or retrieve cached) for blob perimeter regions from 
     * segmented images of type type and binned or not binned.
     * @param type
     * @param applyToBinnedImage 
     * @return  
     */
    public List<List<BlobPerimeterRegion>> getPerimeterRegions(
        SegmentationType type, boolean applyToBinnedImage) {
         
        return generatePerimeterRegions(type, applyToBinnedImage);
    }
    
    protected List<List<BlobPerimeterRegion>> 
        generatePerimeterRegionsUnbinned(SegmentationType type) {
        
        List<List<BlobPerimeterRegion>> bprs = segBPRMap.get(type);

        if (bprs != null) {
            return bprs;
        }
        
        GreyscaleImage segImg = imgHelper.getSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException("segmented image hasn't been created yet.  error?");
        }
        
        boolean useBinned = false;
        
        List<List<CurvatureScaleSpaceContour>> contours = getPerimeterContours(
            type, useBinned);
        
        bprs = BlobsAndContours.populateRegions(imgHelper, contours, type, 
            useBinned);
        
        segBPRMap.put(type, bprs);
        
        return bprs;
    }
    
    protected List<List<BlobPerimeterRegion>> 
        generatePerimeterRegionsForBinned(SegmentationType type) {
        
        List<List<BlobPerimeterRegion>> bprs = segBinnedBPRMap.get(type);

        if (bprs != null) {
            return bprs;
        }
        
        GreyscaleImage segImg = imgHelper.getBinnedSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException(
            "segmented image hasn't been created yet.  error?");
        }
        
        boolean useBinned = true;
        
        List<List<CurvatureScaleSpaceContour>> contours = getPerimeterContours(
            type, useBinned);
        
        bprs = BlobsAndContours.populateRegions(imgHelper, contours, type, 
            useBinned);
        
        segBinnedBPRMap.put(type, bprs);
        
        return bprs;
    }
}
