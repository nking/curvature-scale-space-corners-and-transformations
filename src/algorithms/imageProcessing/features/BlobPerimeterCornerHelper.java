package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceCornerDetector;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class BlobPerimeterCornerHelper {
    
    private final SegmentedImageHelper imgHelper;
    
    private final Map<SegmentationType, List<Set<PairInt>>> segBlobsMap =      
        new HashMap<SegmentationType, List<Set<PairInt>>>();
   
    private final Map<SegmentationType, List<PairIntArray>> 
        segBlobPerimetersMap = new HashMap<SegmentationType, List<PairIntArray>>();
    
    private final Map<SegmentationType, List<Set<PairInt>>> segBinnedBlobsMap =      
        new HashMap<SegmentationType, List<Set<PairInt>>>();
   
    private final Map<SegmentationType, List<PairIntArray>> 
        segBinnedBlobPerimetersMap = new HashMap<SegmentationType, List<PairIntArray>>();
    
    protected Map<SegmentationType, List<List<CornerRegion>>>
        segCornersMap = new HashMap<SegmentationType, 
        List<List<CornerRegion>>>();
    
    protected Map<SegmentationType, List<List<CornerRegion>>>
        segBinnedCornersMap = new HashMap<SegmentationType, 
        List<List<CornerRegion>>>();
    
    // ---------- these 2 maps are specifically for segmentation type NONE --
    protected Map<SegmentationType, List<CornerRegion>>
        segAllCornersMap = new HashMap<SegmentationType, List<CornerRegion>>();
    
    protected Map<SegmentationType, List<CornerRegion>>
        segBinnedAllCornersMap = new HashMap<SegmentationType, List<CornerRegion>>();
    //--------------------------------------------------------------------------
    
    public BlobPerimeterCornerHelper(final ImageExt img) {
        imgHelper = new SegmentedImageHelper(img);
    }
    
    public BlobPerimeterCornerHelper(final ImageExt img, final String debugTag) {
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
    
    public boolean didApplyHistEq() {
        return imgHelper.didApplyHistEq();
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
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored image");
        }
        return imgHelper.getSegmentationImage(type);
    }

    public GreyscaleImage getBinnedSegmentationImage(SegmentationType type) {
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored image");
        }
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
        
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored perimeter list");
        }
        
        boolean filterOutImageBoundaryBlobs = false;
        
        boolean filterOutZeroPixels = true;
                
        return getBlobPerimeters(type, useBinnedImage, 
            filterOutImageBoundaryBlobs, filterOutZeroPixels);
    }
    
    public List<PairIntArray> getBlobPerimeters(SegmentationType type, boolean 
        useBinnedImage, boolean filterOutImageBoundaryBlobs, 
        boolean filterOutZeroPixels) {
                
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored perimeter list");
        }
        
        if (useBinnedImage) {
            return getBinnedBlobPerimeters(type, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
        } else {
            return getUnbinnedBlobPerimeters(type, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
        }
       
    }
    
    private List<PairIntArray> getUnbinnedBlobPerimeters(SegmentationType type,
        boolean filterOutImageBoundaryBlobs, boolean filterOutZeroPixels) {
                
        List<PairIntArray> blobPerimeters = blobPerimeters = segBlobPerimetersMap.get(type);
        
        if (blobPerimeters != null) {
            return blobPerimeters;
        }
        
        List<Set<PairInt>> blobs = segBlobsMap.get(type);
        
        boolean useBinned = false;
        
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        if (blobs == null) {
            
            long t0 = System.currentTimeMillis();
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
            
            long t1 = System.currentTimeMillis();
            long t1Sec = (t1 - t0)/1000;
            Logger.getLogger(this.getClass().getName()).info("blobs(sec)=" 
                + t1Sec);
            
            segBlobsMap.put(type, blobs);
        }
        
        long t0 = System.currentTimeMillis();
        
        blobPerimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
            imgHelper, type, blobs, useBinned,
            discardWhenCavityIsSmallerThanBorder);
        
        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        Logger.getLogger(this.getClass().getName()).info("perimeters(sec)=" 
            + t1Sec);
        
        segBlobPerimetersMap.put(type, blobPerimeters);
        
        // blobs may have been removed, so re-put them
        segBlobsMap.put(type, blobs);
        
        return blobPerimeters;
    }
    
    private List<PairIntArray> getBinnedBlobPerimeters(SegmentationType type,
        boolean filterOutImageBoundaryBlobs, boolean filterOutZeroPixels) {
        
        List<PairIntArray> blobPerimeters = segBinnedBlobPerimetersMap.get(type);
        
        if (blobPerimeters != null) {
            return blobPerimeters;
        }
        
        List<Set<PairInt>> blobs = segBinnedBlobsMap.get(type);
        
        boolean useBinned = true;
        
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
            
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
        
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored blobs list");
        }
        
        boolean filterOutImageBoundaryBlobs = false;
        
        return getBlobs(type, useBinnedImage, filterOutImageBoundaryBlobs);
    }
    
    public List<Set<PairInt>> getBlobs(SegmentationType type, boolean 
        useBinnedImage, boolean filterOutImageBoundaryBlobs) {
         
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored blobs list");
        }
        
        boolean filterOutZeroPixels = true;
        
        if (useBinnedImage) {
            return getBinnedBlobs(type, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
        } else {
            return getUnbinnedBlobs(type, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
        }
    }
    
    public List<Set<PairInt>> getBlobs(SegmentationType type, boolean 
        useBinnedImage, boolean filterOutImageBoundaryBlobs,
        boolean filterOutZeroPixels) {
              
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE does not have a stored blobs list");
        }
        
        if (useBinnedImage) {
            return getBinnedBlobs(type, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
        } else {
            return getUnbinnedBlobs(type, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
        }
    }

    private List<Set<PairInt>> getBinnedBlobs(SegmentationType type,
        boolean filterOutImageBoundaryBlobs) {
        
        boolean filterOutZeroPixels = true;
        
        return getBinnedBlobs(type, filterOutImageBoundaryBlobs, 
            filterOutZeroPixels);
    }
    
    private List<Set<PairInt>> getBinnedBlobs(SegmentationType type,
        boolean filterOutImageBoundaryBlobs, boolean filterOutZeroPixels) {
        
        List<Set<PairInt>> blobs = segBinnedBlobsMap.get(type);
        
        boolean useBinned = true;
                
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
            
            segBinnedBlobsMap.put(type, blobs);
        }
        
        return blobs;            
    }

    private List<Set<PairInt>> getUnbinnedBlobs(SegmentationType type,
        boolean filterOutImageBoundaryBlobs, boolean filterOutZeroPixels) {
        
        List<Set<PairInt>> blobs = null;
        
        blobs = segBlobsMap.get(type);
        
        boolean useBinned = false;
                
        if (blobs == null) {
            
            blobs = BlobsAndPerimeters.extractBlobsFromSegmentedImage(
                imgHelper, type, useBinned, filterOutImageBoundaryBlobs,
                filterOutZeroPixels);
            
            segBlobsMap.put(type, blobs);
        }
        
        return blobs;
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
        
        if (type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type NONE needs to use getAllCorners(...)");
        }
        
        return generatePerimeterCorners(type, applyToBinnedImage);
    }
    
    public List<CornerRegion> getAllCorners(SegmentationType type, 
        boolean applyToBinnedImage) {
        
        if (applyToBinnedImage) {
            return getAllCornersForBinned(type);
        } else {
            return getAllCornersUnbinned(type);
        }
    }
    
    private List<CornerRegion> getAllCornersForBinned(SegmentationType type) {
        
        List<CornerRegion> corners = segBinnedAllCornersMap.get(type);
        
        if (corners != null) {
            return corners;
        }
        
        throw new IllegalArgumentException("corners have not been created" +
            " for type " + type);
    }
    
    private List<CornerRegion> getAllCornersUnbinned(SegmentationType type) {
        
        List<CornerRegion> corners = segAllCornersMap.get(type);
        
        if (corners != null) {
            return corners;
        }
        
        throw new IllegalArgumentException("corners have not been created" +
            " for type " + type);
    }
    
    protected List<List<CornerRegion>> 
        generatePerimeterCornersUnbinned(SegmentationType type) {
        
        List<List<CornerRegion>> corners = segCornersMap.get(type);

        if (corners != null) {
            return corners;
        }
        
        GreyscaleImage segImg = imgHelper.getSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException(
            "segmented image hasn't been created yet.  error?");
        }
        
        boolean useBinned = false;
                        
        final float factorIncreaseForCurvatureMinimum = 1.f;
   
        corners = BlobsAndCorners.populateCorners(this, type, useBinned,
            factorIncreaseForCurvatureMinimum);
       
        segCornersMap.put(type, corners);
        
        return corners;
    }
    
    protected List<List<CornerRegion>> 
        generatePerimeterCornersForBinned(SegmentationType type) {
         
        List<List<CornerRegion>> corners = segBinnedCornersMap.get(type);

        if (corners != null) {
            return corners;
        }
        
        GreyscaleImage segImg = imgHelper.getBinnedSegmentationImage(type);
        
        if (segImg == null) {
            //TODO: consider changing logic to perform this if needed
            throw new IllegalArgumentException(
            "segmented image hasn't been created yet.  error?");
        }
       
        boolean useBinned = true;
                        
        final float factorIncreaseForCurvatureMinimum = 1.f;

        corners = BlobsAndCorners.populateCorners(this, type, useBinned,
            factorIncreaseForCurvatureMinimum);
        
        segBinnedCornersMap.put(type, corners);
        
        return corners;        
    }
        
    /**
     * pre-prepare the corners as every point in the perimeter of the blobs.
     * Note that if perimeters were not already extracted, empty lists for 
     * perimeters are cached.  Also note that the corner regions are not 
     * ordered.
     * 
     * @param type 
     * @param applyToBinnedImage 
    */
    public void extractBlobPerimeterAsCornerRegions(SegmentationType type, 
        boolean applyToBinnedImage) {
        
        if (applyToBinnedImage) {
            extractBlobPerimeterAsCornerRegionsForBinned(type);
        } else {
            extractBlobPerimeterAsCornerRegionsForUnbinned(type);
        }
    }
    
    /**
     * pre-prepare the corners as points in the perimeter of the blobs.
     * Note that if perimeters were not already extracted, empty lists for 
     * perimeters are cached.  Also note that the corner regions are not 
     * ordered.
     * 
     * @param type 
     * @param applyToBinnedImage 
    */
    public void extractBlobPerimeterAsCornerRegions(SegmentationType type, 
        boolean applyToBinnedImage, boolean doNotAddPoints) {
        
        if (applyToBinnedImage) {
            extractBlobPerimeterAsCornerRegionsForBinned(type, doNotAddPoints);
        } else {
            extractBlobPerimeterAsCornerRegionsForUnbinned(type, doNotAddPoints);
        }
    }

    /**
     * extract the corners as second derivative gaussian high value points
     * and keep those associated with blobs created for type and binning.
     * 
     * @param type 
     * @param applyToBinnedImage 
    */
    public void extractSecondDerivativeCorners(SegmentationType type, 
        boolean applyToBinnedImage) throws IOException, NoSuchAlgorithmException {
                
        if (applyToBinnedImage) {
            extractSecondDerivativeCornersForBinned2(type);
        } else {
            extractSecondDerivativeCornersForUnbinned2(type);
        }
    }
    
    /**
     * extract the corners as second derivative gaussian high value points
     * and keep those associated with blobs created for type and binning.
     * 
     * @param type 
     * @param applyToBinnedImage 
    */
    public void extractSecondDerivativeCornersToColor(SegmentationType type, 
        boolean applyToBinnedImage) throws IOException, NoSuchAlgorithmException {
                
        if (applyToBinnedImage) {
            extractSecondDerivativeCornersForBinned3(type);
        } else {
            extractSecondDerivativeCornersForUnbinned3(type);
        }
    }
    
    /**
     * extract the canny edge corners and do not filter for association with blobs.
     * Note, these corners are stored in a data structure that is only accessed
     * if type is NONE and using getAllCorners()
     * 
     * @param type 
     * @param applyToBinnedImage 
    */
    public void extractCannyCornersWithoutBlobs(SegmentationType type, 
        boolean applyToBinnedImage) {
        
        if (!type.equals(SegmentationType.GREYSCALE_CANNY)) {
            throw new IllegalArgumentException(
            "segmentation type must be GREYSCALE_CANNY to store and retrieve the points"
            + " correctly.  Any blobs needed must be created separately with a type that"
            + " is not NONE");
        }
                
        if (applyToBinnedImage) {
            extractCannyCornersForBinnedWithoutBlobs(type);
        } else {
            extractCannyCornersForUnbinnedWithoutBlobs(type);
        }
    }
    
    /**
     * extract the corners as second derivative gaussian high value points
     * and filter or do not filter for association with blobs.
     * Note, these corners are stored in a data structure that is only accessed
     * if type is NONE and using getAllCorners()
     * 
     * @param type 
     * @param applyToBinnedImage 
    */
    public void extractSecondDerivativeCornersWithoutBlobs(SegmentationType type, 
        boolean applyToBinnedImage) {
        
        if (!type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type must be NONE to store and retrieve the points"
            + " correctly.  Blobs must be created separately with a type that"
            + " is not NONE");
        }
                
        if (applyToBinnedImage) {
            extractSecondDerivativeCornersForBinnedWithoutBlobs(type);
        } else {
            extractSecondDerivativeCornersForUnbinnedWithoutBlobs(type);
        }
    }

    /**
     * pre-prepare the corners as every point in the perimeter of the blobs.
     * 
     * @param type 
     */
    private void extractBlobPerimeterAsCornerRegionsForBinned(SegmentationType 
        type) {
        
        boolean doNotAddPoints = false;
        
        extractBlobPerimeterAsCornerRegionsForBinned(type, doNotAddPoints);
    }
    
    /**
     * pre-prepare the corners as every point in the perimeter of the blobs.
     * 
     * @param type 
     */
    private void extractBlobPerimeterAsCornerRegionsForBinned(
        SegmentationType type, boolean doNotAddPoints) {
        
        List<List<CornerRegion>> corners = segBinnedCornersMap.get(type);

        if (corners != null) {
            return;
        }
        
        List<Set<PairInt>> blobs = getBlobs(type, true);
        
        GreyscaleImage gsImg = this.getGreyscaleImageBinned();
        int imageWidth = gsImg.getWidth(); 
        int imageHeight = gsImg.getHeight();
        
        boolean useBinned = true;
        boolean discardWhenCavityIsSmallerThanBorder = false;
        
        List<PairIntArray> perimeterLists = segBinnedBlobPerimetersMap.get(type);
        
        if (perimeterLists == null) {
            
            perimeterLists = BlobsAndPerimeters.extractBoundsOfBlobs(
                imgHelper, type, blobs, useBinned,
                discardWhenCavityIsSmallerThanBorder);
            
            segBinnedBlobPerimetersMap.put(type, perimeterLists);
        }
        
        if (doNotAddPoints) {
            corners = extractCorners1(perimeterLists, type, useBinned);
        } else {
            corners = extractCorners2(perimeterLists, type, useBinned);
        }
            
        segBinnedCornersMap.put(type, corners);
        
        if (imgHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, corners,
                1, gsImg, "blob_corners_" + imgHelper.getDebugTag() + "_" 
                + MiscDebug.getCurrentTimeFormatted());
        }
        
        Logger log = Logger.getLogger(this.getClass().getName());
        log.info("perimeterLists.size()=" + perimeterLists.size() +
            " corners.size()=" + corners.size());
        
        assert(perimeterLists.size() == corners.size());
    }

    private void extractBlobPerimeterAsCornerRegionsForUnbinned(SegmentationType type) {
        
        boolean doNotAddPoints = false;
        
        extractBlobPerimeterAsCornerRegionsForUnbinned(type, doNotAddPoints);
    }
    
    /**
     * pre-prepare the corners as every point in the perimeter of the blobs.
     * 
     * @param type 
     */
    private void extractBlobPerimeterAsCornerRegionsForUnbinned(
        SegmentationType type, boolean doNotAddPoints) {
        
        List<List<CornerRegion>> corners = segCornersMap.get(type);

        if (corners != null) {
            return;
        }
        
        List<Set<PairInt>> blobs = getBlobs(type, false);
        
        GreyscaleImage gsImg = this.getGreyscaleImage();
        int imageWidth = gsImg.getWidth(); 
        int imageHeight = gsImg.getHeight();
        
        boolean useBinned = false;
        boolean discardWhenCavityIsSmallerThanBorder = true;
        
        List<PairIntArray> perimeterLists = segBlobPerimetersMap.get(type);
        
        if (perimeterLists == null) {
            
            perimeterLists = BlobsAndPerimeters.extractBoundsOfBlobs(
                imgHelper, type, blobs, useBinned,
                discardWhenCavityIsSmallerThanBorder);
            
            segBlobPerimetersMap.put(type, perimeterLists);
        }
        
        if (doNotAddPoints) {
            corners = extractCorners1(perimeterLists, type, useBinned);
        } else {
            corners = extractCorners2(perimeterLists, type, useBinned);
        }
        
        segCornersMap.put(type, corners);
        
        if (imgHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, corners,
                1, gsImg, "blob_corners_" + imgHelper.getDebugTag() + "_" 
                + MiscDebug.getCurrentTimeFormatted());
        }
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
    
    public boolean isInDebugMode() {
        return imgHelper.isInDebugMode();
    }

    public String getDebugTag() {
        return imgHelper.getDebugTag();
    }

    private List<List<CornerRegion>> extractCorners2(List<PairIntArray> 
        perimeterLists, SegmentationType type, boolean useBinned) {
        
        boolean useMorePoints = false;
        
        List<List<CornerRegion>> corners = new ArrayList<List<CornerRegion>>();
        
        if (useMorePoints) {
            int skip = DitherDefault.dither + 1;
            corners = new ArrayList<List<CornerRegion>>();
            for (int i = 0; i < perimeterLists.size(); ++i) {
                PairIntArray edge = perimeterLists.get(i);
                List<CornerRegion> bprList = new ArrayList<CornerRegion>();
                for (int j = 0; j < edge.getN(); j += skip) {
                    CornerRegion cr = new CornerRegion(i, 1, 0);
                    cr.setFlagThatNeighborsHoldDummyValues();
                    cr.set(0, Float.MIN_VALUE, edge.getX(j), edge.getY(j));
                    cr.setIndexWithinCurve(j);
                    bprList.add(cr);
                }
                corners.add(bprList);
            }
        } else {        
            final float factorIncreaseForCurvatureMinimum = 1.f;
                        
            corners = BlobsAndCorners.populateCorners(this, type, useBinned,
                factorIncreaseForCurvatureMinimum);
            assert(corners.size() == perimeterLists.size());
            int k = 5;
            for (int i = 0; i < corners.size(); ++i) {
                PairIntArray edge = perimeterLists.get(i);
                int n = edge.getN();
                if (n < k) {
                    continue;
                }
                List<CornerRegion> list = corners.get(i);
                if (list.isEmpty()) {
                    // add a starter point
                    int idx = (n - 1)/2;
                    CornerRegion cr = new CornerRegion(i, 1, 0);
                    cr.setFlagThatNeighborsHoldDummyValues();
                    cr.set(0, Float.MIN_VALUE, edge.getX(idx), edge.getY(idx));
                    cr.setIndexWithinCurve(idx);
                    list.add(cr);
                }                    
                // bisect largest gap until have k points or largest gap is '1'
                while (list.size() < k) {
                    int maxGap = list.get(0).getIndexWithinCurve();
                    if (maxGap == -1) {
                        // corners cannot be evaluated this way due to missing info
                        break;
                    }
                    // maxGapIdx + 1 is the point where to insert into list
                    int maxGapIdx = -1;
                    for (int j = 0; j < list.size(); ++j) {
                        int gap;
                        if (j == (list.size() - 1)) {
                            gap = (n - 1) - list.get(j).getIndexWithinCurve();
                        } else {
                            gap = list.get(j + 1).getIndexWithinCurve() - list.get(j).getIndexWithinCurve();
                        }
                        if (gap > maxGap) {
                            maxGap = gap;
                            maxGapIdx = j;
                        }
                    }
                    if (list.size() == 1) {
                        int idx = ((n - 1) + list.get(0).getIndexWithinCurve())/2;
                        CornerRegion cr = new CornerRegion(i, 1, 0);
                        cr.setFlagThatNeighborsHoldDummyValues();
                        cr.set(0, Float.MIN_VALUE, edge.getX(idx), edge.getY(idx));
                        cr.setIndexWithinCurve(idx);
                        list.add(1, cr);
                        continue;
                    }
                    if (maxGap == 1) {
                        break;
                    }
                    int idx;
                    if (maxGapIdx == -1) {
                        idx = list.get(0).getIndexWithinCurve()/2;
                    } else if (maxGapIdx >= (list.size() - 1)) {
                        idx = ((n - 1) + list.get(maxGapIdx).getIndexWithinCurve())/2;
                    } else {
                        idx = (list.get(maxGapIdx + 1).getIndexWithinCurve() 
                            + list.get(maxGapIdx).getIndexWithinCurve())/2;
                    }
                    CornerRegion cr = new CornerRegion(i, 1, 0);
                    cr.setFlagThatNeighborsHoldDummyValues();
                    cr.set(0, Float.MIN_VALUE, edge.getX(idx), edge.getY(idx));
                    cr.setIndexWithinCurve(idx);
                    list.add(maxGapIdx + 1, cr);
                }
            }
        }
        return corners;
    }

    private List<List<CornerRegion>> extractCorners1(List<PairIntArray> 
        perimeterLists, SegmentationType type, boolean useBinned) {
                         
        final float factorIncreaseForCurvatureMinimum = 1.f;

        List<List<CornerRegion>> corners = BlobsAndCorners.populateCorners(this, 
            type, useBinned, factorIncreaseForCurvatureMinimum);
        
        assert(corners.size() == perimeterLists.size());
            
        return corners;
    }

    private void extractSecondDerivativeCornersForBinned(SegmentationType type) {
        
        List<List<CornerRegion>> corners = segBinnedCornersMap.get(type);

        if (corners != null) {
            return;
        }

        // we keep only the points associated with blobs
        List<Set<PairInt>> blobs = getBlobs(type, true);
        
        GreyscaleImage gsImg = this.getGreyscaleImageBinned();
        
        List<PairIntArray> perimeterLists = segBinnedBlobPerimetersMap.get(type);
        
        if (perimeterLists == null) {
            
            perimeterLists = new ArrayList<PairIntArray>();
            
            // keep the number of items in list consistent even though empty
            for (int i = 0; i < blobs.size(); ++i) {
                perimeterLists.add(new PairIntArray());
            }
            
            segBinnedBlobPerimetersMap.put(type, perimeterLists);
        }
        
        Logger log = Logger.getLogger(this.getClass().getName());
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(gsImg,
            //200, true);
            5000, true);
        
        corners = new ArrayList<List<CornerRegion>>();
           
        List<Set<PairInt>> pixels2 = new ArrayList<Set<PairInt>>();
        
        for (int i = 0; i < blobs.size(); ++i) {
            List<CornerRegion> crList = new ArrayList<CornerRegion>();
            corners.add(crList);
            pixels2.add(new HashSet<PairInt>());
        }
        
        for (PairInt p : pixels) {
            boolean found = false;
            for (int i = 0; i < blobs.size(); ++i) {
                Set<PairInt> blob = blobs.get(i);
                if (blob.contains(p)) {
                    pixels2.get(i).add(p);
                    found = true;
                    break;
                }
            }
            if (!found) {
                log.info("discarding because not in blob: " + p.toString());
            }
        }
        
        // if uncomment this, have to increase dither to larger than 1
        //8NeighborCentroids(pixels2);
        
        for (int i = 0; i < pixels2.size(); ++i) {
            Set<PairInt> p2 = pixels2.get(i);
            for (PairInt p : p2) {
                CornerRegion cr = new CornerRegion(i, 1, 0);
                cr.setFlagThatNeighborsHoldDummyValues();
                cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
                cr.setIndexWithinCurve(-1);
                corners.get(i).add(cr);
            }
        }
        
        segBinnedCornersMap.put(type, corners);
        
        //TODO: corner making will be refactored soon. this addition here allows later logic
        // using getAllCorners to fetch the 2nd deriv corners as expected.
        // (haven't created settings for end deriv + blob assoc filter yet)
        List<CornerRegion> corners2 = segBinnedAllCornersMap.get(SegmentationType.NONE);
        if (corners2 == null) {
            corners2 = new ArrayList<CornerRegion>();
            for (List<CornerRegion> list : corners) {
                corners2.addAll(list);
            }
            segBinnedAllCornersMap.put(SegmentationType.NONE, corners2);
        }

        if (imgHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, corners,
                1, gsImg, "blob_corners_" + imgHelper.getDebugTag() + "_" 
                + MiscDebug.getCurrentTimeFormatted());
        }
        
        log.info("perimeterLists.size()=" + perimeterLists.size() +
            " corners.size()=" + corners.size());
        
        assert(perimeterLists.size() == corners.size());
        
    }
    
    /**
     * places corners in the "all corners" map obtainable by getAllCorners(SegmentationType.NONE).
     * It extracts the 2nd deriv pts and filters to keep only those that are
     * coincident w/ pix value > 0 from the segmentation image.
     * @param type 
     */
    private void extractSecondDerivativeCornersForBinned2(SegmentationType type) 
        throws IOException, NoSuchAlgorithmException {
        
        List<CornerRegion> corners = segBinnedAllCornersMap.get(type);

        if (corners != null) {
            return;
        }

        GreyscaleImage gsImg = this.getGreyscaleImageBinned();
                        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(
            gsImg,
            //200, true);
            5000, true);
        
        this.imgHelper.applySegmentation(SegmentationType.GREYSCALE_WAVELET, true);
        GreyscaleImage segImg = this.imgHelper.getBinnedSegmentationImage(type);
        Set<PairInt> pixSet = new HashSet<PairInt>();
        for (int i = 0; i < segImg.getWidth(); ++i) {
            for (int j = 0; j < segImg.getHeight(); ++j) {
                int v = segImg.getValue(i, j);
                if (v > 0) {
                    pixSet.add(new PairInt(i, j));
                }
            }
        }
        
        corners = new ArrayList<CornerRegion>();
        
        for (PairInt p : pixels) {
            if (pixSet.contains(p)) {
                CornerRegion cr = new CornerRegion(0, 1, 0);
                cr.setFlagThatNeighborsHoldDummyValues();
                cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
                cr.setIndexWithinCurve(-1);
                corners.add(cr);
            }
        }
        
        segBinnedAllCornersMap.put(SegmentationType.NONE, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(pixels, gsImg.copyToColorGreyscale(),
                    "all_2ndderiv_" + imgHelper.getDebugTag() + "_" + ts);
                MiscDebug.writeImage(corners, gsImg.copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
        
    }
    
    /**
     * places corners in the "all corners" map obtainable by getAllCorners(SegmentationType.NONE).
     * It extracts the 2nd deriv pts and filters to keep only those that are
     * coincident w/ pix value > 0 from the segmentation image.
     * @param type 
     */
    private void extractSecondDerivativeCornersForBinned3(SegmentationType type) 
        throws IOException, NoSuchAlgorithmException {
        
        List<CornerRegion> corners = segBinnedAllCornersMap.get(type);

        if (corners != null) {
            return;
        }

        int binFactor = imgHelper.getBinFactor();
        Image clrImg = this.getImage();
        ImageProcessor imageProcessor = new ImageProcessor();
        clrImg = imageProcessor.binImage(clrImg, binFactor);
        GreyscaleImage rImg = clrImg.copyRedToGreyscale();
        GreyscaleImage gImg = clrImg.copyGreenToGreyscale();
        GreyscaleImage bImg = clrImg.copyBlueToGreyscale();
                                        
        Set<PairInt> pR = imageProcessor.extract2ndDerivPoints(rImg,
            //200, true);
            5000, true);
        Set<PairInt> pG = imageProcessor.extract2ndDerivPoints(gImg,
            //200, true);
            5000, true);
        Set<PairInt> pB = imageProcessor.extract2ndDerivPoints(bImg,
            //200, true);
            5000, true);
        Set<PairInt> pixels = new HashSet<PairInt>();
        pixels.addAll(pR);
        pixels.addAll(pG);
        pixels.addAll(pB);
                
        boolean use1D = true;
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        GreyscaleImage rSegImg = imageSegmentation.createGreyscale5(rImg, use1D);
        GreyscaleImage gSegImg = imageSegmentation.createGreyscale5(gImg, use1D);
        GreyscaleImage bSegImg = imageSegmentation.createGreyscale5(bImg, use1D);
        
        GreyscaleImage tmpSegImgForPrinting = rSegImg.copyImage();
        
        Set<PairInt> pixSet = new HashSet<PairInt>();
        for (int i = 0; i < rSegImg.getWidth(); ++i) {
            for (int j = 0; j < rSegImg.getHeight(); ++j) {
                if (rSegImg.getValue(i, j) > 0) {
                    pixSet.add(new PairInt(i, j));
                } else if (gSegImg.getValue(i, j) > 0) {
                    pixSet.add(new PairInt(i, j));
                    tmpSegImgForPrinting.setValue(i, j, gSegImg.getValue(i, j));
                } else if (bSegImg.getValue(i, j) > 0) {
                    pixSet.add(new PairInt(i, j));
                    tmpSegImgForPrinting.setValue(i, j, bSegImg.getValue(i, j));
                }
            }
        }
        
        corners = new ArrayList<CornerRegion>();
                   
        for (PairInt p : pixels) {
            if (pixSet.contains(p)) {
                CornerRegion cr = new CornerRegion(0, 1, 0);
                cr.setFlagThatNeighborsHoldDummyValues();
                cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
                cr.setIndexWithinCurve(-1);
                corners.add(cr);
            }
        }
        
        segBinnedAllCornersMap.put(SegmentationType.NONE, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(pixels, this.getGreyscaleImage(true).copyToColorGreyscale(),
                    "all_2ndderiv_" + imgHelper.getDebugTag() + "_" + ts);
                MiscDebug.writeImage(corners, this.getGreyscaleImage(true).copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_" + ts);
                MiscDebug.writeImage(tmpSegImgForPrinting, 
                    "segmented_" + imgHelper.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
    }
    
    /**
     * places corners in the "all corners" map obtainable by getAllCorners.
     * It extracts the 2nd deriv pts and filters to keep only those that are
     * coincident w/ pix value > 0 from the segmentation image.
     * @param type 
     */
    private void extractSecondDerivativeCornersForUnbinned2(SegmentationType type) throws IOException, NoSuchAlgorithmException {
        
        List<CornerRegion> corners = segAllCornersMap.get(type);

        if (corners != null) {
            return;
        }

        GreyscaleImage gsImg = this.getGreyscaleImageBinned();
                        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(gsImg,
            //200, true);
            5000, true);
        
        this.imgHelper.applySegmentation(SegmentationType.GREYSCALE_WAVELET, false);
        GreyscaleImage segImg = this.imgHelper.getSegmentationImage(type);
        Set<PairInt> pixSet = new HashSet<PairInt>();
        for (int i = 0; i < segImg.getWidth(); ++i) {
            for (int j = 0; j < segImg.getHeight(); ++j) {
                int v = segImg.getValue(i, j);
                if (v > 0) {
                    pixSet.add(new PairInt(i, j));
                }
            }
        }
        
        corners = new ArrayList<CornerRegion>();
                   
        for (PairInt p : pixels) {
            if (pixSet.contains(p)) {
                CornerRegion cr = new CornerRegion(0, 1, 0);
                cr.setFlagThatNeighborsHoldDummyValues();
                cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
                cr.setIndexWithinCurve(-1);
                corners.add(cr);
            }
        }
        
        segAllCornersMap.put(SegmentationType.NONE, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(pixels, gsImg.copyToColorGreyscale(),
                    "all_2ndderiv_" + imgHelper.getDebugTag() + "_" + ts);
                MiscDebug.writeImage(corners, gsImg.copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
        
    }
    
    /**
     * places corners in the "all corners" map obtainable by getAllCorners.
     * It extracts the 2nd deriv pts and filters to keep only those that are
     * coincident w/ pix value > 0 from the segmentation image.
     * @param type 
     */
    private void extractSecondDerivativeCornersForUnbinned3(SegmentationType type) 
        throws IOException, NoSuchAlgorithmException {
        
        List<CornerRegion> corners = segAllCornersMap.get(type);

        if (corners != null) {
            return;
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        
        Image clrImg = this.getImage();        
        GreyscaleImage rImg = clrImg.copyRedToGreyscale();
        GreyscaleImage gImg = clrImg.copyGreenToGreyscale();
        GreyscaleImage bImg = clrImg.copyBlueToGreyscale();
                                        
        Set<PairInt> pR = imageProcessor.extract2ndDerivPoints(rImg,
            //200, true);
            5000, true);
        Set<PairInt> pG = imageProcessor.extract2ndDerivPoints(gImg,
            //200, true);
            5000, true);
        Set<PairInt> pB = imageProcessor.extract2ndDerivPoints(bImg,
            //200, true);
            5000, true);
        Set<PairInt> pixels = new HashSet<PairInt>();
        pixels.addAll(pR);
        pixels.addAll(pG);
        pixels.addAll(pB);
                
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        GreyscaleImage rSegImg = imageSegmentation.createGreyscale5(rImg);
        GreyscaleImage gSegImg = imageSegmentation.createGreyscale5(gImg);
        GreyscaleImage bSegImg = imageSegmentation.createGreyscale5(bImg);
        
        Set<PairInt> pixSet = new HashSet<PairInt>();
        for (int i = 0; i < rSegImg.getWidth(); ++i) {
            for (int j = 0; j < rSegImg.getHeight(); ++j) {
                if (rSegImg.getValue(i, j) > 0) {
                    pixSet.add(new PairInt(i, j));
                } else if (gSegImg.getValue(i, j) > 0) {
                    pixSet.add(new PairInt(i, j));
                } else if (bSegImg.getValue(i, j) > 0) {
                    pixSet.add(new PairInt(i, j));
                }
            }
        }
        
        corners = new ArrayList<CornerRegion>();
                   
        for (PairInt p : pixels) {
            if (pixSet.contains(p)) {
                CornerRegion cr = new CornerRegion(0, 1, 0);
                cr.setFlagThatNeighborsHoldDummyValues();
                cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
                cr.setIndexWithinCurve(-1);
                corners.add(cr);
            }
        }
        
        segAllCornersMap.put(SegmentationType.NONE, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeImage(pixels, this.getGreyscaleImage(false).copyToColorGreyscale(),
                    "all_2ndderiv_" + imgHelper.getDebugTag() + "_" + ts);
                MiscDebug.writeImage(corners, this.getGreyscaleImage(false).copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_" + ts);
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
        
    }
    
    private void extractCannyCornersForBinnedWithoutBlobs(SegmentationType type) {
        
        if (!type.equals(SegmentationType.GREYSCALE_CANNY)) {
            throw new IllegalArgumentException(
            "segmentation type must be GREYSCALE_CANNY to store and retrieve the points"
            + " correctly.  Any blobs needed must be created separately with a type that"
            + " is not NONE");
        }
        
        List<CornerRegion> corners = segBinnedAllCornersMap.get(type);

        if (corners != null) {
            return;
        }
                        
        GreyscaleImage gsImg = this.getGreyscaleImageBinned().copyImage();
        
        CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(gsImg);
        detector.doNotFindJunctions();
        detector.doNotStoreCornerRegions();
        
        detector.findCorners();
                       
        PairIntArray c1 = detector.getCorners();
        //Set<CornerRegion> cornersA = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        
        if (c1 == null) {
            Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).severe("no cornerw were found");
            return;
        }
        
        corners = new ArrayList<CornerRegion>();
        
        for (int i = 0; i < c1.getN(); ++i) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, c1.getX(i), c1.getY(i));
            cr.setIndexWithinCurve(-1);
            corners.add(cr);
        }
                
        segBinnedAllCornersMap.put(type, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                MiscDebug.writeImage(corners, gsImg.copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_"
                        + MiscDebug.getCurrentTimeFormatted());
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
    }
    
    private void extractCannyCornersForUnbinnedWithoutBlobs(SegmentationType type) {
        
        if (!type.equals(SegmentationType.GREYSCALE_CANNY)) {
            throw new IllegalArgumentException(
            "segmentation type must be GREYSCALE_CANNY to store and retrieve the points"
            + " correctly.  Any blobs needed must be created separately with a type that"
            + " is not NONE");
        }
        
        List<CornerRegion> corners = segAllCornersMap.get(type);

        if (corners != null) {
            return;
        }
                        
        GreyscaleImage gsImg = this.getGreyscaleImageBinned().copyImage();
        
        CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(gsImg);
        detector.doNotFindJunctions();
        detector.doNotStoreCornerRegions();
        
        detector.findCorners();
        
        PairIntArray c1 = detector.getCorners();
        //Set<CornerRegion> cornersA = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        
        if (c1 == null) {
            Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).severe("no cornerw were found");
            return;
        }
        
        corners = new ArrayList<CornerRegion>();
        
        for (int i = 0; i < c1.getN(); ++i) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, c1.getX(i), c1.getY(i));
            cr.setIndexWithinCurve(-1);
            corners.add(cr);
        }
                
        segAllCornersMap.put(type, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                MiscDebug.writeImage(corners, gsImg.copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_"
                        + MiscDebug.getCurrentTimeFormatted());
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
    }
    
    private void extractSecondDerivativeCornersForBinnedWithoutBlobs(SegmentationType type) {
        
        if (!type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type must be NONE to store and retrieve the points"
            + " correctly.  Blobs must be created separately with a type that"
            + " is not NONE");
        }
        
        List<CornerRegion> corners = segBinnedAllCornersMap.get(type);

        if (corners != null) {
            return;
        }
                
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage gsImg = this.getGreyscaleImageBinned();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(gsImg,
            //500, true);
             //1000, true);
             5000, true);
        
        corners = new ArrayList<CornerRegion>();
        
        for (PairInt p : pixels) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
            cr.setIndexWithinCurve(-1);
            corners.add(cr);
        }
        
        segBinnedAllCornersMap.put(type, corners);
        
        if (imgHelper.isInDebugMode()) {
            try {
                MiscDebug.writeImage(corners, gsImg.copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_"
                        + MiscDebug.getCurrentTimeFormatted());
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
                
    }
    
    private void extractSecondDerivativeCornersForUnbinned(SegmentationType type) {
        
        List<List<CornerRegion>> corners = segCornersMap.get(type);

        if (corners != null) {
            return;
        }

        // we keep only the points associated with blobs
        List<Set<PairInt>> blobs = getBlobs(type, false);
        
        GreyscaleImage gsImg = this.getGreyscaleImage();
        
        List<PairIntArray> perimeterLists = segBlobPerimetersMap.get(type);
        
        if (perimeterLists == null) {
            
            perimeterLists = new ArrayList<PairIntArray>();
            
            // keep the number of items in list consistent even though empty
            for (int i = 0; i < blobs.size(); ++i) {
                perimeterLists.add(new PairIntArray());
            }
            
            segBlobPerimetersMap.put(type, perimeterLists);
        }
        
        Logger log = Logger.getLogger(this.getClass().getName());
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(gsImg,
            //200, true);
            5000, true);
        
        corners = new ArrayList<List<CornerRegion>>();
        
        List<Set<PairInt>> pixels2 = new ArrayList<Set<PairInt>>();
        
        for (int i = 0; i < blobs.size(); ++i) {
            List<CornerRegion> crList = new ArrayList<CornerRegion>();
            corners.add(crList);
            pixels2.add(new HashSet<PairInt>());
        }
        
        for (PairInt p : pixels) {
            boolean found = false;
            for (int i = 0; i < blobs.size(); ++i) {
                Set<PairInt> blob = blobs.get(i);
                if (blob.contains(p)) {
                    pixels2.get(i).add(p);
                    found = true;
                    break;
                }
            }
            if (!found) {
                log.info("discarding because not in blob: " + p.toString());
            }
        }
        
        // if uncomment this, have to increase dither to larger than 1
        //reduceTo8NeighborCentroids(pixels2);
        
        for (int i = 0; i < pixels2.size(); ++i) {
            Set<PairInt> p2 = pixels2.get(i);
            for (PairInt p : p2) {
                CornerRegion cr = new CornerRegion(i, 1, 0);
                cr.setFlagThatNeighborsHoldDummyValues();
                cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
                cr.setIndexWithinCurve(-1);
                corners.get(i).add(cr);
            }
        }
        
        segCornersMap.put(type, corners);

        //TODO: corner making will be refactored soon. this addition here allows later logic
        // using getAllCorners to fetch the 2nd deriv corners as expected.
        // (haven't created settings for end deriv + blob assoc filter yet)
        List<CornerRegion> corners2 = segAllCornersMap.get(SegmentationType.NONE);
        if (corners2 == null) {
            corners2 = new ArrayList<CornerRegion>();
            for (List<CornerRegion> list : corners) {
                corners2.addAll(list);
            }
            segAllCornersMap.put(SegmentationType.NONE, corners2);
        }
        
        if (imgHelper.isInDebugMode()) {
            MiscDebug.writeEdgesAndCorners(perimeterLists, corners,
                1, gsImg, "blob_corners_" + imgHelper.getDebugTag() + "_"
                + MiscDebug.getCurrentTimeFormatted());
        }
        
        log.info("perimeterLists.size()=" + perimeterLists.size() +
            " corners.size()=" + corners.size());
        
        assert(perimeterLists.size() == corners.size());
        
    }

    private void extractSecondDerivativeCornersForUnbinnedWithoutBlobs(SegmentationType type) {
        
        if (!type.equals(SegmentationType.NONE)) {
            throw new IllegalArgumentException(
            "segmentation type must be NONE to store and retrieve the points"
            + " correctly.  Blobs must be created separately with a type that"
            + " is not NONE");
        }
        
        List<CornerRegion> corners = segAllCornersMap.get(type);

        if (corners != null) {
            return;
        }
                
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage gsImg = this.getGreyscaleImage();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(gsImg,
            //500, true);
            //1000, true);
            5000, true);
        
        corners = new ArrayList<CornerRegion>();
        
        for (PairInt p : pixels) {
            CornerRegion cr = new CornerRegion(0, 1, 0);
            cr.setFlagThatNeighborsHoldDummyValues();
            cr.set(0, Float.MIN_VALUE, p.getX(), p.getY());
            cr.setIndexWithinCurve(-1);
            corners.add(cr);
        }
        
        segAllCornersMap.put(type, corners);

        if (imgHelper.isInDebugMode()) {
            try {
                MiscDebug.writeImage(corners, gsImg.copyToColorGreyscale(),
                    "corners_" + imgHelper.getDebugTag() + "_"
                        + MiscDebug.getCurrentTimeFormatted());
            } catch (IOException ex) {
                Logger.getLogger(BlobPerimeterCornerHelper.class.getName()).
                    log(Level.SEVERE, null, ex);
            }
        }
                
    }
    
    private void reduceTo8NeighborCentroids(List<Set<PairInt>> pixelLists) {
        
        Set<PairInt> processed = new HashSet<PairInt>();
                
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int nBefore = 0;
        int nAfter = 0;
        
        Set<PairInt> neighbors = new HashSet<PairInt>();
        
        for (int i = 0; i < pixelLists.size(); ++i) {
            
            Set<PairInt> pixels = pixelLists.get(i);
            
            Set<PairInt> output = new HashSet<PairInt>();
            
            nBefore += pixels.size();

            for (PairInt p : pixels) {

                if (processed.contains(p)) {
                    continue;
                }

                curveHelper.findNeighbors(p.getX(), p.getY(), pixels, processed, 
                    dxs, dys, neighbors);

                processed.add(p);
                processed.addAll(neighbors);

                if (neighbors.size() == 0) {
                    output.add(p);
                } else {
                    double[] xyCen = curveHelper.calculateXYCentroids(neighbors);
                    int x = (int)Math.round(xyCen[0]);
                    int y = (int)Math.round(xyCen[1]);
                    assert(Math.abs(x - p.getX()) <= 3);
                    assert(Math.abs(y - p.getY()) <= 3);
                    output.add(new PairInt(x, y));
                }
            }
        
            pixels.clear();
            pixels.addAll(output);
            
            nAfter += pixels.size();
        }
        
        Logger.getLogger(this.getClass().getName()).info("nBefore=" + nBefore 
            + " nAfter=" + nAfter);
    }

}
