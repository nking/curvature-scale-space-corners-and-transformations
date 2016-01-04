package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
    
    public void replaceSegmentationWithCanny(SegmentationType segmentationType, 
        boolean applyToBinnedImage) {
        
        imgHelper.replaceSegmentationWithCanny(segmentationType, applyToBinnedImage);
        
        if (applyToBinnedImage) {
            segBinnedBlobsMap.remove(segmentationType);
            segBinnedBlobPerimetersMap.remove(segmentationType);
            assert(!segBinnedBlobsMap.containsKey(segmentationType));
            assert(!segBinnedBlobPerimetersMap.containsKey(segmentationType));
        } else {
            segBlobsMap.remove(segmentationType);
            segBlobPerimetersMap.remove(segmentationType);
            assert(!segBlobsMap.containsKey(segmentationType));
            assert(!segBlobPerimetersMap.containsKey(segmentationType));
        }
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
        
        boolean filterOutImageBoundaryBlobs = false;
        
        boolean filterOutZeroPixels = true;
                
        return getBlobPerimeters(type, useBinnedImage, 
            filterOutImageBoundaryBlobs, filterOutZeroPixels);
    }
    
    public List<PairIntArray> getBlobPerimeters(SegmentationType type, boolean 
        useBinnedImage, boolean filterOutImageBoundaryBlobs, 
        boolean filterOutZeroPixels) {
                
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
        
        boolean filterOutImageBoundaryBlobs = false;
        
        return getBlobs(type, useBinnedImage, filterOutImageBoundaryBlobs);
    }
    
    public List<Set<PairInt>> getBlobs(SegmentationType type, boolean 
        useBinnedImage, boolean filterOutImageBoundaryBlobs) {
                
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
         
        return generatePerimeterCorners(type, applyToBinnedImage);
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
        
        final boolean outdoorMode = false;
        
        final boolean enableJaggedLineCorrections = false;
        
        final float factorIncreaseForCurvatureMinimum = 1.f;
   
        corners = BlobsAndCorners.populateCorners(this, type, useBinned,
            outdoorMode, enableJaggedLineCorrections, 
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
        
        final boolean outdoorMode = false;
        
        final boolean enableJaggedLineCorrections = false;
        
        final float factorIncreaseForCurvatureMinimum = 1.f;

        corners = BlobsAndCorners.populateCorners(this, type, useBinned,
            outdoorMode, enableJaggedLineCorrections, 
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
     * pre-prepare the corners as every point in the perimeter of the blobs.
     * 
     * @param type 
     */
    private void extractBlobPerimeterAsCornerRegionsForBinned(SegmentationType type) {
        
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
        
        corners = extractCorners2(perimeterLists, type, useBinned);     
            
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

    /**
     * pre-prepare the corners as every point in the perimeter of the blobs.
     * 
     * @param type 
     */
    private void extractBlobPerimeterAsCornerRegionsForUnbinned(SegmentationType type) {
        
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
        
        corners = extractCorners2(perimeterLists, type, useBinned);
        
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
            final boolean outdoorMode = false;
            final boolean enableJaggedLineCorrections = false;
            final float factorIncreaseForCurvatureMinimum = 1.f;
                        
            corners = BlobsAndCorners.populateCorners(this, type, useBinned,
                outdoorMode, enableJaggedLineCorrections, 
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

}
