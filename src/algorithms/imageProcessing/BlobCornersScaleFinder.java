package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.Set;

/**
 * class to invoke methods needed to solve for euclidean scale between
 * image1 and image2 using methods specific to corners on closed curves.
 * 
 * @author nichole
 */
public class BlobCornersScaleFinder extends AbstractBlobScaleFinder {

    @Override
    public TransformationParameters solveForScale(
        ISegmentedImageHelper img1Helper, SegmentationType type1,
        boolean useBinned1,
        ISegmentedImageHelper img2Helper, SegmentationType type2,
        boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {
        
        if (!(img1Helper instanceof SegmentedImageBlobContourHelper) ||
            !(img2Helper instanceof SegmentedImageBlobContourHelper)) {
            throw new IllegalArgumentException("img1Helper and img2Helper must"
            + " be instances of SegmentedImageBlobContourHelper");
        }
        
        BlobsAndCorners bc1 = ((SegmentedImageBlobCornerHelper)img1Helper)
            .getBlobsAndCorners(type1, useBinned1);
        
        GreyscaleImage img1 = ((SegmentedImageBlobCornerHelper)img1Helper)
            .getGreyscaleImage(useBinned1);

        BlobsAndCorners bc2 = ((SegmentedImageBlobCornerHelper)img2Helper)
            .getBlobsAndCorners(type2, useBinned2);
        
        GreyscaleImage img2 = img2Helper.getGreyscaleImage(useBinned2);
        
        List<List<CornerRegion>> contours1List = bc1.getCorners();
        List<List<CornerRegion>> contours2List = bc2.getCorners();
        List<Set<PairInt>> blobs1 = bc1.getBlobs();
        List<Set<PairInt>> blobs2 = bc2.getBlobs();
        List<PairIntArray> perimeters1 = bc1.getBlobOrderedPerimeters();
        List<PairIntArray> perimeters2 = bc2.getBlobOrderedPerimeters();
        
        throw new UnsupportedOperationException("Not supported yet."); 
    }
    
}
