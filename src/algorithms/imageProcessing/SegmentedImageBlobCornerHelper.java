package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class SegmentedImageBlobCornerHelper extends AbstractSegmentedImageHelper {

    public SegmentedImageBlobCornerHelper(ImageExt img) {
        super(img);
    }
    
    public SegmentedImageBlobCornerHelper(ImageExt img, String debugTag) {
        super(img, debugTag);
    }

    @Override
    protected void clearBinnedPointsOfInterestMaps() {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    protected void clearUnbinnedPointsOfInterestMaps() {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    public void generatePerimeterPointsOfInterest(SegmentationType type, 
        boolean applyToBinnedImage) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    protected void generatePerimeterPointsOfInterestForUnbinned(SegmentationType type) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    protected void generatePerimeterPointsOfInterestForBinned(SegmentationType type) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    public int sumPointsOfInterest(SegmentationType segmentationType1, boolean useBinned1) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

}
