package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public interface IBlobScaleFinder {

    /**
     * set the instance to debug mode for additional logging and figures
     */
    void setToDebug();

    /**
     * solve for euclidean transformation scale between image1 and image2
     * given the images contained in img1Helper and img2Helper.
     * @param img1Helper
     * @param type1
     * @param useBinned1
     * @param img2Helper
     * @param type2
     * @param useBinned2
     * @param outputScaleRotTransXYStDev
     * @return 
     */
    public TransformationParameters solveForScale(
        ISegmentedImageHelper img1Helper, 
        SegmentationType type1, boolean useBinned1, 
        ISegmentedImageHelper img2Helper, 
        SegmentationType type2, boolean useBinned2, 
        float[] outputScaleRotTransXYStDev);
    
}
