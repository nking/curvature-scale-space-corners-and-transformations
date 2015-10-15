package algorithms.imageProcessing;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;

/**
 *
 * @author nichole
 */
public interface ISegmentedImageHelper {

    /**
     * apply histogram equalization to both image sets if statistics suggest 
     * they need it.
     * @param otherImgHelper
     * @return 
     */
    public boolean applyEqualizationIfNeededByComparison(ISegmentedImageHelper 
        otherImgHelper);

    /**
     * apply histogram equalization to this set of images.
     */
    public void applyEqualization();
    
    /**
     * apply segmentation type type to the greyscale or color binned image
     * depending upon the segmentation type.
     *
     * @param type
     * @param applyToBinnedImage
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    public void applySegmentation(SegmentationType type, boolean applyToBinnedImage) 
        throws IOException, NoSuchAlgorithmException;

    /**
     * create a binned greyscale image.
     * @param maxDimension 
     */
    public void createBinnedGreyscaleImage(int maxDimension);

    /**
     * generate the points of interest for each closed curve, 
     * needed for the matching algorithms.
     * 
     * @param type
     * @param applyToBinnedImage 
     */
    public void generatePerimeterPointsOfInterest(SegmentationType type, 
        boolean applyToBinnedImage);

    /**
     * convenience method that returns '1' if getForBinned is false else returns
     * binFactor
     *
     * @param getForBinned
     * @return
     */
    public int getBinFactor(boolean getForBinned);

    /**
     * get the binned segmented image of kind type.
     * @param type
     * @return 
     */
    public GreyscaleImage getBinnedSegmentationImage(SegmentationType type);

    /**
     * get the binned or full size greyscale image.
     * 
     * @param getTheBinned
     * @return 
     */
    public GreyscaleImage getGreyscaleImage(boolean getTheBinned);

    /**
     * get the full size greyscale image.
     * @return 
     */
    public GreyscaleImage getGreyscaleImage();

    /**
     * get the binned greyscale image.
     * @return 
     */
    public GreyscaleImage getGreyscaleImageBinned();
    
    /**
     * get the original image given to this instance.
     * @return 
     */
    public ImageExt getImage();

    /**
     * get the upper group size limit given to this instance for the full size
     * image
     * @return 
     */
    public int getLargestGroupLimit();

    /**
     * get the upper group size limit given to this instance for the binned size
     * image
     * @return 
     */
    public int getLargestGroupLimitBinned();

    /**
     * get the segmented image of kind type.
     * @param type
     * @return 
     */
    public GreyscaleImage getSegmentationImage(SegmentationType type);

    /**
     * get the lower group size limit given to this instance for the full size
     * image
     * @return 
     */
    public int getSmallestGroupLimit();

    /**
     * get the lower group size limit given to this instance for the binned size
     * image
     * @return 
     */
    public int getSmallestGroupLimitBinned();

    /**
     * sum the number of points of interest over all curves from segmentation
     * type segmentationType
     * @param segmentationType
     * @param useBinned1
     * @return 
     */
    public int sumPointsOfInterest(SegmentationType segmentationType, 
        boolean useBinned1);
    
}
