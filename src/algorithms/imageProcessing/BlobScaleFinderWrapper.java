package algorithms.imageProcessing;

import algorithms.imageProcessing.SegmentedImageHelper.SegmentationType;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * determine scale between 2 images using blob contours.
 * NOT READY FOR USE YET.
 *
 * @author nichole
 */
public class BlobScaleFinderWrapper {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = true;

    protected final int binnedImageMaxDimension = 300;

    protected final SegmentedImageHelper img1Helper;

    protected final SegmentedImageHelper img2Helper;

    /**
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     */
    public BlobScaleFinderWrapper(ImageExt img1, ImageExt img2) {

        img1Helper = new SegmentedImageHelper(img1, "1");

        img2Helper = new SegmentedImageHelper(img2, "2");

    }

    public void setToDebug() {
        debug = true;
    }

    /**
     * NOT READY FOR USE YET.
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.
     *
     * This method does not require pre-processing such as sky subtraction
     * because it uses adaptive mean thresholding, but if sky subtraction is
     * already performed, you might want to use the alternate method
     * calculateScale0().
     *
     * Note that it is expected that this transformation result will be followed
     * by a more rigorous solver such as the FeatureMatcher for a correspondence
     * list (and a better Euclidean transform) to be used in.

     <pre>
     The blobs are found through two different ways depending upon the image
     statistics.
     If the image appears to be very bright, a method which is better at finding
     dark blobs is used:
         img0 = img.copyToGreyscale();
         img0 = imageProcessor.binImage(img0, binFactor);
         imageProcessor.applyImageSegmentation(img0, 2);
         imageProcessor.applyAdaptiveMeanThresholding(img0, 20/binFactor);
     else:
         img0 = imageProcessor.binImage(img0, binFactor);
         img0 = imageProcessor.createGreyscaleFromColorSegmentation(img, 4)
         imageProcessor.applyAdaptiveMeanThresholding(img0, 2); 2 is for unbinned so may need tuning
     </pre>

     * @return Euclidean scale to be applied to image1 to place it in the same
     * scale reference frame as image2.  Rotation and transformation are also
     * roughly solved for.
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public TransformationParameters calculateScale() throws IOException,
        NoSuchAlgorithmException {

        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(
            img1Helper.getGreyscaleImage(), true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(
            img2Helper.getGreyscaleImage(), true);
        
        if (debug) {
            log.info(stats1.toString());
            log.info(stats2.toString());
        }

        /*
        depending on image statistics, different combinations of segmentation
        and binning are tried.
        
        If a segmentation type involves using random, the method might be
        tried again if specified.  
        A different algorithm for the clustering may be needed.
        */

        SegmentationType[] orderOfSeg1;
        boolean[] orderOfBinning1;
        int[] numTriesAllowed1;
        int[] numTries1;
        SegmentationType[] orderOfSeg2;
        boolean[] orderOfBinning2;
        int[] numTriesAllowed2;
        int[] numTries2;
        
        orderOfSeg1 = new SegmentationType[]{
            SegmentationType.BINARY, SegmentationType.COLOR_POLARCIEXY_ADAPT};
        orderOfBinning1 = new boolean[] {false, false};
        numTriesAllowed1 = new int[]{2, 2};
        numTries1 = new int[]{0, 0};
        
        orderOfSeg2 = new SegmentationType[]{
            SegmentationType.BINARY, SegmentationType.COLOR_POLARCIEXY_ADAPT};
        orderOfBinning2 = new boolean[] {false, false};
        numTriesAllowed2 = new int[]{2, 2};
        numTries2 = new int[]{0, 0};
        
        /*
        orderOfSeg1 = new SegmentationType[]{SegmentationType.COLOR_POLARCIEXY_ADAPT};
        orderOfBinning1 = new boolean[] {false};
        numTriesAllowed1 = new int[]{2};
        numTries1 = new int[]{0};
        orderOfSeg2 = new SegmentationType[]{SegmentationType.COLOR_POLARCIEXY_ADAPT};
        orderOfBinning2 = new boolean[] {false};
        numTriesAllowed2 = new int[]{2};
        numTries2 = new int[]{0};
        */
        
        assert(orderOfSeg1.length == orderOfBinning1.length);
        assert(orderOfSeg1.length == numTries1.length);
        assert(orderOfSeg1.length == numTriesAllowed1.length);
        assert(orderOfSeg2.length == orderOfBinning2.length);
        assert(orderOfSeg2.length == numTries2.length);
        assert(orderOfSeg2.length == numTriesAllowed2.length);
        
        int ordered1Idx = 0;
        int ordered2Idx = 0;

        while ((ordered1Idx < orderOfSeg1.length) && (ordered2Idx < orderOfSeg2.length)) {

            boolean useBinned1 = orderOfBinning1[ordered1Idx];

            boolean useBinned2 = orderOfBinning2[ordered2Idx];

            SegmentationType segmentationType1 = orderOfSeg1[ordered1Idx];

            SegmentationType segmentationType2 = orderOfSeg2[ordered2Idx];
            
            log.info("for 1: " + segmentationType1.name() + " binned=" + useBinned1);
            log.info("for 2: " + segmentationType2.name() + " binned=" + useBinned2);

            if (useBinned1) {
                img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            }

            if (useBinned2) {
                img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            }

            img1Helper.applySegmentationToBinned(segmentationType1, useBinned1);
            
            img2Helper.applySegmentationToBinned(segmentationType2, useBinned2);

            img1Helper.extractBlobsAndContours(segmentationType1, useBinned1);

            img2Helper.extractBlobsAndContours(segmentationType2, useBinned2);

            BlobScaleFinder bsFinder = new BlobScaleFinder();

            float[] outputScaleRotTransXYStDev = new float[4];
            TransformationParameters params = bsFinder.solveForScale(
                img1Helper, segmentationType1, useBinned1,
                img2Helper, segmentationType2, useBinned2,
                outputScaleRotTransXYStDev);

            if (params != null) {

                log.info("params for type"
                    + " (" + segmentationType1.name() + ", binned=" + useBinned1 + ")"
                    + " (" + segmentationType2.name() + ", binned=" + useBinned2 + ")"
                    + " : " + params.toString());

                log.info(String.format(
                    "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                    outputScaleRotTransXYStDev[0], outputScaleRotTransXYStDev[1],
                    outputScaleRotTransXYStDev[2], outputScaleRotTransXYStDev[3]));

                //TODO: review this limit
                if (
                    ((outputScaleRotTransXYStDev[0]/params.getScale()) < 0.2)
                    ) {
                    return params;
                }
            }

            // if arrive here, have to decide to keep current segmentation and binning or increment.
            // at least one index has to change
            
            int nContours1 = sumContours(img1Helper, segmentationType1, useBinned1);

            int nContours2 = sumContours(img2Helper, segmentationType2, useBinned2);
            
            // new logic to repeat if process uses random and the num tries allowed is high enough:
            if (((numTries1[ordered1Idx] < numTriesAllowed1[ordered1Idx])
                && (nContours1 > 3)) ||
                ((numTries2[ordered2Idx] < numTriesAllowed2[ordered2Idx])
                && (nContours2 > 3))
                ) {
                                
                numTries1[ordered1Idx]++;
                numTries2[ordered2Idx]++;
                continue;
            }
            
            if (nContours1 > 10) {
                if (nContours2 > 10) {
                    if (nContours1 > nContours2) {
                        ordered2Idx++;
                    } else {
                        ordered1Idx++;
                    }
                } else {
                    ordered1Idx++;
                }
                continue;
            }
            
            if (nContours2 > 10) {
                ordered1Idx++;
                continue;
            }

            ordered1Idx++;
            ordered2Idx++;
        }
        
        return null;
    }

    /**
     * NOT READY FOR USE YET.
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.  Note that image processing
     * such as sky masks should be applied before using this method.
     * Also note that it is expected that it will be followed by a more rigorous
     * solver such as the FeatureMatcher for a correspondence list
     * (and a better Euclidean transform) to be used in.

     This method does not use adaptive mean thresholding and was originally 
     created for use on images where background processing such as sky masking 
     has already occurred.

     * @return Euclidean scale to be applied to image1 to place it in the same
     * scale reference frame as image2.  Rotation and transformation are also
     * roughly solved for.
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public TransformationParameters calculateScale0() throws IOException,
        NoSuchAlgorithmException {

        /*
        Tries to solve using various combinations of binning and segmentation,
        starting with binned images because the solution if possible with them
        is faster.
        (1) Greyscale binned, k=2
        (2) Clr binned, k=2
            -- possibly watershed on blobs if no solution
        (3) Greyscale not binned, k=2
        (4) Clr binned, k=2
            -- possibly watershed on blobs if no solution
        (5) Clr not binned, k=3 or 4
        */

        SegmentationType[] orderOfSeg1 = new SegmentationType[]{
            SegmentationType.GREYSCALE_KMPP, 
            SegmentationType.COLOR_POLARCIEXY,
            SegmentationType.GREYSCALE_KMPP, 
            SegmentationType.COLOR_POLARCIEXY};
        boolean[] orderOfBinning1 = new boolean[] {true, true, false, false};
        
        boolean didApplyHistEq =
            img1Helper.applyEqualizationIfNeededByComparison(img2Helper);
        log.info("didApplyHistEq=" + didApplyHistEq);
        
        SegmentationType[] orderOfSeg2 = Arrays.copyOf(orderOfSeg1, orderOfSeg1.length);
        boolean[] orderOfBinning2 = Arrays.copyOf(orderOfBinning1, orderOfBinning1.length);

        int ordered1Idx = 0;
        int ordered2Idx = 0;

        while ((ordered1Idx < orderOfSeg1.length) && (ordered2Idx < orderOfSeg2.length)) {

            boolean useBinned1 = orderOfBinning1[ordered1Idx];

            boolean useBinned2 = orderOfBinning2[ordered2Idx];

            SegmentationType segmentationType1 = orderOfSeg1[ordered1Idx];

            SegmentationType segmentationType2 = orderOfSeg2[ordered2Idx];

            if (useBinned1) {
                img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            }

            if (useBinned2) {
                img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            }

            img1Helper.applySegmentationToBinned(segmentationType1, useBinned1);

            img2Helper.applySegmentationToBinned(segmentationType2, useBinned2);

            img1Helper.extractBlobsAndContours(segmentationType1, useBinned1);

            img2Helper.extractBlobsAndContours(segmentationType2, useBinned2);

            BlobScaleFinder bsFinder = new BlobScaleFinder();

            float[] outputScaleRotTransXYStDev = new float[4];
            TransformationParameters params = bsFinder.solveForScale(
                img1Helper, segmentationType1, useBinned1,
                img2Helper, segmentationType2, useBinned2,
                outputScaleRotTransXYStDev);

            if (params != null) {

                log.info("params for type"
                    + " (" + segmentationType1.name() + ", binned=" + useBinned1 + ")"
                    + " (" + segmentationType2.name() + ", binned=" + useBinned2 + ")"
                    + " : " + params.toString());

                log.info(String.format(
                    "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                    outputScaleRotTransXYStDev[0], outputScaleRotTransXYStDev[1],
                    outputScaleRotTransXYStDev[2], outputScaleRotTransXYStDev[3]));

                //TODO: review this limit
                if (
                    ((outputScaleRotTransXYStDev[0]/params.getScale()) < 0.2)
                    ) {
                    return params;
                }
            }

            // if arrive here, have to decide to keep current segmentation and binning or increment.
            // at least one index has to change
            
            int nContours1 = sumContours(img1Helper, segmentationType1, useBinned1);

            int nContours2 = sumContours(img2Helper, segmentationType2, useBinned2);
            
            if (nContours1 > 10) {
                if (nContours2 > 10) {
                    if (nContours1 > nContours2) {
                        ordered2Idx++;
                    } else {
                        ordered1Idx++;
                    }
                } else {
                    ordered1Idx++;
                }
                continue;
            }
            
            if (nContours2 > 10) {
                ordered1Idx++;
                continue;
            }

            ordered1Idx++;
            ordered2Idx++;
        }

        return null;
    }

    private int sumContours(SegmentedImageHelper imgHelper, 
        SegmentationType segmentationType, boolean useBinned) {
        
        BlobsAndContours bc = imgHelper.getBlobsAndContours(segmentationType, useBinned);
        
        int n = 0;
        
        for (List<CurvatureScaleSpaceContour> list : bc.getContours()) {
            n += list.size();
        }
        
        return n;
    }

}
