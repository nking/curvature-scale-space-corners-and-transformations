package algorithms.imageProcessing;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.logging.Logger;

/**
 * determine scale between 2 images using blob contours or corners.
 * NOT READY FOR USE YET.
 *
 * @author nichole
 */
public class BlobScaleFinderWrapper {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = true;

    protected final int binnedImageMaxDimension = 300;

    /*
    choices for solutions:

    (1) contours and curvature scale space matches followed by features to refine
        and validate
    (2) corners and feature matches w/ combinations
    (3) corners and feature matches in simplest ordered pairings
    (4) contours and curvature scale space matches using simplest ordered
        pairings

    when finished, should try (3) and/or (4) first then (2)
    */
    private enum AlgType {
        CONTOURS_ORDERED,
        CORNERS_UNORDERED, CONTOURS_UNORDERED
    }
    protected AlgType algType = AlgType.CONTOURS_ORDERED;

    protected final BlobPerimeterHelper img1Helper;

    protected final BlobPerimeterHelper img2Helper;

    protected BlobCornerHelper blobCornerHelper1 = null;
    protected BlobCornerHelper blobCornerHelper2 = null;

    protected BlobContourHelper blobContourHelper1 = null;
    protected BlobContourHelper blobContourHelper2 = null;

    protected final IntensityFeatures features1;
    protected final IntensityFeatures featuresBinned1;
    protected final IntensityFeatures features2;
    protected final IntensityFeatures featuresBinned2;
    
    private boolean useSameSegmentation = false;

    /**
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     */
    public BlobScaleFinderWrapper(ImageExt img1, ImageExt img2) {

        img1Helper = new BlobPerimeterHelper(img1, "1");

        img2Helper = new BlobPerimeterHelper(img2, "2");

        img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);

        img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);

        features1 = new IntensityFeatures(img1, 5, true);

        featuresBinned1 = new IntensityFeatures(
            img1Helper.getGreyscaleImageBinned(), 5, true);

        features2 = new IntensityFeatures(img2, 5, true);

        featuresBinned2 = new IntensityFeatures(
            img2Helper.getGreyscaleImageBinned(), 5, true);
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
         imageSegmentation.applyUsingKMPP(img0, 2);
         imageProcessor.applyAdaptiveMeanThresholding(img0, 20/binFactor);
     else:
         img0 = imageProcessor.binImage(img0, binFactor);
         img0 = imageSegmentation.applyUsingCIEXYPolarTheta(img, 4)
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

        /*
        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(
            img1Helper.getGreyscaleImage(), true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(
            img2Helper.getGreyscaleImage(), true);
        if (debug) {
            log.info(stats1.toString());
            log.info(stats2.toString());
        }
        */

        ImageStatistics statsR1 = ImageStatisticsHelper.examine(
            img1Helper.getImage().r, true);
        ImageStatistics statsB1 = ImageStatisticsHelper.examine(
            img1Helper.getImage().b, true);
        ImageStatistics statsG1 = ImageStatisticsHelper.examine(
            img1Helper.getImage().g, true);

        ImageStatistics statsR2 = ImageStatisticsHelper.examine(
            img2Helper.getImage().r, true);
        ImageStatistics statsB2 = ImageStatisticsHelper.examine(
            img2Helper.getImage().b, true);
        ImageStatistics statsG2 = ImageStatisticsHelper.examine(
            img2Helper.getImage().g, true);

        log.info("stats R1=" + statsR1.toString());
        log.info("stats G1=" + statsG1.toString());
        log.info("stats B1=" + statsB1.toString());

        log.info("stats R2=" + statsR2.toString());
        log.info("stats G2=" + statsG2.toString());
        log.info("stats B2=" + statsB2.toString());

        int limit = 20;
        useSameSegmentation = false;
        if ((Math.abs(statsR1.getMode() - statsR2.getMode()) < limit) &&
            (Math.abs(statsG1.getMode() - statsG2.getMode()) < limit) &&
            (Math.abs(statsB1.getMode() - statsB2.getMode()) < limit) &&
            (Math.abs(statsR1.getMedian() - statsR2.getMedian()) < limit) &&
            (Math.abs(statsG1.getMedian() - statsG2.getMedian()) < limit) &&
            (Math.abs(statsB1.getMedian() - statsB2.getMedian()) < limit)) {
            useSameSegmentation = true;
        }
        
        TransformationParameters params = null;
        
        params = calculateScaleImpl();
        
        if (params == null) {
            algType = AlgType.CORNERS_UNORDERED;
            params = calculateScaleImpl();
        }
        
        if (params == null) {
            algType = AlgType.CONTOURS_UNORDERED;
            params = calculateScaleImpl();
        }
        
        return params;
    }

    private TransformationParameters calculateScaleImpl() throws IOException,
        NoSuchAlgorithmException {
        
        /*
        depending on image statistics, different combinations of segmentation
        and binning are tried.

        If a segmentation type involves using random, the method might be
        tried again if specified.
        A different algorithm for the clustering may be needed.
        */

        /*
        SegmentationOrder(SegmentationType sType, int numExtraBinnedAllowed,
            int numExtraUnbinnedAllowed)
        */

        SegmentationOrder[] seg1 = new SegmentationOrder[]{
            new SegmentationOrder(SegmentationType.COLOR_POLARCIEXY, 0, 0),
            new SegmentationOrder(SegmentationType.GREYSCALE_KMPP, 0, 0),
            //new SegmentationOrder(SegmentationType.BINARY, 1, 1)
        };
        SegmentationOrder[] seg2 = new SegmentationOrder[]{
            new SegmentationOrder(SegmentationType.COLOR_POLARCIEXY, 0, 0),
            new SegmentationOrder(SegmentationType.GREYSCALE_KMPP, 0, 0),
            //new SegmentationOrder(SegmentationType.BINARY, 1, 1)
        };

        //TODO: currently caching segmentation, but because of the random
        // portion of some of the methods, need ability to redo segmentation
        // sometimes
        
        int ordered1Idx = 0;
        int ordered2Idx = 0;

        while ((ordered1Idx < seg1.length) && (ordered2Idx < seg2.length)) {

            boolean useBinned1 = seg1[ordered1Idx].currentIsBinned();

            boolean useBinned2 = seg2[ordered2Idx].currentIsBinned();

            SegmentationType segmentationType1 = seg1[ordered1Idx].geSegmentationType();

            SegmentationType segmentationType2 = seg2[ordered2Idx].geSegmentationType();

            log.info("for 1: " + segmentationType1.name() + " binned=" + useBinned1
                + " useSameSegmentation=" + useSameSegmentation
                + " ordered1Idx=" + ordered1Idx);
            log.info("for 2: " + segmentationType2.name() + " binned=" + useBinned2
                + " ordered2Idx=" + ordered2Idx);

            IntensityFeatures f1;
            IntensityFeatures f2;

            if (useBinned1) {
                img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
                f1 = featuresBinned1;
            } else {
                f1 = features1;
            }

            if (useBinned2) {
                img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
                f2 = featuresBinned2;
            } else {
                f2 = features2;
            }

            img1Helper.applySegmentation(segmentationType1, useBinned1);

            img2Helper.applySegmentation(segmentationType2, useBinned2);

            float[] outputScaleRotTransXYStDev = new float[4];

            TransformationParameters params = null;

            int n1 = 0;
            int n2 = 0;

            if (algType.equals(AlgType.CONTOURS_ORDERED) || 
                algType.equals(AlgType.CONTOURS_UNORDERED)) {
                
                if (blobContourHelper1 == null) {
                    blobContourHelper1 = new BlobContourHelper(img1Helper);
                    blobContourHelper2 = new BlobContourHelper(img2Helper);
                }

                blobContourHelper1.generatePerimeterContours(
                    segmentationType1, useBinned1);
                blobContourHelper2.generatePerimeterContours(
                    segmentationType2, useBinned2);
                
            } else if (algType.equals(AlgType.CORNERS_UNORDERED)) {

                if (blobCornerHelper1 == null) {
                    blobCornerHelper1 = new BlobCornerHelper(img1Helper);
                    blobCornerHelper2 = new BlobCornerHelper(img2Helper);
                }

                blobCornerHelper1.generatePerimeterCorners(
                    segmentationType1, useBinned1);
                blobCornerHelper2.generatePerimeterCorners(
                    segmentationType2, useBinned2);

            }
            
            if (algType.equals(AlgType.CONTOURS_ORDERED)) {
                
                BlobContoursScaleFinder0 bsFinder = new BlobContoursScaleFinder0();

                params = bsFinder.solveForScale(blobContourHelper1, f1,
                    segmentationType1, useBinned1, blobContourHelper2, f2,
                    segmentationType2, useBinned2, outputScaleRotTransXYStDev);

                n1 = blobContourHelper1.sumPointsOfInterest(segmentationType1, useBinned1);
                n2 = blobContourHelper2.sumPointsOfInterest(segmentationType2, useBinned2);
                
            } else if (algType.equals(AlgType.CORNERS_UNORDERED)) {
                
                BlobCornersScaleFinder bsFinder = new BlobCornersScaleFinder();

                params = bsFinder.solveForScale(blobCornerHelper1, f1,
                    segmentationType1, useBinned1, blobCornerHelper2, f2,
                    segmentationType2, useBinned2, outputScaleRotTransXYStDev);

                n1 = blobCornerHelper1.sumPointsOfInterest(segmentationType1, useBinned1);
                n2 = blobCornerHelper2.sumPointsOfInterest(segmentationType2, useBinned2);

            } else if (algType.equals(AlgType.CONTOURS_UNORDERED)) {

                BlobContoursScaleFinder bsFinder = new BlobContoursScaleFinder();

                params = bsFinder.solveForScale(blobContourHelper1, f1,
                    segmentationType1, useBinned1, blobContourHelper2, f2,
                    segmentationType2, useBinned2, outputScaleRotTransXYStDev);

                n1 = blobContourHelper1.sumPointsOfInterest(segmentationType1, useBinned1);
                n2 = blobContourHelper2.sumPointsOfInterest(segmentationType2, useBinned2);
                
            }

            if (params != null) {

                log.info("params for type"
                    + " (" + segmentationType1.name() + ", binned=" + useBinned1 + ")"
                    + " (" + segmentationType2.name() + ", binned=" + useBinned2 + ")"
                    + " : " + params.toString());

                log.info(String.format(
                    "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                    outputScaleRotTransXYStDev[0], outputScaleRotTransXYStDev[1],
                    outputScaleRotTransXYStDev[2], outputScaleRotTransXYStDev[3]));

                // TODO: consider returning the number of points used in the
                // calculation
                float factor0 = (outputScaleRotTransXYStDev[0]/params.getScale());
                float factor1 = (float)(2.*Math.PI/outputScaleRotTransXYStDev[1]);

                // consider comparing stdev in translations to a fraction of the image
                float factor2 = outputScaleRotTransXYStDev[2];
                float factor3 = outputScaleRotTransXYStDev[3];

                //TODO: review these limits
                if ((factor0 < 0.2) && (factor1 >= 18.) && (factor2 < 100)
                    && (factor3 < 100)) {

                    return params;
                }
            }

            // if arrive here, have to decide to keep current segmentation and
            // binning or increment.  at least one index has to change

            log.info("for 1: " + segmentationType1.name() + " binned=" + useBinned1
                + " nC1=" + n1);
            log.info("for 2: " + segmentationType2.name() + " binned=" + useBinned2
                + " nC2=" + n2);

            if (useSameSegmentation) {
                /*if (
                    ((useBinned1 && (nContours1 < 3)) && (useBinned2 && (nContours2 < 10))) ||
                    ((useBinned1 && (nContours1 < 10)) && (useBinned2 && (nContours2 < 3)))
                    ) {
                    // if this is the last segmentation to try, do not skip...
                    // TODO: improve segmentation types
                    if (ordered1Idx < (seg1.length - 1)) {
                        seg1[ordered1Idx].setToSkip();
                        seg2[ordered2Idx].setToSkip();
                        ordered1Idx++;
                        ordered2Idx++;
                        continue;
                    }
                }*/
                boolean incr1 = !seg1[ordered1Idx].incrementAndHasNext();
                boolean incr2 = !seg2[ordered2Idx].incrementAndHasNext();
                if (incr1 || incr2) {
                    ordered1Idx++;
                    ordered2Idx++;
                }
                continue;
            }

            if (n1 > 10) {
                if (n2 > 10) {
                    if (n1 > n2) {
                        if (!seg2[ordered2Idx].incrementAndHasNext()) {
                            ordered2Idx++;
                        }
                    } else {
                        if (!seg1[ordered1Idx].incrementAndHasNext()) {
                            ordered1Idx++;
                        }
                    }
                } else {
                    if (!seg1[ordered1Idx].incrementAndHasNext()) {
                        ordered1Idx++;
                    }
                }
                continue;
            }

            if (n2 > 10) {
                if (!seg1[ordered1Idx].incrementAndHasNext()) {
                    ordered1Idx++;
                }
                continue;
            }

            if (!seg1[ordered1Idx].incrementAndHasNext()) {
                ordered1Idx++;
            }
            if (!seg2[ordered2Idx].incrementAndHasNext()) {
                ordered2Idx++;
            }
        }

        return null;
    }

    private class SegmentationOrder {
        private final SegmentationType type;
        private int numBinnedAttempts = 0;
        private int numUnbinnedAttempts = 0;
        private final int numExtraBinnedAllowed;
        private final int numExtraUnbinnedAllowed;
        private boolean currentIsBinned = true;

        /** flag to use when a quick binned attempt shows that this segmentation
        is not the best choice for image and the full segmentation should be skipped*/
        private boolean skip = false;

        public SegmentationOrder(SegmentationType sType, int numExtraBinnedAllowed,
            int numExtraUnbinnedAllowed) {
            this.type = sType;
            this.numExtraBinnedAllowed = numExtraBinnedAllowed;
            this.numExtraUnbinnedAllowed = numExtraUnbinnedAllowed;
        }

        public boolean incrementAndHasNext() {
            if (skip) {
                return false;
            }
            if (currentIsBinned && (numBinnedAttempts < numExtraBinnedAllowed)) {
                numBinnedAttempts++;
                return true;
            } else if (currentIsBinned) {
                numBinnedAttempts++;
                currentIsBinned = false;
                return true;
            } else if (!currentIsBinned && (numUnbinnedAttempts < numExtraUnbinnedAllowed)) {
                numUnbinnedAttempts++;
                return true;
            }
            return false;
        }

        public void setToSkip() {
            skip = true;
        }

        public boolean currentIsBinned() {
            return currentIsBinned;
        }

        public SegmentationType geSegmentationType() {
            return type;
        }
    }

}
