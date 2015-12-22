package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.HoughTransform.HoughTransformLines;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * determine scale between 2 images using blob contours or corners.
 * NOT READY FOR USE YET.
 *
 * @author nichole
 */
public class BlobScaleFinderWrapper {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected final int binnedImageMaxDimension = 512;

    /*
    choices for solutions:

    (1) contours and curvature scale space matches followed by features to refine
        and validate
    (2) corners and feature matches w/ combinations
    (3) corners and feature matches in simplest ordered pairings
    (4) contours and curvature scale space matches using simplest ordered
        pairings

    should try (3) and/or (4) first then (2) and (1)
    */
    public static enum AlgType {
        CORNERS_COMBINATIONS, CONTOURS_COMBINATIONS
    }
    protected AlgType algType = AlgType.CORNERS_COMBINATIONS;

    protected final BlobPerimeterHelper img1Helper;
    protected final BlobPerimeterHelper img2Helper;

    protected BlobCornerHelper blobCornerHelper1 = null;
    protected BlobCornerHelper blobCornerHelper2 = null;

    protected BlobContourHelper blobContourHelper1 = null;
    protected BlobContourHelper blobContourHelper2 = null;

    // use with img1Helper.getImage() or getGreyscaleImage(), but not both
    protected final IntensityFeatures features1;
    // use img2Helper.getGreyscaleImageBinned(), 
    protected final IntensityFeatures2 featuresBinned1;
    protected final IntensityFeatures features2;
    protected final IntensityFeatures2 featuresBinned2;

    private final FeatureMatcherSettings settings;
    
    private int binFactor1 = 1;
    private int binFactor2 = 1;
    
    private boolean useSameSegmentation = false;

    private MatchingSolution solution = null;
    private AlgType solutionAlgType = null;
    private SegmentationType solutionSegmentationType1 = null;
    private SegmentationType solutionSegmentationType2 = null;
    private boolean solutionUsedBinned1 = false;
    private boolean solutionUsedBinned2 = false;
    
    /**
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     * @param settings
     * @param rotatedOffsets
     */
    public BlobScaleFinderWrapper(ImageExt img1, ImageExt img2, 
        FeatureMatcherSettings settings, RotatedOffsets rotatedOffsets) {

        this.settings = settings.copy();
        
        if (settings.debug()) {
            
            img1Helper = new BlobPerimeterHelper(img1, settings.getDebugTag() + "_1");

            img2Helper = new BlobPerimeterHelper(img2, settings.getDebugTag() + "_2");
            
        } else {
            
            img1Helper = new BlobPerimeterHelper(img1);

            img2Helper = new BlobPerimeterHelper(img2);
        }
                
        features1 = new IntensityFeatures(5, settings.useNormalizedFeatures(),
            rotatedOffsets);

        features2 = new IntensityFeatures(5, settings.useNormalizedFeatures(),
            rotatedOffsets);
                        
        if (settings.startWithBinnedImages()) {
            
            img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            
            binFactor1 = img1Helper.getBinFactor();

            img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            
            binFactor2 = img2Helper.getBinFactor();
            
            featuresBinned1 = new IntensityFeatures2(5, 
                settings.useNormalizedFeatures(), rotatedOffsets);
            featuresBinned1.calculateGradientWithGreyscale(img1Helper.getGreyscaleImageBinned());
            
            featuresBinned2 = new IntensityFeatures2(5, 
                settings.useNormalizedFeatures(), rotatedOffsets);
            featuresBinned2.calculateGradientWithGreyscale(img2Helper.getGreyscaleImageBinned());
            
        } else {
            
            featuresBinned1 = null;
            
            featuresBinned2 = null;  
        }
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
            img1Helper.getImage().getRValues(), true);
        ImageStatistics statsB1 = ImageStatisticsHelper.examine(
            img1Helper.getImage().getBValues(), true);
        ImageStatistics statsG1 = ImageStatisticsHelper.examine(
            img1Helper.getImage().getGValues(), true);

        ImageStatistics statsR2 = ImageStatisticsHelper.examine(
            img2Helper.getImage().getRValues(), true);
        ImageStatistics statsB2 = ImageStatisticsHelper.examine(
            img2Helper.getImage().getBValues(), true);
        ImageStatistics statsG2 = ImageStatisticsHelper.examine(
            img2Helper.getImage().getGValues(), true);

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
        
        boolean[] useBinned = settings.startWithBinnedImages()? 
            new boolean[]{true, false} : new boolean[]{false};
        
        for (boolean ub : useBinned) {
            
            if (params == null) {
                algType = AlgType.CORNERS_COMBINATIONS;
                params = calculateScaleImpl(ub);
            }
            
            /*
            if (params == null) {
                algType = AlgType.CONTOURS_COMBINATIONS;
                params = calculateScaleImpl();
            }*/
        }
        
        return params;
    }

    private TransformationParameters calculateScaleImpl(boolean useBinned) 
        throws IOException, NoSuchAlgorithmException {
        
        /*
        depending on image statistics, different combinations of segmentation
        and binning are tried.
 
        ADAPTIVE_MEAN is a good quick segmentation algorithm (O(N)), but it produces
        many blobs, so the total calculation takes twice as long
        as some of the other methods (empirically derived...).
        
        TODO:  Need to reduce the space complexity of the images to be able to more
        easily cache all of these images and products.
        */

        SegmentationType[] seg1 = new SegmentationType[]{
            ////SegmentationType.COLOR_POLARCIEXY,
        //      SegmentationType.DT_CLUSTERING,
              //SegmentationType.GREYSCALE_KMPP,
                SegmentationType.GREYSCALE_WAVELET,
        //      SegmentationType.GREYSCALE_CANNY,
        //      SegmentationType.GREYSCALE_HIST,
        //      SegmentationType.COLOR_POLARCIEXY_LARGE,
            ////SegmentationType.ADAPTIVE_MEAN
        };
        SegmentationType[] seg2 = new SegmentationType[]{
            ////SegmentationType.COLOR_POLARCIEXY,
       //     SegmentationType.DT_CLUSTERING,
            //SegmentationType.GREYSCALE_KMPP,
            SegmentationType.GREYSCALE_WAVELET,
        //      SegmentationType.GREYSCALE_CANNY,
        //      SegmentationType.GREYSCALE_HIST,
       //     SegmentationType.COLOR_POLARCIEXY_LARGE,
            ////SegmentationType.ADAPTIVE_MEAN
        };
        
        int ordered1Idx = 0;
        int ordered2Idx = 0;

        while ((ordered1Idx < seg1.length) && (ordered2Idx < seg2.length)) {

            SegmentationType segmentationType1 = seg1[ordered1Idx];

            SegmentationType segmentationType2 = seg2[ordered2Idx];

            log.info("for 1: " + segmentationType1.name() + " alg=" + algType.name()
                + " binned=" + useBinned + " useSameSegmentation=" + useSameSegmentation
                + " ordered1Idx=" + ordered1Idx);
            log.info("for 2: " + segmentationType2.name()
                + " binned=" + useBinned + " ordered2Idx=" + ordered2Idx);

            IntensityFeatures f1;
            IntensityFeatures f2;

            if (useBinned) {
                img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
                f1 = featuresBinned1;
            } else {
                f1 = features1;
            }

            if (useBinned) {
                img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
                f2 = featuresBinned2;
            } else {
                f2 = features2;
            }

            long t0 = System.currentTimeMillis();
            
            img1Helper.applySegmentation(segmentationType1, useBinned);

            long t1 = System.currentTimeMillis();
            
            img2Helper.applySegmentation(segmentationType2, useBinned);

            long t2 = System.currentTimeMillis();
            
            long t1Sec = (t1 - t0)/1000;
            long t2Sec = (t2 - t1)/1000;
            log.info("segmentation1(sec)=" + t1Sec 
                + " segmentation2(sec)=" + t2Sec);
            
            MatchingSolution soln = null;

            int n1 = 0;
            int n2 = 0;

            t0 = System.currentTimeMillis();
            
            if (algType.equals(AlgType.CORNERS_COMBINATIONS)) {

                if (blobCornerHelper1 == null) {
                    if (settings.debug()) {
                        blobCornerHelper1 = new BlobCornerHelper(img1Helper, 
                            settings.getDebugTag() + "_!");
                        blobCornerHelper2 = new BlobCornerHelper(img2Helper, 
                            settings.getDebugTag() + "_2");
                    } else {
                        blobCornerHelper1 = new BlobCornerHelper(img1Helper);
                        blobCornerHelper2 = new BlobCornerHelper(img2Helper);
                    }
                }

                List<HoughTransformLines> houghTransformLines1 = 
                    findLinesUsingHoughTransform(blobCornerHelper1,
                    segmentationType1, useBinned);
                
                boolean hasManyIntersectingLines1 = 
                    hasManyIntersectingLines(houghTransformLines1);
                
                if (hasManyIntersectingLines1) {
                    // use canny edges segmentation to replace segmentationType1
                    img1Helper.replaceSegmentationWithCanny(segmentationType1, 
                        useBinned);
                    houghTransformLines1 = 
                        findLinesUsingHoughTransform(blobCornerHelper1,
                        segmentationType1, useBinned);
                }
                
                List<HoughTransformLines> houghTransformLines2 = 
                    findLinesUsingHoughTransform(blobCornerHelper2,
                    segmentationType2, useBinned);
                
                boolean hasManyIntersectingLines2 = 
                    hasManyIntersectingLines(houghTransformLines2);
                
                if (hasManyIntersectingLines2) {
                    // use canny edges segmentation to replace segmentationType2
                    img2Helper.replaceSegmentationWithCanny(segmentationType2, 
                        useBinned);
                    houghTransformLines2 = 
                        findLinesUsingHoughTransform(blobCornerHelper2,
                        segmentationType2, useBinned);
                }
                
                t0 = System.currentTimeMillis();
                
                blobCornerHelper1.generatePerimeterCorners(segmentationType1, 
                    useBinned);

                t1 = System.currentTimeMillis();
                
                blobCornerHelper2.generatePerimeterCorners(segmentationType2, 
                    useBinned);
                
                //if (hasManyIntersectingLines1) {
                    removeLineArtifactCorners(houghTransformLines1, blobCornerHelper1,
                        segmentationType1, useBinned);
                //}
                //if (hasManyIntersectingLines2) {
                    removeLineArtifactCorners(houghTransformLines2, blobCornerHelper2,
                        segmentationType2, useBinned);
                //}
        
                t2 = System.currentTimeMillis();
                t1Sec = (t1 - t0)/1000;
                t2Sec = (t2 - t1)/1000;
                 Logger.getLogger(this.getClass().getName()).info("corners1(sec)=" 
                    + t1Sec + " sec corners1(sec)=" + t2Sec + " sec");
            
                BlobCornersScaleFinder bsFinder = new BlobCornersScaleFinder();

                if (settings.debug()) {
                    bsFinder.setToDebug();
                }
        
                soln = bsFinder.solveForScale(blobCornerHelper1, f1,
                    segmentationType1, useBinned, blobCornerHelper2, f2,
                    segmentationType2, useBinned);

                n1 = blobCornerHelper1.sumPointsOfInterest(segmentationType1, 
                    useBinned);
                n2 = blobCornerHelper2.sumPointsOfInterest(segmentationType2, 
                    useBinned);

            } else if (algType.equals(AlgType.CONTOURS_COMBINATIONS)) {

                if (blobContourHelper1 == null) {
                    if (settings.debug()) {
                        blobContourHelper1 = new BlobContourHelper(img1Helper, 
                            settings.getDebugTag() + "_1");
                        blobContourHelper2 = new BlobContourHelper(img2Helper,
                            settings.getDebugTag() + "_2");
                    } else {
                        blobContourHelper1 = new BlobContourHelper(img1Helper);
                        blobContourHelper2 = new BlobContourHelper(img2Helper);
                    }
                }

                blobContourHelper1.generatePerimeterContours(
                    segmentationType1, useBinned);
                blobContourHelper2.generatePerimeterContours(
                    segmentationType2, useBinned);
                
                BlobContoursScaleFinder bsFinder = new BlobContoursScaleFinder();

                if (settings.debug()) {
                    bsFinder.setToDebug();
                }
                
                soln = bsFinder.solveForScale(blobContourHelper1, f1,
                    segmentationType1, useBinned, blobContourHelper2, f2,
                    segmentationType2, useBinned);

                n1 = blobContourHelper1.sumPointsOfInterest(segmentationType1, 
                    useBinned);
                n2 = blobContourHelper2.sumPointsOfInterest(segmentationType2, 
                    useBinned);
                
            }
       
            t1 = System.currentTimeMillis();
            t1Sec = (t1 - t0)/1000;
            Logger.getLogger(this.getClass().getName()).info("matching(sec)=" 
                + t1Sec);
                
            if (soln != null) {
                
                TransformationParameters params = soln.getParams();
                
                log.info("params for type"
                    + " (" + segmentationType1.name() + ", binned=" + useBinned + ")"
                    + " (" + segmentationType2.name() + ", binned=" + useBinned + ")"
                    + " : " + params.toString());

                log.info(String.format(
                    "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                    params.getStandardDeviations()[0], 
                    params.getStandardDeviations()[1],
                    params.getStandardDeviations()[2], 
                    params.getStandardDeviations()[3]));

                boolean small = MiscStats.standardDeviationsAreSmall(params);
                
                if (small) {

                    solutionAlgType = algType;
                    solutionSegmentationType1 = segmentationType1;
                    solutionSegmentationType2 = segmentationType2;
                    solutionUsedBinned1 = useBinned;
                    solutionUsedBinned2 = useBinned;
                    solution = soln;
                    
                    return params;
                }
            }

            // if arrive here, have to decide to keep current segmentation and
            // binning or increment.  at least one index has to change

            log.info("for 1: " + segmentationType1.name() + " binned=" + useBinned
                + " nC1=" + n1);
            log.info("for 2: " + segmentationType2.name() + " binned=" + useBinned
                + " nC2=" + n2);

            if (useSameSegmentation) {
                ordered1Idx++;
                ordered2Idx++;
                continue;
            }

            if (n1 > 10) {
                if (n2 > 10) {
                    if (n1 > n2) {
                        ordered2Idx++;
                    } else {
                        ordered1Idx++;
                    }
                } else {
                    ordered1Idx++;
                }
                continue;
            }

            if (n2 > 10) {
                ordered1Idx++;
                continue;
            }

            ordered1Idx++;
            ordered2Idx++;
        }

        return null;
    }
    
    private List<HoughTransformLines> findLinesUsingHoughTransform(
        BlobCornerHelper blobCornerHelper,
        SegmentationType segmentationType, boolean useBinnedImage) {
     
        List<PairIntArray> perimeterLists = blobCornerHelper.imgHelper.
            getBlobPerimeters(segmentationType, useBinnedImage);
                
        int imageWidth = useBinnedImage ? 
            blobCornerHelper.imgHelper.getGreyscaleImageBinned().getWidth() :
            blobCornerHelper.imgHelper.getGreyscaleImage().getWidth();
        
        int imageHeight = useBinnedImage ? 
            blobCornerHelper.imgHelper.getGreyscaleImageBinned().getHeight() :
            blobCornerHelper.imgHelper.getGreyscaleImage().getHeight();
        
        int thetaTol = 1;
        int radiusTol = 7;
        
        HoughTransform ht = new HoughTransform();
        
        List<HoughTransformLines> lineList = new ArrayList<HoughTransformLines>();
        
        for (int ii = 0; ii < perimeterLists.size(); ++ii) {

            // NOTE: in testable method for this, should allow ability to
            // pass in junctions and not delete corners that are in
            // junctions.
            // For these blob perimeters, there are not junctions.

            PairIntArray edge = perimeterLists.get(ii);

            Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap =
               ht.calculateLineGivenEdge(edge, imageWidth, imageHeight);

            List<PairInt> outSortedKeys = ht.sortByVotes(outputPolarCoordsPixMap);

            // === find indiv lines within the edge ====

            HoughTransformLines htl = ht.createPixTRMapsFromSorted(
                outSortedKeys, outputPolarCoordsPixMap,
                thetaTol, radiusTol);
            
            lineList.add(htl);
        }
        
        return lineList;
    }
    
    private boolean hasManyIntersectingLines(List<HoughTransformLines> lineList) {
        
        int sizeLimit = 15;
        
        // key={polar theta, radius}, value=number of lines with key
        Map<PairInt, Integer> lineMap = new HashMap<PairInt, Integer>();
                
        for (int ii = 0; ii < lineList.size(); ++ii) {

            HoughTransformLines htl = lineList.get(ii);
            
            Map<PairInt, PairInt> pixToTRMap = htl.getPixelToPolarCoordMap();
            List<Set<PairInt>> outputSortedGroups = htl.getSortedLineGroups();
            
            for (Set<PairInt> line : outputSortedGroups) {
                
                if (line.size() < sizeLimit) {
                    break;
                }
                
                PairInt aPoint = line.iterator().next();
                
                PairInt tr = pixToTRMap.get(aPoint);                
                
                Integer count = lineMap.get(tr);
                
                if (count == null) {
                    count = Integer.valueOf(1);
                } else {
                    count = Integer.valueOf(count.intValue() + 1);
                }
                
                lineMap.put(tr, count);
            }
        }
        
        double avg = 0;
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (Entry<PairInt, Integer> entry : lineMap.entrySet()) {
            int c = entry.getValue().intValue();
            avg += c;
            if (c < min) {
                min = c;
            } 
            if (c > max) {
                max = c;
            }
        }
        avg /= (double)lineMap.size();
        
        if ((avg > 1.5) && (max > 1)) {
            return true;
        }
        //log.info("avg=" + avg + " min=" + min + " max=" + max);
        
        return false;
    }
    
    private void removeLineArtifactCorners(List<HoughTransformLines> 
        houghTransformLines, BlobCornerHelper blobCornerHelper, 
        SegmentationType segmentationType, boolean useBinnedImage) {
        
        List<PairIntArray> perimeterLists = blobCornerHelper.imgHelper.
            getBlobPerimeters(segmentationType, useBinnedImage);
                
        List<List<CornerRegion>> cornerRegionLists = 
            blobCornerHelper.getPerimeterCorners(segmentationType, useBinnedImage);
            
        int imageWidth = useBinnedImage ? 
            blobCornerHelper.imgHelper.getGreyscaleImageBinned().getWidth() :
            blobCornerHelper.imgHelper.getGreyscaleImage().getWidth();
        
        int imageHeight = useBinnedImage ? 
            blobCornerHelper.imgHelper.getGreyscaleImageBinned().getHeight() :
            blobCornerHelper.imgHelper.getGreyscaleImage().getHeight();
        
        int thetaTol = 1;
        int radiusTol = 7;

        //use hough transform for lines to remove corners from line artifacts
        CornerCorrector.removeCornersFromLineArtifacts(houghTransformLines,
            perimeterLists, cornerRegionLists, thetaTol, radiusTol, imageWidth, 
            imageHeight);
    }

    public MatchingSolution getSolution() {
        return solution;
    }
    public AlgType getSolutionAlgType() {
        return solutionAlgType;
    }
    public SegmentationType getSolutionSegmentationType1() {
        return solutionSegmentationType1;
    }
    public SegmentationType getSolutionSegmentationType2() {
        return solutionSegmentationType2;
    }
    public boolean getSolutionUsedBinned1() {
        return solutionUsedBinned1;
    }
    public boolean getSolutionUsedBinned2() {
        return solutionUsedBinned2;
    }
    
    public int getBinFactor1() {
        return binFactor1;
    }
    
    public int getBinFactor2() {
        return binFactor2;
    }
    
    public List<List<CornerRegion>> getAllCornerRegions1OfSolution() {
        
        if (algType.equals(AlgType.CONTOURS_COMBINATIONS)) {
            return null;
        }
        
        return blobCornerHelper1.generatePerimeterCorners(
            solutionSegmentationType1, solutionUsedBinned1);
    }
    
    public List<List<CornerRegion>> getAllCornerRegions2OfSolution() {
        
        if (algType.equals(AlgType.CONTOURS_COMBINATIONS)) {
            return null;
        }
        
        return blobCornerHelper2.generatePerimeterCorners(
            solutionSegmentationType2, solutionUsedBinned2);
    }
    
    public List<List<BlobPerimeterRegion>> getAllBlobRegions1OfSolution() {
        
        if (algType.equals(AlgType.CORNERS_COMBINATIONS)) {
            return null;
        }
        
        return blobContourHelper1.generatePerimeterRegions(
            solutionSegmentationType1, solutionUsedBinned1);
    }
    
    public List<List<BlobPerimeterRegion>> getAllBlobRegions2OfSolution() {
        
        if (algType.equals(AlgType.CORNERS_COMBINATIONS)) {
            return null;
        }
        
        return blobContourHelper2.generatePerimeterRegions(
            solutionSegmentationType2, solutionUsedBinned2);
    }
    
    public List<List<CurvatureScaleSpaceContour>> getAllContours1OfSolution() {
        
        if (algType.equals(AlgType.CORNERS_COMBINATIONS)) {
            return null;
        }
        
        return blobContourHelper1.generatePerimeterContours(
            solutionSegmentationType1, solutionUsedBinned1);
    }
    
    public List<List<CurvatureScaleSpaceContour>> getAllContours2OfSolution() {
        
        if (algType.equals(AlgType.CORNERS_COMBINATIONS)) {
            return null;
        }
        
        return blobContourHelper2.generatePerimeterContours(
            solutionSegmentationType2, solutionUsedBinned2);
    }
    
    public IntensityFeatures getSolutionFeatures1() {
        if (solutionUsedBinned1) {
            return featuresBinned1;
        }
        return features1;
    }
    
    public IntensityFeatures getSolutionFeatures2() {
        if (solutionUsedBinned2) {
            return featuresBinned2;
        }
        return features2;
    }
    
}
