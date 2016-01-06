package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.HoughTransform.HoughTransformLines;
import algorithms.compGeometry.RotatedOffsets;
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
 * create lists of matched points present in 2 images.
 * 
 * @author nichole
 */
public class NonEuclideanSegmentFeatureMatcher {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected final int binnedImageMaxDimension = 512;

    protected final BlobPerimeterCornerHelper img1Helper;
    protected final BlobPerimeterCornerHelper img2Helper;

    // use with img1Helper.getImage() or getGreyscaleImage(), but not both
    protected final IntensityFeatures features1;
    // use img2Helper.getGreyscaleImageBinned(), 
    protected final IntensityFeatures featuresBinned1;
    protected final IntensityFeatures features2;
    protected final IntensityFeatures featuresBinned2;

    private final FeatureMatcherSettings settings;
    
    private boolean useSameSegmentation = false;
    
    private List<FeatureComparisonStat> solutionStats = null;
    
    private List<PairInt> solutionMatched1 = null;
    
    private List<PairInt> solutionMatched2 = null;
    
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
    public NonEuclideanSegmentFeatureMatcher(ImageExt img1, ImageExt img2, 
        FeatureMatcherSettings settings, RotatedOffsets rotatedOffsets) {

        this.settings = settings.copy();
                
        if (settings.debug()) {
            
            img1Helper = new BlobPerimeterCornerHelper(img1, settings.getDebugTag() + "_1");

            img2Helper = new BlobPerimeterCornerHelper(img2, settings.getDebugTag() + "_2");
            
        } else {
            
            img1Helper = new BlobPerimeterCornerHelper(img1);

            img2Helper = new BlobPerimeterCornerHelper(img2);
        }
              
        // delaying creation of gradient images for full images until needed:
        features1 = new IntensityFeatures(5, settings.useNormalizedFeatures(),
            rotatedOffsets);

        features2 = new IntensityFeatures(5, settings.useNormalizedFeatures(),
            rotatedOffsets);
                        
        if (settings.startWithBinnedImages()) {
            
            img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            
            img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
                        
            featuresBinned1 = new IntensityFeatures(5, 
                settings.useNormalizedFeatures(), rotatedOffsets);
            featuresBinned1.calculateGradientWithGreyscale(img1Helper.getGreyscaleImageBinned());
            
            featuresBinned2 = new IntensityFeatures(5, 
                settings.useNormalizedFeatures(), rotatedOffsets);
            featuresBinned2.calculateGradientWithGreyscale(img2Helper.getGreyscaleImageBinned());
            
        } else {
            
            featuresBinned1 = null;
            
            featuresBinned2 = null;  
        }
    }

    public void match() throws IOException, NoSuchAlgorithmException {

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
        
        SegmentationType type = SegmentationType.GREYSCALE_WAVELET;
        if (settings.doOverrideWithCannySegmentation()) {
            type = SegmentationType.GREYSCALE_CANNY;
        }
                
        boolean[] useBinned = settings.startWithBinnedImages()? 
            new boolean[]{true, false} : new boolean[]{false};
        
        for (boolean ub : useBinned) {
            
            if (!ub) {
                if (!features1.gradientWasCreated()) {
                    features1.calculateGradientWithGreyscale(
                        img1Helper.getGreyscaleImage());
                }
                if (!features2.gradientWasCreated()) {
                    features2.calculateGradientWithGreyscale(
                        img2Helper.getGreyscaleImage());
                }
            }
            
            boolean solved = generateAndMatchCornerRegions(type, ub);
            
            if (solved) {
                break;
            }
        }        
    }

    private void prepareCorners(SegmentationType type,
        boolean useBinned) throws IOException, NoSuchAlgorithmException {
        
        log.info(type.name() + " binned=" + useBinned 
            + " useSameSegmentation=" + useSameSegmentation);

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

        img1Helper.applySegmentation(type, useBinned);

        img2Helper.applySegmentation(type, useBinned);

        List<HoughTransformLines> houghTransformLines1, houghTransformLines2;
        
        boolean useCanny = type.equals(SegmentationType.GREYSCALE_CANNY);
        
        if (useCanny) {
            
            boolean filterOutImageBoundaryBlobs = true;
            boolean filterOutZeroPixels = false;

            // pre-make the blobs using non-default variables:
            img1Helper.getBlobs(type, useBinned,
                filterOutImageBoundaryBlobs, filterOutZeroPixels);
            img1Helper.extractBlobPerimeterAsCornerRegions(
                type, useBinned);

            // pre-make the blobs using non-default variables:
            img2Helper.getBlobs(type, useBinned,
                filterOutImageBoundaryBlobs, filterOutZeroPixels);
            img2Helper.extractBlobPerimeterAsCornerRegions(type,
                useBinned);

        } else {
            img1Helper.extractBlobPerimeterAsCornerRegions(type, useBinned);
            img2Helper.extractBlobPerimeterAsCornerRegions(type, useBinned);
            
            img1Helper.generatePerimeterCorners(type, useBinned);
                
            img2Helper.generatePerimeterCorners(type, useBinned);
        }
        
        houghTransformLines1 = findLinesUsingHoughTransform(img1Helper, type, 
            useBinned);

        houghTransformLines2 = findLinesUsingHoughTransform(img2Helper, type, 
            useBinned);
        
        removeLineArtifactCorners(houghTransformLines1, img1Helper, type, 
            useBinned);
        
        removeLineArtifactCorners(houghTransformLines2, img2Helper, type,
            useBinned);
        
        
    }
    
    private boolean generateAndMatchCornerRegions(SegmentationType type,
        boolean useBinned) throws IOException, NoSuchAlgorithmException {
        
        log.info(type.name() + " binned=" + useBinned 
            + " useSameSegmentation=" + useSameSegmentation);
        
        prepareCorners(type, useBinned);
        
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        
        boolean matchCurveToCurve = false;
        
        if (matchCurveToCurve) {
            return matchCurveByCurve(type, useBinned, rotatedOffsets);
        } else {
            return matchAllPoints(type, useBinned, rotatedOffsets);
        }
        
    }
    
    private List<HoughTransformLines> findLinesUsingHoughTransform(
        BlobPerimeterCornerHelper blobCornerHelper,
        SegmentationType segmentationType, boolean useBinnedImage) {
     
        List<PairIntArray> perimeterLists = blobCornerHelper.
            getBlobPerimeters(segmentationType, useBinnedImage);
                
        int imageWidth = useBinnedImage ? 
            blobCornerHelper.getGreyscaleImageBinned().getWidth() :
            blobCornerHelper.getGreyscaleImage().getWidth();
        
        int imageHeight = useBinnedImage ? 
            blobCornerHelper.getGreyscaleImageBinned().getHeight() :
            blobCornerHelper.getGreyscaleImage().getHeight();
        
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
            
            if (edge.getN() == 0) {
                HoughTransformLines htl = ht.new HoughTransformLines(
                    new HashMap<PairInt, PairInt>(), new ArrayList<Set<PairInt>>());
                lineList.add(htl);
                continue;
            }

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
    
    private void removeLineArtifactCorners(List<HoughTransformLines> 
        houghTransformLines, BlobPerimeterCornerHelper blobCornerHelper, 
        SegmentationType segmentationType, boolean useBinnedImage) {
        
        List<PairIntArray> perimeterLists = blobCornerHelper.
            getBlobPerimeters(segmentationType, useBinnedImage);
                
        List<List<CornerRegion>> cornerRegionLists = 
            blobCornerHelper.getPerimeterCorners(segmentationType, useBinnedImage);
            
        int imageWidth = useBinnedImage ? 
            blobCornerHelper.getGreyscaleImageBinned().getWidth() :
            blobCornerHelper.getGreyscaleImage().getWidth();
        
        int imageHeight = useBinnedImage ? 
            blobCornerHelper.getGreyscaleImageBinned().getHeight() :
            blobCornerHelper.getGreyscaleImage().getHeight();
        
        int thetaTol = 1;
        int radiusTol = 7;

        //use hough transform for lines to remove corners from line artifacts
        CornerCorrector.removeCornersFromLineArtifacts(houghTransformLines,
            perimeterLists, cornerRegionLists, thetaTol, radiusTol, imageWidth, 
            imageHeight);
    }

    private boolean matchAllPoints(SegmentationType type, boolean useBinned,
        RotatedOffsets rotatedOffsets) {
        
        List<List<CornerRegion>> corners1List = img1Helper.getPerimeterCorners(
            type, useBinned);
        
        List<List<CornerRegion>> corners2List = img2Helper.getPerimeterCorners(
            type, useBinned);
        
        final int blockHalfWidth = 5;
        final boolean useNormalizedIntensities = true;
        
        int binFactor1, binFactor2;
        GreyscaleImage img1, img2;
        if (useBinned) {
            binFactor1 = img1Helper.getBinFactor(useBinned);
            binFactor2 = img2Helper.getBinFactor(useBinned);
            img1 = img1Helper.getGreyscaleImageBinned();
            img2 = img2Helper.getGreyscaleImageBinned();
        } else {
            binFactor1 = 1;
            binFactor2 = 1;
            img1 = img1Helper.getGreyscaleImage();
            img2 = img2Helper.getGreyscaleImage();
        }
        
        IntensityFeatures features1 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities, rotatedOffsets);
        features1.calculateGradientWithGreyscale(img1.copyImage());
        
        IntensityFeatures features2 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities, rotatedOffsets);
        features2.calculateGradientWithGreyscale(img2.copyImage());
        
        List<CornerRegion> corners1 = new ArrayList<CornerRegion>();
        List<CornerRegion> corners2 = new ArrayList<CornerRegion>();
        for (int i = 0; i < corners1List.size(); ++i) {
            corners1.addAll(corners1List.get(i));
        }
        for (int i = 0; i < corners2List.size(); ++i) {
            corners2.addAll(corners2List.get(i));
        }
        
        List<FeatureComparisonStat> stats = match(img1, img2, features1,
            features2, corners1, corners2, binFactor1, binFactor2);
        
        if (stats.isEmpty()) {
            return false;
        }
        
//TODO: here, need to remove points inconsistent with homology
        
        if (useBinned) {
            stats = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), stats, 
                binFactor1, binFactor2, rotatedOffsets);
        }
        
        this.solutionStats = stats;
        
        solutionMatched1 = new ArrayList<PairInt>();
        solutionMatched2 = new ArrayList<PairInt>();
        for (FeatureComparisonStat stat : stats) {
            solutionMatched1.add(stat.getImg1Point().copy());
            solutionMatched2.add(stat.getImg2Point().copy());
        }
        
        return true;
    }
    
    private  List<FeatureComparisonStat> match(GreyscaleImage img1, GreyscaleImage img2,
        IntensityFeatures features1, IntensityFeatures features2,
        List<CornerRegion> corners1, List<CornerRegion> corners2,
        int binFactor1, int binFactor2) {
        
        int dither = 1;
        
        CornerMatcher<CornerRegion> matcher = new CornerMatcher<CornerRegion>(dither);
        
        boolean matched = matcher.matchCorners(features1, features2, 
            corners1, corners2, img1, img2,
            binFactor1, binFactor2);
                
        List<FeatureComparisonStat> stats = matcher.getSolutionStats();
        
        return stats;
    }

    private List<FeatureComparisonStat> reviseStatsForFullImages(
        GreyscaleImage gsImg1, GreyscaleImage gsImg2,
        List<FeatureComparisonStat> stats, int prevBinFactor1, 
        int prevBinFactor2, RotatedOffsets rotatedOffsets) {

        log.info("refine stats for full image reference frames");

        List<FeatureComparisonStat> revised = new ArrayList<FeatureComparisonStat>();

        FeatureMatcher featureMatcher = new FeatureMatcher();

        IntensityFeatures features1 = new IntensityFeatures(5,
            settings.useNormalizedFeatures(), rotatedOffsets);
        features1.calculateGradientWithGreyscale(gsImg1);

        IntensityFeatures features2 = new IntensityFeatures(5,
            settings.useNormalizedFeatures(), rotatedOffsets);
        features2.calculateGradientWithGreyscale(gsImg2);

        int dither = 2;

        for (int i = 0; i < stats.size(); ++i) {

            FeatureComparisonStat stat = stats.get(i);

            int x1 = stat.getImg1Point().getX() * prevBinFactor1;
            int y1 = stat.getImg1Point().getY() * prevBinFactor1;
            int x2 = stat.getImg2Point().getX() * prevBinFactor2;
            int y2 = stat.getImg2Point().getY() * prevBinFactor2;

            // have to discard the best angles found in stat and derive new
            // for these higher resolution images
            FeatureComparisonStat compStat =
                featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, x1, y1, x2, y2, dither, 
                    gsImg1, gsImg2);

            if (compStat == null ||
                (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())) {
                continue;
            }

            revised.add(compStat);
        }
        
        return revised;
    }

    public List<FeatureComparisonStat> getSolutionStats() {
        return solutionStats;
    }

    private boolean matchCurveByCurve(SegmentationType type, boolean useBinned,
        RotatedOffsets rotatedOffsets) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    /**
     * @return the solutionMatched1
     */
    public List<PairInt> getSolutionMatched1() {
        return solutionMatched1;
    }

    /**
     * @return the solutionMatched2
     */
    public List<PairInt> getSolutionMatched2() {
        return solutionMatched2;
    }
    
}
