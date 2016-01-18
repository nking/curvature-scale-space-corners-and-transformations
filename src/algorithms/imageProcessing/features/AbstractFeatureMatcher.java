package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.HoughTransform.HoughTransformLines;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceCornerDetector;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.util.MiscStats;
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
public abstract class AbstractFeatureMatcher {
    
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
    
    protected final FeatureMatcherSettings settings;
    
    protected boolean useSameSegmentation = false;
    
    protected List<FeatureComparisonStat> solutionStats = null;
    
    protected List<PairInt> solutionMatched1 = null;
    
    protected List<PairInt> solutionMatched2 = null;
    
    public AbstractFeatureMatcher(ImageExt img1, ImageExt img2, 
        FeatureMatcherSettings settings) {
        
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        
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

    public boolean match() throws IOException, NoSuchAlgorithmException {
        
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
        
        //TODO: document this logic for the user
        SegmentationType type = SegmentationType.GREYSCALE_WAVELET;
        if (settings.doOverrideWithCannySegmentation()) {
            type = SegmentationType.GREYSCALE_CANNY;
        } else if (settings.doUse2ndDerivCorners()) {
            type = SegmentationType.NONE;
        }
        
        boolean[] useBinned = settings.startWithBinnedImages() ? 
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
                return true;
            }
        }
        
        return false;
    }

    protected void prepareCorners(SegmentationType type, boolean useBinned) 
        throws IOException, NoSuchAlgorithmException {
        
        log.info(type.name() + " binned=" + useBinned + " useSameSegmentation=" 
            + useSameSegmentation);
        
        IntensityFeatures f1, f2;
        
        if (useBinned) {
            img1Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            img2Helper.createBinnedGreyscaleImage(binnedImageMaxDimension);
            f1 = featuresBinned1;
            f2 = featuresBinned2;
        } else {
            f1 = features1;
            f2 = features2;
        }
        
        img1Helper.applySegmentation(type, useBinned);
        img2Helper.applySegmentation(type, useBinned);
        
        List<HoughTransformLines> houghTransformLines1;
        List<HoughTransformLines> houghTransformLines2;
                
        if (type.equals(SegmentationType.GREYSCALE_CANNY)) {
            
            boolean filterOutImageBoundaryBlobs = true;
            boolean filterOutZeroPixels = false;
            boolean doNotAddPoints = true;
            
            // pre-make the blobs using non-default variables:
            img1Helper.getBlobs(type, useBinned, filterOutImageBoundaryBlobs, 
                filterOutZeroPixels);
            
            img2Helper.getBlobs(type, useBinned, filterOutImageBoundaryBlobs, 
                filterOutZeroPixels); 
            
            if (settings.doUse2ndDerivCorners()) {
                //TODO: consider using extractSecondDerivativeCornersithoutBlobs instead
                img1Helper.extractSecondDerivativeCorners(type, useBinned);
                img2Helper.extractSecondDerivativeCorners(type, useBinned);
                
                Logger.getLogger(this.getClass().getName()).info(
                    "nPts 2nd deriv img1 = " + img1Helper.getAllCorners(
                    SegmentationType.NONE, useBinned).size());
                
                Logger.getLogger(this.getClass().getName()).info(
                    "nPts 2nd deriv img2 = " + img2Helper.getAllCorners(
                    SegmentationType.NONE, useBinned).size());
                
            } else {
                img1Helper.extractBlobPerimeterAsCornerRegions(type, useBinned, 
                    doNotAddPoints);
                img2Helper.extractBlobPerimeterAsCornerRegions(type, useBinned, 
                    doNotAddPoints);
            }
        }
        
        if (settings.doUse2ndDerivCorners()) {
            img1Helper.extractSecondDerivativeCornersWithoutBlobs(type, useBinned);
            img2Helper.extractSecondDerivativeCornersWithoutBlobs(type, useBinned);
            
            Logger.getLogger(this.getClass().getName()).info(
                "nPts 2nd deriv img1 = " + img1Helper.getAllCorners(
                    SegmentationType.NONE, useBinned).size());

            Logger.getLogger(this.getClass().getName()).info(
                "nPts 2nd deriv img2 = " + img2Helper.getAllCorners(
                    SegmentationType.NONE, useBinned).size());
        } else {
            img1Helper.extractBlobPerimeterAsCornerRegions(type, useBinned);
            img2Helper.extractBlobPerimeterAsCornerRegions(type, useBinned);
        
            houghTransformLines1 = findLinesUsingHoughTransform(img1Helper, 
                type, useBinned);
            houghTransformLines2 = findLinesUsingHoughTransform(img2Helper, 
                type, useBinned);
            removeLineArtifactCorners(houghTransformLines1, img1Helper, type, 
                useBinned);
            removeLineArtifactCorners(houghTransformLines2, img2Helper, type, 
                useBinned);
        }
    }

    protected boolean generateAndMatchCornerRegions(SegmentationType type, 
        boolean useBinned) throws IOException, NoSuchAlgorithmException {
        
        log.info(type.name() + " binned=" + useBinned + " useSameSegmentation=" 
            + useSameSegmentation);
        
        prepareCorners(type, useBinned);
                        
        return match(type, useBinned);
    }
    
    protected abstract boolean match(SegmentationType type, boolean useBinned);

    protected List<HoughTransform.HoughTransformLines> 
        findLinesUsingHoughTransform(BlobPerimeterCornerHelper blobCornerHelper, 
        SegmentationType segmentationType, boolean useBinnedImage) {
        
        List<PairIntArray> perimeterLists = blobCornerHelper.getBlobPerimeters(
            segmentationType, useBinnedImage);
        
        int imageWidth = useBinnedImage ? 
            blobCornerHelper.getGreyscaleImageBinned().getWidth() : 
            blobCornerHelper.getGreyscaleImage().getWidth();
        
        int imageHeight = useBinnedImage ? 
            blobCornerHelper.getGreyscaleImageBinned().getHeight() : 
            blobCornerHelper.getGreyscaleImage().getHeight();
        
        int thetaTol = 1;
        int radiusTol = 7;
        
        HoughTransform ht = new HoughTransform();
        
        List<HoughTransform.HoughTransformLines> lineList = 
            new ArrayList<HoughTransform.HoughTransformLines>();
        
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
            HoughTransform.HoughTransformLines htl = 
                ht.createPixTRMapsFromSorted(outSortedKeys, 
                outputPolarCoordsPixMap, thetaTol, radiusTol);
            
            lineList.add(htl);
        }
        return lineList;
    }

    protected void removeLineArtifactCorners(List<HoughTransformLines> 
        houghTransformLines, BlobPerimeterCornerHelper blobCornerHelper, 
        SegmentationType segmentationType, boolean useBinnedImage) {
        
        List<PairIntArray> perimeterLists = blobCornerHelper.getBlobPerimeters(
            segmentationType, useBinnedImage);
        
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

    protected List<FeatureComparisonStat> reviseStatsForFullImages(
        GreyscaleImage gsImg1, GreyscaleImage gsImg2, 
        List<FeatureComparisonStat> stats, int prevBinFactor1, 
        int prevBinFactor2, RotatedOffsets rotatedOffsets) {
        
        log.info("refine stats for full image reference frames");
        
        if (stats.isEmpty()) {
            return stats;
        }
        
        List<FeatureComparisonStat> revised = new ArrayList<FeatureComparisonStat>();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        IntensityFeatures f1 = new IntensityFeatures(5, 
            settings.useNormalizedFeatures(), rotatedOffsets);
        f1.calculateGradientWithGreyscale(gsImg1);
        
        IntensityFeatures f2 = new IntensityFeatures(5, 
            settings.useNormalizedFeatures(), rotatedOffsets);
        f2.calculateGradientWithGreyscale(gsImg2);
        
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
                featureMatcher.ditherAndRotateForBestLocation2(f1, 
                    f2, x1, y1, x2, y2, dither, gsImg1, gsImg2);
            
            if (compStat == null || (compStat.getSumIntensitySqDiff() > 
                compStat.getImg2PointIntensityErr())) {
                continue;
            }
            
            revised.add(compStat);
        }
        
        return revised;
    }

    public void copyToInstanceVars(List<FeatureComparisonStat> stats) {
        
        this.solutionStats = new ArrayList<FeatureComparisonStat>();
        this.solutionMatched1 = new ArrayList<PairInt>();
        this.solutionMatched2 = new ArrayList<PairInt>();
        
        for (FeatureComparisonStat stat : stats) {
            solutionStats.add(stat.copy());
            solutionMatched1.add(stat.getImg1Point().copy());
            solutionMatched2.add(stat.getImg2Point().copy());
        }
    }
    
    protected void extractCannyCornerRegions(GreyscaleImage img1,
        GreyscaleImage img2, Set<CornerRegion> outputCornerRegions1,
        Set<CornerRegion> outputCornerRegions2) {

        GreyscaleImage gsImg1 = img1.copyImage();
        GreyscaleImage gsImg2 = img2.copyImage();

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(gsImg1, SIGMA.ONE);

        CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(gsImg1);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();

        Set<CornerRegion> cornerRegions1 = detector.getEdgeCornerRegions(true);
        //cornerRegions1 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);

        outputCornerRegions1.addAll(cornerRegions1);

        imageProcessor.blur(gsImg2, SIGMA.ONE);

        detector = new CurvatureScaleSpaceCornerDetector(gsImg2);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
        Set<CornerRegion> cornerRegions2 = detector.getEdgeCornerRegions(true);
        //cornerRegions2 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);

        outputCornerRegions2.addAll(cornerRegions2);
    }
    
    protected List<FeatureComparisonStat> findCorrespondence(GreyscaleImage img1,
        GreyscaleImage img2, Set<CornerRegion> cornerRegions1,
        Set<CornerRegion> cornerRegions2, TransformationParameters parameters,
        RotatedOffsets rotatedOffsets, int dither2, int transXYTol,
        float scaleTol, float rotationInRadiansTol) {

        FeatureMatcher featureMatcher = new FeatureMatcher();

        List<FeatureComparisonStat> stats = 
            featureMatcher.findSimilarFeaturesAsStats(img1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            img2, cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]),
            parameters, scaleTol, rotationInRadiansTol, transXYTol,
            dither2, rotatedOffsets);
        
        return stats;
    }

    protected MatchingSolution transformSolutionToFullFrames(MatchingSolution 
        soln, BlobPerimeterCornerHelper img1Helper, 
        BlobPerimeterCornerHelper img2Helper, int binFactor1, int binFactor2) {
        
        if (binFactor1 == 1 && binFactor2 == 1) {
            return soln;
        }
       
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        
        assert(rotatedOffsets.containsData());
        
        List<FeatureComparisonStat> stats = soln.getComparisonStats();
        
        for (int i = 0; i < stats.size(); ++i) {
            FeatureComparisonStat stat = stats.get(i);
            stat.setBinFactor1(binFactor1);
            stat.setBinFactor2(binFactor2);
        }
        
        if (settings.debug()) {
            GreyscaleImage im1 = (binFactor1 != 1) ? 
                img1Helper.getGreyscaleImageBinned() :
                img1Helper.getGreyscaleImage();
            GreyscaleImage im2 = (binFactor2 != 1) ? 
                img2Helper.getGreyscaleImageBinned() :
                img2Helper.getGreyscaleImage();
            MiscDebug.writeImages(im1, im2, stats, 
                "_matched_binned_" + settings.getDebugTag() 
                + MiscDebug.getCurrentTimeFormatted(), 1);
        }
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        List<FeatureComparisonStat> fullStats = matcher.reviseStatsForFullImages(
            img1Helper.getGreyscaleImage(),
            img2Helper.getGreyscaleImage(),
            settings,
            soln.getParams(), soln.getComparisonStats(),
            binFactor1, binFactor2, rotatedOffsets);
            
        if ((fullStats == null) || fullStats.isEmpty()) {
            return null;
        }
        
        TransformationParameters revisedParams = 
            MiscStats.calculateTransformation(1, 1, fullStats,
                new float[4], false);
        
        if (revisedParams == null) {
            return null;
        }
        
        if (settings.debug()) {
            GreyscaleImage im1 = img1Helper.getGreyscaleImage();
            GreyscaleImage im2 = img2Helper.getGreyscaleImage();
            MiscDebug.writeImages(im1, im2, fullStats, 
                "_matched_" + settings.getDebugTag() +
                MiscDebug.getCurrentTimeFormatted(), 2);
        }
        
        MatchingSolution fullSoln = new MatchingSolution(revisedParams, 
            fullStats, 1, 1);
        
        return fullSoln;
    }

    protected void filterForLocalization(GreyscaleImage img, 
        IntensityFeatures f, List<CornerRegion> corners) {
                
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < corners.size(); ++i) {
            CornerRegion cr = corners.get(i);
            
            try {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                if (f.removeDueToLocalization(img, x, y,
                    f.calculateOrientation(x, y))) {
                    remove.add(Integer.valueOf(i));
                }
            } catch (CornerRegion.CornerRegionDegneracyException ex) {
            }
        }
        
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i);
            corners.remove(idx);
        }
    }
    
    protected void filterForLocalization2(GreyscaleImage img, 
        IntensityFeatures f, List<List<CornerRegion>> corners) {
                        
        for (int i = 0; i < corners.size(); ++i) {
            List<CornerRegion> cr = corners.get(i);
            filterForLocalization(img, f, cr);
        }
    }
    
    /**
     * get a copy of the solution's feature stats.
     * @return 
     */
    public List<FeatureComparisonStat> getSolutionStats() {
        return solutionStats;
    }

    /**
     * get the matched points from image 1
     * @return the solutionMatched1
     */
    public List<PairInt> getSolutionMatched1() {
        return solutionMatched1;
    }

    /**
     * get the matched points from image 2
     * @return the solutionMatched2
     */
    public List<PairInt> getSolutionMatched2() {
        return solutionMatched2;
    }
    
}
