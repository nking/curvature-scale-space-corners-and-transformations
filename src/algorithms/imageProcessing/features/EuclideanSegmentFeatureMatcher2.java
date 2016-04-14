package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.SegmentationType;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 class whose goal is to find best single euclidean transformation for 
 image1 to image2.  It uses segmentation to create blobs.
 The "points of interest" are either made from the blob perimeters or
 from 2nd derivatives of gaussian convolution with greyscale image.
 It matches the points from image1 against points from image 2
 using features and those with a close 2nd best SSD match are discarded.
 The points are then associated with the enclosing blobs and
 then combinations of the blobs are made to create euclidean transformations.
 The transformations are evaluated against all points.
 (Note if that is successful, will try to evaluate against only the
 matched points which are fewer in number and will change class to use that
 logic if tests pass).
 
 The checkerboard tests do not match well with this one because the points are
 not unique enough, but EuclideanSegmentFeatureMatcher does solve it.
 SO, when this method returns very few or no matched points, the invoker
 should then follow with EuclideanSegmentFeatureMatcher.
 
 * @author nichole
 */
public class EuclideanSegmentFeatureMatcher2 extends AbstractFeatureMatcher {

    protected TransformationParameters solutionTransformation = null;
    
    public EuclideanSegmentFeatureMatcher2(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {
        
          super(image1, image2, settings);
    }
    
    @Override
    protected void prepareCorners(SegmentationType type, boolean useBinned) 
        throws IOException, NoSuchAlgorithmException {
        
        // overriding to add logic for settings for 2nd deriv points
        
        super.prepareCorners(type, useBinned);
        
        if (settings.doUse2ndDerivCorners() && type.equals(SegmentationType.NONE)) {
            
            // create blobs, and associate 2nd deriv pts with the blobs.
            // Note that this is repeating some work.
            // these feature matchers will be refactored soon.
            
            SegmentationType type2 = SegmentationType.GREYSCALE_WAVELET;

            img1Helper.applySegmentation(type2, useBinned);
            img2Helper.applySegmentation(type2, useBinned);

            img1Helper.extractSecondDerivativeCorners(type2, useBinned);
            img2Helper.extractSecondDerivativeCorners(type2, useBinned);
        }
        
    }

    @Override
    protected boolean match(SegmentationType type, boolean useBinned) {
        
        if (settings.doUse2ndDerivCorners() && type.equals(SegmentationType.NONE)) {
            // for 2nd pt derivs, need to override the type to wavelet.
            // this should be handled better after refacoring.
            type = SegmentationType.GREYSCALE_WAVELET;
        }
        
        int binFactor1, binFactor2;
        GreyscaleImage img1, img2;
        IntensityFeatures f1, f2;
        
        if (useBinned) {
            binFactor1 = img1Helper.getBinFactor(useBinned);
            binFactor2 = img2Helper.getBinFactor(useBinned);
            img1 = img1Helper.getGreyscaleImageBinned();
            img2 = img2Helper.getGreyscaleImageBinned();
            f1 = featuresBinned1;
            f2 = featuresBinned2;
        } else {
            binFactor1 = 1;
            binFactor2 = 1;
            img1 = img1Helper.getGreyscaleImage();
            img2 = img2Helper.getGreyscaleImage();
            f1 = features1;
            f2 = features2;
            if (!f1.gradientWasCreated()) {
                f1.calculateGradientWithGreyscale(
                    img1Helper.getGreyscaleImage().copyImage());
            }
            if (!f2.gradientWasCreated()) {
                f2.calculateGradientWithGreyscale(
                    img2Helper.getGreyscaleImage().copyImage());
            }
        }
        
        List<List<CornerRegion>> corners1List = img1Helper.getPerimeterCorners(
            type, useBinned);
        
        List<List<CornerRegion>> corners2List = img2Helper.getPerimeterCorners(
            type, useBinned);
                 
        int dither = 1;        
        
        List<CornerRegion> corners1 = new ArrayList<CornerRegion>();
        List<CornerRegion> corners2 = new ArrayList<CornerRegion>();
        for (int i = 0; i < corners1List.size(); ++i) {
            corners1.addAll(corners1List.get(i));
        }
        for (int i = 0; i < corners2List.size(); ++i) {
            corners2.addAll(corners2List.get(i));
        }

        CornerMatcher<CornerRegion> matcher = new CornerMatcher<CornerRegion>(dither);

        boolean matched = matcher.matchCorners(f1, f2, corners1, corners2, 
            img1, img2, binFactor1, binFactor2);

        if (!matched) {
            return false;
        }
        
        /*MiscDebug.writeImagesInAlternatingColor(img1.copyToColorGreyscaleExt(), 
            img2.copyToColorGreyscaleExt(), matcher.getSolutionStats(), 
            "_matched_non_euclid_" + MiscDebug.getCurrentTimeFormatted(), 2);
        */
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>(); 
        for (FeatureComparisonStat stat :  matcher.getSolutionStats()) {
            stats.add(stat.copy());
        }
        
        matcher = null;
        
        if (stats.isEmpty()) {
            return false;
        }
        
        /*
        matched1 and matched2 are both associated with blobs.
        
        key = index1, index2
        value = all matches for those 2 blobs
        
        Map<PairInt, List<FeatureComparisonStat>> 
        */
        
        Map<PairInt, List<FeatureComparisonStat>> index1Index2Matches 
            = associateMatchesWithBlobs(stats, 
            img1Helper.getBlobs(type, useBinned),
            img2Helper.getBlobs(type, useBinned));
        
        filterForLocalization(img1Helper.getGreyscaleImage(useBinned),
            img2Helper.getGreyscaleImage(useBinned), f1, f2,
            index1Index2Matches);
        /*
        if (true) {
            GreyscaleImage im1 = img1Helper.getGreyscaleImage(useBinned);
            GreyscaleImage im2 = img2Helper.getGreyscaleImage(useBinned);
            List<FeatureComparisonStat> blobAssoc = new ArrayList<FeatureComparisonStat>();
            for (java.util.Map.Entry<PairInt, List<FeatureComparisonStat>> entry : index1Index2Matches.entrySet()) {
                blobAssoc.addAll(entry.getValue());
            }
            MiscDebug.writeImagesInAlternatingColor(img1.copyToColorGreyscaleExt(), 
                img2.copyToColorGreyscaleExt(), blobAssoc, 
                "_matched_assoc_blobs_" + MiscDebug.getCurrentTimeFormatted(), 2);
            log.info(index1Index2Matches.size() + " blob pairs to match from");
            log.info(blobAssoc.size() + " matched points associated w/ blobs");
        }*/
                
        BlobCornersEuclideanCalculator2 calculator = 
            new BlobCornersEuclideanCalculator2();
        
        MatchingSolution soln = calculator.solveTransformation(
            img1Helper.getGreyscaleImage(useBinned),
            img2Helper.getGreyscaleImage(useBinned),
            f1, f2, dither, index1Index2Matches,
            corners1, corners2, binFactor1, binFactor2);
        
        if (soln == null) {
            return false;
        }
        
        log.info(soln.getParams().toString());
            
        // transform images to full size
        soln = transformSolutionToFullFrames(soln, img1Helper, img2Helper, 
            binFactor1, binFactor2);

        if (soln == null) {
            return false;
        }
            
        log.info("full frame soln: " + soln.getParams().toString());
        
        copyToInstanceVars(soln.getComparisonStats());
        
        solutionTransformation = soln.getParams().copy();
                
        return true;
    }

    public TransformationParameters getSolutionTransformationParameters() {
        return solutionTransformation;
    }

    private Map<PairInt, List<FeatureComparisonStat>> associateMatchesWithBlobs(
        List<FeatureComparisonStat> stats, List<Set<PairInt>> blobs1, 
        List<Set<PairInt>> blobs2) {
        
        /*
        to make the lookups for stats points O(1), could make a large map for each
        blobs collection w/ key = coord and value= list index.  iterating over
        all blob points is nBlobs * avgBlobSize * 2 where 2 is once for image1,
        and then image2.
        
        if there are 100 stats and blobs1 size = blobs2 = 50,
        then finding indexes for all stats is at most 100*50 + 100*50 = 10,000
        
        if make large map, and if each blob has n points,
          50*n*2 + 100*2  is number of lookups, so it's only a better runtime if
        avg size of blobs is < 100 for this example.
        */
        
        Map<PairInt, List<FeatureComparisonStat>> matchedBlobs = 
            new HashMap<PairInt, List<FeatureComparisonStat>>();

        for (int i = 0; i < stats.size(); ++i) {
            
            FeatureComparisonStat stat = stats.get(i);
            
            PairInt p1 = stat.getImg1Point();
            PairInt p2 = stat.getImg2Point();
            
            PairInt indexes = new PairInt(0, 0);
            
            boolean found = false;
            for (int j = 0; j < blobs1.size(); ++j) {
                Set<PairInt> blob = blobs1.get(j);
                if (blob.contains(p1)) {
                    indexes.setX(j);
                    found = true;
                    break;
                }
            }
            if (!found) {
                continue;
            }
            found = false;
            for (int j = 0; j < blobs2.size(); ++j) {
                Set<PairInt> blob = blobs2.get(j);
                if (blob.contains(p2)) {
                    indexes.setY(j);
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                continue;
            }
            
            List<FeatureComparisonStat> list = matchedBlobs.get(indexes);
            if (list == null) {
                list = new ArrayList<FeatureComparisonStat>();
                matchedBlobs.put(indexes, list);
            }
            list.add(stat);
        }
        
        return matchedBlobs;
    }

    private void filterForLocalization(GreyscaleImage img1, 
        GreyscaleImage img2, IntensityFeatures features1, 
        IntensityFeatures features2, Map<PairInt, 
            List<FeatureComparisonStat>> index1Index2Matches) {
        
        Set<PairInt> remove = new HashSet<PairInt>();
        
        for (Entry<PairInt, List<FeatureComparisonStat>> entry : index1Index2Matches.entrySet()) {
            List<FeatureComparisonStat> stats = entry.getValue();
            for (int i = (stats.size() - 1); i > -1; --i) {
                FeatureComparisonStat stat = stats.get(i);
                PairInt p1 = stat.getImg1Point();
                if (features1.removeDueToLocalization(img1, p1.getX(), p1.getY(),
                    features1.calculateOrientation(p1.getX(), p1.getY()))) {
                    stats.remove(i);
                    continue;
                }
                PairInt p2 = stat.getImg2Point();
                if (features2.removeDueToLocalization(img2, p2.getX(), p2.getY(),
                    features2.calculateOrientation(p2.getX(), p2.getY()))) {
                    stats.remove(i);
                    continue;
                }
            }
            if (stats.isEmpty()) {
                remove.add(entry.getKey());
            }
        }
        
        for (PairInt key : remove) {
            index1Index2Matches.remove(key);
        }
    }
    
}
