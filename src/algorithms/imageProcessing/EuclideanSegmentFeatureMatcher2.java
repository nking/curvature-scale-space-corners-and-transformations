package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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
 
 * @author nichole
 */
public class EuclideanSegmentFeatureMatcher2 extends AbstractFeatureMatcher {

    protected TransformationParameters solutionTransformation = null;
    
    public EuclideanSegmentFeatureMatcher2(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {
        
          super(image1, image2, settings);
    }

    @Override
    protected boolean match(SegmentationType type, boolean useBinned) {
        
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
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>(); 
        //List<PairInt> matched1 = new ArrayList<PairInt>();
        //List<PairInt> matched2 = new ArrayList<PairInt>();
        for (FeatureComparisonStat stat :  matcher.getSolutionStats()) {
            stats.add(stat.copy());
            //matched1.add(stat.getImg1Point().copy());
            //matched2.add(stat.getImg2Point().copy());
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
        
        throw new UnsupportedOperationException("not yet implemented");
        /*
        Map<PairInt, List<FeatureComparisonStat>> index1Index2Matches 
            = associateMatchesWithBlobs(stats, 
            img1Helper.getBlobs(type, useBinned),
            img2Helper.getBlobs(type, useBinned));
                
        BlobCornersEuclideanCalculator2 calculator = 
            new BlobCornersEuclideanCalculator2();
        
        MatchingSolution soln = calculator.solveTransformation(
            img1Helper.getGreyscaleImage(useBinned),
            img2Helper.getGreyscaleImage(useBinned),
            f1, f2, dither, index1Index2Matches,
            corners1, corners2);
        
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
        
        if (useBinned) {
            stats = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), stats, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
        }
        
        copyToInstanceVars(stats);
        
        solutionTransformation = bestParams.copy();
                
        return true;*/
    }
    
    public TransformationParameters getSolutionTransformationParameters() {
        return solutionTransformation;
    }
    
}
