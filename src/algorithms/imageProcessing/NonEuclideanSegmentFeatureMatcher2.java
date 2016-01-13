package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * create lists of possibly degenerately matched points between 2 images.
 * It uses the cosine similarity of the descriptors along with a threshold
 * of 0.95 to keep matches above the threshold even if there are more than
 * one.
 * 
 * @author nichole
 */
public class NonEuclideanSegmentFeatureMatcher2 extends AbstractFeatureMatcher {
    
    public NonEuclideanSegmentFeatureMatcher2(ImageExt img1, ImageExt img2, 
        FeatureMatcherSettings settings) {

        super(img1, img2, settings);
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
        
        List<CornerRegion> corners1 = new ArrayList<CornerRegion>();
        List<CornerRegion> corners2 = new ArrayList<CornerRegion>();
        for (int i = 0; i < corners1List.size(); ++i) {
            corners1.addAll(corners1List.get(i));
        }
        for (int i = 0; i < corners2List.size(); ++i) {
            corners2.addAll(corners2List.get(i));
        }
         
        int dither = 1;
        
        List<FeatureComparisonStat> stats;
        
        CosineSimilarityCornerMatcher<CornerRegion> matcher = 
            new CosineSimilarityCornerMatcher<CornerRegion>(dither);
        
        boolean matched = matcher.matchCorners(f1, f2, 
            corners1, corners2, img1, img2);

        if (!matched) {
            return false;
        }
        
        Map<PairInt, List<PairInt>> matched1Matched2 = matcher.getMatched1Matched2();
        Map<PairInt, List<FeatureComparisonStat>> solutionStats0
            = matcher.getSolutionStats();
        
        throw new UnsupportedOperationException("not yet implemented");
        
        /*
    use a ransac solver that chooses from a degenerate matched list
        
        if (stats.isEmpty()) {
            return false;
        }
                
        if (useBinned) {
            stats = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), stats, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
        }
        
        copyToInstanceVars(stats);
                
        return true;*/
    }
    
}
