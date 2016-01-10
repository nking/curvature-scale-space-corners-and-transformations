package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;

/**
 * create lists of matched points present in 2 images.
 * 
 * @author nichole
 */
public class NonEuclideanSegmentFeatureMatcher extends AbstractFeatureMatcher {
    
    /**
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     * @param settings
     */
    public NonEuclideanSegmentFeatureMatcher(ImageExt img1, ImageExt img2, 
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
        
        boolean matchCurveToCurve = true;
         
        int dither = 1;
        
        List<FeatureComparisonStat> stats;
        
        if (matchCurveToCurve) {
        
            CurveToCurveCornerMatcher<CornerRegion> matcher = 
                new CurveToCurveCornerMatcher<CornerRegion>(dither);
        
            boolean matched = matcher.matchCorners(f1, f2, 
                corners1List, corners2List, img1, img2, binFactor1, binFactor2);

            if (!matched) {
                return false;
            }
            
            stats = matcher.getSolutionStats();

        } else {
            
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
                
            stats = matcher.getSolutionStats();
        }
        
        if (stats.isEmpty()) {
            return false;
        }
                
        if (useBinned) {
            stats = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), stats, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
        }
        
        copyToInstanceVars(stats);
                
        return true;
    }
    
}
