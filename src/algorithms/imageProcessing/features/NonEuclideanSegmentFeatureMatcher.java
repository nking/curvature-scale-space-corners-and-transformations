package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.SegmentationType;
import java.util.ArrayList;
import java.util.List;

/**
 * create lists of singly matched points between 2 images.
 * It uses the criteria that matches are discarded if a point has a second
 * best match whose SSD is within 0.8*SSD of best match.
 * 
 * @author nichole
 */
public class NonEuclideanSegmentFeatureMatcher extends AbstractFeatureMatcher {
    
    protected List<FeatureComparisonStat> rejectedBy2ndBest = new ArrayList<FeatureComparisonStat>();
        
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

        if (settings.doUse2ndDerivCorners() || settings.doOverrideWithCannySegmentation()) {
            return matchFor2ndDerivPts(type, useBinned);
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
        
        boolean filterForLocalization = true;
        if (filterForLocalization) {
            filterForLocalization2(img1, f1, corners1List);
            filterForLocalization2(img2, f2, corners2List);
        }
        
        boolean matchCurveToCurve = true;
         
        int dither = 1;
        
        List<FeatureComparisonStat> rb2j = new ArrayList<FeatureComparisonStat>();
            
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
            
            rb2j.addAll(matcher.getRejectedBy2ndBest());
        }
        
        if (stats.isEmpty()) {
            return false;
        }
                
        if (useBinned) {
            stats = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), stats, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
            
            rb2j = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), rb2j, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
        }
        
        copyToInstanceVars(stats, rb2j);
                
        return true;
    }
    
    protected boolean matchFor2ndDerivPts(SegmentationType type, boolean useBinned) {
        
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
                
        List<CornerRegion> corners1 = new ArrayList<CornerRegion>(
            img1Helper.getAllCorners(type, useBinned));
        
        List<CornerRegion> corners2 = new ArrayList<CornerRegion>(
            img2Helper.getAllCorners(type, useBinned));
        
        boolean filterForLocalization = true;
        if (filterForLocalization) {
            filterForLocalization(img1, f1, corners1);
            filterForLocalization(img2, f2, corners2);
        }
        
        log.info("nPts after localization filter img1 = " + corners1.size());

        log.info("nPts after localization filter img2 = " + corners2.size());
               
        int dither = 1;
                
        CornerMatcher<CornerRegion> matcher = new CornerMatcher<CornerRegion>(dither);

        boolean matched = matcher.matchCorners(f1, f2, corners1, corners2, 
            img1, img2, binFactor1, binFactor2);

        if (!matched) {
            return false;
        }
                
        List<FeatureComparisonStat> stats = matcher.getSolutionStats();
        
        log.info("nPts after SSD match (incl filter for 2nd best) = " + stats.size());
        
        log.info("nPts in 2nd best rejection list = " + matcher.getRejectedBy2ndBest().size());
        
        if (stats.isEmpty()) {
            return false;
        }
        
        List<FeatureComparisonStat> rb2j = 
            new ArrayList<FeatureComparisonStat>(matcher.getRejectedBy2ndBest());
                
        if (useBinned) {
            stats = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), stats, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
            
            rb2j = reviseStatsForFullImages(img1Helper.getGreyscaleImage(), 
                img2Helper.getGreyscaleImage(), rb2j, 
                binFactor1, binFactor2, f1.getRotatedOffsets());
        }
        
        copyToInstanceVars(stats, rb2j);
                
        return true;
    }

    private void copyToInstanceVars(List<FeatureComparisonStat> stats, 
        List<FeatureComparisonStat> rejBy2ndBest) {
        
        super.copyToInstanceVars(stats);
        
        rejectedBy2ndBest.clear();
        
        rejectedBy2ndBest.addAll(rejBy2ndBest);        
    }
    
    public List<FeatureComparisonStat> getRejectedBy2ndBest() {
        return rejectedBy2ndBest;
    }
    
}
