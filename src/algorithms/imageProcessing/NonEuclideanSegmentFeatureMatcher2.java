package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * create lists of possibly degenerately matched points between 2 images.
 * It uses the cosine similarity of the descriptors along with a threshold
 * of 0.95 to keep matches above the threshold even if there are more than
 * one.
 * 
 * @author nichole
 */
public class NonEuclideanSegmentFeatureMatcher2 extends AbstractFeatureMatcher {
    
    protected Map<PairInt, List<FeatureComparisonStat>> solutionStatsMultiplicity = null;
        
    protected List<List<PairInt>> solutionMatched2Multiplicity = null;
    
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
        List<CornerRegion> corners1 = new ArrayList<CornerRegion>();
        List<CornerRegion> corners2 = new ArrayList<CornerRegion>();
        
        if (settings.doUse2ndDerivCorners()) {
         
            corners1.addAll(img1Helper.getAllCorners(type, useBinned));
            corners2.addAll(img2Helper.getAllCorners(type, useBinned));
            
        } else {
            
            List<List<CornerRegion>> corners1List = img1Helper.getPerimeterCorners(
                type, useBinned);

            List<List<CornerRegion>> corners2List = img2Helper.getPerimeterCorners(
                type, useBinned);


            for (int i = 0; i < corners1List.size(); ++i) {
                corners1.addAll(corners1List.get(i));
            }
            for (int i = 0; i < corners2List.size(); ++i) {
                corners2.addAll(corners2List.get(i));
            }
        }
        
        boolean filterForLocalization = true;
        if (filterForLocalization) {
            filterForLocalization(img1, f1, corners1);
            filterForLocalization(img2, f2, corners2);
        }
        
        int dither = 1;
                
        CosineSimilarityCornerMatcher<CornerRegion> matcher = 
            new CosineSimilarityCornerMatcher<CornerRegion>(dither);
        
        boolean matched = matcher.matchCorners(f1, f2, corners1, corners2, 
            img1, img2);

        if (!matched) {
            return false;
        }
        
        Map<PairInt, List<FeatureComparisonStat>> solutionStats0
            = matcher.getSolutionStats();

        if (solutionStats0 == null || solutionStats0.isEmpty()) {
            return false;
        }
        
        if (useBinned) {
            solutionStats0 = reviseStatsForFullImages(
                img1Helper.getGreyscaleImage(), img2Helper.getGreyscaleImage(), 
                solutionStats0, binFactor1, binFactor2, f1.getRotatedOffsets());
        }
        
        if (settings.debug()) {
            List<FeatureComparisonStat> tmpDbg = new ArrayList<FeatureComparisonStat>();
            for (Entry<PairInt, List<FeatureComparisonStat>> entry : solutionStats0.entrySet()) {
                tmpDbg.addAll(entry.getValue());
            }
            MiscDebug.plotImages(tmpDbg, img1Helper.getGreyscaleImage(), 
            img2Helper.getGreyscaleImage(), 2, "_non_euclid_degen_" + settings.getDebugTag());        
        }
        
        copyToInstanceVars(solutionStats0);
          
        return true;
    }

    private Map<PairInt, List<FeatureComparisonStat>> reviseStatsForFullImages(
        GreyscaleImage img1, GreyscaleImage img2, 
        Map<PairInt, List<FeatureComparisonStat>> solutionStats, 
        int prevBinFactor1, int prevBinFactor2, RotatedOffsets rotatedOffsets) {
      
        log.info("refine stats for full image reference frames");
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        IntensityFeatures f1 = new IntensityFeatures(5, 
            settings.useNormalizedFeatures(), rotatedOffsets);
        f1.calculateGradientWithGreyscale(img1);
        
        IntensityFeatures f2 = new IntensityFeatures(5, 
            settings.useNormalizedFeatures(), rotatedOffsets);
        f2.calculateGradientWithGreyscale(img2);
        
        int dither = 2;
        
        Map<PairInt, List<FeatureComparisonStat>> revised = 
            new HashMap<PairInt, List<FeatureComparisonStat>>();
        
        for (Entry<PairInt, List<FeatureComparisonStat>> entry : solutionStats.entrySet()) {
            
            List<FeatureComparisonStat> stats = entry.getValue();
            
            List<FeatureComparisonStat> output = new ArrayList<FeatureComparisonStat>();
            
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
                        f2, x1, y1, x2, y2, dither, img1, img2);

                if (compStat == null || (compStat.getSumIntensitySqDiff() > 
                    compStat.getImg2PointIntensityErr())) {
                    continue;
                }
            
                output.add(compStat);
            }
            
            PairInt revisedKey = entry.getKey().copy();
            revisedKey.setX(revisedKey.getX() * prevBinFactor1);
            revisedKey.setY(revisedKey.getY() * prevBinFactor1);
            
            revised.put(revisedKey, output);
        }
        
        return revised;
    }

    private void copyToInstanceVars(Map<PairInt, List<FeatureComparisonStat>> 
        degenerateSolutionStats) {   
        
        this.solutionStats = new ArrayList<FeatureComparisonStat>();
        
        this.solutionMatched1 = new ArrayList<PairInt>();
        
        this.solutionMatched2 = new ArrayList<PairInt>();
        
        this.solutionStatsMultiplicity = new HashMap<PairInt, List<FeatureComparisonStat>>();
        
        this.solutionMatched2Multiplicity = new ArrayList<List<PairInt>>();
  
        for (Entry<PairInt, List<FeatureComparisonStat>> entry :
            degenerateSolutionStats.entrySet()) {
            
            if (entry.getValue().size() == 0) {
                continue;
            }
            
            PairInt p1 = entry.getKey();
            
            PairInt p2 = null;
                        
            double bestP2SSD = Double.MAX_VALUE;
            FeatureComparisonStat bestStat = null;
            
            List<PairInt> p2List = new ArrayList<PairInt>();
            for (FeatureComparisonStat stat : entry.getValue()) {
                if (stat.getSumIntensitySqDiff() < bestP2SSD) {
                    p2 = stat.getImg2Point();
                    bestStat = stat;
                }
                p2List.add(stat.getImg2Point());
            }
            
            solutionStats.add(bestStat);
            solutionMatched1.add(p1);
            solutionMatched2.add(p2);
            
            this.solutionStatsMultiplicity.put(p1, entry.getValue());
            this.solutionMatched2Multiplicity.add(p2List);
        }
    }
    
    public Map<PairInt, List<FeatureComparisonStat>> getSolutionStatsMultiplicity() {
        return solutionStatsMultiplicity;
    }
            
    public List<List<PairInt>> getSolutionMatched2Multiplicity() {
        return solutionMatched2Multiplicity;
    }

}
