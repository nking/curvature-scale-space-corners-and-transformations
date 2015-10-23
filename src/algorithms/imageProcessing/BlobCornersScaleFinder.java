package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * class to invoke methods needed to solve for euclidean scale between
 * image1 and image2 using methods specific to corners on closed curves.
 * 
 * @author nichole
 */
public class BlobCornersScaleFinder extends AbstractBlobScaleFinder {

    public TransformationParameters solveForScale(
        BlobCornerHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1, 
        BlobCornerHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {
        
        List<List<CornerRegion>> corners1List = img1Helper.getPerimeterCorners(
            type1, useBinned1);
        List<List<CornerRegion>> corners2List = img2Helper.getPerimeterCorners(
            type2, useBinned2);
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(
            type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(
            type2, useBinned2);
                
        assert(blobs1.size() == perimeters1.size());
        assert(blobs1.size() == corners1List.size());
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
         
        Map<Integer, IntensityFeatureComparisonStats> 
            index1BestMap = new HashMap<Integer, IntensityFeatureComparisonStats>();
        
        Map<Integer, TransformationParameters> 
            index1BestParamsMap = new HashMap<Integer, TransformationParameters>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            if (corners1List.get(idx1).size() < 3) {
                continue;
            }

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = perimeters1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);
            
            List<CornerRegion> corners1 = corners1List.get(idx1);
            Collections.sort(corners1, new DescendingKComparator());

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            IntensityFeatureComparisonStats bestStats = null;
            
            TransformationParameters bestParams = null;
            
            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                if (corners2List.get(idx2).size() < 3) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = perimeters2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);
                
                List<CornerRegion> corners2 = corners2List.get(idx2);
                Collections.sort(corners2, new DescendingKComparator());

double[] xyCen2 = curveHelper.calculateXYCentroids(curve2);
log.info("index1=" + index1.toString() + " index2=" + index2.toString()
+ " xyCen1=" + Arrays.toString(xyCen1) + " xyCen2=" + Arrays.toString(xyCen2));

log.info(
String.format("[%d](%d,%d) [%d](%d,%d)  nCurvePoints=%d, %d", 
idx1, (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
idx2, (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
curve1.getN(), curve2.getN()));               

                ClosedCurveCornerMatcherWrapper mapper =
                    new ClosedCurveCornerMatcherWrapper(features1, features2, 
                    corners1, corners2, true);
                
                boolean matched = mapper.matchCorners();
              
                if (!matched) {
                    continue;
                }
                
                TransformationParameters params = mapper.getSolvedParameters();
                
                if (params == null) {
                    continue;
                }
                
                List<FeatureComparisonStat> compStats = mapper.getSolutionMatchedCompStats();
                    
                log.info("theta diff filtered: " + printToString(compStats) 
                    + " combinedStat=" + calculateCombinedIntensityStat(compStats));
            
                removeIntensityOutliers(compStats);
                
                if (compStats.size() < 3) {
                    continue;
                }
                
                IntensityFeatureComparisonStats stats = new 
                    IntensityFeatureComparisonStats(index1.intValue(), 
                    index2.intValue(), mapper.getSolvedCost(), 
                        mapper.getSolvedParameters().getScale());
                
                stats.addAll(compStats);
                
                int comp = -1;
                if (bestStats != null) {
//TODO: revise the comparison
                    comp = stats.compareTo(bestStats);
                }
                if (comp == -1) {
                    bestStats = stats;
                    bestParams = params;                
                    log.info("  added to best for [" + index1.toString() + "] ["
                        + index2.toString() + "] cost=" + stats.getCost()
                        + " with n=" + stats.getComparisonStats().size());
                }
            }
            
            if (bestStats == null) {
                continue;
            }
            
            index1BestMap.put(index1, bestStats);
            
            index1BestParamsMap.put(index1, bestParams);
        }

        List<FeatureComparisonStat> bestOverall = null;
        if (!index1BestMap.isEmpty()) {
            bestOverall = filterToBestConsistent(index1BestMap, 
                index1BestParamsMap, corners1List, corners2List);
        }

        if (bestOverall == null) {
            return null;
        }

        TransformationParameters params = calculateTransformation(
            img1Helper.imgHelper.getBinFactor(useBinned1),
            img2Helper.imgHelper.getBinFactor(useBinned2),
            bestOverall, outputScaleRotTransXYStDev);

        if (params == null) {
            return null;
        }

        return params;
    }

    private List<FeatureComparisonStat> filterToBestConsistent(
        Map<Integer, IntensityFeatureComparisonStats> index1BestMap, 
        Map<Integer, TransformationParameters> index1BestParamsMap, 
        List<List<CornerRegion>> corners1List, 
        List<List<CornerRegion>> corners2List) {
        
        int bestCostIdx = -1;
        double bestCost = Double.MAX_VALUE;
               
        for (int i = 0; i < corners1List.size(); ++i) {
            
            Integer key = Integer.valueOf(i);
        
            TransformationParameters params = index1BestParamsMap.get(key);
            
            if (params == null) {
                continue;
            }
            
            IntensityFeatureComparisonStats ifs = index1BestMap.get(key);
            
            log.info("params=" + params.toString() + " cost=" + ifs.getCost());
            
            if (ifs.getCost() < bestCost) {
                bestCostIdx = i;
                bestCost = ifs.getCost();
            }
        }
        List<FeatureComparisonStat> compStats = new ArrayList<FeatureComparisonStat>();
        if (bestCostIdx == -1) {
            return compStats;
        }
        
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        
        Transformer transformer = new Transformer();
        
        /*
        params similar (by rot, scale, tx and ty) to the bestCost params can be 
        aggregated
        */
        TransformationParameters bestCostParams 
            = index1BestParamsMap.get(Integer.valueOf(bestCostIdx));
        IntensityFeatureComparisonStats bestIFS 
            = index1BestMap.get(Integer.valueOf(bestCostIdx));
        
        compStats.addAll(bestIFS.getComparisonStats());
        
        if (bestCostParams.getOriginX() != 0 || bestCostParams.getOriginY() != 0) {
            transformer.transformToOrigin(0, 0, bestCostParams);
        }
        
        for (int i = 0; i < corners1List.size(); ++i) {
            
            if (bestCostIdx == i) {
                continue;
            }
            
            Integer key = Integer.valueOf(i);
        
            TransformationParameters params = index1BestParamsMap.get(key);
            
            if (params == null) {
                continue;
            }
            
            if (!tc.areSimilarByScaleAndRotation(bestCostParams, params)) {
                continue;
            }
            
            IntensityFeatureComparisonStats ifs = index1BestMap.get(key);
            
            if (params.getOriginX() != 0 || params.getOriginY() != 0) {
                transformer.transformToOrigin(0, 0, params);
            }
            
            // if transX and transY are similar, add, else check difference with transformed
            boolean areSimilar = (Math.abs(bestCostParams.getTranslationX() - params.getTranslationX()) < 10)
                && (Math.abs(bestCostParams.getTranslationY() - params.getTranslationY()) < 10);
            
            if (areSimilar) {
                compStats.addAll(ifs.getComparisonStats());
            }            
        }
        
        removeDiscrepantThetaDiff(compStats);
        
        return compStats;
    }
    
}
