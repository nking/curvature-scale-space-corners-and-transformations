package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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

    public TransformationParameters solveForScale(BlobCornerHelper img1Helper, 
        SegmentationType type1, boolean useBinned1, BlobCornerHelper img2Helper, 
        SegmentationType type2, boolean useBinned2,
        float[] outputScaleRotTransXYStDev) {

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
            
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);
        
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
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);
        
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            index1BestMap = new HashMap<Integer, 
            FixedSizeSortedVector<IntensityFeatureComparisonStats>>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            if (corners1List.get(idx1).isEmpty()) {
                continue;
            }

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = perimeters1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);
            
            List<CornerRegion> corners1 = corners1List.get(idx1);
            Collections.sort(corners1, new DescendingKComparator());

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            // keeping the top '1' for each index1.  comparison is by cost.
            // choosing more than one because later bipartite matching attempts
            // to match best for all index1 matchings
            FixedSizeSortedVector<IntensityFeatureComparisonStats> bestStats 
                = new FixedSizeSortedVector<IntensityFeatureComparisonStats>(1, 
                IntensityFeatureComparisonStats.class);
            
            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                if (corners2List.get(idx2).isEmpty()) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = perimeters2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);
                
                List<CornerRegion> corners2 = corners1List.get(idx2);
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
                                
                removeDiscrepantThetaDiff(compStats);

                log.info("theta diff filtered: " + printToString(compStats) 
                    + " combinedStat=" + calculateCombinedIntensityStat(compStats));
            
                removeIntensityOutliers(compStats);
                
                IntensityFeatureComparisonStats stats = new 
                    IntensityFeatureComparisonStats(index1.intValue(), 
                    index2.intValue(), mapper.getSolvedCost(), 
                        mapper.getSolvedParameters().getScale());
                
                stats.addAll(compStats);
                
                // bestStats keeps the top '1' smallest cost solutions added to it
                // (though combinedStats are used when nMatched is 2 or less)
                boolean added = bestStats.add(stats);
                
                if (added) {
                    log.info("  added to best for [" + index1.toString() + "] ["
                        + index2.toString() + "] cost=" + stats.getCost()
                        + " with n=" + stats.getComparisonStats().size());
                }
            }
            
            if (bestStats.getNumberOfItems() == 0) {
                continue;
            }
            
            index1BestMap.put(index1, bestStats);
        }

        List<FeatureComparisonStat> bestOverall = null;
        if (!index1BestMap.isEmpty()) {
            bestOverall = filterToBestConsistent(index1BestMap, corners1List,
                corners2List);
        }

        if (bestOverall == null) {
            return null;
        }

        TransformationParameters params = calculateTransformation(
            img1Helper.imgHelper.getBinFactor(useBinned1), type1, useBinned1,
            img2Helper.imgHelper.getBinFactor(useBinned2), type2, useBinned2,
            bestOverall, outputScaleRotTransXYStDev);

        if (params == null) {
            return null;
        }

        return params;
    }

    private List<FeatureComparisonStat> filterToBestConsistent(
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
        index1BestMap, List<List<CornerRegion>> corners1List, 
        List<List<CornerRegion>> corners2List) {
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
}
