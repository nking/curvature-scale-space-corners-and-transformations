package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

/**
 * class to invoke methods needed to solve for euclidean scale between
 * image1 and image2 using methods specific to corners on closed curves.
 * 
 * @author nichole
 */
public class BlobCornersScaleFinder extends AbstractBlobScaleFinder {

    public MatchingSolution solveForScale(
        BlobCornerHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1, 
        BlobCornerHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2) {
        
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
        
        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);
                
        assert(blobs1.size() == perimeters1.size());
        assert(blobs1.size() == corners1List.size());
                 
        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
        
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            mMap = match(features1, features2, img1, img2, perimeters1, 
            perimeters2, corners1List, corners2List);

        if (mMap.isEmpty()) {
            return null;
        }
        
        MatchingSolution soln = checkForNonDegenerateSolution(mMap, binFactor1, 
            binFactor2);
       
        if (soln != null) {
            return soln;
        }
        
        List<IntensityFeatureComparisonStats> ifsList = new ArrayList<IntensityFeatureComparisonStats>();
        List<TransformationParameters> paramsList = new ArrayList<TransformationParameters>();
        Map<PairInt, Float> indexScore = new HashMap<PairInt, Float>();
        
        //boolean useBipartite = true;
        
        //if (useBipartite) {
            resolveUsingBipartite(mMap, corners1List, corners2List, binFactor1,
                binFactor2, ifsList, paramsList, indexScore);
        //} else {
        //    resolveWithoutBipartite(mMap, corners1List, corners2List, binFactor1,
        //        binFactor2, ifsList, paramsList, indexScore);
        //}
        
        // to correct for wrap around from 360 to 0, repeating same calc with shifted values
       
        int[] indexesToKeep = MiscStats.filterForRotationUsingHist(paramsList, 0);
        
        int[] indexesToKeepShifted = MiscStats.filterForRotationUsingHist(paramsList, 30);
        
        if (indexesToKeepShifted.length > indexesToKeep.length) {
            indexesToKeep = indexesToKeepShifted;
        }
        
        filter(ifsList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForScaleUsingHist(paramsList);
        
        filter(ifsList, paramsList, indexesToKeep);
                
        indexesToKeep = MiscStats.filterForTranslationXUsingHist(paramsList);
        
        filter(ifsList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForTranslationYUsingHist(paramsList);
        
        filter(ifsList, paramsList, indexesToKeep);
        
        if (paramsList.size() == 0) {
            return null;
        }
        
        List<FeatureComparisonStat> combined = new ArrayList<FeatureComparisonStat>();
        for (int i = 0; i < ifsList.size(); ++i) {
            IntensityFeatureComparisonStats ifs = ifsList.get(i);
            combined.addAll(ifs.getComparisonStats());
        }
        
        TransformationParameters combinedParams = calculateTransformation(
            binFactor1, binFactor2, combined, new float[4]);
        
        if (combinedParams == null) {
            return null;
        }
        
        // pre-check for delta tx, delta ty essentially.  the parameter limits
        // are hard wired and the same as used elsewhere.
        boolean check = true;
        while (check && (paramsList.size() > 1)) {
            boolean small = MiscStats.standardDeviationsAreSmall(combinedParams);
            if (small) {
                check = false;
            } else {
                // --- either keep only smallest SSD or remove highest SSD ---
                double maxCost = Double.MIN_VALUE;
                int maxCostIdx = -1;
                for (int ii = 0; ii < ifsList.size(); ++ii) {
                    IntensityFeatureComparisonStats ifs = ifsList.get(ii);
                    PairInt p = new PairInt(ifs.getIndex1(), ifs.getIndex2());
                    float score1 = indexScore.get(p).floatValue();
                    if (score1 > maxCost) {
                        maxCost = score1;
                        maxCostIdx = ii;
                    }
                }
                ifsList.remove(maxCostIdx);
                paramsList.remove(maxCostIdx);
                combined.clear();
                for (int i = 0; i < ifsList.size(); ++i) {
                    IntensityFeatureComparisonStats ifs = ifsList.get(i);
                    combined.addAll(ifs.getComparisonStats());
                }
                combinedParams = calculateTransformation(
                    binFactor1, binFactor2, combined, new float[4]);
                if (combinedParams == null) {
                    return null;
                }
            }
        }
        
        soln = new MatchingSolution(combinedParams, combined);
            
        return soln;      
    }

    private <T extends CornerRegion> Map<Integer, 
        FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
        match(IntensityFeatures features1, IntensityFeatures features2, 
        GreyscaleImage img1, GreyscaleImage img2, 
        List<PairIntArray> perimeters1, List<PairIntArray> perimeters2, 
        List<List<T>> corners1List, List<List<T>> corners2List) {
/*        
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
float[] xPoints1 = new float[perimeters1.size()];
float[] yPoints1 = new float[perimeters1.size()];
double[][] xy1 = new double[perimeters1.size()][2];
for (int i = 0; i < perimeters1.size(); ++i) {
xy1[i] = curveHelper.calculateXYCentroids(perimeters1.get(i));
xPoints1[i] = (float)xy1[i][0];
yPoints1[i] = (float)xy1[i][1];
}
float[] xPoints2 = new float[perimeters2.size()];
float[] yPoints2 = new float[perimeters2.size()];
double[][] xy2 = new double[perimeters2.size()][2];
for (int i = 0; i < perimeters2.size(); ++i) {
xy2[i] = curveHelper.calculateXYCentroids(perimeters2.get(i));
xPoints2[i] = (float)xy2[i][0];
yPoints2[i] = (float)xy2[i][1];
}
ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
plotter.plotLabeledPoints(0, img1.getWidth(), 0, img1.getHeight(), xPoints1, yPoints1,
"img1", "X", "Y");
ScatterPointPlotterPNG plotter2 = new ScatterPointPlotterPNG();
plotter2.plotLabeledPoints(0, img2.getWidth(), 0, img2.getHeight(), xPoints2, yPoints2,
"img2", "X", "Y");
try {
    plotter.writeToFile("img1_labelled.png");
    plotter2.writeToFile("img2_labelled.png");
} catch (IOException ex) {
    Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
} 
*/

        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> trMap 
            = new HashMap<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>>();

        int n1 = corners1List.size();
        int n2 = corners2List.size();
        
        if (n1 == 0 || n2 == 0) {
            return trMap;
        }
        
        for (int idx1 = 0; idx1 < n1; ++idx1) {

            if (corners1List.get(idx1).size() < 3) {
                continue;
            }

            Integer index1 = Integer.valueOf(idx1);

            List<T> corners1 = corners1List.get(idx1);
            Collections.sort(corners1, new DescendingKComparator());
            
            //TODO: for extreme case of image full of repeated patterns, may need nTop=n2
            int nTop = n2/2;
            if (nTop == 0) {
                nTop = n2;
            }
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> bestMatches = 
                new FixedSizeSortedVector<IntensityFeatureComparisonStats>(nTop, 
                IntensityFeatureComparisonStats.class);
                                    
            for (int idx2 = 0; idx2 < n2; ++idx2) {

                if (corners2List.get(idx2).size() < 3) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                List<T> corners2 = corners2List.get(idx2);
                
                Collections.sort(corners2, new DescendingKComparator());

                ClosedCurveCornerMatcherWrapper<T> mapper =
                    new ClosedCurveCornerMatcherWrapper<T>();

                boolean matched = mapper.matchCorners(features1, features2, 
                    corners1, corners2, true, img1, img2);

                if (!matched) {
                    continue;
                }
                
                TransformationPair2<T> transformationPair = mapper.getSolution();
                transformationPair.setCornerListIndex1(idx1);
                transformationPair.setCornerListIndex2(idx2);
                
                TransformationParameters params = 
                    transformationPair.getTransformationParameters();
                
                if (params == null) {
                    continue;
                }
                
                @SuppressWarnings({"unchecked"})
                List<FeatureComparisonStat> compStats = 
                    transformationPair.getNextCorner().getMatchedFeatureComparisonStats();
                                    
                log.fine("theta diff filtered: " + printToString(compStats) 
                    + " combinedStat=" + calculateCombinedIntensityStat(compStats));
            
                FeatureMatcher.removeIntensityOutliers(compStats);
                
                if (compStats.size() < 3) {
                    continue;
                }
                
    FeatureMatcher.removeDiscrepantThetaDiff(compStats);

    if (compStats.size() < 3) {
        continue;
    }
    
                if (mapper.getSolvedCost() >= 1500) {
                    continue;
                }

                IntensityFeatureComparisonStats stats = new 
                    IntensityFeatureComparisonStats(index1.intValue(), 
                    index2.intValue(), mapper.getSolvedCost(), 
                    params.getScale());
              
                stats.addAll(compStats);
                         
                boolean added = bestMatches.add(stats);
            }
            
            if (bestMatches.getNumberOfItems() > 0) {
                trMap.put(Integer.valueOf(idx1), bestMatches);
            }            
        }

        return trMap;
    }

    private MatchingSolution checkForNonDegenerateSolution(Map<Integer, 
        FixedSizeSortedVector<IntensityFeatureComparisonStats>> mMap, 
        int binFactor1, int binFactor2) {
        
         List<IntensityFeatureComparisonStats> ifcsList = new ArrayList<IntensityFeatureComparisonStats>();
        
        Set<Integer> indexes2 = new HashSet<Integer>();
        
        for (Map.Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> entry :
            mMap.entrySet()) {
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> vector = entry.getValue();
            
            if (vector.getArray().length > 1) {
                return null;
            }
            assert(vector.getArray().length == 1);
            
            IntensityFeatureComparisonStats ifs = vector.getArray()[0];
            
            Integer key2 = Integer.valueOf(ifs.getIndex2());
            
            if (indexes2.contains(key2)) {
                return null;
            }
            
            indexes2.add(key2);
            
            ifcsList.add(ifs);
        }
        
        // the matches are unique, so will look for the smallest cost
        // and then add similar to it
        
        double minCost = Double.MAX_VALUE;
        int minCostIdx = -1;
        for (int i = 0; i < ifcsList.size(); ++i) {
            IntensityFeatureComparisonStats ifs = ifcsList.get(i);
            if (ifs.getCost() < minCost) {
                minCost = ifs.getCost();
                minCostIdx = i;
            }
        }
        
        TransformationParameters minCostParams = calculateTransformation(
            binFactor1, binFactor2, ifcsList.get(minCostIdx).getComparisonStats(),
            new float[4]);
        
        if (minCostParams == null) {
            return null;
        }
        
        List<FeatureComparisonStat> combined = new ArrayList<FeatureComparisonStat>();
        combined.addAll(ifcsList.get(minCostIdx).getComparisonStats());
        
        for (int i = 0; i < ifcsList.size(); ++i) {
            if (i == minCostIdx) {
                continue;
            }
            IntensityFeatureComparisonStats ifs = ifcsList.get(i);
            TransformationParameters params = calculateTransformation(
               binFactor1, binFactor2, ifs.getComparisonStats(),
                new float[4]);
            
            if (params == null) {
                continue;
            }
            
            if (Math.abs(params.getScale() - minCostParams.getScale()) < 0.05) {
                float angleDiff = AngleUtil.getAngleAverageInDegrees(
                    params.getRotationInDegrees(), minCostParams.getRotationInDegrees());
                if (Math.abs(angleDiff) < 10) {
                    if (Math.abs(params.getTranslationX() - minCostParams.getTranslationX()) < 10) {
                        if (Math.abs(params.getTranslationY() - minCostParams.getTranslationY()) < 10) {
                            combined.addAll(ifs.getComparisonStats());
                        }
                    }
                }
            }
        }
        
        TransformationParameters combinedParams = calculateTransformation(
            binFactor1, binFactor2, combined, new float[4]);
        
        if (combinedParams == null) {
            return null;
        }
        
        MatchingSolution soln = new MatchingSolution(combinedParams, combined);
        
        return soln;
    }

    private void filter(List<IntensityFeatureComparisonStats> ifsList,
        List<TransformationParameters> paramsList, int[] indexesToKeep) {
        
        if (indexesToKeep.length < 2) {
            return;
        }
        
        List<IntensityFeatureComparisonStats> ifsList2 = 
            new ArrayList<IntensityFeatureComparisonStats>(indexesToKeep.length);
        
        List<TransformationParameters> paramsList2 = 
            new ArrayList<TransformationParameters>();
        
        for (int i = 0; i < indexesToKeep.length;++i) {
            int idx = indexesToKeep[i];
            ifsList2.add(ifsList.get(idx));
            paramsList2.add(paramsList.get(idx));
        }
        
        ifsList.clear();
        ifsList.addAll(ifsList2);
        
        paramsList.clear();
        paramsList.addAll(paramsList2);
    }

    private void resolveUsingBipartite(Map<Integer, 
        FixedSizeSortedVector<IntensityFeatureComparisonStats>> mMap, 
        List<List<CornerRegion>> corners1List, 
        List<List<CornerRegion>> corners2List, 
        int binFactor1, int binFactor2,
        List<IntensityFeatureComparisonStats> ifsList, 
        List<TransformationParameters> paramsList,
        Map<PairInt, Float> indexScore) {
        
        int n1 = corners1List.size();
        int n2 = corners2List.size();
        
        float[][] cost = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new float[n2];
            Arrays.fill(cost[i], Float.MAX_VALUE);
        }
        
        int nMaxMatchable = countMaxMatchable(corners1List, corners2List);
        
        Set<PairInt> present = new HashSet<PairInt>();
        
        int tolTransXY = 5;//10;
        
        for (Map.Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> entry 
            : mMap.entrySet()) {
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> vector = entry.getValue();
            
            IntensityFeatureComparisonStats[] ind1To2Pairs = vector.getArray();            
            
            /* for the cost, need to consider the evaluation of the parameters,
            and the SSD of the point.
            The score for the evaluation is nMaxMatchable/nEval and its range is
            1 to nMaxMatchable.
            The SSD is filtered above to a max of 1500.  Will add '1' to it
            to avoid a zero for a perfect match, then the score for SSD ranges
            from 1 to 1500.
            Can make the cost the multiplication of the two scores as long as
            the value (nMaxMatchable/1500) stays below ((1<<31)-1).
            
            have normalized both scores by their maximum values so that their
            contributions to the cost are equal.
            */
            
            for (int i = 0; i < vector.getNumberOfItems(); ++i) {
                
                IntensityFeatureComparisonStats ifs = ind1To2Pairs[i];              
                
                TransformationParameters params = calculateTransformation(
                    binFactor1, binFactor2, ifs.getComparisonStats(),
                    new float[4]);
                if (params == null) {
                    continue;
                }
                int idx1 = ifs.getIndex1();
                int idx2 = ifs.getIndex2();
                PairInt p = new PairInt(idx1, idx2);
                int nEval = evaluate(params, corners1List, corners2List, tolTransXY);
                if (nEval == 0) {
                    continue;
                }
                //float score1 = (float)nMaxMatchable/(float)nEval;
                float score2 = (float)ifs.getCost() + 1;
                //float score = score1 * score2;
                float normalizedScore = (score2/(float)nEval)/1500.f;
                cost[idx1][idx2] = normalizedScore;
                indexScore.put(p, Float.valueOf(normalizedScore));
                present.add(p);
                ifsList.add(ifs);
                paramsList.add(params); 
            }
        }

        boolean transposed = false;
        if (n1 > n2) {
            cost = MatrixUtil.transpose(cost);
            transposed = true;
        }

        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(cost);
        
        Set<PairInt> matched = new HashSet<PairInt>();
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            PairInt p = new PairInt(idx1, idx2);
            if (present.contains(p)) {
                 matched.add(p);
            }
        }
        
        int n = ifsList.size();
        int i = 0;
        while (i < n) {
            IntensityFeatureComparisonStats ifs = ifsList.get(i);
            PairInt p = new PairInt(ifs.getIndex1(), ifs.getIndex2());
            if (matched.contains(p)) {
                ++i;
                continue;
            }
            ifsList.remove(i);
            paramsList.remove(i);
            n = ifsList.size();
        }
        
    }

    /*
    private void resolveWithoutBipartite(Map<Integer, 
        FixedSizeSortedVector<IntensityFeatureComparisonStats>> mMap, 
        List<List<CornerRegion>> corners1List, 
        List<List<CornerRegion>> corners2List, 
        int binFactor1, int binFactor2, 
        List<IntensityFeatureComparisonStats> ifsList, 
        List<TransformationParameters> paramsList, 
        Map<PairInt, Float> indexScore) {
                
        Map<Integer, Integer> index1Map = new HashMap<Integer, Integer>();
        Map<Integer, Set<Integer>> index2Map = new HashMap<Integer, Set<Integer>>();
        Map<Integer, IntensityFeatureComparisonStats> index1StatMap = new HashMap<Integer, IntensityFeatureComparisonStats>();
        Map<Integer, Double> index1ScoreMap = new HashMap<Integer, Double>();
        
        Map<PairInt, TransformationParameters> index12ParamsMap = new HashMap<PairInt, TransformationParameters>();
                
        int tolTransXY = 5;//10;
        
        for (Map.Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> entry 
            : mMap.entrySet()) {
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> vector = entry.getValue();
            
            IntensityFeatureComparisonStats[] ind1To2Pairs = vector.getArray();            
          
            for (int i = 0; i < vector.getNumberOfItems(); ++i) {
                
                IntensityFeatureComparisonStats ifs = ind1To2Pairs[i];              
                
                TransformationParameters params = calculateTransformation(
                    binFactor1, binFactor2, ifs.getComparisonStats(),
                    new float[4]);
                if (params == null) {
                    continue;
                }
                int idx1 = ifs.getIndex1();
                int idx2 = ifs.getIndex2();
                PairInt p = new PairInt(idx1, idx2);
                
                int nEval = evaluate(params, corners1List, corners2List, tolTransXY);
                if (nEval == 0) {
                    continue;
                }
                //double score1 = (double)nMaxMatchable/(double)nEval;
                double score2 = (double)ifs.getCost() + 1.;
                //float score = score1 * score2;
                double normalizedScore = (score2/(double)nEval)/1500.;
                
                Integer key1 = Integer.valueOf(idx1);
                Double bestCost = index1ScoreMap.get(key1);
                
                if ((bestCost == null) || (bestCost.doubleValue()> normalizedScore)) {
                    index1StatMap.put(key1, ifs);
                    index1ScoreMap.put(key1, Double.valueOf(normalizedScore));
                    index12ParamsMap.put(p, params);
                }
            }
        }
        
        for (Entry<Integer, IntensityFeatureComparisonStats> entry :
            index1StatMap.entrySet()) {
            
            IntensityFeatureComparisonStats ifs = entry.getValue();
            
            Integer index1 = entry.getKey();
            Integer index2 = Integer.valueOf(ifs.getIndex2());
            assert(ifs.getIndex1() == index1.intValue());
            
            index1Map.put(index1, index2);
            
            Set<Integer> indexes1 = index2Map.get(index2);
            if (indexes1 == null) {
                indexes1 = new HashSet<Integer>();
                index2Map.put(index2, indexes1);
            }
            indexes1.add(index1);
        }
        
        // resolve any double matchings, but discard the higher cost matches
        //   from conflicted matches rather than re-trying a solution for them
        
        Set<Integer> resolved = new HashSet<Integer>();
        for (Entry<Integer, Set<Integer>> entry : index2Map.entrySet()) {
            if (resolved.contains(entry.getKey())) {
                continue;
            }
            Set<Integer> set = entry.getValue();
            if (set.size() > 1) {
                double bestCost = Double.MAX_VALUE;
                Integer bestIndex1 = -1;
                for (Integer index1 : set) {
                    Double cost = index1ScoreMap.get(index1);
                    assert(cost != null);
                    if (cost < bestCost) {
                        bestCost = cost;
                        bestIndex1 = index1;
                    }
                    resolved.add(index1);
                }
                for (Integer index1 : set) {
                    if (index1.equals(bestIndex1)) {
                        continue;
                    }
                    index1Map.remove(index1);
                    index1StatMap.remove(index1);
                }
            }
        }
        
        for (Entry<Integer, IntensityFeatureComparisonStats> entry :
            index1StatMap.entrySet()) {
            
            IntensityFeatureComparisonStats ifs = entry.getValue();
            PairInt p = new PairInt(ifs.getIndex1(), ifs.getIndex2());
            
            TransformationParameters params = index12ParamsMap.get(p);
            
            Double score = index1ScoreMap.get(entry.getKey());
            
            ifsList.add(ifs);
            paramsList.add(params);
            indexScore.put(p, Float.valueOf(score.floatValue()));
        }
    }*/
    
}
