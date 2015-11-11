package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ScatterPointPlotterPNG;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
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
        
        int n1 = perimeters1.size();
        int n2 = perimeters2.size();
        
        float[][] cost = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new float[n2];
            Arrays.fill(cost[i], Float.MAX_VALUE);
        }
        
        int nMaxMatchable = countMaxMatchable(corners1List, corners2List);
        
        Set<PairInt> present = new HashSet<PairInt>();
        List<IntensityFeatureComparisonStats> ifsList = new ArrayList<IntensityFeatureComparisonStats>();
        List<TransformationParameters> paramsList = new ArrayList<TransformationParameters>();
        
        Map<PairInt, Float> indexScore = new HashMap<PairInt, Float>();
        
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
                float score1 = (float)nMaxMatchable/(float)nEval;
                float score2 = (float)ifs.getCost() + 1;
                float score = score1 * score2;
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
        
        // to correct for wrap around from 360 to 0, repeating same calc with shifted values
        
        //int[] indexesToKeep = MiscStats.filterForScaleAndRotationUsingHist(paramsList, 0);
        //int[] indexesToKeep2 = MiscStats.filterForScaleAndRotation(paramsList, 0);
        //int[] indexesToKeepShifted = MiscStats.filterForScaleAndRotationUsingHist(paramsList, 30);
        
        int[] indexesToKeep = MiscStats.filterForRotationUsingHist(paramsList, 0);
        
        int[] indexesToKeepShifted = MiscStats.filterForRotationUsingHist(paramsList, 30);
        
        if (indexesToKeepShifted.length > indexesToKeep.length) {
            indexesToKeep = indexesToKeepShifted;
        }
        
        filter(ifsList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForScaleUsingHist(paramsList);
        
        filter(ifsList, paramsList, indexesToKeep);
        
        //indexesToKeep = MiscStats.filterForTranslation(paramsList);
        
        indexesToKeep = MiscStats.filterForTranslationXUsingHist(paramsList);
        
        filter(ifsList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForTranslationYUsingHist(paramsList);
        
        filter(ifsList, paramsList, indexesToKeep);
        
        if (paramsList.size() == 0) {
            return null;
        }
        
        List<FeatureComparisonStat> combined = new ArrayList<FeatureComparisonStat>();
        for (i = 0; i < ifsList.size(); ++i) {
            IntensityFeatureComparisonStats ifs = ifsList.get(i);
            combined.addAll(ifs.getComparisonStats());
        }
        
        TransformationParameters combinedParams = calculateTransformation(
            binFactor1, binFactor2, combined, new float[4]);
        
        if (combinedParams == null) {
            return null;
        }
        
        // pre-check for delta tx, deltaty essentially
        boolean check = true;
        while (check && (paramsList.size() > 1)) {
            float tS = (combinedParams.getStandardDeviations()[0]/combinedParams.getScale());
            float tR = (float)(2.*Math.PI/combinedParams.getStandardDeviations()[1]);
            float tTx = combinedParams.getStandardDeviations()[2];
            float tTy = combinedParams.getStandardDeviations()[3];
            float tXConstraint = 20;
            float tYConstraint = 20;
            if (combinedParams.getNumberOfPointsUsed() < 3) {
                tXConstraint = 10;
                tYConstraint = 10;
            }
            if ((tS < 0.2) && (tR >= 18.) && (tTx < tXConstraint)
                && (tTy < tYConstraint)) {
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
                for (i = 0; i < ifsList.size(); ++i) {
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

    private <T extends CornerRegion> Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
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
                    new ClosedCurveCornerMatcherWrapper<>();
                
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
    
}
