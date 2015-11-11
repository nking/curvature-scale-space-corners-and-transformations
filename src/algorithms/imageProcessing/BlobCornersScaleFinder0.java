package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.clustering.FixedDistanceGroupFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.MiscStats;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

/**
 * class to attempt to find Euclidean scale transformation between image1
 * and image2 using corners of a blob and an assumption of
 * ordered matching between them.
 * The method is fast, but cannot be used to compare the same object boundaries
 * if the boundaries have different numbers of corners.
 * 
 * @author nichole
 */
public class BlobCornersScaleFinder0 extends AbstractBlobScaleFinder {

    protected void filterToSameNumberIfPossible(
        PairIntArray curve1, List<CornerRegion> regions1,
        PairIntArray curve2, List<CornerRegion> regions2) {
        
        int dist = 2;
        
        if (regions1.size() != regions2.size()) {
            filterCorners(curve1, regions1, dist);
            filterCorners(curve2, regions2, dist);
        }
    }
    
    protected void filterCorners(PairIntArray curve, 
        List<CornerRegion> regions, int dist) {
        
        /*
        if there are more than 1 corner within dist of 2 or so of on another,
        remove all except strongest corner.
        */
        List<Set<Integer>> closeCornerIndexes = findCloseCorners(dist, regions);
        
        if (closeCornerIndexes.isEmpty()) {
            return;
        }
        
        List<Integer> remove = new ArrayList<Integer>();
        for (Set<Integer> set : closeCornerIndexes) {
            float maxK = Float.MIN_VALUE;
            Integer maxKIndex = null;
            for (Integer index : set) {
                CornerRegion cr = regions.get(index.intValue());
                float k = cr.getK()[cr.getKMaxIdx()];
                if (k > maxK) {
                    maxK = k;
                    maxKIndex = index;
                }
            }
            for (Integer index : set) {
                if (!index.equals(maxKIndex)) {
                    remove.add(index);
                }
            }
        }
        
        if (remove.size() > 1) {
            Collections.sort(remove);
        }
        for (int i = (remove.size() - 1); i > -1; --i) {
            regions.remove(remove.get(i).intValue());
        }
    }
    
    private static List<Set<Integer>> findCloseCorners(int tolD, 
        List<CornerRegion> regions) {
       
        float[] x = new float[regions.size()];
        float[] y = new float[regions.size()];
        for (int i = 0; i < regions.size(); ++i) {
            CornerRegion cr = regions.get(i);            
            x[i] = cr.getX()[cr.getKMaxIdx()];
            y[i] = cr.getY()[cr.getKMaxIdx()];
        }
        
        List<Set<Integer>> close = new ArrayList<Set<Integer>>();
        
        FixedDistanceGroupFinder groupFinder = new FixedDistanceGroupFinder(x, y);

        groupFinder.findGroupsOfPoints(tolD);
        
        int nGroups = groupFinder.getNumberOfGroups();
        
        for (int i = 0; i < nGroups; ++i) {
            Set<Integer> group = groupFinder.getGroupIndexes(i);
            if (group.size() > 1) {
                close.add(group);
            }
        }
        
        return close;
    }
  
    public MatchingSolution solveForScale(
        BlobCornerHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobCornerHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2) {

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
            
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);

        List<List<CornerRegion>> corners1List = 
            img1Helper.getPerimeterCorners(type1, useBinned1);
        
        List<List<CornerRegion>> corners2List = 
            img2Helper.getPerimeterCorners(type2, useBinned2);
        
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(type2, useBinned2);

        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
                
        int n1 = corners1List.size();
        int n2 = corners2List.size();
                   
        /*
        for some images, there is a repetitive pattern.  projection effects
        such as shadows or occlusion may lead to a false match with a lower
        cost than the true match.
        
        so the frequency of similar transformation solutions has to be considered.

        The matches are kept in a sorted vector until all matches are tried.
        */

        Map<Integer, FixedSizeSortedVector<TransformationPair4>> mMap =
            match(features1, features2, img1, img2, 
                perimeters1, perimeters2, corners1List, corners2List);
                        
        if (mMap.isEmpty()) {
            return null;
        }
        
        MatchingSolution soln = checkForNonDegenerateSolution(mMap, binFactor1, 
            binFactor2);
       
        if (soln != null) {
            return soln;
        }
        
        float[][] cost = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new float[n2];
            Arrays.fill(cost[i], Float.MAX_VALUE);
        }
        
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
        
        int nMaxMatchable = countMaxMatchable(corners1List, corners2List);
                
        Map<PairInt, Float> indexScore = new HashMap<PairInt, Float>();

        Set<PairInt> present = new HashSet<PairInt>();
        List<TransformationPair4> tpList = new ArrayList<TransformationPair4>();
        List<TransformationParameters> paramsList = new ArrayList<TransformationParameters>();
        int tolTransXY = 10;
        for (Entry<Integer, FixedSizeSortedVector<TransformationPair4>> entry 
            : mMap.entrySet()) {
            FixedSizeSortedVector<TransformationPair4> vector = entry.getValue();
            TransformationPair4[] ind1To2Pairs = vector.getArray();            
            for (int i = 0; i < vector.getNumberOfItems(); ++i) {
                TransformationPair4 tp4 = ind1To2Pairs[i];              
                
                TransformationParameters params = calculateTransformation(
                    binFactor1, binFactor2, tp4.getMatchedCompStats(),
                    new float[4]);
                if (params == null) {
                    continue;
                }
                int idx1 = tp4.getCornerListIndex1();
                int idx2 = tp4.getCornerListIndex2();
                PairInt p = new PairInt(idx1, idx2);
                
                int nEval = evaluate(params, corners1List, corners2List, tolTransXY);
                if (nEval == 0) {
                    continue;
                }
                float score1 = (float)nMaxMatchable/(float)nEval;
                float score2 = (float)tp4.getCost() + 1;
                float score = score1 * score2;
                float normalizedScore = (score2/(float)nEval)/1500.f;
                cost[idx1][idx2] = normalizedScore;
                indexScore.put(p, Float.valueOf(normalizedScore));
                
                present.add(p);
                
                tpList.add(tp4);
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
        
        int n = tpList.size();
        int i = 0;
        while (i < n) {
            TransformationPair4 tp4 = tpList.get(i);
            PairInt p = new PairInt(tp4.getCornerListIndex1(), tp4.getCornerListIndex2());
            if (matched.contains(p)) {
                ++i;
                continue;
            }
            tpList.remove(i);
            paramsList.remove(i);
            n = tpList.size();
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
        
        filter(tpList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForScaleUsingHist(paramsList);
        
        filter(tpList, paramsList, indexesToKeep);
        
        //indexesToKeep = MiscStats.filterForTranslation(paramsList);
        
        indexesToKeep = MiscStats.filterForTranslationXUsingHist(paramsList);
        
        filter(tpList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForTranslationYUsingHist(paramsList);
        
        filter(tpList, paramsList, indexesToKeep);
        
        if (paramsList.size() == 0) {
            return null;
        }
        
        List<FeatureComparisonStat> combined = new ArrayList<FeatureComparisonStat>();
        for (i = 0; i < tpList.size(); ++i) {
            TransformationPair4 ifs = tpList.get(i);
            combined.addAll(ifs.getMatchedCompStats());
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
                for (int ii = 0; ii < tpList.size(); ++ii) {
                    TransformationPair4 ifs = tpList.get(ii);
                    PairInt p = new PairInt(ifs.getCornerListIndex1(), ifs.getCornerListIndex2());
                    float score1 = indexScore.get(p).floatValue();
                    if (score1 > maxCost) {
                        maxCost = score1;
                        maxCostIdx = ii;
                    }
                }
                tpList.remove(maxCostIdx);
                paramsList.remove(maxCostIdx);
                combined.clear();
                for (i = 0; i < tpList.size(); ++i) {
                    TransformationPair4 ifs = tpList.get(i);
                    combined.addAll(ifs.getMatchedCompStats());
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

    private Map<Integer, FixedSizeSortedVector<TransformationPair4>> match(
        IntensityFeatures features1, IntensityFeatures features2, 
        GreyscaleImage img1, GreyscaleImage img2, 
        List<PairIntArray> perimeters1, List<PairIntArray> perimeters2, 
        List<List<CornerRegion>> corners1List, 
        List<List<CornerRegion>> corners2List) {
        
        Map<Integer, FixedSizeSortedVector<TransformationPair4>> trMap 
            = new HashMap<Integer, FixedSizeSortedVector<TransformationPair4>>();

        int n1 = corners1List.size();
        int n2 = corners2List.size();
        
        if (n1 == 0 || n2 == 0) {
            return trMap;
        }
/*        
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
double[][] xy1 = new double[perimeters1.size()][2];
for (int i = 0; i < perimeters1.size(); ++i) {
xy1[i] = curveHelper.calculateXYCentroids(perimeters1.get(i));
}
double[][] xy2 = new double[perimeters2.size()][2];
for (int i = 0; i < perimeters2.size(); ++i) {
xy2[i] = curveHelper.calculateXYCentroids(perimeters2.get(i));
}
*/
        for (int i = 0; i < n1; ++i) {
                                    
            List<CornerRegion> corners1 = corners1List.get(i);
            int nc1 = corners1.size();
            if (nc1 < 3) {
                continue;
            }
            
            //TODO: for extreme case of image full of repeated patterns, may need nTop=n2
            int nTop = n2/2;
            if (nTop == 0) {
                nTop = n2;
            }
            
            FixedSizeSortedVector<TransformationPair4> bestMatches = 
                new FixedSizeSortedVector<TransformationPair4>(nTop, 
                TransformationPair4.class);
            
            for (int j = 0; j < n2; ++j) {
                                
                List<CornerRegion> corners2 = corners2List.get(j);
                int nc2 = corners2.size();
                if (nc2 < 3) {
                    continue;
                }

                int minN = Math.min(nc1, nc2);
                int diffNC = Math.abs(nc1 - nc2);
                
                if (minN >= 9) {
                    if (diffNC > 4) {
                        continue;
                    }
                } else if ((minN < 9) && (diffNC > 2)){
                    continue;
                }
                
                List<CornerRegion> c1;
                List<CornerRegion> c2;
                
                if (nc1 == nc2) {
                    c1 = corners1;
                    c2 = corners2;
                } else {
                    c1 = new ArrayList<CornerRegion>(corners1);
                    c2 = new ArrayList<CornerRegion>(corners2);                  
                    filterToSameNumberIfPossible(perimeters1.get(i), c1,
                        perimeters2.get(j), c2);
                    
                    if (c1.size() != c2.size()) {
                        continue;
                    }
                }
                
                ClosedCurveCornerMatcher0 matcher = new ClosedCurveCornerMatcher0();
                
                // discrepant thetas are removed from these results:
                boolean solved = matcher.matchCorners(features1, features2, c1, 
                    c2, true, img1, img2);
                
                if (solved) {
                    
                    TransformationPair4 tr = matcher.getTransformationPair();
                    tr.setCornerListIndex1(i);
                    tr.setCornerListIndex2(j);
                    
                    if (tr.getCost() < 1500) {
                        boolean added = bestMatches.add(tr); 
                    }
                }
            }
            
            if (bestMatches.getNumberOfItems() > 0) {
                trMap.put(Integer.valueOf(i), bestMatches);
            }
        }
                
        return trMap;
    }

    private void filter(List<TransformationPair4> tpList, 
        List<TransformationParameters> paramsList, int[] indexesToKeep) {
        
        if (indexesToKeep.length < 2) {
            return;
        }
        
        List<TransformationPair4> tpList2 = 
            new ArrayList<TransformationPair4>(indexesToKeep.length);
        
        List<TransformationParameters> paramsList2 = 
            new ArrayList<TransformationParameters>();
        
        for (int i = 0; i < indexesToKeep.length;++i) {
            int idx = indexesToKeep[i];
            tpList2.add(tpList.get(idx));
            paramsList2.add(paramsList.get(idx));
        }
        
        tpList.clear();
        tpList.addAll(tpList2);
        
        paramsList.clear();
        paramsList.addAll(paramsList2);
    }

    /**
     * checking for case where each key has only 1 value in the vector,
     * and each index1 has a unique mapping to index2, and if found,
     * choosing the lowest cost solution and adding similar solutions
     * to it.
     * @param mMap
     * @return 
     */
    private MatchingSolution checkForNonDegenerateSolution(
        Map<Integer, FixedSizeSortedVector<TransformationPair4>> mMap,
        int binFactor1, int binFactor2) {
        
        List<TransformationPair4> tpList = new ArrayList<TransformationPair4>();
        
        Set<Integer> indexes2 = new HashSet<Integer>();
        
        for (Entry<Integer, FixedSizeSortedVector<TransformationPair4>> entry :
            mMap.entrySet()) {
            
            FixedSizeSortedVector<TransformationPair4> vector = entry.getValue();
            
            if (vector.getArray().length > 1) {
                return null;
            }
            assert(vector.getArray().length == 1);
            
            TransformationPair4 tp4 = vector.getArray()[0];
            
            Integer key2 = Integer.valueOf(tp4.getCornerListIndex2());
            
            if (indexes2.contains(key2)) {
                return null;
            }
            
            indexes2.add(key2);
            
            tpList.add(tp4);
        }
        
        // the matches are unique, so will look for the smallest cost
        // and then add similar to it
        
        double minCost = Double.MAX_VALUE;
        int minCostIdx = -1;
        for (int i = 0; i < tpList.size(); ++i) {
            TransformationPair4 tp4 = tpList.get(i);
            if (tp4.getCost() < minCost) {
                minCost = tp4.getCost();
                minCostIdx = i;
            }
        }
        
        TransformationParameters minCostParams = calculateTransformation(
            binFactor1, binFactor2, tpList.get(minCostIdx).getMatchedCompStats(),
            new float[4]);
        
        if (minCostParams == null) {
            return null;
        }
        
        List<FeatureComparisonStat> combined = new ArrayList<FeatureComparisonStat>();
        combined.addAll(tpList.get(minCostIdx).getMatchedCompStats());
        
        for (int i = 0; i < tpList.size(); ++i) {
            if (i == minCostIdx) {
                continue;
            }
            TransformationPair4 tp4 = tpList.get(i);
            TransformationParameters params = calculateTransformation(
               binFactor1, binFactor2, tp4.getMatchedCompStats(),
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
                            combined.addAll(tp4.getMatchedCompStats());
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

}
