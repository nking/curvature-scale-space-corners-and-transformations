package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscMath;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class BlobContoursScaleFinder extends AbstractBlobScaleFinder {

    public MatchingSolution solveForScale(
        BlobContourHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobContourHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2) {

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
            
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);
        
        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);
        
        List<List<CurvatureScaleSpaceContour>> contours1List = 
            img1Helper.getPerimeterContours(type1, useBinned1);
        
        List<List<CurvatureScaleSpaceContour>> contours2List = 
            img2Helper.getPerimeterContours(type2, useBinned2);
        
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, useBinned1);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, useBinned2);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(type1, useBinned1);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(type2, useBinned2);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        Map<PairInt, CSSContourMatcherWrapper> singleSolnMap =
            new HashMap<PairInt,  CSSContourMatcherWrapper>();

        Map<Integer, TransformationParameters> pMap = new HashMap<Integer,
            TransformationParameters>();
        
        /*
        Note that the matching contours are improved and filtered by using 
        feature descriptors.
        
        Finding the best match if any for each index1.
        If a subsequent match to the that index1's index2 is present, it will
        only be considered if the cost is less than the previous.
        It's a greedy best approach.
        After all of the index1 solutions are gathered, the best is found
        by adjusting the costs by the largest sigma from the first peak contours 
        from edges1 (i.e. applying the penalty from the paper to smaller peak
        sigma's costs).
        The best agreeing solutions are kept.
        
        If instead wanted a bipartite matching of best costs, would want to
        keep the top 2 or so of each index1.
        Would need to be careful to not include large peak contour false matches.
        All costs in the top would need to be adjusted by the penalty mentioned
        above.
        Then the best among more than one would be defined by cost.  
        Then more than one solutions points would be kept if similarity
        in difference of theta and scale were within a limit.
        */
       
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            index1BestMap = new HashMap<Integer, 
            FixedSizeSortedVector<IntensityFeatureComparisonStats>>();

        for (int idx1 = 0; idx1 < blobs1.size(); ++idx1) {

            if (contours1List.get(idx1).isEmpty()) {
                continue;
            }

            Integer index1 = Integer.valueOf(idx1);

            PairIntArray curve1 = perimeters1.get(idx1);

            Set<PairInt> blob1 = blobs1.get(idx1);

            double[] xyCen1 = curveHelper.calculateXYCentroids(blob1);

            // keeping the top '2' for each index1.  comparison is by cost.
            // choosing more than one because later bipartite matching attempts
            // to match best for all index1 matchings
            FixedSizeSortedVector<IntensityFeatureComparisonStats> bestStats 
                = new FixedSizeSortedVector<IntensityFeatureComparisonStats>(2, 
                IntensityFeatureComparisonStats.class);

            for (int idx2 = 0; idx2 < blobs2.size(); ++idx2) {

                if (contours2List.get(idx2).isEmpty()) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                PairIntArray curve2 = perimeters2.get(idx2);

                Set<PairInt> blob2 = blobs2.get(idx2);
                
                CurvatureScaleSpaceInflectionSingleEdgeMapper mapper =
                    new CurvatureScaleSpaceInflectionSingleEdgeMapper(
                    0, 0, 0, 0);

                TransformationParameters params = mapper.matchContours(
                    contours1List.get(idx1), contours2List.get(idx2));

                if ((params == null) ||
                    (mapper.getMatcher().getSolutionMatchedContours1().size() < 2)) {

                    if ((mapper.getMatcher() != null) &&
                        (mapper.getMatcher().getSolutionMatchedContours1().size() == 1)) {

                        singleSolnMap.put(new PairInt(idx1, idx2), mapper.getMatcher());
                    }

                    continue;
                }

                // edit points using feature descriptors and remove outliers:
                List<FeatureComparisonStat> compStats =
                    filterContourPointsByFeatures(img1, img2, index1, index2,
                    blob1, blob2, curve1, curve2, features1, features2,
                    mapper.getMatcher());

                if (compStats.size() < 2) {
                    continue;
                }

                //TODO: consider moving this type of statistic into the
                //cost during contour matching.  wanting to avoid accepting
                //solutions which are a small number of spurious matches due
                //to one contour having many points to match to.
                int nc1 = contours1List.get(index1.intValue()).size();
                int nc2 = contours2List.get(index2.intValue()).size();
                float frac = (float)nc1/(float)nc2;
                boolean lgDiffN = ((nc1 > nc2) && frac > 2)
                    || ((nc1 < nc2) && frac < 0.5);
                
                // also, discard if fraction of matched/maxmatchable is too low
                if (!lgDiffN) {
                    float nmm = Math.min(nc1, nc2);
                    float nm = compStats.size();
                    if ((nm/nmm) < 0.5) {
                        lgDiffN = true;
                    }
                }
                
                log.info(String.format(
                    "   nc1=%d nc2=%d (frac=%.2f) nMCs=%d",
                    nc1, nc2, frac, 
                    mapper.getMatcher().getSolutionMatchedContours1().size()));
                
                if (lgDiffN) {                    
                    log.info(
                        "discarding a good match because frac of maxMatchable is low.");
                    continue;
                }
                
                //double combinedStat = calculateCombinedIntensityStat(compStats);

                IntensityFeatureComparisonStats stats = new 
                    IntensityFeatureComparisonStats(index1.intValue(), 
                    index2.intValue(), mapper.getMatcher().getSolvedCost(), 
                        mapper.getMatcher().getSolvedScale());
                stats.addAll(compStats);
                
                // bestStats keeps the top '2' smallest cost solutions added to it
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
        
        if (index1BestMap.isEmpty()) {
            return null;
        }

        MatchingSolution soln = checkForNonDegenerateSolution(index1BestMap, binFactor1, 
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
     
        int nMaxMatchable = countMaxMatchable(contours1List, contours2List);

        Set<PairInt> present = new HashSet<PairInt>();
        List<IntensityFeatureComparisonStats> ifcsList = new ArrayList<IntensityFeatureComparisonStats>();
        List<TransformationParameters> paramsList = new ArrayList<TransformationParameters>();

        int tolTransXY = 10;

        for (Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> entry 
            : index1BestMap.entrySet()) {
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> vector = entry.getValue();
            IntensityFeatureComparisonStats[] ind1To2Pairs = vector.getArray();            
            
            for (int i = 0; i < vector.getNumberOfItems(); ++i) {
                
                IntensityFeatureComparisonStats ifcs = ind1To2Pairs[i];   
                                
                TransformationParameters params = calculateTransformation(
                    binFactor1, binFactor2, ifcs.getComparisonStats(),
                    new float[4]);
                if (params == null) {
                    continue;
                }
                int idx1 = ifcs.getIndex1();
                int idx2 = ifcs.getIndex2();
                
                PairInt p = new PairInt(idx1, idx2);
                int nEval = evaluate(params, contours1List, contours2List, tolTransXY);
                cost[idx1][idx2] = (float)nMaxMatchable/(float)nEval;
                present.add(p);
                
                ifcsList.add(ifcs);
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
        
        int n = ifcsList.size();
        int i = 0;
        while (i < n) {
            IntensityFeatureComparisonStats ifcs = ifcsList.get(i);
            PairInt p = new PairInt(ifcs.getIndex1(), ifcs.getIndex2());
            if (matched.contains(p)) {
                ++i;
                continue;
            }
            ifcsList.remove(i);
            paramsList.remove(i);
            n = ifcsList.size();
        }
        
        /*
        several ways to determine the most frequent transformation parameters.
        
        (1) (a) use my clustering code on (x,y) as (scale, rotation).
            want scale and rotation to occupy the same
            range of values, so scale should be altered by a factor and an offset
            to have the same min and max range as the rotations in degrees do.
            -- get the largest cluster and from that filter the         
               tpList and paramsList by it.
               (b) use similar methods on the filtered transX and transY.
                   offsets to make them all positive are necessary for the clustering
                   code.
                   -- get the largest cluster and from that, filter the remaining
                      tpList and paramsList.
        
        OR
        (2) use histograms to perform same functions as in (1)
        */
        
        // to correct for wrap around from 360 to 0, repeating same calc with shifterd values
        
        int[] indexesToKeep = MiscStats.filterForScaleAndRotation(paramsList, 0);
        
        int[] indexesToKeepShifted = MiscStats.filterForScaleAndRotation(paramsList, 30);
        
        if (indexesToKeepShifted.length > indexesToKeep.length) {
            indexesToKeep = indexesToKeepShifted;
        }
        
        filter(ifcsList, paramsList, indexesToKeep);
        
        indexesToKeep = MiscStats.filterForTranslation(paramsList);
        
        filter(ifcsList, paramsList, indexesToKeep);
        
        if (paramsList.size() == 0) {
            return null;
        }
        
        List<FeatureComparisonStat> combined = new ArrayList<FeatureComparisonStat>();
        for (i = 0; i < ifcsList.size(); ++i) {
            IntensityFeatureComparisonStats ifcs = ifcsList.get(i);
            combined.addAll(ifcs.getComparisonStats());
        }
        
        TransformationParameters combinedParams = calculateTransformation(
            binFactor1, binFactor2, combined, new float[4]);
        
        if (combinedParams == null) {
            return null;
        }
        
        soln = new MatchingSolution(combinedParams, combined);
            
        return soln;
    }

    protected List<FeatureComparisonStat> filterContourPointsByFeatures(
        GreyscaleImage img1, GreyscaleImage img2,
        Integer index1, Integer index2,
        Set<PairInt> blob1, Set<PairInt> blob2,
        PairIntArray curve1, PairIntArray curve2,
        IntensityFeatures features1, IntensityFeatures features2,
        CSSContourMatcherWrapper matcher) {

        FeatureMatcher featureMatcher = new FeatureMatcher();

        int dither = 1;

        List<FeatureComparisonStat> compStats = new ArrayList<FeatureComparisonStat>();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        double[] xyCen1 = curveHelper.calculateXYCentroids(curve1);
        double[] xyCen2 = curveHelper.calculateXYCentroids(curve2);

        double statSqSum = 0;

        int nMaxMatchable = matcher.getNMaxMatchable();
        int nMaxStats = 0;
        int nStats = 0;

        double cost = matcher.getSolvedCost();

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(
            "[%d](%d,%d) [%d](%d,%d) cost=%.1f scale=%.2f  nMatched=%d ",
            index1.intValue(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
            index2.intValue(), (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
            (float)cost, (float)matcher.getSolvedScale(),
            matcher.getSolutionMatchedContours1().size()));

        FeatureComparisonStat bestCompStat = null;

//TODO: refactor this to be able to reuse blob perimeter regions
//  from the invoker, the blob scale finer wrapper
        
        for (int j = 0; j < matcher.getSolutionMatchedContours1().size(); ++j) {

            CurvatureScaleSpaceContour c1 = matcher.getSolutionMatchedContours1().get(j);
            CurvatureScaleSpaceContour c2 = matcher.getSolutionMatchedContours2().get(j);

            // the sizes of the peak details will be the same
            CurvatureScaleSpaceImagePoint[] details1 = c1.getPeakDetails();
            CurvatureScaleSpaceImagePoint[] details2 = c2.getPeakDetails();

            nMaxStats += details1.length;

            for (int jj = 0; jj < details1.length; ++jj) {

                BlobPerimeterRegion region1 = 
                    BlobsAndContours.extractBlobPerimeterRegion(
                    index1.intValue(), details1[jj], curve1, blob1
                );
                region1.setIndexWithinCurve(details1[jj].getCoordIdx());

                BlobPerimeterRegion region2 = 
                    BlobsAndContours.extractBlobPerimeterRegion(
                    index2.intValue(), details2[jj], curve2, blob2
                );
                region1.setIndexWithinCurve(details2[jj].getCoordIdx());
   
                FeatureComparisonStat compStat = null;
                
                try {
                    
                    compStat =
                        featureMatcher.ditherAndRotateForBestLocation(
                        features1, features2, region1, region2, dither, img1, img2);

sb.append(
    String.format(" (%d,%d) theta1=%d   (%d,%d) theta2=%d",
    region1.getX()[1], region1.getY()[1],
    Math.round(region1.getRelativeOrientationInDegrees()),
    region2.getX()[1], region2.getY()[1],
    Math.round(region2.getRelativeOrientationInDegrees())));

                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                }
                
                if (compStat != null) {
                    float sumIntSqDiff = compStat.getSumIntensitySqDiff();
                    float intErrDiff = compStat.getImg2PointIntensityErr();
                    if (sumIntSqDiff < intErrDiff) {
                        if (bestCompStat == null) {
                            bestCompStat = compStat;
                        } else {
                            if (sumIntSqDiff < bestCompStat.getSumIntensitySqDiff()) {
                                bestCompStat = compStat;
                            }
                        }
                        statSqSum += (sumIntSqDiff*sumIntSqDiff);
                        sb.append(String.format("  %.1f(%.1f), ", sumIntSqDiff, intErrDiff));
                        compStats.add(compStat);
                        nStats++;
                    }
                }
            } // end details
        }// end matching contours for index1, index2

        if (bestCompStat == null) {
            return compStats;
        }

        log.info(sb.toString());

        // if bestCompStat's difference in orientation is different than the
        // others', re-do the others to see if have an improved calculation.
        // the "re-do" should try a dither of 1 or 2
        float bestDiffTheta = AngleUtil.getAngleDifference(
            bestCompStat.getImg2PointRotInDegrees(),
            bestCompStat.getImg1PointRotInDegrees());
        if (bestDiffTheta < 0) {
            bestDiffTheta += 360;
        }

        boolean redoStats = false;
        for (FeatureComparisonStat cStat : compStats) {
            float diffTheta = AngleUtil.getAngleDifference(
                cStat.getImg2PointRotInDegrees(),
                cStat.getImg1PointRotInDegrees());
            if (diffTheta < 0) {
                diffTheta += 360;
            }
            if (Math.abs(bestDiffTheta - diffTheta) > 25) {
                redoStats = true;
                break;
            }
        }

        //TODO: may need to consider re-doing if compStats.size() is << nMaxStats too
redoStats = true;
        if (redoStats) {

            compStats = redoFilterContourPointsByFeatures(img1, img2, index1,
                index2, blob1, blob2, curve1, curve2,
                bestCompStat.getImg1PointRotInDegrees(),
                bestCompStat.getImg2PointRotInDegrees(), matcher);
            
            log.fine("redone: " + printToString(compStats) + " combinedStat="
                + calculateCombinedIntensityStat(compStats));
       
            FeatureMatcher.removeDiscrepantThetaDiff(compStats);

            log.fine("theta diff filtered: " + printToString(compStats) + " combinedStat="
                + calculateCombinedIntensityStat(compStats));
            
            FeatureMatcher.removeIntensityOutliers(compStats);
        }

        return compStats;
    }

    List<FeatureComparisonStat> redoFilterContourPointsByFeatures(
        GreyscaleImage img1, GreyscaleImage img2,
        Integer index1, Integer index2,
        Set<PairInt> blob1, Set<PairInt> blob2,
        PairIntArray curve1, PairIntArray curve2,
        float theta1, float theta2,
        CSSContourMatcherWrapper matcher) {

        log.info("redoFilterContourPointsByFeatures");
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        // for the redo, because the orientations are set rather than found,
        // not going to re-use the invoker's instance of Features
        IntensityFeatures features1 = new IntensityFeatures(img1, 5, true);
        IntensityFeatures features2 = new IntensityFeatures(img2, 5, true);

        int dither = 2;

        List<FeatureComparisonStat> compStats = new ArrayList<FeatureComparisonStat>();

        int nMaxMatchable = matcher.getNMaxMatchable();
        int nMaxStats = 0;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("[%d] [%d] ", index1.intValue(), index2.intValue()));

        float theta1Radians = (float)(theta1 * Math.PI/180.);
        float theta2Radians = (float)(theta2 * Math.PI/180.);

        for (int j = 0; j < matcher.getSolutionMatchedContours1().size(); ++j) {

            CurvatureScaleSpaceContour c1 = matcher.getSolutionMatchedContours1().get(j);
            CurvatureScaleSpaceContour c2 = matcher.getSolutionMatchedContours2().get(j);

            // the sizes of the peak details will be the same
            CurvatureScaleSpaceImagePoint[] details1 = c1.getPeakDetails();
            CurvatureScaleSpaceImagePoint[] details2 = c2.getPeakDetails();

            nMaxStats += details1.length;

            for (int jj = 0; jj < details1.length; ++jj) {

                BlobPerimeterRegion region1 = 
                    BlobsAndContours.extractBlobPerimeterRegion(
                    index1.intValue(), details1[jj], curve1, blob1
                );
                region1.setIndexWithinCurve(details1[jj].getCoordIdx());
                region1.overrideRelativeOrientation(theta1Radians);

                BlobPerimeterRegion region2 = 
                    BlobsAndContours.extractBlobPerimeterRegion(
                    index2.intValue(), details2[jj], curve2, blob2
                );
                region2.setIndexWithinCurve(details2[jj].getCoordIdx());
                region2.overrideRelativeOrientation(theta2Radians);

                FeatureComparisonStat compStat = null;
                
                try {
                    compStat =
                        featureMatcher.ditherAndRotateForBestLocation(
                        features1, features2, region1, region2, dither);
                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                }
                
                if (compStat != null) {

                    float sumIntSqDiff = compStat.getSumIntensitySqDiff();
                    float intErrDiff = compStat.getImg2PointIntensityErr();

                    if (sumIntSqDiff < intErrDiff) {
                        compStats.add(compStat);
                    }

                }
            } // end details
        }// end matching contours for index1, index2

        if (debug) {
            sb.append(printToString(compStats)).append(" nMaxStats=").append(nMaxStats);

            log.info(sb.toString());
        }
        
        //removeOutliers(compStats);

        return compStats;
    }

    private MatchingSolution checkForNonDegenerateSolution(
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            index1BestMap, int binFactor1, int binFactor2) {
        
        List<IntensityFeatureComparisonStats> ifcsList = new ArrayList<IntensityFeatureComparisonStats>();
        
        Set<Integer> indexes2 = new HashSet<Integer>();
        
        for (Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> entry :
            index1BestMap.entrySet()) {
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> vector = entry.getValue();
            
            if (vector.getArray().length > 1) {
                return null;
            }
            assert(vector.getArray().length == 1);
            
            IntensityFeatureComparisonStats ifcs = vector.getArray()[0];
            
            Integer key2 = Integer.valueOf(ifcs.getIndex2());
            
            if (indexes2.contains(key2)) {
                return null;
            }
            
            indexes2.add(key2);
            
            ifcsList.add(ifcs);
        }
        
        // the matches are unique, so will look for the smallest cost
        // and then add similar to it
        
        double minCost = Double.MAX_VALUE;
        int minCostIdx = -1;
        for (int i = 0; i < ifcsList.size(); ++i) {
            IntensityFeatureComparisonStats tp4 = ifcsList.get(i);
            if (tp4.getCost() < minCost) {
                minCost = tp4.getCost();
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
            IntensityFeatureComparisonStats ifcs = ifcsList.get(i);
            TransformationParameters params = calculateTransformation(
               binFactor1, binFactor2, ifcs.getComparisonStats(),
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
                            combined.addAll(ifcs.getComparisonStats());
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

    private int countMaxMatchable(List<List<CurvatureScaleSpaceContour>> c1List, 
        List<List<CurvatureScaleSpaceContour>> c2List) {
        
        int n1 = 0;
        int n2 = 0;
        
        for (List<CurvatureScaleSpaceContour> list : c1List) {
            n1 += list.size();
        }
        
        for (List<CurvatureScaleSpaceContour> list : c2List) {
            n2 += list.size();
        }
        
        return Math.max(n1, n2);
    }
    
    private int evaluate(TransformationParameters params, 
        List<List<CurvatureScaleSpaceContour>> c1List, 
        List<List<CurvatureScaleSpaceContour>> c2List, int tolTransXY) {
        
        //TODO: move this method to a class for utility methods
        
        int nMatched = 0;
        
        int[] xPoints = convertToXPoints(c2List);
        int[] yPoints = convertToYPoints(c2List);
        
        Transformer transformer = new Transformer();
        
        NearestPoints np = new NearestPoints(xPoints, yPoints);
        
        for (int i = 0; i < c1List.size(); ++i) {
            
            List<CurvatureScaleSpaceContour> contours1 = c1List.get(i);
            
            for (int ii = 0; ii < contours1.size(); ++ii) {
                
                CurvatureScaleSpaceContour cr = contours1.get(ii);
                
                double[] xyTr = transformer.applyTransformation(params, 
                    cr.getPeakDetails()[0].getXCoord(), 
                    cr.getPeakDetails()[0].getYCoord());
                
                Set<Integer> indexes = np.findNeighborIndexes(
                    (int)Math.round(xyTr[0]), (int)Math.round(xyTr[1]), 
                    tolTransXY);
                
                if (indexes != null && indexes.size() > 0) {
                    nMatched++;
                }
            }
        }
        
        return nMatched;
    }

    private int[] convertToXPoints(List<List<CurvatureScaleSpaceContour>> cList) {
        
        int n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            n += list.size();
        }
        
        int[] x = new int[n];
        n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            for (CurvatureScaleSpaceContour cr : list) {
                x[n] = cr.getPeakDetails()[0].getXCoord();
                n++;
            }
        }
        
        return x;
    }
    
    private int[] convertToYPoints(List<List<CurvatureScaleSpaceContour>> cList) {
        
        int n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            n += list.size();
        }
        
        int[] y = new int[n];
        n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            for (CurvatureScaleSpaceContour cr : list) {
                y[n] = cr.getPeakDetails()[0].getYCoord();
                n++;
            }
        }
        
        return y;
    }

    private void filter(List<IntensityFeatureComparisonStats> ifcsList, 
        List<TransformationParameters> paramsList, int[] indexesToKeep) {
        
        if (indexesToKeep.length < 2) {
            return;
        }
        
        List<IntensityFeatureComparisonStats> ifcsList2 = 
            new ArrayList<IntensityFeatureComparisonStats>(indexesToKeep.length);
        
        List<TransformationParameters> paramsList2 = 
            new ArrayList<TransformationParameters>();
        
        for (int i = 0; i < indexesToKeep.length;++i) {
            int idx = indexesToKeep[i];
            ifcsList2.add(ifcsList.get(idx));
            paramsList2.add(paramsList.get(idx));
        }
        
        ifcsList.clear();
        ifcsList.addAll(ifcsList2);
        
        paramsList.clear();
        paramsList.addAll(paramsList2);
    }
}
