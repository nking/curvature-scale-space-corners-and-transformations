package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.util.AngleUtil;
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

        Map<PairInt, CSSContourMatcherWrapper> singleSolnMap =
            new HashMap<PairInt,  CSSContourMatcherWrapper>();

        Map<Integer, TransformationParameters> pMap = new HashMap<Integer,
            TransformationParameters>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

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

            if (debug) {
                for (int k = 0; k < bestStats.getNumberOfItems(); ++k) {
                    IntensityFeatureComparisonStats stats = bestStats.getArray()[k];
                    double[] xyCen2 = curveHelper.calculateXYCentroids(
                        blobs2.get(stats.getIndex2()));
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format(
                        "==>[%d](%d,%d) [%d](%d,%d) cost=%.2f scale=%.2f nMatched=(%d,%d) intSqDiff=%.1f",
                        stats.getIndex1(), (int)Math.round(xyCen1[0]), (int)Math.round(xyCen1[1]),
                        stats.getIndex2(), (int)Math.round(xyCen2[0]), (int)Math.round(xyCen2[1]),
                        (float)stats.getCost(), (float)stats.getScale(), 
                        stats.getComparisonStats().size(),
                        contours1List.get(stats.getIndex1()).size(), 
                        (float)calculateCombinedIntensityStat(stats.getComparisonStats())));
                    log.info(sb.toString());
                }
            }

            index1BestMap.put(index1, bestStats);
            
            IntensityFeatureComparisonStats bestMatches = bestStats.getArray()[0];
            
            float[] scaleRotTransXYStDev00 = new float[4];
            TransformationParameters params = calculateTransformation(
                binFactor1, binFactor2, bestMatches.getComparisonStats(), 
                scaleRotTransXYStDev00);
            
            if (params == null) {
                continue;
            }
            
            if ((bestMatches.getComparisonStats().size() > 3) 
                && (idx1 < (blobs1.size() - 1))) {
                
                if (stDevsAreSmall(params, scaleRotTransXYStDev00)) {
                    
                    double c = calculateCombinedIntensityStat(
                        bestMatches.getComparisonStats());
                    log.info("MATCHED EARLY: combined compStat=" + c);
                        
                    MatchingSolution soln = new MatchingSolution(params,
                        bestMatches.getComparisonStats());
                        
                    return soln;
                }
            }
            
            pMap.put(index1, params);
        }
        
        List<FeatureComparisonStat> bestOverall = null;
        if (index1BestMap.isEmpty()) {
            
            // -------- process the single solution compStats ------------
            
            if (singleSolnMap.size() > 1) {
                bestOverall = processSingleSolutionsIfNoBest(img1, img2,
                    singleSolnMap, blobs1, blobs2, perimeters1, perimeters2,
                    features1, features2);
            }

            if (index1BestMap.isEmpty()) {
                return null;
            }
        } else {
            bestOverall = filterToBestConsistent(index1BestMap, contours1List,
                contours2List);
        }

        if (bestOverall == null) {
            return null;
        }

        TransformationParameters params = calculateTransformation(
            img1Helper.imgHelper.getBinFactor(useBinned1),
            img2Helper.imgHelper.getBinFactor(useBinned2),
            bestOverall, new float[4]);

        if (params == null) {
            return null;
        }

        MatchingSolution soln = new MatchingSolution(params, bestOverall);
            
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

        for (int j = 0; j < matcher.getSolutionMatchedContours1().size(); ++j) {

            CurvatureScaleSpaceContour c1 = matcher.getSolutionMatchedContours1().get(j);
            CurvatureScaleSpaceContour c2 = matcher.getSolutionMatchedContours2().get(j);

            // the sizes of the peak details will be the same
            CurvatureScaleSpaceImagePoint[] details1 = c1.getPeakDetails();
            CurvatureScaleSpaceImagePoint[] details2 = c2.getPeakDetails();

            nMaxStats += details1.length;

            for (int jj = 0; jj < details1.length; ++jj) {

                BlobPerimeterRegion region1 = extractBlobPerimeterRegion(
                    index1.intValue(), details1[jj], curve1, blob1
                );

                BlobPerimeterRegion region2 = extractBlobPerimeterRegion(
                    index2.intValue(), details2[jj], curve2, blob2
                );
   
                FeatureComparisonStat compStat =
                    featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, region1, region2, dither, img1, img2);

sb.append(
    String.format(" (%d,%d) theta1=%d   (%d,%d) theta2=%d",
    region1.getX(), region1.getY(),
    Math.round(region1.getRelativeOrientationInDegrees()),
    region2.getX(), region2.getY(),
    Math.round(region2.getRelativeOrientationInDegrees())));

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

                BlobPerimeterRegion region1 = extractBlobPerimeterRegion(
                    index1.intValue(), details1[jj], curve1, blob1
                );
                region1.overrideRelativeOrientation(theta1Radians);

                BlobPerimeterRegion region2 = extractBlobPerimeterRegion(
                    index2.intValue(), details2[jj], curve2, blob2
                );
                region2.overrideRelativeOrientation(theta2Radians);

                FeatureComparisonStat compStat =
                    featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, region1, region2, dither);

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

    protected List<FeatureComparisonStat> filterToBestConsistent(
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
        index1StatsMap, List<List<CurvatureScaleSpaceContour>> contours1Lists,
        List<List<CurvatureScaleSpaceContour>> contours2Lists) {

        if (index1StatsMap == null || index1StatsMap.isEmpty()) {
            return null;
        }
        
        if (index1StatsMap.size() > 1) {
            // make corrections for cost between different edges            
            correctCostsUsingMaxSigma(index1StatsMap, contours1Lists);            
        }
       
        TreeMap<Double, List<IntensityFeatureComparisonStats>> bestMatches = 
            findBestMatchesUsingBipartite(index1StatsMap, 
            contours1Lists.size(), contours2Lists.size());

        int nMaxMatchable = 2 * Math.max(contours1Lists.size(), contours2Lists.size());
     
        Map<Integer, IntensityFeatureComparisonStats> index1Map = new
            HashMap<Integer, IntensityFeatureComparisonStats>();
            
        Map<PairInt, PairFloat> indexesDiffTheta = new HashMap<PairInt, PairFloat>();
        
        /* calculate the highest number of similar transformations and the
        lowest cost from those.  store nSimilar, indexes, cost for each iteration*/
        int[] nSimilarSummary = new int[nMaxMatchable];
        Integer[][] indexesSummary = new Integer[nMaxMatchable][];
        double[] costsSummary = new double[nMaxMatchable];
        int[] mainIndexSummary = new int[nMaxMatchable];

        int count = 0;
        
        for (Map.Entry<Double, List<IntensityFeatureComparisonStats>> entry : bestMatches.entrySet()) {

            Double adjustedCost = entry.getKey();
            
            List<IntensityFeatureComparisonStats> stats = entry.getValue();
            
            for (IntensityFeatureComparisonStats stat : stats) {
                int idx1 = stat.getIndex1();
                int idx2 = stat.getIndex2();
                assert(index1Map.get(Integer.valueOf(idx1)) == null);
                index1Map.put(Integer.valueOf(idx1), stat);
                PairInt p = new PairInt(idx1, idx2);
                PairFloat diffThetaMnStdv = indexesDiffTheta.get(p);
                if (diffThetaMnStdv == null) {
                    float[] diffThetas = new float[stat.getComparisonStats().size()];
                    for (int i = 0; i < stat.getComparisonStats().size(); ++i) {
                        FeatureComparisonStat fcs = stat.getComparisonStats().get(i);                
                        float diff = AngleUtil.getAngleDifference(
                           fcs.getImg1PointRotInDegrees(), fcs.getImg2PointRotInDegrees());
                        diffThetas[i] = diff;
                    }
                    float[] msv = MiscMath.getAvgAndStDev(diffThetas);
                    
                    // TODO: consider whether need to exclude points further from
                    // the mean than a difference in thetaDiff of 20
                    // unless edits have changed it, this was run previously,
                    // so what is present here should already be consistent
                    // diffThetas:
                    //removeDiscrepantThetaDiff(stat.getComparisonStats());
                    
                    diffThetaMnStdv = new PairFloat(msv[0], msv[1]);
                }
                
                double scale = stat.getScale();
                
                Set<Integer> similarParamsIndexes1 = new HashSet<Integer>();
                similarParamsIndexes1.add(Integer.valueOf(idx1));
                
                // count the number of solutions similar in diffTheta and scale
                for (Map.Entry<Double, List<IntensityFeatureComparisonStats>> entry2 : bestMatches.entrySet()) {
                    Double adjustedCost2 = entry2.getKey();
                    if (adjustedCost2.equals(adjustedCost)) {
                        continue;
                    }
                    for (IntensityFeatureComparisonStats stat2 : entry2.getValue()) {
                        int idx1P = stat2.getIndex1();
                        int idx2P = stat2.getIndex2();
                        // bipartite matching should have made unique idx1 already:
                        assert(idx1 != idx1P);
                        PairInt p2 = new PairInt(idx1P, idx2P);
                        PairFloat diffThetaMnStdv2 = indexesDiffTheta.get(p2);
                        if (diffThetaMnStdv2 == null) {
                            float[] diffThetas = new float[stat2.getComparisonStats().size()];
                            for (int i = 0; i < stat2.getComparisonStats().size(); ++i) {
                                FeatureComparisonStat fcs = stat2.getComparisonStats().get(i);                
                                float diff = AngleUtil.getAngleDifference(
                                    fcs.getImg1PointRotInDegrees(), fcs.getImg2PointRotInDegrees());
                                diffThetas[i] = diff;
                            }
                            float[] msv = MiscMath.getAvgAndStDev(diffThetas);
                            diffThetaMnStdv2 = new PairFloat(msv[0], msv[1]);
                        }
                        
                        //--- compare diffTheta and scale with ---
                        if ((Math.abs(diffThetaMnStdv.getX() - diffThetaMnStdv2.getX()) < 10)
                            && (Math.abs(scale - stat2.getScale()) < 0.1)) {
                            similarParamsIndexes1.add(Integer.valueOf(idx1P));
                        }
                    }
                }
                
                // store for sort and combine
                nSimilarSummary[count] = similarParamsIndexes1.size();
                indexesSummary[count] = similarParamsIndexes1.toArray(new Integer[similarParamsIndexes1.size()]);
                costsSummary[count] = adjustedCost;
                mainIndexSummary[count] = idx1;
                
                count++;
            }
        }
        
        nSimilarSummary = Arrays.copyOf(nSimilarSummary, count);
        indexesSummary = Arrays.copyOf(indexesSummary, count);
        costsSummary = Arrays.copyOf(costsSummary, count);
        mainIndexSummary = Arrays.copyOf(mainIndexSummary, count);
        
        // sort to prefer the solution w/ largest number of similar solutions:
        //MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary, indexesSummary, mainIndexSummary);

        // OR sort to refer the solution w/ best cost and any it is similar to:
        //--- these are still sorted by costs already, so no need to resort again ---
        // MultiArrayMergeSort.sortBy1stAscThen2ndDesc(costsSummary, nSimilarSummary, indexesSummary, mainIndexSummary, 0, costsSummary.length - 1);

        Integer[] indexes = indexesSummary[0];
       
        List<FeatureComparisonStat> output = new ArrayList<FeatureComparisonStat>();
        for (Integer index1 : indexes) {
            IntensityFeatureComparisonStats stats = index1Map.get(index1);
            output.addAll(stats.getComparisonStats());
        }
        return output;
    }
    
    private List<FeatureComparisonStat> processSingleSolutionsIfNoBest(
        GreyscaleImage img1, GreyscaleImage img2,
        Map<PairInt, CSSContourMatcherWrapper> singleSolnMap,
        List<Set<PairInt>> blobs1, List<Set<PairInt>> blobs2,
        List<PairIntArray> perimeters1, List<PairIntArray> perimeters2,
        IntensityFeatures features1, IntensityFeatures features2) {

log.info("WARNING: processing single solutions... may remove these in future");

        List<FeatureComparisonStat> csList = new ArrayList<FeatureComparisonStat>();

        int nm = singleSolnMap.size();
        
        double[] costs = new double[nm];
        int[] idx1s = new int[nm];
        int[] idx2s = new int[nm];
        FeatureComparisonStat[] ifs = new FeatureComparisonStat[nm];
        
        int count = 0;
        for (Entry<PairInt, CSSContourMatcherWrapper> entry : singleSolnMap.entrySet()) {
            
            PairInt index1Index2 = entry.getKey();
            int idx1 = index1Index2.getX();
            int idx2 = index1Index2.getY();
            
            CSSContourMatcherWrapper matcher = entry.getValue();
            
            PairIntArray curve1 = perimeters1.get(idx1);
            Set<PairInt> blob1 = blobs1.get(idx1);
            PairIntArray curve2 = perimeters2.get(idx2);
            Set<PairInt> blob2 = blobs2.get(idx2);

            List<FeatureComparisonStat> compStats =
                filterContourPointsByFeatures(img1, img2,
                Integer.valueOf(idx1), Integer.valueOf(idx2),
                blob1, blob2, curve1, curve2, features1, features2,
                matcher);

            if (compStats.isEmpty()) {
                continue;
            }
            
            assert(compStats.size() == 1);

            double combStat = calculateCombinedIntensityStat(compStats);

            idx1s[count] = idx1;
            idx2s[count] = idx2;
            ifs[count] = compStats.get(0);

            //ssd of intensity is a better selector for one dataset. this may change w/ more testing
            //costs[count] = stats.getAdjustedCost();
            costs[count] = combStat;
            
            count++;
        }
        if (count == 0) {
            return null;
        }
        costs = Arrays.copyOf(costs, count);
        idx1s = Arrays.copyOf(idx1s, count);
        idx2s = Arrays.copyOf(idx2s, count);
        ifs = Arrays.copyOf(ifs, count);
        int[] lookupIndexes = new int[count];
        for (int i = 0; i < count; ++i) {
            lookupIndexes[i] = i;
        }
        // best will be at bottom of list:
        MultiArrayMergeSort.sortByDecr(costs, lookupIndexes);
        
        float ssdBest = ifs[lookupIndexes[count - 1]].getSumIntensitySqDiff();
        
        float thetaDiffBest = AngleUtil.getAngleDifference(
            ifs[lookupIndexes[count - 1]].getImg1PointRotInDegrees(),
            ifs[lookupIndexes[count - 1]].getImg2PointRotInDegrees());
        
        int count2 = 0;
        double[] costs2 = new double[count];
        int[] idx1s2 = new int[count];
        int[] idx2s2 = new int[count];
        FeatureComparisonStat[] ifs2 = new FeatureComparisonStat[count];
        for (int i = 0; i < count; ++i) {
            int idx0 = lookupIndexes[i];
            FeatureComparisonStat fcs = ifs[idx0];
            double ssd = fcs.getSumIntensitySqDiff();

            float dtm = AngleUtil.getAngleDifference(
                fcs.getImg1PointRotInDegrees(), fcs.getImg2PointRotInDegrees());

            if (Math.abs(dtm - thetaDiffBest) > 20) {
                continue;
            }
            
            if (ssd > (3 * ssdBest)) {
                continue;
            }

            costs2[count2] = costs[i];
            idx1s2[count2] = idx1s[idx0];
            idx2s2[count2] = idx2s[idx0];
            ifs2[count2] = fcs;
            
            count2++;
        }
        costs2 = Arrays.copyOf(costs2, count2);
        idx1s2 = Arrays.copyOf(idx1s2, count2);
        idx2s2 = Arrays.copyOf(idx2s2, count2);
        ifs2 = Arrays.copyOf(ifs2, count2);
        
        Set<Integer> chosen1 = new HashSet<Integer>();
        Set<Integer> chosen2 = new HashSet<Integer>();

        for (int i = 0; i < count2; ++i) {

            Integer index1 = Integer.valueOf(idx1s2[i]);
            Integer index2 = Integer.valueOf(idx2s2[i]);

            if (chosen1.contains(index1)) {
                continue;
            }
            if (chosen2.contains(index2)) {
                continue;
            }
            
            Double cost = Double.valueOf(costs2[i]);
            FeatureComparisonStat ifcs = ifs2[i];

            // only true if still using costs rather than SSD:
            //assert(Math.abs(cost.doubleValue() - ifcs.getAdjustedCost()) < 0.01);

            csList.add(ifcs);

            chosen1.add(index1);
            chosen2.add(index2);
        }
        
        return csList;
    }

    private void correctCostsUsingMaxSigma(
        Map<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>> 
            index1StatsMap, List<List<CurvatureScaleSpaceContour>> contours1Lists) {
                
        float maxPeakSigma = Float.MIN_VALUE;
        
        Map<Integer, Float> index1PeakSigma = new HashMap<Integer, Float>();
        
        for (Entry<Integer, FixedSizeSortedVector<
            IntensityFeatureComparisonStats>> entry : index1StatsMap.entrySet()) {
            
            Integer index1 = entry.getKey();
            
            List<CurvatureScaleSpaceContour> contours = contours1Lists.get(index1.intValue());
           
            // contours are ordered, so just need the first.
            if (contours.isEmpty()) {
                continue;
            }
            float peakSigma = contours.get(0).getPeakSigma();
            if (peakSigma > maxPeakSigma) {
                maxPeakSigma = peakSigma;
            }
            // contour 6 is not a strong blob by has maxsigma=166 which is >4*expected from 0
            index1PeakSigma.put(index1, Float.valueOf(peakSigma));
        }
        
        for (Entry<Integer, Float> entry : index1PeakSigma.entrySet()) {
            
            Integer index1 = entry.getKey();
            double penalty = maxPeakSigma - entry.getValue().floatValue();
            
            FixedSizeSortedVector<IntensityFeatureComparisonStats> index1Stats =
                index1StatsMap.get(index1);
            
            for (int i = 0; i < index1Stats.getNumberOfItems(); ++i) {
                
                IntensityFeatureComparisonStats stats = index1Stats.getArray()[i];
                                
                stats.setAdjustedCost(stats.getCost() + penalty);
            }
        }
    }

    /**
     * bipartite matching after sorting by largest number of matches and lowest
     * cost then greedily choosing from that order.
     * @param index1StatsMap
     * @param n1
     * @param n2
     * @return 
     */
    private TreeMap<Double, List<IntensityFeatureComparisonStats>> 
    findBestMatchesUsingBipartite(Map<Integer, 
        FixedSizeSortedVector<IntensityFeatureComparisonStats>> index1StatsMap,
        int n1, int n2) {
        
        TreeMap<Double, List<IntensityFeatureComparisonStats>> matched = new
            TreeMap<Double, List<IntensityFeatureComparisonStats>>();
        
        int maxMatchable = 2 * Math.max(n1, n2);
                
        double[] costs = new double[maxMatchable];
        int[] idx1s = new int[maxMatchable];
        int[] idx2s = new int[maxMatchable];
        int[] nMatches = new int[maxMatchable];
        IntensityFeatureComparisonStats[] ics = new IntensityFeatureComparisonStats[maxMatchable];
        
        int count = 0;
        for (Entry<Integer, FixedSizeSortedVector<IntensityFeatureComparisonStats>>
            entry : index1StatsMap.entrySet()) {
            Integer index1 = entry.getKey();            
            for (IntensityFeatureComparisonStats stats : entry.getValue().getArray()) {
                if (stats == null) {
                    continue;
                }
                assert(index1.intValue() == stats.getIndex1());
                
                idx1s[count] = stats.getIndex1();
                idx2s[count] = stats.getIndex2();
                ics[count] = stats;
                nMatches[count] = stats.getComparisonStats().size();
                
                //ssd of intensity is a better selector for one dataset. this may change w/ more testing
                //costs[count] = stats.getAdjustedCost();
                costs[count] = calculateCombinedIntensityStat(stats.getComparisonStats());
                
                count++;
            }
        }
        costs = Arrays.copyOf(costs, count);
        idx1s = Arrays.copyOf(idx1s, count);
        idx2s = Arrays.copyOf(idx2s, count);
        ics = Arrays.copyOf(ics, count);
        nMatches = Arrays.copyOf(nMatches, count);
        int[] lookupIndexes = new int[count];
        for (int i = 0; i < count; ++i) {
            lookupIndexes[i] = i;
        }
        
        //sort for the highest number of matches having the lowest costs.
        // decr nMatches, asc costs
        MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nMatches, costs, lookupIndexes);
        
        /*
        adding another filter based upon the SSD of intensity features for the
        top number of matches and lowest cost.
        This uses the average and standard deviation from it if there are
        more than 2 SSD stats and removes all other matches in the list
        where the SSD is larger than 2 sigma or so from that difference.
        It's not the most stable filter considering that some regions surrounding
        a contour may have extremely high variability, but it should usually be
        better to use this filter.
        
        Also, looks like it's necessary to filter for the difference in theta
        when it's much larger than the top difference in theta.
        */
        float[] ssdMeanStDv = FeatureMatcher.calcIntensitySSDMeanAndStDev(
            ics[lookupIndexes[0]].getComparisonStats());
        
        float diffInThetaMean = calculateDiffThetaMean(ics[lookupIndexes[0]].getComparisonStats());
        
        int count2 = 0;
        double[] costs2 = new double[count];
        int[] idx1s2 = new int[count];
        int[] idx2s2 = new int[count];
        int[] nMatches2 = new int[count];
        IntensityFeatureComparisonStats[] ics2 = new IntensityFeatureComparisonStats[count];
        for (int i = 0; i < count; ++i) {
            int idx0 = lookupIndexes[i];
            IntensityFeatureComparisonStats ifcs = ics[idx0];
            double ssd = calculateCombinedIntensityStat(ifcs.getComparisonStats());
            
            if (ssd > (ssdMeanStDv[0] + (2 * ssdMeanStDv[1]))) {
                continue;
            }
                            
            float dtm = calculateDiffThetaMean(ifcs.getComparisonStats());
                
            if (Math.abs(dtm - diffInThetaMean) > 20) {
                continue;
            }
                
            costs2[count2] = costs[i];
            idx1s2[count2] = idx1s[idx0];
            idx2s2[count2] = idx2s[idx0];
            ics2[count2] = ics[idx0];
            nMatches2[count2] = nMatches[i];
            count2++;
        }
        costs2 = Arrays.copyOf(costs2, count2);
        idx1s2 = Arrays.copyOf(idx1s2, count2);
        idx2s2 = Arrays.copyOf(idx2s2, count2);
        ics2 = Arrays.copyOf(ics2, count2);
        nMatches2 = Arrays.copyOf(nMatches2, count2);
        
        Set<Integer> chosen1 = new HashSet<Integer>();
        Set<Integer> chosen2 = new HashSet<Integer>();
        
        for (int i = 0; i < count2; ++i) {
                        
            Integer index1 = Integer.valueOf(idx1s2[i]);
            Integer index2 = Integer.valueOf(idx2s2[i]);
            
            if (chosen1.contains(index1)) {
                continue;
            }
            if (chosen2.contains(index2)) {
                continue;
            }
            
            Double cost = Double.valueOf(costs2[i]);
            IntensityFeatureComparisonStats ifcs = ics2[i];
              
            // only true if still using costs rather than SSD:
            //assert(Math.abs(cost.doubleValue() - ifcs.getAdjustedCost()) < 0.01);
            
            /*
            TreeMap<Double, List<IntensityFeatureComparisonStats>> matched
            */
            List<IntensityFeatureComparisonStats> ifcsList = matched.get(cost);
            if (ifcsList == null) {
                ifcsList = new ArrayList<IntensityFeatureComparisonStats>();
                matched.put(cost, ifcsList);
            }
            ifcsList.add(ifcs);
                
            chosen1.add(index1);
            chosen2.add(index2);
        }
        
        return matched;
    }

}
