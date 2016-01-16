package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.SubsetChooser;
import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * given a map of lists of blob matched points, calculate euclidean
 * transformations and evaluate them against all points.
 *
 * @author nichole
 */
public class BlobCornersEuclideanCalculator2 {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * returns a solution for given images and their binFactors.  Note that
     * corrections to make the solution full frame transformation have not
     * occurred, the solution is in the given binned reference frames.
     * 
     * @param image1
     * @param image2
     * @param features1
     * @param features2
     * @param dither
     * @param matchedLists a map with keys=pairint of index1, index2 of matched
     * blob lists and values=lists of features for that index1 index2 match.
     * @param allCorners1
     * @param allCorners2
     * @param binFactor1
     * @param binFactor2
     * @return 
     */
    public MatchingSolution solveTransformation(
        GreyscaleImage image1, GreyscaleImage image2,
        IntensityFeatures features1, IntensityFeatures features2, int dither,
        Map<PairInt, List<FeatureComparisonStat>> matchedLists, 
        List<CornerRegion> allCorners1, List<CornerRegion> allCorners2,
        int binFactor1, int binFactor2) {
        
        // calculate all transformations, returns transformations from all
        // combinations of 2 if enough points per blob, else combinations of 3
        List<TransformationParameters> params = calculateTransformations(
            matchedLists);
        
        params = MiscStats.filterToSimilarParamSets2(params, binFactor1, 
            binFactor2);
        
        boolean evalWithAllCorners = true;
        
        if (evalWithAllCorners) {
            
            return evaluateWithAllCorners(image1, image2, features1, features2, 
                dither, params, allCorners1, allCorners2, 
                binFactor1, binFactor2);
            
        } else {
            
            return evaluateWithMatchedLists(image1, image2, features1, 
                features2, dither, params, matchedLists, binFactor1, binFactor2);
        }
    }
    
    protected MatchingSolution evaluateWithMatchedLists(
        GreyscaleImage image1, GreyscaleImage image2,
        IntensityFeatures features1, IntensityFeatures features2, int dither,
        List<TransformationParameters> paramsList, 
        Map<PairInt, List<FeatureComparisonStat>> matchedLists,
        int binFactor1, int binFactor2) {
     
        throw new UnsupportedOperationException("not yet implemented");
    }

    private List<TransformationParameters> calculateTransformations(
        Map<PairInt, List<FeatureComparisonStat>> matchedLists) {
        
        /*
        ways to calculate:
        
        (1) if there are usually at least two points per matchedLists values,
            then can use the points within to calculate transformations.
            use combinations of each 2 keys in matchedLists to create
            transformations.
        (2) if there is only one point per values of matchedLists, then can
            use combinations of 3 keys in matchedLists to create transformations.        
        */
        
        List<List<FeatureComparisonStat>> listsOfMatchedBlobStats = 
            new ArrayList<List<FeatureComparisonStat>>();
        for (Entry<PairInt, List<FeatureComparisonStat>> entry : matchedLists.entrySet()) {
            listsOfMatchedBlobStats.add(entry.getValue());
        }

        List<PairInt> descendingNFreq = new ArrayList<PairInt>();
        List<PairInt> descendingNMembers = new ArrayList<PairInt>();
        countFrequencyThenSortDesc(listsOfMatchedBlobStats, descendingNFreq,
            descendingNMembers);
        
        List<TransformationParameters> parametersList = 
            new ArrayList<TransformationParameters>();
        
        // if there are 2 or more blobs with 3 or more members, make a sublist
        // and calculate transformations for those.  
        // testing whether can use only those instead of all combinations below it
        List<List<FeatureComparisonStat>> matchedBlobStatsSublist = 
            new ArrayList<List<FeatureComparisonStat>>();
        for (Entry<PairInt, List<FeatureComparisonStat>> entry : matchedLists.entrySet()) {
            List<FeatureComparisonStat> list = entry.getValue();
            if (list.size() > 2) {
                matchedBlobStatsSublist.add(list);
            }
        }
        if (matchedBlobStatsSublist.size() == 1) {
            TransformationParameters params = MiscStats.calculateTransformation(
                1, 1, matchedBlobStatsSublist.get(0), new float[4], true);
            if (params != null) {
                parametersList.add(params);
                log.info("single high matching params=" + params.toString());
            }
        } else if (matchedBlobStatsSublist.size() > 1) {
            calculateWithCombinations(matchedBlobStatsSublist, 2, parametersList);
            for (TransformationParameters params : parametersList) {
                log.info("blob high n matches params=" + params.toString());
            }
        }
        
        if (parametersList.size() < 10) {
            int k = 3;
            if ((descendingNFreq.get(0).getX() > 1) && 
                (descendingNFreq.get(0).getY() > (matchedLists.size()/2))) {
                k = 2;
            }

            calculateWithCombinations(listsOfMatchedBlobStats, k, parametersList);
        }
        
        return parametersList;
    }

    /**
     * count the number of keys have size of values and return in descending
     * order, pairints holding the sizes of the values in matchedList and the
     * numbers of keys for those sizes.
     * runtime complexity is at most O(N) where N is matchedLists.size().
     * @param listsOfMatchedBlobStats
     * @param outputDescendingNFreq - list of pairints with 
     *     x=size of lists of stats, y=number of blobs with that size,
     *     sorted by decreasing value of y
     * @param outputDescendingNMembers - list of pairints with 
     *     x=size of lists of stats, y=number of blobs with that size,
     *     sorted by decreasing value of x
     */
    protected void countFrequencyThenSortDesc(
        List<List<FeatureComparisonStat>> listsOfMatchedBlobStats,
        List<PairInt> outputDescendingNFreq, List<PairInt> outputDescendingNMembers) {
        
        Map<Integer, Integer> valueCounts = new HashMap<Integer, Integer>();
        
        for (List<FeatureComparisonStat> list : listsOfMatchedBlobStats) {
            
            int n = list.size();
            
            Integer key = Integer.valueOf(n);
            
            Integer c = valueCounts.get(key);
            
            if (c == null) {
                valueCounts.put(key, Integer.valueOf(1));
            } else {
                valueCounts.put(key, Integer.valueOf(c.intValue() + 1));
            }
        }
        
        int[] v = new int[valueCounts.size()];
        int[] c = new int[v.length];
        int count = 0;
        int maxC = Integer.MIN_VALUE;
        int maxV = Integer.MIN_VALUE;
        for (Entry<Integer, Integer> entry : valueCounts.entrySet()) {
            v[count] = entry.getKey().intValue();
            c[count] = entry.getValue().intValue();
            if (c[count] > maxC) {
                maxC = c[count];
            }
            if (v[count] > maxV) {
                maxV = v[count];
            }
            count++;
        }
        
        CountingSort.sortByDecr(c, v, maxC);
        
        for (int i = 0; i < c.length; ++i) {
            outputDescendingNFreq.add(new PairInt(v[i], c[i]));
        }
        
        CountingSort.sortByDecr(v, c, maxV);
        for (int i = 0; i < v.length; ++i) {
            outputDescendingNMembers.add(new PairInt(v[i], c[i]));
        }
    }

    private List<TransformationParameters> calculateWithCombinations(
        List<List<FeatureComparisonStat>> listsOfMatchedBlobStats, 
        int k, List<TransformationParameters> output) {
                
        boolean removeIntensityOutliers = false;
        float[] outputScaleRotTransXYStDev = new float[4];
        
        int nPoints = listsOfMatchedBlobStats.size();
        
        int[] selectedIndexes = new int[k];
        
        SubsetChooser subsetChooser = new SubsetChooser(nPoints, k);
            
        int nV = subsetChooser.getNextSubset(selectedIndexes);

        while (nV != -1) {
            
            List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();

            for (int bitIndex : selectedIndexes) {

                int idx = bitIndex;

                List<FeatureComparisonStat> list = listsOfMatchedBlobStats.get(idx);

                stats.addAll(list);                
            }
            
            TransformationParameters params = MiscStats.calculateTransformation(
                1, 1, stats, outputScaleRotTransXYStDev, 
                removeIntensityOutliers);
            
            if (params != null) {
                output.add(params);
            }

            nV = subsetChooser.getNextSubset(selectedIndexes);
        }
        
        return output;
    }
      
    protected MatchingSolution evaluateWithAllCorners(
        GreyscaleImage image1, GreyscaleImage image2,
        IntensityFeatures features1, IntensityFeatures features2, int dither,
        List<TransformationParameters> parameterList, 
        List<CornerRegion> allCorners1, List<CornerRegion> allCorners2,
        int binFactor1, int binFactor2) {
        
        int n2c = allCorners2.size();
        int[] xC2 = new int[n2c];
        int[] yC2 = new int[n2c];
        for (int i = 0; i < n2c; ++i) {
            CornerRegion cr2 = allCorners2.get(i);
            xC2[i] = cr2.getX()[cr2.getKMaxIdx()];
            yC2[i] = cr2.getY()[cr2.getKMaxIdx()];
        }
               
        NearestPoints np2 = new NearestPoints(xC2, yC2);
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        Transformer transformer = new Transformer();

        final int rotationTolerance = 20;
        int tolTransXY = 5;

        TransformationParameters bestParams = null;
        float bestCost = Float.MAX_VALUE;
        float bestCost1Norm = Float.MAX_VALUE;
        List<FeatureComparisonStat> bestStats = null;
        int bestTolTransXY2 = -1;

        for (TransformationParameters params : parameterList) {

            if (!paramsAreValid(params)) {
                continue;
            }
            
            double rotInRadians = params.getRotationInRadians();
            double cos = Math.cos(rotInRadians);
            double sin = Math.sin(rotInRadians);

            int rotD = Math.round(params.getRotationInDegrees());

            int tolTransXY2 = tolTransXY;
            if (params.getScale() < 1) {
                tolTransXY2 = Math.round(tolTransXY * params.getScale());
            }
            if (tolTransXY2 == 0) {
                tolTransXY2 = 1;
            }
            int dither2 = Math.round(dither * params.getScale());
            if (dither2 == 0) {
                dither2 = 1;
            } else if (dither2 > dither) {
                // large dither makes runtime larger
                dither2 = dither;
            }

            int nEval = 0;
            double sumSSD = 0;
            double sumDist = 0;
            List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
            List<Double> distances = new ArrayList<Double>();

            for (CornerRegion cr : allCorners1) {

                CornerRegion crTr = transformer.applyTransformation(params, cr, cos, sin);

                Set<Integer> indexes2 = np2.findNeighborIndexes(
                    crTr.getX()[crTr.getKMaxIdx()],
                    crTr.getY()[crTr.getKMaxIdx()], tolTransXY2);

                double bestCostPerIndex = Double.MAX_VALUE;
                Integer bestCostPerIndexIndex = null;
                FeatureComparisonStat bestCostPerIndexStat = null;
                double bestCostPerIndexDist = Double.MAX_VALUE;
                    
                for (Integer index : indexes2) {

                    int idx2 = index.intValue();

                    CornerRegion corner2 = allCorners2.get(idx2);

                    FeatureComparisonStat compStat =
                        featureMatcher.ditherAndRotateForBestLocation2(
                        features1, features2, cr, corner2, dither2,
                        rotD, rotationTolerance, image1, image2);

                    if ((compStat == null) ||
                        (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                        ) {
                        continue;
                    }

                    double xTr = (compStat.getImg1Point().getX() *
                        params.getScale() * cos) +
                        (compStat.getImg1Point().getY() *
                        params.getScale() * sin);
                    xTr += params.getTranslationX();

                    double yTr = (-compStat.getImg1Point().getX() *
                        params.getScale() * sin) +
                        (compStat.getImg1Point().getY() *
                        params.getScale()* cos);
                    yTr += params.getTranslationY();

                    double dist = distance(xTr, yTr,
                        compStat.getImg2Point().getX(),
                        compStat.getImg2Point().getY());

                    double cost = 
                        (((float)dist + 0.01f)/(float)tolTransXY2) *
                        (compStat.getSumIntensitySqDiff() + 1);

                    if (cost < bestCostPerIndex) {
                        bestCostPerIndex = cost;
                        bestCostPerIndexIndex = index;
                        bestCostPerIndexStat = compStat;
                        bestCostPerIndexDist = dist;
                    }
                }

                if (bestCostPerIndexIndex != null) {
                    // cost is 
                    stats.add(bestCostPerIndexStat);
                    distances.add(Double.valueOf(bestCostPerIndexDist));
                    sumDist += bestCostPerIndexDist;
                    sumSSD += bestCostPerIndexStat.getSumIntensitySqDiff();
                    nEval++;
                }
            }
            
            List<Integer> removedIndexes = MiscStats.filterForDegeneracy(stats);
            for (int i = (removedIndexes.size() - 1); i > -1; --i) {
                int rmIdx = removedIndexes.get(i);
                distances.remove(rmIdx);
            }
            
            removedIndexes = FeatureMatcher.removeIntensityOutliers(stats, 1.25f); 
            if (removedIndexes.size() < distances.size()) {
                for (int i = (removedIndexes.size() - 1); i > -1; --i) {
                    int rmIdx = removedIndexes.get(i);
                    distances.remove(rmIdx);
                }
                sumSSD = 0;
                sumDist = 0;
                for (int i = 0; i < stats.size(); ++i) {
                    sumSSD += stats.get(i).getSumIntensitySqDiff();
                    sumDist += distances.get(i).doubleValue();
                }
                nEval = stats.size();
            }
            
            if (nEval == 0) {
                continue;
            }

            // distance needs to be adjusted by scale, else the cost prefers
            // small scale solutions
            sumDist /= params.getScale();

            sumSSD /= (double)nEval;
            sumDist /= (double)nEval;

            // add eps to sums so a zero doesn't cancel out the result of the other cost components
            float cost1Norm = 1.f/(float)nEval;
            float cost2Norm = (float)sumSSD + 1;
            float cost3Norm = ((float)sumDist + 0.01f)/(float)tolTransXY2;
            float normalizedCost = cost1Norm * cost2Norm * cost3Norm;

            boolean t1 = (normalizedCost < bestCost);
/*      
float diffRot = AngleUtil.getAngleDifference(params.getRotationInDegrees(), 350);
StringBuilder sb = new StringBuilder();
if ((Math.abs(diffRot) < 20) && (Math.abs(params.getScale() - 1) < 0.15) && 
(Math.abs(params.getTranslationX() - -125) < 30) && 
(Math.abs(params.getTranslationY() - -45) < 30)) {
    sb.append("*** ");
}
sb.append(String.format(" nEval=%d  normCost=%.1f  %s", nEval, normalizedCost, params.toString()));
log.info(sb.toString());
*/
            if (t1 && (nEval > 2)) {
                bestCost = normalizedCost;
                bestParams = params;
                bestStats = stats;
                bestCost1Norm = cost1Norm;
                params.setNumberOfPointsUsed(stats.size());
                bestTolTransXY2 = tolTransXY2;
            }
        }
            
        // calculate the quality array
        if (bestParams != null) {
            int n = bestStats.size();
            
            double[] sumDistSSD = null;
            float sigmaFactor = 1.5f;            
            int nIter = 0;
            int nMaxIter = 5;
            while ((nIter == 0) || (nIter < nMaxIter)) {   
                log.info("before bestStats.size()=" + bestStats.size());
                sumDistSSD = MiscStats.filterStatsForTranslation(bestParams, 
                    bestStats, sigmaFactor);
                log.info("after bestStats.size()=" + bestStats.size());
                if (sumDistSSD != null && !bestStats.isEmpty()) {
                    break;
                }                
                sigmaFactor += 1;
                nIter++;
            }            
            
            if (sumDistSSD != null) {
                                
                TransformationParameters combinedParams =
                    MiscStats.calculateTransformation(1, 1, bestStats,
                        new float[4], false);
                
                if (combinedParams != null) {
                    bestParams = combinedParams;
                    float cost1Norm = 1.f/(float)bestStats.size();
                    float cost2Norm = (float)sumDistSSD[1] + 1;
                    float cost3Norm = ((float)sumDistSSD[0] + 0.01f)/(float)bestTolTransXY2;
                    bestCost = cost1Norm * cost2Norm * cost3Norm;
                }
            }
        }
        
        if (bestParams != null) {
            
            if (binFactor1 != 1 || binFactor2 != 1) {
                for (int i = 0; i < bestStats.size(); ++i) {
                    FeatureComparisonStat stat = bestStats.get(i);
                    stat.setBinFactor1(binFactor1);
                    stat.setBinFactor2(binFactor2);
                }
            }

            MatchingSolution soln = new MatchingSolution(bestParams, bestStats,
                binFactor1, binFactor2);
            return soln;
        }

        return null;
    }
    
    private boolean paramsAreValid(TransformationParameters params) {
        if (params == null) {
            return false;
        }
        if (Float.isNaN(params.getScale())  || Float.isNaN(params.getRotationInRadians())) {
            return false;
        }
        return true;
    }
    
    private double distance(double x1, double y1, double x2, double y2) {

        double diffX = x1 - x2;
        double diffY = y1 - y2;

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }
}
