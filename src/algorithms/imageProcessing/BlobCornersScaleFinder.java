package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.clustering.FixedDistanceGroupFinder;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ScatterPointPlotterPNG;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

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
        SegmentationType type2, boolean useBinned2, int dither) {
        
        boolean filterOutImageBoundaryBlobs = false;
        boolean filterOutZeroPixels = true;
            
        return solveForScale(img1Helper, features1, type1, useBinned1, 
            img2Helper, features2, type2, useBinned2, dither, 
            filterOutImageBoundaryBlobs, filterOutZeroPixels);
    }
    
    public MatchingSolution solveForScale(
        BlobCornerHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobCornerHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2, int dither,
        boolean filterOutImageBoundaryBlobs, 
        boolean filterOutZeroPixels) {

        List<List<CornerRegion>> corners1List = img1Helper.getPerimeterCorners(
            type1, useBinned1);
        List<List<CornerRegion>> corners2List = img2Helper.getPerimeterCorners(
            type2, useBinned2);
        List<Set<PairInt>> blobs1 = img1Helper.imgHelper.getBlobs(type1, 
            useBinned1, filterOutImageBoundaryBlobs);
        List<Set<PairInt>> blobs2 = img2Helper.imgHelper.getBlobs(type2, 
            useBinned2, filterOutImageBoundaryBlobs);
        List<PairIntArray> perimeters1 = img1Helper.imgHelper.getBlobPerimeters(
            type1, useBinned1, filterOutImageBoundaryBlobs, filterOutZeroPixels);
        List<PairIntArray> perimeters2 = img2Helper.imgHelper.getBlobPerimeters(
            type2, useBinned2, filterOutImageBoundaryBlobs, filterOutZeroPixels);

        GreyscaleImage img1 = img1Helper.imgHelper.getGreyscaleImage(useBinned1);
        GreyscaleImage img2 = img2Helper.imgHelper.getGreyscaleImage(useBinned2);

        assert(blobs1.size() == perimeters1.size());
        assert(blobs1.size() == corners1List.size());
        assert(blobs2.size() == perimeters2.size());
        assert(blobs2.size() == corners2List.size());

        float dist = 2.5f;
        for (int i = 0; i < perimeters1.size(); ++i) {
            filterCorners(perimeters1.get(i), corners1List.get(i), dist);
        }
        for (int i = 0; i < perimeters2.size(); ++i) {
            filterCorners(perimeters2.get(i), corners2List.get(i), dist);
        }

        if (true) {
            for (List<CornerRegion> list : corners1List) {
                sortCornersToCCW(list);
            }
            for (List<CornerRegion> list : corners2List) {
                sortCornersToCCW(list);
            }
        }
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

StringBuilder sb = new StringBuilder("xy1:\n");
for (int i = 0; i < xy1.length; ++i) {
    sb.append(String.format("[%2d] (%3d, %3d)\n", i,
        (int)Math.round(xy1[i][0]), (int)Math.round(xy1[i][1])));
}
sb.append("xy2:\n");
for (int i = 0; i < xy2.length; ++i) {
    sb.append(String.format("[%2d] (%3d, %3d)\n", i,
        (int)Math.round(xy2[i][0]), (int)Math.round(xy2[i][1])));
}
System.out.println(sb.toString());
*/
/*
PairInt[] im1Chk = new PairInt[]{
    new PairInt(59, 178), new PairInt(42, 110), new PairInt(27, 105),
    new PairInt(68,  80), new PairInt(25,  55), new PairInt(80, 144)
};
PairInt[] im2Chk = new PairInt[]{
    new PairInt(189, 179), new PairInt(164, 109), new PairInt(164, 109),
    new PairInt(189, 179), new PairInt(154,  59), new PairInt(193, 146)
};
int[] im1ChkIdxs = new int[im1Chk.length];
int[] im2ChkIdxs = new int[im2Chk.length];
for (int i = 0; i < im1Chk.length; ++i) {
    PairInt p = im1Chk[i];
    for (int j = 0; j < xy1.length; ++j) {
        double diffX = p.getX() - xy1[j][0];
        double diffY = p.getY() - xy1[j][1];
        if (Math.abs(diffX) < 20 && Math.abs(diffY) < 20) {
            im1ChkIdxs[i] = j;
            break;
        }
    }
}
for (int i = 0; i < im2Chk.length; ++i) {
    PairInt p = im2Chk[i];
    for (int j = 0; j < xy2.length; ++j) {
        double diffX = p.getX() - xy2[j][0];
        double diffY = p.getY() - xy2[j][1];
        if (Math.abs(diffX) < 20 && Math.abs(diffY) < 20) {
            im2ChkIdxs[i] = j;
            break;
        }
    }
}
sb = new StringBuilder("expeected matches:\n");
for (int i = 0; i < im1ChkIdxs.length; ++i) {
    sb.append(String.format("[%d] to [%d]", im1ChkIdxs[i], im2ChkIdxs[i]));
    sb.append("\n");
}
System.out.println(sb.toString());
*/

        MatchingSolution soln = match(img1Helper, img2Helper,
            features1, features2, img1, img2, corners1List, corners2List,
            useBinned1, useBinned2, dither, filterOutImageBoundaryBlobs);

        return soln;
    }

    private <T extends CornerRegion> MatchingSolution match(
        BlobCornerHelper img1Helper, BlobCornerHelper img2Helper,
        IntensityFeatures features1, IntensityFeatures features2,
        GreyscaleImage img1, GreyscaleImage img2,
        List<List<T>> corners1List, List<List<T>> corners2List,
        boolean useBinned1, boolean useBinned2, int dither,
        boolean filterOutImageBoundaryBlobs) {

        int binFactor1 = img1Helper.imgHelper.getBinFactor(useBinned1);
        int binFactor2 = img2Helper.imgHelper.getBinFactor(useBinned2);

        // filter these corners to remove featureless patches when possible
        List<List<T>> filteredCorners1List = corners1List;
            //filterByLowLimitError(corners1List,
            //img1, features1, binFactor1, useBinned1, img1Helper.debugTag);

        List<List<T>> filteredCorners2List = corners2List;
            //filterByLowLimitError(corners2List,
            //img2, features2, binFactor2, useBinned2, img2Helper.debugTag);

        Map<PairInt, TransformationParameters> trMap
            = new HashMap<PairInt, TransformationParameters>();

        int n1 = filteredCorners1List.size();
        int n2 = filteredCorners2List.size();

        if (n1 == 0 || n2 == 0) {
            return null;
        }

        /*
        -- get best TransformationParameters for each idx1
           (this is at most 40 of them)
        -- consider combining similar parameters to reduce the eval step
        -- evaluate each param against all points.
            -- make 2 eval methods, hopefully the fastest is enough
               -- (1) eval by finding an existing point within tolerance of predicted position
               -- (2) eval by finding best SSD of points within tolerance of predicted and
                      use nEval, SSD and dist to return a normalized cost
        */

        for (int idx1 = 0; idx1 < n1; ++idx1) {

            List<T> corners1 = filteredCorners1List.get(idx1);

            if (corners1.size() < 2) {
                continue;
            }

            /*
            first, see if nEval alone finds the true matches for curve to curve
            */
            int maxNEval = Integer.MIN_VALUE;
            TransformationParameters maxNEvalParams = null;
            Integer maxNEvalIndex2 = null;
            double minCost = Double.MAX_VALUE;
            List<FeatureComparisonStat> minCostStats = null;

            for (int idx2 = 0; idx2 < n2; ++idx2) {

                List<T> corners2 = filteredCorners2List.get(idx2);

                if (corners2.size() < 2) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);

                ClosedCurveCornerMatcher2<T> mapper =
                    new ClosedCurveCornerMatcher2<T>(dither);

                boolean matched = mapper.matchCorners(features1, features2,
                    corners1, corners2, img1, img2, binFactor1, binFactor2);

                if (!matched) {
                    continue;
                }

                // NOTE: if solving for binned, the params are in binned reference frame
                TransformationParameters params = mapper.getSolution();
                int nEval = mapper.getNEval();
                double cost = mapper.getSolutionCost();

                if (nEval < 2) {
                    continue;
                }

                /*
                TODO: consider whether need to further use the SSD error as
                a component in cost as the ability to distinguish a match for
                the region.  A filter was applied above to remove the smallest
                SSD errors because they are nearly featureless patches that
                might be more easily degenerate matches if similar regions are
                present.
                */

                if (cost < minCost) {
                    if ((nEval < 3) && (maxNEval > 4)) {
                        // do not accept if nEval is much lower than maxNEval
                        continue;
                    }
                    maxNEval = nEval;
                    maxNEvalParams = params;
                    maxNEvalIndex2 = index2;
                    minCost = cost;
                    minCostStats = mapper.getSolutionStats();

                } else if ((maxNEval == 2) && (nEval > 3)) {
                    //TODO: may need to revise this
                    double avgCost = (cost + minCost)/2.;
                    if ((Math.abs(cost - avgCost)/(0.1*avgCost)) < 2) {
                        maxNEval = nEval;
                        maxNEvalParams = params;
                        maxNEvalIndex2 = index2;
                        minCost = cost;
                        minCostStats = mapper.getSolutionStats();
                    }
                }
            }

            if (maxNEvalParams != null) {
                trMap.put(new PairInt(idx1, maxNEvalIndex2.intValue()),
                    maxNEvalParams);
            }
        }

        int n2c = 0;
        for (List<T> corners2 : corners2List) {
            n2c += corners2.size();
        }

        List<T> cr2List = new ArrayList<T>();
        int[] xC2 = new int[n2c];
        int[] yC2 = new int[n2c];
        n2c = 0;
        for (List<T> corners2 : corners2List) {
            for (T cr2: corners2) {
                xC2[n2c] = cr2.getX()[cr2.getKMaxIdx()];
                yC2[n2c] = cr2.getY()[cr2.getKMaxIdx()];
                cr2List.add(cr2);
                n2c++;
            }
        }

        NearestPoints np2 = new NearestPoints(xC2, yC2);

        MatchingSolution soln = evaluateForBestUsingFeatures(
            img1Helper, img2Helper, features1, features2, img1, img2,
            trMap, corners1List, cr2List, np2, binFactor1, binFactor2, dither,
            filterOutImageBoundaryBlobs);

        return soln;
    }

    /**
     * given parameters map, evaluate the corner lists and return the solution.
     * Note that if binFactors are not equal to '1', it's assumed that the
     * paramsMap are in the binned reference frame and corrections to that
     * are made for the returned solution.
     * @param <T>
     * @param img1Helper
     * @param img2Helper
     * @param features1
     * @param features2
     * @param img1
     * @param img2
     * @param paramsMap
     * @param corners1List
     * @param corners2List
     * @param np2
     * @param binFactor1
     * @param binFactor2
     * @return
     */
    private <T extends CornerRegion> MatchingSolution evaluateForBestUsingFeatures(
        BlobCornerHelper img1Helper, BlobCornerHelper img2Helper,
        IntensityFeatures features1, IntensityFeatures features2,
        GreyscaleImage img1, GreyscaleImage img2,
        Map<PairInt, TransformationParameters> paramsMap,
        List<List<T>> corners1List, List<T> corners2List,
        NearestPoints np2, int binFactor1, int binFactor2, int dither,
        boolean filterOutImageBoundaryBlobs) {

        int tolTransXY = 5;//10;

        double maxDistance = Math.sqrt(2) * tolTransXY;

        /*
        -- evaluate each param against all points.
           -- (2) eval by finding best SSD of points within tolerance of predicted and
              use nEval, SSD and dist to return a normalized cost
        */

        //List<TransformationParameters> parameterList =
        //    MiscStats.filterToSimilarParamSets(paramsMap, binFactor1, binFactor2);

        List<TransformationParameters> parameterList =
            MiscStats.filterToSimilarParamSets2(paramsMap, binFactor1, binFactor2);

        int n1 = 0;
        for (List<T> corners1 : corners1List) {
            n1 += corners1.size();
        }
        int n2 = corners2List.size();
        int nMaxMatchable = Math.min(n1, n2);

        FeatureMatcher featureMatcher = new FeatureMatcher();
        Transformer transformer = new Transformer();

        final int rotationTolerance = 20;

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

            for (List<T> corners1 : corners1List) {

                for (T cr : corners1) {

                    T crTr = transformer.applyTransformation(params, cr, cos, sin);

                    Set<Integer> indexes2 = np2.findNeighborIndexes(
                        crTr.getX()[crTr.getKMaxIdx()],
                        crTr.getY()[crTr.getKMaxIdx()], tolTransXY2);

                    double bestCostPerIndex = Double.MAX_VALUE;
                    Integer bestCostPerIndexIndex = null;
                    FeatureComparisonStat bestCostPerIndexStat = null;
                    double bestCostPerIndexDist = Double.MAX_VALUE;
                    
                    for (Integer index : indexes2) {

                        int idx2 = index.intValue();

                        T corner2 = corners2List.get(idx2);

                        FeatureComparisonStat compStat =
                            featureMatcher.ditherAndRotateForBestLocation2(
                            features1, features2, cr, corner2, dither2,
                            rotD, rotationTolerance, img1, img2);

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
                            ((float)compStat.getSumIntensitySqDiff() + 1);
                        
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

            // add '1' to sums so a zero doesn't cancel out the result of the other cost components
            float cost1Norm = 1.f/(float)nEval;
            float cost2Norm = (float)sumSSD + 1;
            float cost3Norm = ((float)sumDist + 0.01f)/(float)tolTransXY2;
            float normalizedCost = cost1Norm * cost2Norm * cost3Norm;

            //TODO: cost1Norm's proportion in normalizedCost should be higher

            boolean t1 = (normalizedCost < bestCost);

            if (t1 && (nEval > 2)) {
                bestCost = normalizedCost;
                bestParams = params;
                bestStats = stats;
                bestCost1Norm = cost1Norm;
                params.setNumberOfPointsUsed(stats.size());
                bestTolTransXY2 = tolTransXY2;
            }
        }

        boolean extractMore = true;
            
        if (bestParams != null) {
            
            int n = bestStats.size();
            
            double[] sumDistSSD = null;
            float sigmaFactor = 1.5f;            
            int nIter = 0;
            int nMaxIter = 5;
            while ((nIter == 0) || (nIter < nMaxIter)) {            
                sumDistSSD = MiscStats.filterStatsForTranslation(bestParams, 
                    bestStats, sigmaFactor);
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
                    
                } else {
                    extractMore = false;
                }
            } else {
                extractMore = false;
            }
        }
        
        if (extractMore) {

            extractMore = false;
            
            int img1Width = img1.getWidth();
            int img1Height = img1.getWidth();
            int img2Width = img2.getWidth();
            int img2Height = img2.getHeight();

            int[] qCounts =
                ImageStatisticsHelper.getQuadrantCountsForIntersection(
                bestParams, bestStats,
                img1Width, img1Height, img2Width, img2Height);

            for (int count : qCounts) {
                if (count < 5) {
                    extractMore = true;
                    break;
                }
            }

            if (extractMore) {

                ImageExt imgExt1 = img1Helper.imgHelper.getImage().copyToImageExt();
                ImageExt imgExt2 = img2Helper.imgHelper.getImage().copyToImageExt();
                ImageProcessor imageProcessor = new ImageProcessor();
                int smallestGroupLimit, largestGroupLimit;
                if (binFactor1 == 1) {
                    smallestGroupLimit = img1Helper.imgHelper.getSmallestGroupLimit();
                    largestGroupLimit = img1Helper.imgHelper.getLargestGroupLimit();
                } else {
                    imgExt1 = imageProcessor.binImage(imgExt1, binFactor1);
                    smallestGroupLimit = img1Helper.imgHelper.getSmallestGroupLimitBinned();
                    largestGroupLimit = img1Helper.imgHelper.getLargestGroupLimitBinned();
                }
                if (binFactor2 != 1) {
                    imgExt2 = imageProcessor.binImage(imgExt2, binFactor2);
                }
                ImageSegmentation imageSegmentation = new ImageSegmentation();
                GreyscaleImage imgSeg1Tmp = imageSegmentation.createGreyscale7(imgExt1);
                GreyscaleImage imgSeg2Tmp = imageSegmentation.createGreyscale7(imgExt2);

                BlobCornerFinderForParameters finder = new BlobCornerFinderForParameters();

                //TODO: possible problem here using same group limit size on both images
                List<FeatureComparisonStat> stats2 = finder.extractFeatures(
                    bestParams, img1, img2,
                    imgSeg1Tmp, imgSeg2Tmp,
                    binFactor1, binFactor2,
                    smallestGroupLimit, largestGroupLimit,
                    features1.getRotatedOffsets(), filterOutImageBoundaryBlobs,
                    img1Helper.imgHelper.isInDebugMode(),
                    img1Helper.imgHelper.getDebugTag());
                
                /*
                if there are degenerate matches for either pt1 or pt2, keep
                the one with smallest SSD 
                */
                
                List<Integer> removedIndexes = MiscStats.filterForDegeneracy(bestStats);
                List<Integer> removedIndexes2 = MiscStats.filterForDegeneracy(stats2);
                
                log.info("filtered points for degeneracies");
                
                if (removedIndexes.size() > 1 || removedIndexes2.size() > 1) {

                    double[] sumDistSSD = null;
                    float sigmaFactor = 1.5f;
                    int nIter = 0;
                    int nMaxIter = 5;
                    while ((nIter == 0) || (nIter < nMaxIter)) {
                        sumDistSSD = MiscStats.filterStatsForTranslation(bestParams,
                            bestStats, sigmaFactor);
                        if (sumDistSSD != null && !bestStats.isEmpty()) {
                            break;
                        }
                        sigmaFactor += 1;
                        nIter++;
                    }

                    log.info("filtered points for translation outliers");
                    
                    log.info("calculateTransformation for " + bestStats.size() 
                        + " stats");
                    
                    bestParams = MiscStats.calculateTransformation(binFactor1,
                        binFactor2, bestStats, new float[4], false);
        
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

            MatchingSolution soln = new MatchingSolution(bestParams, bestStats);
            return soln;
        }

        return null;
    }

    protected void filterCorners(PairIntArray curve,
        List<CornerRegion> regions, float dist) {

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

    private static List<Set<Integer>> findCloseCorners(float tolD,
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

    private double distance(double x1, double y1, double x2, double y2) {

        double diffX = x1 - x2;
        double diffY = y1 - y2;

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }

    private void sortCornersToCCW(List<CornerRegion> corners) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        PairIntArray cornerXY = new PairIntArray();
        for (int ii = 0; ii < corners.size(); ++ii) {
            CornerRegion cr = corners.get(ii);
            cornerXY.add(cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]);
        }
        boolean isCW = curveHelper.curveIsOrderedClockwise(cornerXY);
        if (isCW) {
            int n = corners.size();
            if (n > 1) {
                int end = n >> 1;
                for (int ii = 0; ii < end; ii++) {
                    int idx2 = n - ii - 1;
                    CornerRegion swap = corners.get(ii);
                    corners.set(ii, corners.get(idx2));
                    corners.set(idx2, swap);
                }
            }
        }
    }

    /**
     * calculate a low limit for the SSD errors from autocorrelation to use as
     * a high pass filter for the corners.  Note that binFactor is not currently
     * used but may be in the future, allowing adjusted descriptor sizes.
     * @param <T>
     * @param img
     * @param features
     * @param binFactor
     * @param useBinned
     * @return
     */
    private <T extends CornerRegion> List<List<T>> filterByLowLimitError(
        List<List<T>> cornerLists,
        GreyscaleImage img, IntensityFeatures features, int binFactor,
        boolean useBinned, String debugTag) {

        /*
        for normalized descriptors, using near 1.25%
        math.pow((255.* 0.012), 2) * 36 = 340
        math.pow((255.* 0.0125), 2) * 36 = 370

        for unnormalized, may need to determine it per image.  it should usually
        be a higher limit.
        10% is the very high limit of 23,409
        */

 //histograms... lowLimit = 0.6 * 1st peak if y > 1 ?
        float lowLimit = 0;

        List<List<T>> filteredCornerLists = new ArrayList<List<T>>();

        int nTot = 0;

        List<List<Float>> ssdErrors = new ArrayList<List<Float>>();
        for (int i = 0; i < cornerLists.size(); ++i) {

            List<Float> normList = new ArrayList<Float>();

            List<T> corners = cornerLists.get(i);

            for (T cr : corners) {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                int rot0;
                try {
                    rot0 = features.calculate45DegreeOrientation(img, x, y);
                    IntensityDescriptor desc0 = features.extractIntensity(img, x, y, rot0);
                    if (desc0 != null) {
                        float e0 = desc0.sumSquaredError();
                        normList.add(Float.valueOf(e0));
                    }
                } catch (CornerRegion.CornerRegionDegneracyException e) {
                }

                nTot++;
            }
            ssdErrors.add(normList);
        }

        float[] values = new float[nTot];
        int count = 0;
        for (List<Float> list : ssdErrors) {
            for (Float err : list) {
                values[count] = err.floatValue();
                count++;
            }
        }
        float binWidth = 150.f;
        HistogramHolder hist = Histogram.createSimpleHistogram(
            0.f, 4000.f, binWidth, values,
            Errors.populateYErrorsBySqrt(values));
        try {
            hist.plotHistogram("norm SSDErr " + debugTag, debugTag + "_norm_ssd_errors");
        } catch (IOException ex) {
            Logger.getLogger(BlobCornersScaleFinder.class.getName()).log(Level.SEVERE, null, ex);
        }
        for (int i = 0; i < cornerLists.size(); ++i) {

            List<T> filtered = new ArrayList<T>();

            List<T> corners = cornerLists.get(i);

            for (T cr : corners) {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                int rot0;
                try {
                    rot0 = features.calculate45DegreeOrientation(img, x, y);
                    IntensityDescriptor desc0 = features.extractIntensity(img, x, y, rot0);
                    if (desc0 != null) {
                        float e0 = desc0.sumSquaredError();
                        if (e0 > lowLimit) {
                            filtered.add(cr);
                            nTot++;
                        }
                    }
                } catch (CornerRegion.CornerRegionDegneracyException e) {
                }
            }
            filteredCornerLists.add(filtered);
        }

        return filteredCornerLists;
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
    
}
